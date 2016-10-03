#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdint.h>
#include <linux/tcp.h>
#include <linux/ip.h>
#include <linux/if_ether.h>
/* DFC */
#include "dfc.h"
#include "dfc_framework.h"
/*-------------------------------------------------------------------------------------*/
#define BUF_LEN				4096
/*-------------------------------------------------------------------------------------*/
#define MAX_PACKET_SIZE		1514
#define FIXED_PACKET_SIZE	1
#define TOTAL_PACKETS		10000000
int packet_size;
unsigned char *packets[TOTAL_PACKETS];
/*-------------------------------------------------------------------------------------*/
# define timersub(a, b, result)									  \
		do {													  \
				(result)->tv_sec = (a)->tv_sec - (b)->tv_sec;		  \
				(result)->tv_usec = (a)->tv_usec - (b)->tv_usec;	  \
				if ((result)->tv_usec < 0) {					      \
						--(result)->tv_sec;						      \
						(result)->tv_usec += 1000000;				  \
				}													  \
		} while (0)
/*-------------------------------------------------------------------------------------*/
int
found (void *id, int index, unsigned char *text)
{
		//printf ("found\n");
		return 0;
}
/*-------------------------------------------------------------------------------------*/
__inline__ uint64_t rdtsc(void) {
		uint32_t lo, hi;
		__asm__ __volatile__ (      // serialize
							  "xorl %%eax,%%eax \n        cpuid"
							  ::: "%rax", "%rbx", "%rcx", "%rdx");
		/* 
		 * We cannot use "=A", since this would use %rax on x86_64 
		 * and return only the lower 32bits of the TSC
		 */
		__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
		return (uint64_t)hi << 32 | lo;
}
/*-------------------------------------------------------------------------------------*/
void
createStrings()
{
		uint64_t seed64;
		uint32_t seed32_1;
		struct ethhdr *eth;
        struct iphdr *ip;
        struct tcphdr *tcp;
        uint32_t rand_val;
		int size, i, j;
		
		for (i = 0; i < TOTAL_PACKETS; i++) {
				packets[i] = calloc(1, packet_size);
				if (packets[i] == NULL) {
						fprintf(stderr, "%s: Can't allocate memory for packets!\n",
								__FUNCTION__);
						exit(EXIT_FAILURE);
				}
		}
		
		for (j = 0; j < TOTAL_PACKETS; j++) {
				seed64 = rdtsc();
				seed32_1 = seed64 & 0xFFFFFFFF;
				size = MAX_PACKET_SIZE;
				
				srand(seed32_1);
				rand_val = lrand48();
				
#ifdef FIXED_PACKET_SIZE
				size = packet_size;
#else
				size = rand_val % MAX_PACKET_SIZE;
#endif
				//fprintf(stderr, "Packet %d is of size: %u\n", j, size);
				for (i = 0; i < packet_size; i += sizeof(rand_val)) {
						memcpy(&packets[j][i], &rand_val, sizeof(rand_val));
				}
				
				/* build an ethernet header */
				eth = (struct ethhdr *)packets[j];
				
				eth->h_dest[0] = 0x00;
				eth->h_dest[1] = 0x1b;
				eth->h_dest[2] = 0x21;
				eth->h_dest[3] = 0x6b;
				eth->h_dest[4] = 0x1d;
				eth->h_dest[5] = 0xd4;
				
				eth->h_source[0] = 0x00;
				eth->h_source[1] = 0x00;
				eth->h_source[2] = 0x00;
				eth->h_source[3] = 0x00;
				eth->h_source[4] = 0x00;
				eth->h_source[5] = 0x01;
				
				eth->h_proto = 0x0800;
				
				/* build an IP header */
				ip = (struct iphdr *)(packets[j] + sizeof(*eth));
				
				ip->version = 4;
				ip->ihl = 5;
				ip->tos = 0;
				ip->tot_len = size - sizeof(*eth);
				ip->id = 0;
				ip->frag_off = 0;
				ip->ttl = 32;
				ip->protocol = 6/*IPPROTO_TCP*/;
				ip->saddr = 0x0A000001;
				ip->daddr = 0x0A000003;
				ip->check = 0;
				
				tcp = (struct tcphdr *)((char *)ip + sizeof(*ip));
				
				rand_val = 80;
				tcp->source = rand_val & 0xFFFF;
				tcp->dest = (rand_val >> 16) & 0xFFFF;
				
				tcp->doff = 5;
				tcp->check = 0;
				if (j % 10000 == 0) printf ("created: %d\r", j);
		}
		printf("created: %d", j);
		printf ("\n");
		
}
/*-------------------------------------------------------------------------------------*/
static int
dfc_rule_match()
{
		/* Do nothing */
		return 0;
}
/*-------------------------------------------------------------------------------------*/
int
main (int argc, char **argv)
{
		FILE *fp;
		char buf[BUF_LEN];
		char *content_start, *content_end;
		DFC_STRUCTURE *dfc;
		int len, id;
		int num_patterns;
		
		if (argc < 3) {
				printf ("Usage: %s [packet_size] [# patterns]\n", argv[0]);
				exit(EXIT_FAILURE);
		}
		packet_size = atoi (argv[1]);
		num_patterns = atoi (argv[2]);
		
		dfc = DFC_New();
		
		/* Read the rule file. */
		fp = fopen ("snort_valid_rules", "r");
		if (fp == NULL) {
				perror ("fopen");
				exit(EXIT_FAILURE);
		}
		
		/* Parse the first 'content' from the rule file. */
		id = 0;
		while (fgets (buf, BUF_LEN, fp) != NULL) {
				if (buf[0] == '#')
						continue;
				content_start = strstr (buf, "content");
				if (content_start == NULL)
						continue;
				content_start += 9;
				content_end = strstr (content_start, "\"; ");
				if (content_end == NULL) {
						printf ("no end of content\n");
						exit(EXIT_FAILURE);
				}
				len = content_end - content_start;
				
				/* Add to acsm state machine. */
				content_start[len] = '\0';
				//printf("Content is: %s\t, Length is %d\n", content_start, len);
				if (id < num_patterns) {
						DFC_AddPattern(dfc, (unsigned char *)content_start,
									   len, 1 /* case-sensitive pattern */,
									   id /* ID used in Snort */,
									   0 /* internal ID */);
				} else
						break;
				id++;
		}
		fclose(fp);
		
		printf ("number of patterns = %d\n", id);
		
		DFC_Compile(dfc);
		
		/*----------------------------------------------------------*/
		printf("Creating Packets\n");
		createStrings();
		printf("Packets Created\n");
		/*----------------------------------------------------------*/
		
		printf ("search start\n");
		int i;
		struct timeval tv_begin, tv_end, tv_result;
		gettimeofday(&tv_begin, NULL);
		
		for (i = 0; i < TOTAL_PACKETS; i++) {
				int rand1;
				struct ethhdr *eth;
				struct iphdr *ip;
				struct tcphdr *tcph;

				rand1 = lrand48() % TOTAL_PACKETS;
				ip = (struct iphdr *)(packets[rand1] + sizeof(*eth));
				tcph = (struct tcphdr *)((unsigned char *)ip + (ip->ihl<<2));

				DFC_Search(dfc, packets[rand1] + 54,
						   ip->tot_len - (ip->ihl<<2) - (tcph->doff<<2),
						   dfc_rule_match, NULL);
#if 0
				fprintf(stderr, "Packet length is: %u\n", ip->tot_len);
#endif
		}
		
		gettimeofday(&tv_end, NULL);
		
		timersub(&tv_end, &tv_begin, &tv_result);
		printf("Time taken to 'content' analyze %u packets: %f secs\n",
			   TOTAL_PACKETS, (float)tv_result.tv_sec + ((float)tv_result.tv_usec)/1000000.0);

		fprintf (stderr, "%d %f\n", num_patterns, (float)tv_result.tv_sec + ((float)tv_result.tv_usec)/1000000.0);

		return 0;
}
