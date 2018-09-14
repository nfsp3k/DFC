#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <dirent.h>
#include "entry.h"
#include "../src/dfc.h"
#include "mpm_engines/ac/acsmx2.h"
#include "match_func.h"
#include <float.h>

#define MAX_BUF 3000
#define INPUT_SIZE 1000000 // byte
#define RANDOM_TEXT 0 // 0: no random, 1: random
#define TEXT_SIZE 1500 // MUST be multiples of 4

struct Entry * loadAndGen(int argc, char **argv, int * num);
struct Entry * loadRules(int argc, char **argv, int *num);
struct Entry * loadRulesDir(char **argv, int *num);
struct Entry * loadRulesRegular(char *filePath, int *num);
void parseContentOptions(struct Entry *lst, int *num);
void genText(struct Entry *lst, int *num, int random);
void runMPM(struct Entry *lst, int num);
void constructMPMStructs(struct Entry *lst, int num);
void detection(struct Entry *lst, int num);
void fast_srand(int seed);

static unsigned int g_seed;

int main(int argc, char **argv){
    // 1. Loads rules and generate test payload.
    int num=0;
    struct Entry *inputList = loadAndGen(argc, argv, &num);

    if (inputList == NULL){
        printf("Failed to generate test payload.\n");
        return -1;
    }

    // 2. Runs mpm over test payloads
    runMPM(inputList, num);

    return 0;
}

void runMPM(struct Entry *lst, int num){
    struct timeval tv_begin, tv_end, tv_result;

    // 1. constructs mpm structures for AC and DFC
    printf("Constructing MPM structures ... \n");
    gettimeofday(&tv_begin, NULL);
    constructMPMStructs(lst, num);
    gettimeofday(&tv_end, NULL);
    timersub(&tv_end, &tv_begin, &tv_result);
    printf("Done. %.2f seconds\n", (float)tv_result.tv_sec + ((float)tv_result.tv_usec)/1000000.0);

    // 2. performs detect processes for AC and DFC
    detection(lst, num);
}

float getGbps(int total_length_bit, struct timeval tv){
    if (((float)tv.tv_sec + ((float)tv.tv_usec)/1000000.0) == 0)
        return 0;
    return total_length_bit/((float)tv.tv_sec + ((float)tv.tv_usec)/1000000.0)/1000000000;
}

void detection(struct Entry *lst, int num){
    int i, j;
    int m_ac, m_dfc;
    float max_improvement = 0;
    float min_improvement = FLT_MAX;
    struct timeval tv_begin, tv_end, tv_ac, tv_dfc;

    printf("\nAC vs. DFC\n");
    for (i = 0; i < num; i++){
        //printf("Text: %s\n\n", lst[i].testPayload[0]);
        if (lst[i].numOfPayload == 0)
            continue;

        m_ac = 0;
        m_dfc = 0;

        // 1. run AC
        gettimeofday(&tv_begin, NULL);
        for (j = 0; j < lst[i].numOfPayload; j++){
        //for (j = 0; j < 1; j++){
            int start_state = 0;
            m_ac += acsmSearch2(lst[i].ac_struct, (unsigned char*)lst[i].testPayload[j], 
                                                lst[i].payloadLen, ac_match, NULL, &start_state);
        }
        gettimeofday(&tv_end, NULL);
        timersub(&tv_end, &tv_begin, &tv_ac);

        // 2. run DFC
        gettimeofday(&tv_begin, NULL);
        for (j = 0; j < lst[i].numOfPayload; j++){
        //for (j = 0; j < 1; j++){
            m_dfc += DFC_Search(lst[i].dfc_struct, (unsigned char*)lst[i].testPayload[j], 
                                                lst[i].payloadLen, dfc_match);
        }
        gettimeofday(&tv_end, NULL);
        timersub(&tv_end, &tv_begin, &tv_dfc);


        // Print result
        int total_length_bit = lst[i].numOfPayload * lst[i].payloadLen * 8;
        float ac_gbps = getGbps(total_length_bit, tv_ac);
        float dfc_gbps = getGbps(total_length_bit, tv_dfc);
        if (ac_gbps != 0 && max_improvement < dfc_gbps/ac_gbps)
            max_improvement = dfc_gbps/ac_gbps;
        if (ac_gbps != 0 && min_improvement > dfc_gbps/ac_gbps)
            min_improvement = dfc_gbps/ac_gbps;

        printf("[%s][%s][%2.1fx] %.3f Gbps - %10d  vs  %.3f Gbps - %10d  (%5d, %s)\n", 
                        (m_ac == m_dfc)? "SUCC" : "FAIL", 
                        (dfc_gbps >= ac_gbps)? "FAST": "SLOW",
                        (ac_gbps != 0)? dfc_gbps/ac_gbps : 0.0,
                        ac_gbps, m_ac, dfc_gbps, m_dfc,
                        lst[i].numOfContents, lst[i].ruleFileName);
    }

    printf("Min-Max: %.1f-%.1f\n", min_improvement, max_improvement);
}

void constructMPMStructs(struct Entry *lst, int num){
    int i, j;
    for (i = 0; i < num; i++){
        printf("(%2d/%2d) %5d %s ... ", i+1, num, lst[i].numOfContents, lst[i].ruleFileName);
        fflush(stdout);

        // 1-1. constructs AC structure
        lst[i].ac_struct = acsmNew2(NULL, NULL, NULL);
        lst[i].ac_struct->acsmFormat = ACF_FULL;
        lst[i].ac_struct->acsmFSA     = FSA_DFA;
        lst[i].ac_struct->acsmAlphabetSize      = 256;
        lst[i].ac_struct->acsmSparseMaxRowNodes = 256;
        lst[i].ac_struct->acsmSparseMaxZcnt     = 10;
        acsmCompressStates((ACSM_STRUCT2*)lst[i].ac_struct, 0);

        // 1-2. constructs DFC structure
        lst[i].dfc_struct = DFC_New();

        //printf("\n");
        for (j = 0; j < lst[i].numOfContents; j++){
            // 2-1. Add pattern to AC structure
            acsmAddPattern2(lst[i].ac_struct, (unsigned char*)lst[i].contents[j], strlen(lst[i].contents[j]), lst[i].noCase[j], 0, 0, 0, NULL, j);
            
            // 2-2. Add pattern to DFC structure
            DFC_AddPattern(lst[i].dfc_struct, (unsigned char*)lst[i].contents[j], strlen(lst[i].contents[j]), lst[i].noCase[j], j);
        }


        // 3-1. Compile AC structure
        acsmCompile2 (lst[i].ac_struct, NULL, NULL);

        // 3-2. Compile DFC structure
        DFC_Compile(lst[i].dfc_struct);
        //DFC_PrintInfo(lst[i].dfc_struct);

        printf("Done. (# of rules: %d)\n", lst[i].numOfRules);
    }

    return;
}


struct Entry * loadAndGen(int argc, char **argv, int *num){
    // 0. sets random seed
    if (RANDOM_TEXT){
        struct timeval tv;
        gettimeofday(&tv, NULL);
        fast_srand(tv.tv_usec % UINT32_MAX);
    }

    // 1. loads rules
    printf("Loading rule files ... ");
    struct Entry *inputList = loadRules(argc, argv, num);
    printf("Done. (# rule files: %d)\n", *num);

    // 2. parses content options
    printf("Parsing rules ... ");
    parseContentOptions(inputList, num);
    printf("Done.\n");

    // 3. generates test payloads
    printf("Generating test payloads ... ");
    genText(inputList, num, RANDOM_TEXT);
    printf("Done.\n");

    return inputList;
}

void fast_srand(int seed){
    g_seed = seed;
}

char fast_rand(void){
    g_seed = (214013*g_seed+2531011);
    return (char)(((g_seed>>16)&0x7FFF) % 256);
}

void genText(struct Entry *lst, int *num, int random){
    int i,j;

    for (i = 0; i < *num; i++){
        if (lst[i].numOfContents == 0)
            continue;

        char *buf;
        int c = 0;

        if (random){
            lst[i].payloadLen = TEXT_SIZE;
            buf = (char*)malloc(sizeof(char) * lst[i].payloadLen + 1);
            for (j = 0; j < TEXT_SIZE; j++){
                buf[j] = fast_rand();
            }
            buf[lst[i].payloadLen] = '\0';
        }else{
            buf = (char*)malloc(sizeof(char) * lst[i].payloadLen + 1);
            for (j = 0; j < lst[i].numOfContents; j++){
                strcpy(buf + c, lst[i].contents[j]);
                c += strlen(lst[i].contents[j]);
            }
            buf[lst[i].payloadLen] = '\0';
        }


        c = INPUT_SIZE / lst[i].payloadLen;
        lst[i].testPayload = (char**)malloc(sizeof(char*) * c);

        for (j = 0; j < c; j++){
            lst[i].testPayload[j] = (char*)malloc(sizeof(char) * lst[i].payloadLen + 1);
            strcpy(lst[i].testPayload[j], buf);
            lst[i].testPayload[j][lst[i].payloadLen] = '\0';
        }
        lst[i].numOfPayload = c;

        free(buf);
    }
}

void parseContentOptions(struct Entry *lst, int *num){
    int i, j;
    
    for (i = 0; i < *num; i++){
        for (j = 0; j < lst[i].numOfRules; j++){
            char *p1 = lst[i].rules[j];
            char *p2, *tmp;

            while(1){
                p2 = strstr(p1, "content:\"");

                if (p2 == NULL)
                    break;

                p2 += 9;
                p1 = strstr(p2, "\";");

                lst[i].contents=(char**)realloc(lst[i].contents, 
                                    sizeof(char*)*(lst[i].numOfContents+1));
                lst[i].noCase=(int*)realloc(lst[i].noCase, 
                                    sizeof(int)*(lst[i].numOfContents+1));

                // copy content
                lst[i].contents[lst[i].numOfContents] = (char*)malloc(sizeof(char) * (p1-p2+1));
                strncpy(lst[i].contents[lst[i].numOfContents], p2, p1-p2);
                lst[i].contents[lst[i].numOfContents][p1-p2] = '\0';

                // set noCase
                if(strncmp("nocase", p1+3, 6) == 0)
                    lst[i].noCase[lst[i].numOfContents] = 1;
                else
                    lst[i].noCase[lst[i].numOfContents] = 0;

                lst[i].numOfContents++;
                lst[i].payloadLen += p1-p2;
            }
        }
    }
}

int isDirectory(const char *path)
{
    struct stat path_stat;
    stat(path, &path_stat);
    return S_ISDIR(path_stat.st_mode);
}

struct Entry * loadRules(int argc, char **argv, int *num){
    if (isDirectory(argv[1]))
        return loadRulesDir(argv, num);
    else
        return loadRulesRegular(argv[1], num);
}

struct Entry * loadRulesDir(char **argv, int *num){
    DIR *dir;
    struct dirent *ent;
    char buf[MAX_BUF];

    if((dir = opendir(argv[1])) != NULL){
        struct Entry * res = NULL;
        while((ent = readdir(dir)) != NULL){
            if (ent->d_name[0] == '.')
                continue;
            if (ent->d_name[0] == '.' && ent->d_name[1] == '.')
                continue;

            sprintf(buf, "%s%s", argv[1], ent->d_name);
            struct Entry *e = loadRulesRegular(buf, num);

            if (e == NULL)
                continue;

            res = (struct Entry *)realloc(res, sizeof(struct Entry) * (*num));
            res[(*num)-1] = *e;
        }
        return res;
    }else
        return NULL;
}

struct Entry * loadRulesRegular(char *filePath, int *num){
    FILE *fp;
    char buf[MAX_BUF];
    struct Entry *e = (struct Entry*)calloc(1, sizeof(struct Entry));

    if ((fp = fopen(filePath, "r")) == NULL){
        perror("No such file.");
        return NULL;
    }

    // sets file name
    e->ruleFileName = (char*)malloc(strlen(filePath)+1);
    strcpy(e->ruleFileName, filePath);
    
    // loads rules
    int numOfRules = 0;
    while(fgets(buf, sizeof(buf), fp) != NULL){
        buf[strlen(buf)-1] = '\0';
        
        if (buf[0] == '#' || buf[0] == '\0')
            continue;

        e->rules = (char**)realloc(e->rules, (numOfRules+1)*sizeof(char*));
        e->rules[numOfRules] = (char*)malloc((strlen(buf) + 1)*sizeof(char));
        strcpy(e->rules[numOfRules], buf);

        numOfRules++;
    }

    e->numOfRules = numOfRules;

    fclose(fp);
    *num = *num + 1;
    return e;
}

void showsResult(int mAC, int mDFC){
}


/* Utils */
void frees(char ** arr, int size){
    int i =0;
    for (i = 0; i < size; i++)
        free(arr[i]);

    free(arr);
    return;
}
