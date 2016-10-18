/* Microbenchmark for DFC */

#include <stdio.h>
#include <string.h>
#include "../src/dfc_framework.h"
#include "../src/dfc.h"

void Print_Result(unsigned char *, uint32_t *, uint32_t);

int main(void){
    char *buf = "This input includes an attack pattern. It might CRASH your machine.";
    char *pat1 = "attack";
    char *pat2 = "crash";
    char *pat3 = "Piolink";
    char *pat4 = "ATTACK";

    // 0. print info
    printf("\n* Text & Patterns Info\n");
    printf(" - (Text) %s\n", buf);
    printf(" - (Pattern) ID: 0, pat: %s, case-sensitive\n", pat1);
    printf(" - (Pattern) ID: 1, pat: %s, case-insensitive\n", pat2);
    printf(" - (Pattern) ID: 2, pat: %s, case-insensitive\n", pat3);
    printf(" - (Pattern) ID: 3, pat: %s, case-insensitive\n", pat4);
    printf("\n");

    // 1. [DFC] Initiate DFC structure
    DFC_STRUCTURE *dfc = DFC_New();

    // 2. [DFC] Add patterns
    DFC_AddPattern(dfc, (unsigned char*)pat1, strlen(pat1), 
                          0/*case-sensitive pattern*/,   0/*Pattern ID*/);
    DFC_AddPattern(dfc, (unsigned char*)pat2, strlen(pat2), 
                          1/*case-insensitive pattern*/, 1/*Pattern ID*/);
    DFC_AddPattern(dfc, (unsigned char*)pat3, strlen(pat3), 
                          1/*case-insensitive pattern*/, 2/*Pattern ID*/);
    DFC_AddPattern(dfc, (unsigned char*)pat4, strlen(pat4), 
                          1/*case-insensitive pattern*/, 3/*Pattern ID*/);

    // 3. [DFC] Build DFC structure
    DFC_Compile(dfc);

    // 4. [DFC] Search
    printf("* Result:\n");
    int res = DFC_Search(dfc, (unsigned char*)buf, strlen(buf), Print_Result);
    printf("\n* Total match count: %d\n", res);

    // 5. [DFC] Free DFC structure
    DFC_FreeStructure(dfc);

    return 0;
}

void Print_Result(unsigned char *pattern, uint32_t *id_list, uint32_t list_size){
    int i;
    printf(" [Matched!] Pattern: %s, IDs:", pattern);

    for(i = 0; i < list_size; i++){
        printf("%u", id_list[i]);
        if(i != list_size - 1)
            printf(", ");
    }

    printf("\n");
}
