#ifndef ENTRY_H
#define ENTRY_H
#include "../src/dfc.h"
#include "mpm_engines/ac/acsmx2.h"

struct Entry{
    char * ruleFileName;
    char ** rules;
    char ** contents;
    int * noCase; // 0: case-sensitive, 1: case-insensitive

    int numOfRules;
    int numOfContents;

    char **testPayload;
    int payloadLen;
    int numOfPayload;

    ACSM_STRUCT2 *ac_struct;
    DFC_STRUCTURE *dfc_struct;
};


#endif
