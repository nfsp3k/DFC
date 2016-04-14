/*********************************/
/* Author  - Byungkwon Choi      */
/* Contact - cbkbrad@kaist.ac.kr */
/*********************************/
#ifndef DFC_H
#define DFC_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>
#include <nmmintrin.h>
#include <emmintrin.h>

#include "dfc_framework.h"

/****************************************************/
/*         Parameters: DF size, CT size             */
/****************************************************/
#define DF_SIZE 0x10000
#define DF_SIZE_REAL 0x2000

#define CT_TYPE1_PID_CNT_MAX 100
#define CT1_TABLE_SIZE 256
#define CT2_TABLE_SIZE 0x1000
#define CT3_TABLE_SIZE 0x1000
#define CT4_TABLE_SIZE 0x20000
#define CT8_TABLE_SIZE 0x20000
#define RECURSIVE_CT_SIZE 4096

#define PID_TYPE        uint32_t
#define PID_CNT_TYPE    uint32_t
#define BUC_CNT_TYPE    uint32_t

#define DTYPE           uint16_t
#define BTYPE           register uint16_t
/****************************************************/



/****************************************************/
/*        You don't need to care about this         */
/****************************************************/
#define BINDEX(x) ((x) >> 3)
#define BMASK(x)  (1 << ((x) & 0x7))

#define DF_MASK (DF_SIZE - 1)

#define CT2_TABLE_SIZE_MASK (CT2_TABLE_SIZE-1)
#define CT3_TABLE_SIZE_MASK (CT3_TABLE_SIZE-1)
#define CT4_TABLE_SIZE_MASK (CT4_TABLE_SIZE-1)
#define CT8_TABLE_SIZE_MASK (CT8_TABLE_SIZE-1)

#ifndef likely
#define likely(expr) __builtin_expect(!!(expr), 1)
#endif
#ifndef unlikely
#define unlikely(expr) __builtin_expect(!!(expr), 0)
#endif

#define MEMASSERT_DFC(p,s) if(!p){printf("DFC-No Memory: %s!\n",s);}
/****************************************************/


/****************************************************/

/* Compact Table Structures */
/* Compact Table (CT1) */
typedef struct pid_list_{
    PID_TYPE pid[CT_TYPE1_PID_CNT_MAX];
    uint16_t cnt;
} CT_Type_1;

/****************************************************/
/*                For New designed CT2              */
/****************************************************/
typedef struct CT_Type_2_2B_Array_{
	uint16_t pat;     // 2B pattern
    PID_CNT_TYPE cnt;     // Number of PIDs
    PID_TYPE *pid;	  // list of PIDs
} CT_Type_2_2B_Array;

/* Compact Table (CT2) */
typedef struct CT_Type_2_2B_{
    BUC_CNT_TYPE cnt;
    CT_Type_2_2B_Array *array;
} CT_Type_2_2B;
/****************************************************/

typedef struct CT_Type_2_Array_{
	uint32_t pat;     // Maximum 4B pattern
    PID_CNT_TYPE cnt;     // Number of PIDs
    PID_TYPE *pid;	  // list of PIDs
	uint8_t *DirectFilter;
	CT_Type_2_2B *CompactTable;
} CT_Type_2_Array;

/* Compact Table (CT2) */
typedef struct CT_Type_2_{
    BUC_CNT_TYPE cnt;
    CT_Type_2_Array *array;
} CT_Type_2;

/****************************************************/
/*                For New designed CT8              */
/****************************************************/
typedef struct CT_Type_2_8B_Array_{
	uint64_t pat;     // 8B pattern
    PID_CNT_TYPE cnt;     // Number of PIDs
    PID_TYPE *pid;	  // list of PIDs
	uint8_t *DirectFilter;
	CT_Type_2_2B *CompactTable;
} CT_Type_2_8B_Array;

/* Compact Table (CT2) */
typedef struct CT_Type_2_8B_{
    BUC_CNT_TYPE cnt;
    CT_Type_2_8B_Array *array;
} CT_Type_2_8B;
/****************************************************/

typedef struct _dfc_pattern
{
    struct _dfc_pattern *next;

    unsigned char       *patrn;     // upper case pattern
    unsigned char       *casepatrn; // original pattern
    int                  n;         // Patternlength
    int                  nocase;    // Flag for case-sensitivity. (0: case-sensitive pattern, 1: opposite)

    uint32_t             sids_size;
    PID_TYPE            *sids;      // external id (unique)
    PID_TYPE             pid;       // external id ()
    PID_TYPE             iid;       // internal id (used in DFC library only)

} DFC_PATTERN;


typedef struct {
    DFC_PATTERN   ** init_hash; // To cull duplicate patterns
    DFC_PATTERN    * dfcPatterns;
    DFC_PATTERN   ** dfcMatchList;

    int          numPatterns;
    PID_TYPE     max_pid;

	/* Direct Filter (DF1) for all patterns */
    uint8_t DirectFilter1[DF_SIZE_REAL];

    uint8_t cDF0[256];
    uint8_t cDF1[DF_SIZE_REAL];
    uint8_t cDF2[DF_SIZE_REAL];

    uint8_t ADD_DF_4_plus[DF_SIZE_REAL];
    uint8_t ADD_DF_4_1[DF_SIZE_REAL];

    uint8_t ADD_DF_8_1[DF_SIZE_REAL];
    uint8_t ADD_DF_8_2[DF_SIZE_REAL];

	

    /* Compact Table (CT1) for 1B patterns */
    CT_Type_1 CompactTable1[CT1_TABLE_SIZE];

    /* Compact Table (CT2) for 2B patterns */
    CT_Type_2 CompactTable2[CT2_TABLE_SIZE];

    /* Compact Table (CT4) for 4B ~ 7B patterns */
    CT_Type_2 CompactTable4[CT4_TABLE_SIZE];

    /* Compact Table (CT8) for 8B ~ patterns */
    CT_Type_2_8B CompactTable8[CT8_TABLE_SIZE];

}DFC_STRUCTURE;

/****************************************************/
typedef enum _dfcMemoryType
{
    DFC_MEMORY_TYPE__NONE = 0,
    DFC_MEMORY_TYPE__DFC,
    DFC_MEMORY_TYPE__PATTERN,
    DFC_MEMORY_TYPE__CT1,
    DFC_MEMORY_TYPE__CT2,
    DFC_MEMORY_TYPE__CT3,
    DFC_MEMORY_TYPE__CT4,
    DFC_MEMORY_TYPE__CT8
} dfcMemoryType;


typedef enum _dfcDataType
{
    DFC_NONE = 0,
    DFC_PID_TYPE,
    DFC_CT_Type_2_Array,
    DFC_CT_Type_2_2B_Array,
    DFC_CT_Type_2_8B_Array
} dfcDataType;
/****************************************************/

/****************************************************/
DFC_STRUCTURE * DFC_New (void);
int DFC_AddPattern (DFC_STRUCTURE *dfc, unsigned char *pat, int n, int nocase, PID_TYPE pid, PID_TYPE sid);
int DFC_Compile(DFC_STRUCTURE *dfc);

int DFC_Search(SEARCH_ARGUMENT);
void DFC_PrintInfo(DFC_STRUCTURE* dfc); // Print info
void DFC_FreeStructure(DFC_STRUCTURE *dfc);
/****************************************************/

#endif
