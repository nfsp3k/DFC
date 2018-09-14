#include "../src/dfc.h"
#include "mpm_engines/ac/acsmx2.h"
#include <stdio.h>

#include <smmintrin.h>

#define SCMEMCMP_BYTES  16
#ifndef likely
#define likely(expr) __builtin_expect(!!(expr), 1)
#endif
#ifndef unlikely
#define unlikely(expr) __builtin_expect(!!(expr), 0)
#endif

static inline int SCMemcmp(const unsigned char *s1, const unsigned char *s2, size_t len)
{
    size_t offset = 0;
    __m128i b1, b2, c;

    do {
        /* apparently we can't just read 16 bytes even though
         * it almost always works fine :) */
        if (likely(len - offset < 16)) {
            return memcmp(s1, s2, len - offset) ? 0 : 1;
        }

        /* do unaligned loads using _mm_loadu_si128. On my Core2 E6600 using
         * _mm_lddqu_si128 was about 2% slower even though it's supposed to
         * be faster. */
        b1 = _mm_loadu_si128((const __m128i *) s1);
        b2 = _mm_loadu_si128((const __m128i *) s2);
        c = _mm_cmpeq_epi8(b1, b2);

        int diff = len - offset;
        if (diff < 16) {
            int rmask = ~(0xFFFFFFFF << diff);

            if ((_mm_movemask_epi8(c) & rmask) != rmask) {
                return 0;
            }
        } else {
            if (_mm_movemask_epi8(c) != 0x0000FFFF) {
                return 0;
            }
        }

        offset += SCMEMCMP_BYTES;
        s1 += SCMEMCMP_BYTES;
        s2 += SCMEMCMP_BYTES;
    } while (len > offset);

    return 1;
}

int dfc_match(unsigned char *pat, uint32_t *sids, uint32_t sids_cnt){
//    int i;
//    for (i = 0; i < sids_cnt; i++){
//        printf("DFC matched: %s (id: %u)\n", pat, sids[i]);
//    }
    return sids_cnt;
}

int ac_match( void *id, void *tree, int index, void * data, void * neg_list)
{
    ACSM_PATTERN2 *mlist = (ACSM_PATTERN2 *)id;
    int matches = 0;

    if(!mlist->nocase){
        int success = SCMemcmp((const unsigned char*)data-index, mlist->casepatrn, index);
        if(success){
            matches++;
            //printf("AC matched: %s (id: %d)\n", mlist->casepatrn, mlist->iid);
        }
    }else{
        matches++;
        //printf("AC matched: %s (id: %d)\n", mlist->casepatrn, mlist->iid);
    }

    if(mlist->next)
        matches += ac_match(mlist->next, tree, mlist->next->n, data, neg_list);

    return matches;
}


