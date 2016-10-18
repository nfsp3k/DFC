/*********************************/
/* Author  - Byungkwon Choi      */
/* Contact - cbkbrad@kaist.ac.kr */
/*********************************/
#ifndef DFC_FRAMEWORK_H
#define DFC_FRAMEWORK_H

/****************************************************/
/*   Func Argument - to define action for matching  */
/****************************************************/
#define ARGUMENT_FOR_MATCH \
                        void (*Match)(unsigned char *, uint32_t *, uint32_t)

#define SEARCH_ARGUMENT \
                        DFC_STRUCTURE *dfc, \
                        unsigned char *buf, \
                        int buflen, \
                        ARGUMENT_FOR_MATCH

#define PROGRE_ARGUMENT \
                        DFC_STRUCTURE *dfc, \
                        unsigned char *buf, \
                        int matches, \
                        BTYPE idx, \
                        BTYPE msk, \
                        ARGUMENT_FOR_MATCH, \
                        const unsigned char *starting_point, \
                        int rest_len

#define VERIFI_ARGUMENT \
                        DFC_STRUCTURE *dfc, \
                        unsigned char *buf, \
                        int matches, \
                        ARGUMENT_FOR_MATCH, \
                        const unsigned char *starting_point

#define PROGRE_PARAMETER \
                        dfc, &buf[i+2], matches, index, mask, Match, buf, buflen-i

#define VERIFI_PARAMETER \
                        dfc, buf, matches, Match, starting_point

#define ACTION_FOR_MATCH \
        Match(mlist->casepatrn, mlist->sids, mlist->sids_size);\
        matches++;

/****************************************************/

#endif

