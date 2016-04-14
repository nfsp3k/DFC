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
                        int (*Match)(void *, void *, int, void *, void *), \
                        void *pkt_data

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
                        dfc, &buf[i+2], matches, index, mask, Match, pkt_data, buf, buflen-i

#define VERIFI_PARAMETER \
                        dfc, buf, matches, Match, pkt_data, starting_point

#define ACTION_FOR_MATCH \
        Match(NULL, NULL, pid, pkt_data, NULL);\
        matches++;

/****************************************************/

#endif

