#ifndef MATCH_FUNC_H
#define MATCH_FUNC_H

int dfc_match(unsigned char *pat, uint32_t *sids, uint32_t sids_cnt);
static inline int SCMemcmp(const unsigned char *s1, const unsigned char *s2, size_t len);
int ac_match( void *id, void *tree, int index, void * data, void * neg_list);

#endif
