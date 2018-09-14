/*
** Copyright (C) 2008 Sourcefire Inc., All Rights Reserved.
**
** Please see the LICENSE file distributed with this source
** code archive for information covering the modification and
** redistribution of this file and binaries built from it.
*/

//#ifdef __cplusplus

//extern "C" {
#ifndef _UTIL_STR_H
#define _UTIL_STR_H
#include <stdio.h>
#include <sys/types.h>

const unsigned char* GetXlateTable(void);
void ConvertCase(unsigned char* dst, unsigned char* src, int len);

#endif /* _UTIL_STR_H */

//}
//#endif
