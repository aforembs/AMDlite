/* amd_lite.c: lightweight AMD ordering, int32_t and int64_t entry points.
 *
 * AMD, Copyright (c) 1996-2024, Timothy A. Davis, Patrick R. Amestoy, and
 * Iain S. Duff.  All Rights Reserved.  SPDX-License-Identifier: BSD-3-clause
 */

#include "amd_lite.h"
#include <stddef.h>
#include <stdint.h>

void amd_lite_defaults (amd_lite_options *options)
{
    if (options == NULL) return ;
    options->dense = AMD_LITE_DEFAULT_DENSE ;
    options->aggressive = AMD_LITE_DEFAULT_AGGRESSIVE ;
    options->format = AMD_LITE_CSR ;
}

#define AMDL_JOIN2(a,b) a ## b
#define AMDL_JOIN(a,b) AMDL_JOIN2(a,b)
#define AMDL_FN(name) AMDL_JOIN(AMDL_PREFIX, name)

/* ------------------------------------------------------------------------- */
/* int32_t version                                                           */
/* ------------------------------------------------------------------------- */

#define AMDL_IDX      int32_t
#define AMDL_UINT     uint32_t
#define AMDL_IDX_MAX  INT32_MAX
#define AMDL_PREFIX   amd_lite_
#define AMDL_PUBLIC   amd_lite_order

#include "amd_lite_impl.h"

#undef AMDL_IDX
#undef AMDL_UINT
#undef AMDL_IDX_MAX
#undef AMDL_PREFIX
#undef AMDL_PUBLIC

/* ------------------------------------------------------------------------- */
/* int64_t version                                                           */
/* ------------------------------------------------------------------------- */

#define AMDL_IDX      int64_t
#define AMDL_UINT     uint64_t
#define AMDL_IDX_MAX  INT64_MAX
#define AMDL_PREFIX   amd_lite_l_
#define AMDL_PUBLIC   amd_lite_l_order

#include "amd_lite_impl.h"

#undef AMDL_IDX
#undef AMDL_UINT
#undef AMDL_IDX_MAX
#undef AMDL_PREFIX
#undef AMDL_PUBLIC

#undef AMDL_FN
#undef AMDL_JOIN
#undef AMDL_JOIN2
