#!/bin/sh
# Build SuiteSparse AMD from a source tree and compare it with AMD-lite.
#
# Usage:
#   ./bench/run_compare.sh /path/to/SuiteSparse/AMD
#   AMD_SRC=/path/to/SuiteSparse/AMD ./bench/run_compare.sh
#
# The path may point either at the AMD directory (containing Source/ and
# Include/) or at a SuiteSparse root (containing AMD/Source/ and AMD/Include/).
set -e
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)

AMD_DIR=${1:-${AMD_SRC:-}}
if [ -z "$AMD_DIR" ]; then
    echo "usage: $0 /path/to/SuiteSparse/AMD" >&2
    echo "   or: AMD_SRC=/path/to/SuiteSparse/AMD $0" >&2
    exit 1
fi

# Also accept a SuiteSparse root instead of the AMD subdirectory.
if [ ! -d "$AMD_DIR/Source" ] && [ -d "$AMD_DIR/AMD/Source" ]; then
    AMD_DIR="$AMD_DIR/AMD"
fi

if [ ! -d "$AMD_DIR/Source" ] || [ ! -d "$AMD_DIR/Include" ]; then
    echo "error: '$AMD_DIR' is not an AMD source directory" >&2
    echo "       expected to find Source/ and Include/ there" >&2
    exit 1
fi

SRC="$AMD_DIR/Source"
INC="$AMD_DIR/Include"
OUT="$ROOT/bench/build-original"
CC=${CC:-cc}
CFLAGS=${CFLAGS:--O3 -std=c99 -DNDEBUG}

mkdir -p "$OUT"

# Minimal SuiteSparse_config replacement: AMD only needs malloc/free.
cat > "$OUT/SuiteSparse_config.h" <<'HDR'
#ifndef SUITESPARSE_CONFIG_H
#define SUITESPARSE_CONFIG_H
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#define SUITESPARSE_VER_CODE(main,sub) ((main)*1000+(sub))
#define SUITESPARSE__VERCODE(main,sub,subsub) ((main)*1000000+(sub)*1000+(subsub))
#define SUITESPARSE__VERSION SUITESPARSE__VERCODE(7,11,0)
#define SUITESPARSE_PRINTF(params) printf params
static inline void *SuiteSparse_malloc(size_t n, size_t size)
{
    if (n == 0) n = 1 ;
    return malloc (n * size) ;
}
static inline void SuiteSparse_free(void *p) { free (p) ; }
#endif
HDR

FILES="amd_1.c amd_2.c amd_aat.c amd_control.c amd_defaults.c amd_dump.c
amd_info.c amd_order.c amd_postorder.c amd_post_tree.c amd_preprocess.c
amd_valid.c amd_version.c amd_l1.c amd_l2.c amd_l_aat.c amd_l_control.c
amd_l_defaults.c amd_l_dump.c amd_l_info.c amd_l_order.c amd_l_postorder.c
amd_l_post_tree.c amd_l_preprocess.c amd_l_valid.c"

for f in $FILES; do
    $CC $CFLAGS -w -I"$OUT" -I"$INC" -c "$SRC/$f" -o "$OUT/${f%.c}.o"
done
ar rcs "$OUT/libamd_orig.a" "$OUT"/*.o

$CC $CFLAGS -Wall -Wextra -I"$ROOT/include" \
    "$ROOT/bench/compare_amd.c" "$ROOT/src/amd_lite.c" \
    "$OUT/libamd_orig.a" -lm -o "$OUT/compare_amd"

"$OUT/compare_amd"
