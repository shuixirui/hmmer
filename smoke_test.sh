#!/bin/bash
# Functional check for jackhmmer after building or changing the code.
# Deliberately NOT a byte-exact diff: optimizations may legitimately change
# output details. We assert only that a working build
#   1. runs to completion (exit 0)
#   2. finds a reasonable number of homologs (>= 5) for the demo query
set -u

# 探索1000: 登录节点(cln01-04)禁止运行程序; 请 ssh test01 等测试节点执行
case $(hostname) in
    cln0[1-4]*|ibcln*)
        echo "ERROR: 禁止在登录节点运行程序。请先 ssh test01 (或 test02-04) 再执行本脚本。" >&2
        exit 1;;
esac

ROOT=$(dirname "$0")
QUERY=${QUERY:-$ROOT/data/query.fasta}
DB=${DB:-$ROOT/data/demodb.fasta}
BIN=$ROOT/src/jackhmmer
CPU=${CPU:-8}
TBL=$(mktemp)
trap 'rm -f "$TBL"' EXIT

FILTERS=${FILTERS:---noali --F1 0.0005 --F2 5e-05 --F3 5e-07}
if ! "$BIN" --cpu "$CPU" $FILTERS -N 3 -E 0.001 --tblout "$TBL" -o /dev/null "$QUERY" "$DB"; then
    echo "FAIL: jackhmmer exited non-zero"; exit 1
fi
hits=$(grep -vc '^#' "$TBL")
if [ "$hits" -lt 5 ]; then
    echo "FAIL: only $hits hits (expected >= 5)"; exit 1
fi
echo "PASS: $hits hits"
