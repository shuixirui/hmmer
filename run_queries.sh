#!/bin/bash
# Batch-run jackhmmer for every query in a release data directory.
#
# Usage:  run_queries.sh [data_dir] [out_dir]
#   [data_dir]  contains query/*.fasta and database.fasta
#               default: ../preliminary next to the repo, i.e. the layout
#                 ~/WORK/hmmer         (this repo)
#                 ~/WORK/preliminary   (released data)
#   [out_dir]   default: <repo>/results  (next to this script, NOT inside the
#               data directory -- keeps results out of the released data tree)
#               may also be set via the OUT env var
#
# Outputs:
#   out_dir/tbl/<sid>.tbl   per-query hit table (--tblout, scoring input)
#   out_dir/msa/<sid>.sto   per-query MSA (Stockholm)
#   out_dir/timing.tsv      sid, wall seconds, hit count (per query, for reference)
#   out_dir/wall_time.tsv   OFFICIAL time: wall clock from the end of warm-up
#                           until the last query finishes
#
# Tunables (env):  JACKHMMER, CPU (default 8), ITERS (3), EVALUE (0.001),
#   FILTERS (AlphaFold3 的三级过滤阈值), WARMUP (默认 5),
#   FORCE=1 to re-run queries whose .tbl already exists (default: skip).
#
# 预热：正式计时前对第一条 query 重复跑 WARMUP 次，把数据库读进 page cache。
# 冷缓存下单条耗时会高 70%（实测 305 s vs 180 s），不预热首条计时不可比。
# 预热结果全部丢弃，不计入 timing.tsv。WARMUP=0 可关闭。
#
# NOTE (探索1000): a full run is heavy (100 queries, hours) — submit it via
# sbatch (see job_baseline.slurm); on test nodes (test01-04) only run small
# subsets. Login nodes (cln01-04) refuse to run. Results are written under the
# repo, so the repo itself MUST live under ~/WORK: home is READ-ONLY on compute
# nodes. (Alternatively point OUT / [out_dir] somewhere writable.)
set -euo pipefail

case $(hostname) in
    cln0[1-4]*|ibcln*)
        echo "ERROR: 禁止在登录节点运行程序。小规模测试请 ssh test01；完整运行请用 sbatch 提交（见 job_baseline.slurm）。" >&2
        exit 1;;
esac

REPO=$(cd "$(dirname "$0")" && pwd)
DATA=${1:-${DATA:-$(dirname "$REPO")/preliminary}}
OUT=${2:-${OUT:-$REPO/results}}
JACKHMMER=${JACKHMMER:-$REPO/src/jackhmmer}
CPU=${CPU:-8}
ITERS=${ITERS:-3}
EVALUE=${EVALUE:-0.001}
FILTERS=${FILTERS:---noali --F1 0.0005 --F2 5e-05 --F3 5e-07}
WARMUP=${WARMUP:-5}
FORCE=${FORCE:-0}

DB=$DATA/database.fasta
[ -f "$DB" ] || { echo "ERROR: $DB not found" >&2; exit 1; }
[ -d "$DATA/query" ] || { echo "ERROR: $DATA/query not found" >&2; exit 1; }
[ -x "$JACKHMMER" ] || { echo "ERROR: jackhmmer not executable: $JACKHMMER" >&2; exit 1; }

if ! mkdir -p "$OUT/tbl" "$OUT/msa" 2>/dev/null || ! touch "$OUT/.write_test" 2>/dev/null; then
    echo "ERROR: 输出目录不可写: $OUT" >&2
    echo "提示: 结果默认写到仓库目录下的 results/。探索1000 的 home 在计算节点上只读，" >&2
    echo "      请把仓库克隆到 ~/WORK 下（数据同样放 ~/WORK），或显式指定可写的输出目录:" >&2
    echo "      bash run_queries.sh <data_dir> ~/WORK/msa/results" >&2
    exit 1
fi
rm -f "$OUT/.write_test"
TIMING=$OUT/timing.tsv
[ -f "$TIMING" ] && [ "$FORCE" = 0 ] || echo -e "sid\tseconds\thits" > "$TIMING"

total=$(ls "$DATA"/query/*.fasta | wc -l)

# ---- 预热：把数据库读进 page cache，结果丢弃 ----------------------------
if [ "$WARMUP" -gt 0 ]; then
    first=$(ls "$DATA"/query/*.fasta | head -1)
    echo "预热: $(basename "$first" .fasta) x $WARMUP 次（不计入计时）"
    for w in $(seq 1 "$WARMUP"); do
        w0=$(date +%s.%N)
        "$JACKHMMER" --cpu "$CPU" $FILTERS -N "$ITERS" -E "$EVALUE" \
            --tblout /dev/null -o /dev/null "$first" "$DB"
        w1=$(date +%s.%N)
        echo "  预热 $w/$WARMUP: $(echo "$w1 $w0" | awk '{printf "%.1f", $1-$2}')s"
    done
    echo "预热完成，开始正式计时"
fi

# ---- 正式计时：从预热结束开始，到最后一条 query 跑完为止（墙钟）----------
i=0; ran=0; skipped=0
t_all0=$(date +%s.%N)
for qf in "$DATA"/query/*.fasta; do
    sid=$(basename "$qf" .fasta)
    i=$((i + 1))
    if [ "$FORCE" = 0 ] && [ -s "$OUT/tbl/$sid.tbl" ]; then
        echo "[$i/$total] $sid: exists, skip"
        skipped=$((skipped + 1))
        continue
    fi
    ran=$((ran + 1))
    t0=$(date +%s.%N)
    "$JACKHMMER" --cpu "$CPU" $FILTERS -N "$ITERS" -E "$EVALUE" \
        --tblout "$OUT/tbl/$sid.tbl" -A "$OUT/msa/$sid.sto" \
        -o /dev/null "$qf" "$DB"
    t1=$(date +%s.%N)
    secs=$(echo "$t1 $t0" | awk '{printf "%.1f", $1-$2}')
    hits=$(grep -vc '^#' "$OUT/tbl/$sid.tbl" || true)
    echo -e "$sid\t$secs\t$hits" >> "$TIMING"
    echo "[$i/$total] $sid: ${secs}s, $hits hits"
done

t_all1=$(date +%s.%N)
wall=$(echo "$t_all1 $t_all0" | awk '{printf "%.1f", $1-$2}')
echo -e "seconds\tqueries_run\tskipped\n$wall\t$ran\t$skipped" > "$OUT/wall_time.tsv"
echo "----------------------------------------"
echo "正式计时（预热结束 → 全部跑完，墙钟）: ${wall} s   [$ran 条运行, $skipped 条跳过]"
if [ "$skipped" -gt 0 ]; then
    echo "WARNING: 有 $skipped 条沿用了已有结果，本次计时不含它们，不能作为正式成绩。" >&2
    echo "         正式测试请清空输出目录或设 FORCE=1。" >&2
fi
awk -F'\t' 'NR>1{s+=$2; n++} END{if(n) printf "逐条耗时之和 %.1f s，平均 %.1f s/条（仅供参考）\n", s, s/n}' "$TIMING"

# 自动评分：score.py 从数据目录里 query 与 database.fasta 的 header 读取 SCOP 标签，
# 输出每条 query 的召回率/精确度/ROC1/ROC10 与平均值（ROC10 为主指标）。
python3 "$REPO/score.py" --data "$DATA" --results "$OUT" || echo "warning: scoring step failed" >&2
