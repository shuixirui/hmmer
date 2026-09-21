#!/bin/bash
#SBATCH -J msa_baseline
#SBATCH -p cnall
#SBATCH -N 1
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=56

# 探索1000 基线作业脚本。目录布局：
#   ~/WORK/hmmer         本仓库（结果写在 hmmer/results/）
#   ~/WORK/preliminary   发放的数据（query/ + database.fasta）
# 用法（在 ~/WORK 下提交，stdout/stderr 写到提交目录）：
#   cd ~/WORK && sbatch hmmer/job_baseline.slurm
# 基线耗时参考：100 条 query 逐条串行约 5–6 小时（另有 5 次预热），
# 若分区有时长上限请用 #SBATCH -t 设足时间，或拆成多个作业。
REPO=${REPO:-$HOME/WORK/hmmer}                  # 代码仓库（必须在 WORK：home 在计算节点只读）
DATA=${DATA:-$HOME/WORK/preliminary}            # 数据目录（必须在 WORK）
OUT=${OUT:-$REPO/results}                       # 结果目录（命中表/MSA/计时/评分，与数据目录分开）

# 若 jackhmmer 是用 module 中的 gcc 编译的，请在此加载同一环境，例如：
# module load compilers/gcc/v12.2.0

# CPU: jackhmmer 的工作线程数（基线口径为 8；节点共 56 核，可自行调整）
CPU=${CPU:-8} bash "$REPO/run_queries.sh" "$DATA" "$OUT"
