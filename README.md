# MSA 搜索加速赛题 —— jackhmmer 基线代码

本仓库是为超算比赛裁剪的 [HMMER 3.4](http://hmmer.org)（Eddy/Rivas 实验室，BSD 协议）。
原版包含 19 个工具，这里只保留比赛用到的 **jackhmmer**：以单条蛋白质序列为查询，
在序列数据库中迭代搜索同源序列并构建多序列比对（MSA）。

你们的任务是在**保证结果正确性**的前提下，让搜索尽可能快。

## 快速开始（探索1000 集群）

集群规则速记：**登录节点（cln01-04）禁止运行任何程序**；编译和小规模调试在
测试节点（`ssh test01`，test01-04 任选空闲的，`top` 查看占用）；完整计算必须
`sbatch` 提交到计算节点（56 核 / 192 GB / CentOS 7.9）；**home 目录在计算节点
上只读**，代码、数据与输出一律放 `~/WORK` 下。

目录布局（仓库和数据并列放在 `~/WORK` 下，脚本默认按这个布局找数据）：

```
~/WORK/
├── hmmer/          本仓库；运行结果写在 hmmer/results/
└── preliminary/    初赛数据：query/（100 条查询）+ database.fasta
```

```bash
# 1. 登录节点：克隆代码到 WORK
#    仓库必须放 WORK：搜索结果默认写在仓库下的 results/，home 在计算节点只读
#    如果集群网络状态差，可以考虑直接下载压缩包并scp上传
mkdir -p ~/WORK && cd ~/WORK
git clone https://github.com/shuixirui/hmmer.git && cd hmmer

# 2. 测试节点：编译 + 冒烟测试
ssh test01
module load compilers/gcc/v12.2.0    # 系统默认 gcc 4.8.5 过老，建议加载新版
cd ~/WORK/hmmer
autoconf            # 由 configure.ac 生成 configure 脚本
./configure
make
./smoke_test.sh     # 预期输出: PASS: 8 hits

# 3. 数据解压到 WORK，与仓库并列
cd ~/WORK
tar xzf preliminary.tar.gz           # 得到 ~/WORK/preliminary/{query/, database.fasta}

# 4. 提交完整基线作业（在 WORK 下提交，输出 stdout.%j 写在提交目录）
sbatch hmmer/job_baseline.slurm      # 路径可在脚本顶部 REPO/DATA/OUT 三行修改
squeue -u $USER                      # 查看状态：PD 排队 / R 运行
```

## 项目结构

```
├── configure.ac        # 构建系统源文件：autoconf 由它生成 configure
├── config.guess/sub    # configure 依赖的平台识别脚本
├── install-sh          # make install 依赖
├── smoke_test.sh       # 冒烟测试（见上）
├── run_queries.sh      # 批量搜索驱动：遍历 query/ 对 database.fasta 逐条搜索
├── job_baseline.slurm  # 探索1000 的 sbatch 作业模板（调用 run_queries.sh）
├── score.py            # 自评打分（从数据 header 读 SCOP 标签）
├── data/
│   ├── query.fasta     # 演示查询：截短血红蛋白 d1dlwa_（116 aa）
│   └── demodb.fasta    # 演示数据库：SCOPe 结构域序列集（15177 条）
├── src/                # HMMER 核心（63 个 .c）
│   ├── jackhmmer.c     # ← 程序入口：读序列、驱动迭代、汇总输出
│   ├── p7_pipeline.c   # ← 单条目标序列的加速流水线（过滤器级联）
│   ├── p7_builder.c    # 由比对构建 profile HMM（迭代的"建模"步骤）
│   ├── p7_tophits.c    # 命中收集、排序与输出
│   ├── impl_sse/       # ← SSE 向量化核心：MSV/Viterbi/Forward 滤波器
│   │                   #    (msvfilter.c, vitfilter.c, fwdback.c ...)
│   ├── impl_vmx/ impl_neon/   # PowerPC/ARM 向量实现（x86 上不参与编译）
│   └── generic_*.c     # 各算法的标量参考实现
└── easel/              # 底层库：序列/比对 I/O、字母表、数学与统计
    └── miniapps/       # easel 自带小工具（构建会编译，比赛不使用）
```

**读代码的建议入口**：`jackhmmer.c` 的主循环 → `p7_Pipeline()`（`p7_pipeline.c`）。
每条目标序列都要过一遍 MSV → Viterbi → Forward 的级联，逐级淘汰。

## 注意事项

- 优化不得改变搜索结果的正确性（评分以平均 ROC10 检验）。
- **不得依据序列 ID 或 header 内容（包括 SCOP 分类码、`annotated_` 前缀）筛选、排序或剔除
  数据库序列，也不得利用这些信息构造命中表。** header 中的标注只供评分使用，
  组委会将审查提交代码，违者取消成绩。
- 原版 HMMER 的完整文档见 [hmmer.org](http://hmmer.org)；本仓库删除了
  文档、测试套件与其余 18 个工具，如需参考请查阅上游仓库
  [EddyRivasLab/hmmer](https://github.com/EddyRivasLab/hmmer)。
