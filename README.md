# HeavyOil

重油体系分子动力学模拟中芳香环聚集体分析工具。

基于 MDAnalysis + RDKit + 周期性 DBSCAN，自动识别芳香环、检测 π-π 堆积聚集体，并生成统计与可视化输出。

## 功能

- 读取 GROMACS 格式的拓扑文件 (.gro) 和轨迹文件 (.xtc)
- 使用 RDKit 自动识别分子中的芳香环并计算环心坐标
- 周期性边界条件下的 DBSCAN 聚类 — 检测芳香环 π-π 堆积形成的聚集体
- 按帧统计聚集体数量、大小分布（均值、标准差、四分位数等）
- 输出 VMD 可视化脚本和聚集体组成数据

## 依赖

| 包 | 用途 |
|---|---|
| MDAnalysis | 分子动力学轨迹读取与原子选择 |
| RDKit | 芳香环识别 |
| OpenBabel | 分子文件格式转换 |
| NumPy / Pandas | 数值计算与统计 |
| PyYAML | 配置文件解析 |
| tqdm | 进度条 |

Python ≥ 3.10，推荐通过 Conda 安装。

## 安装

```bash
# 克隆仓库
git clone https://github.com/lhypop/HeavyOil.git
cd HeavyOil

# 创建 conda 环境（推荐）
conda create -n heavyoil -c conda-forge --file package-list.txt
conda activate heavyoil
```

## 使用方法

基本用法 — 仅做聚类：

```bash
python -m HeavyOil -s system.gro -t traj.xtc
```

指定帧范围和采样间隔：

```bash
python -m HeavyOil -s system.gro -t traj.xtc -b 0 -e 1000 -i 10
```

启用统计信息输出（每帧聚类数量、大小分布）：

```bash
python -m HeavyOil -s system.gro -t traj.xtc --info
```

选择要分析的残基：

```bash
python -m HeavyOil -s system.gro -t traj.xtc --residueselect ASP1,ASN2
```

### 命令行参数

| 参数 | 说明 | 默认值 |
|---|---|---|
| `-s, --gro` | 拓扑文件 (.gro)，**必选** | — |
| `-t, --trj` | 轨迹文件 (.xtc)，可选 | — |
| `-b, --begin` | 起始帧 (0-based) | 0 |
| `-e, --end` | 结束帧 (不含) | 最后一帧 |
| `-i, --interval` | 采样间隔 | 1 |
| `--info` | 输出每帧聚类统计 | 关闭 |
| `--infoselect` | 统计范围：`part`(仅聚集体) / `all`(全体系) | part |
| `--residueselect` | 残基选择，支持通配符 `AS*` | all |

## 配置文件

编辑 `config/config.yaml` 调整 DBSCAN 参数：

```yaml
eps: 0.45          # 邻域半径 (nm)
min_samples: 2     # 核心点最小邻居数
```

## 输出文件

运行后在当前目录生成：

```
<gro文件名>_move.tcl           VMD 可视化脚本
data/cluster_nums_<stem>.txt   各帧聚类统计表
temp/                          中间文件（自动清理）
```

## 工作原理

1. **预处理** — 读取 GRO/XTC 文件，根据残基选择过滤分子列表
2. **芳香环识别** — 用 RDKit 检测每个分子中的芳香环，计算环几何中心
3. **周期性 DBSCAN** — 在考虑周期性边界条件的前提下，对环心坐标做密度聚类
4. **后处理** — 合并/去噪，统计聚集体大小分布
5. **可视化** — 生成 VMD Tcl 脚本，可将不同聚集体分别着色显示

## 许可

GNU General Public License v3.0

## 作者

OCEAN
