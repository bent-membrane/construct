# 脚本索引

50 个数值实验，每个都能独立运行。图号对应 `../figures/` 里的文件名前缀；没有图号的脚本只打印数值。

## 膜与投影

| 脚本 | 内容 |
|---|---|
| `membrane_jet.py` | 图 01 · 膜凹陷 → 喷流 → 投影晕 |
| `mist_layer.py` | 图 02 · 雾化层 |
| `inverse_speed.py` | 图 03 · 从暗物质圈的轮廓反推喷射速度 |
| `anisotropy.py` | 图 07 · 雾化投影的各向异性 |
| `ellipsoid.py` | 图 10 · 三维椭球密度场 + 倾斜投影 |
| `resonance_lock.py` | 图 08 · 共振对齐 |
| `multimode.py` | 图 09 · 连续谱 + 离散模 |
| `dim_cycle.py` | 图 15 · 维度循环 |

## 子弹星团与并合系统

| 脚本 | 内容 |
|---|---|
| `bullet_clock.py` | 图 04 · 用 X 射线激波给碰撞配一块表 |
| `bullet_collision.py` | 图 05 · 两个"中心点 + 雾"相向而行 |
| `bullet_collision_asym.py` | 图 06 · 主团 : 子弹团 = 10 : 1 |
| `bullet_core_fog.py` | 图 16 · 核/雾位移比与分离方向 |
| `bullet_lens_model.py` | 图 22 · 真实双峰透镜模型（Markevitch+2004、Clowe+2006） |
| `four_clusters.py` | 图 18 · 四系统（MACS / Bullet / ZwCl / A520）落到同一条相线 |
| `phi_ceiling.py` | 图 19 · φ 有没有天花板 |
| `lensing.py` | 图 21 · 膜凹陷的引力透镜：光路、偏折角、剪切场 |

## 辫群与拓扑

| 脚本 | 内容 |
|---|---|
| `braid_group.py` | 图 40 · 按群论重算 B2 vs B3，零手调参数 |
| `braid_four_projection.py` | 图 41 · 四股粗粒化 B4 → B2 |
| `cluster_braid_classify.py` | 图 42 · 七系统按辫群粗粒化统一分类 |
| `fourd_topology.py` | 图 38 · 四维矢量场临界点分类（Hofmann+2018 框架） |
| `fourd_k_inverse.py` | 四维投影 → 反推折叠模数 k |
| ⛔ `braid_topology.py` | 图 24 · **已作废**：能量景观含手调参数 λ = 0.35 |
| ⛔ `braid_unravel.py` | 图 25 · **已作废**：三个模型全在机制层用质量/动量 |

## 结构与成像

| 脚本 | 内容 |
|---|---|
| `pit_shape.py` | 图 30 · 同深同宽的坑，坡度差一倍 |
| `vortex_projection.py` | 图 32 · 漩涡投影：黑洞为什么看起来像圆 |
| `noodle_imaging.py` | 图 33 · 面条化：成像被拉长，不是物质被拉长 |
| `radial_profile.py` | 图 39 · 径向剖面：中心薄且没断，外围厚且断裂 |
| `swirl_to_structure.py` | 图 37 · 搅动 → 地形 → 物质会聚（**未做对，已由 38 取代**） |
| `droplet_spiral_break.py` | 图 36 · 油滴被螺旋搅动 → 细丝 → 断成一串 |
| `pitch_power_spectrum.py` | 图 35 · 泡泡尺寸 → 角功率谱（后半段纯几何，前半段是推论） |
| `momentum_return.py` | 图 34 · 回程 = 动量：冲过谷底再荡回来 |

## 量子与波动

| 脚本 | 内容 |
|---|---|
| `slit_polarizer.py` | 图 45 · 双缝的偏振选择度：耦合项的极限与过渡 |
| `phase_visibility.py` | 图 46 · 相位随机化 → 条纹可见度 |
| `screen_modulation.py` | 图 47 · 屏幕上的概率调制 |
| `erase_by_filter.py` | 图 48 · 擦除：不改落点，只分堆 |
| `sigma_measurement.py` | 图 49 · σ 怎么测：Hornberger 2003 的 C70 碰撞退相干 |
| `exp_checks.py` | 图 50 · 两条待办账的数值核对 |

## EEG 与睡眠

| 脚本 | 内容 |
|---|---|
| `sleep_step.py` | 图 51 · 按睡眠周期反推驻波步长 |
| `eeg_harmonics.py` | 图 52 · 驻波谐波序列 vs EEG 频段（纯计算） |
| `eeg_bicoherence.py` | 图 53 · 睡眠 EEG：频谱 vs 谐波（需 Sleep-EDF 数据） |
| `eeg_standing_wave.py` | 图 54 · 64 通道 alpha：驻波还是行波（需 eegmmidb 数据） |

## 通讯几何

| 脚本 | 内容 |
|---|---|
| `comm_schematic.py` | 图 11 · 通讯是怎么走的 |
| `comm_process.py` | 图 12 · 卷曲 → 穿层（全或无）→ 闭环 |
| `comm_torus.py` | 图 13 · 三维卷曲闭环（环面版，(1,4) 螺旋） |
| `comm_ladder.py` | 图 14 · 线 → 面 → 体 → 绳索，绳索上跑闭环 |

## 数据与外部查询

| 脚本 | 内容 | 备注 |
|---|---|---|
| `crossmatch_first.py` | 18 个团的晕椭率/位置角 × FIRST 射电 PA | 需联网查询 Vizier |
| `nvss_match.py` | 18 团 × HEASARC NVSS Position_Angle | 需联网查询 HEASARC |
| `modulus_rate.py` | 模与率：三组数值实验 | 纯计算 |
| `make_story_figures.py` | 为对外文章配图 43、44 | 纯计算 |

## 运行

```bash
python braid_group.py                      # 打印推导与数值
python braid_group.py --plot ../figures/40_braid_group.png
```

依赖：`numpy`、`scipy`、`matplotlib`。EEG 脚本另需 `pyedflib` 与公开数据集，见根目录 README。
