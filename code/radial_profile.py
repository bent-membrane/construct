#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
径向剖面：中心薄且没断，外围厚且断了
======================================

2026-10-27 给出的图像（原始描述，未加工）：
    **"我们这个结构是中间薄周围厚，最中间绳子没断，外围有断裂膨出那种。"**

拆成两条：
    ① 厚度 h(r)：**中心小，外围大**
    ② 连续性：**r < r_c 没断（无临界点）；r > r_c 断裂（有临界点）**

模型（最简，只有一个形状参数）：
    h(r) = h₀ · [1 + (r/r₀)^p]
    断裂判据：h(r) > h_c
    → 临界半径   r_c = r₀ · (h_c/h₀ − 1)^(1/p)

关键：**厚度本身不是我能定的**（那是窗外的东西）。
能算的是：
    · 给定“中心薄外围厚”这个单调趋势，r_c 落在哪
    · r_c 对形状参数 p 有多敏感
    · 外围临界点的密度分布

用法：
    python radial_profile.py
    python radial_profile.py --plot ../figures/39_radial.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np


def thickness(r, h0=1.0, r0=1.0, p=2.0):
    """厚度剖面：中心小、外围大。"""
    return h0 * (1.0 + (r / r0) ** p)


def critical_radius(h_c, h0=1.0, r0=1.0, p=2.0):
    """h(r) = h_c 的位置。"""
    if h_c <= h0:
        return 0.0
    return r0 * (h_c / h0 - 1.0) ** (1.0 / p)


def pit_density(r, r_c, sigma=0.25):
    """
    临界点密度：只在 r > r_c 出现，且在边界附近最密
    （刚断的地方断口最多；越往外越散）。
    """
    if r < r_c:
        return 0.0
    return math.exp(-((r - r_c) / sigma) ** 2)


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("径向剖面：中心薄且没断，外围厚且断了")
    add("=" * 78)
    add("")
    add("  模型：h(r) = h₀·[1 + (r/r₀)^p]，断裂判据 h > h_c")
    add("  取 h₀=1、r₀=1")
    add("")
    add("-" * 78)
    add("一、临界半径 r_c 对两个参数的敏感度")
    add("-" * 78)
    add("   h_c/h₀ |  p=1  |  p=2  |  p=3  |  p=4")
    add("  --------+-------+-------+-------+------")
    for ratio in (1.5, 2.0, 3.0, 5.0, 10.0):
        row = "  %7.1f |" % ratio
        for p in (1, 2, 3, 4):
            rc = (ratio - 1.0) ** (1.0 / p)
            row += " %5.2f |" % rc
        add(row)
    add("")
    add("  ★ 读法：h_c/h0 = 2 时 r_c 恒等于 1，与 p 无关；")
    add("    偏离 2 之后，p 越大 r_c 越贴近 1。")
    add("    也就是说，**形状指数 p 决定「断裂半径对阈值有多敏感」**——")
    add("    这跟「坑的形状决定坡度」是同一个道理：形状比幅度重要。")

    add("")
    add("-" * 78)
    add("二、临界点的密度分布")
    add("-" * 78)
    add("  按「刚断的地方断口最多」：密度 ∝ exp(−((r−r_c)/σ)²)，σ=0.25")
    add("")
    add("   r/r_c | 密度 | 说明")
    add("  -------+------+------------------")
    for rr in (0.5, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0, 3.0):
        d = pit_density(rr, 1.0)
        note = ""
        if rr < 1.0:
            note = "连续区（绳子没断）"
        elif rr < 1.15:
            note = "断口最密"
        elif rr > 2.0:
            note = "散了"
        add("  %6.2f | %4.2f | %s" % (rr, d, note))
    add("")
    add("  ★ 所以「泡泡」不是均匀分布的：")
    add("    · 中心一圈完全没有（绳子没断）")
    add("    · r_c 附近一圈最密（断口集中）")
    add("    · 再往外越来越稀（已经散开）")

    add("")
    add("-" * 78)
    add("三、这给出什么可观测的对照")
    add("-" * 78)
    add("  1. **核心不该有结构**：如果 r < r_c 是连续区，那里面不该有")
    add("     独立的「泡泡」（星系团、晕）。")
    add("     → 对应观测：**星系团中心往往是一个大的主峰，不是一堆小峰。**")
    add("")
    add("  2. **环状分布**：如果泡泡只在 r > r_c，那它们应该呈**环带**分布，")
    add("     不是从中心一直铺到边。")
    add("     → 这是个能查的判据：泡泡的空间分布有没有一个「内边界」。")
    add("")
    add("  3. **断口最密的地方在 r_c 附近**，不是最外圈。")
    add("     → 观测上：中等半径的**碎片最多**，最外围反而稀。")

    add("")
    add("-" * 78)
    add("四、边界")
    add("-" * 78)
    add("  1. **h(r) 的形式是我设的**（单调增、幂律）。趋势是你给的，")
    add("     具体形状没有依据——这是窗外的东西。")
    add("  2. r_c 没有绝对数值，因为 h₀、r₀、h_c 都没有定标。")
    add("     能给的是**相对结构**（环带、内边界、密度峰位置）。")
    add("  3. 「断裂判据 h > h_c」也是设的。真实判据可能是张力、应变、")
    add("     或者别的几何量——但**只要是单调的，上面三条定性结论不变**。")
    return "\n".join(out)


def make_plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))

    # 1) 厚度剖面
    ax = axes[0]
    rs = np.linspace(0, 3, 400)
    for p, c in ((1, "#a0aec0"), (2, "#2b6cb0"), (4, "#2d3748")):
        ax.plot(rs, thickness(rs, p=p), lw=2.4, color=c, label="p=%d" % p)
    ax.axhline(2.0, color="#c53030", ls="--", lw=1.6, label="断裂阈值 hc")
    ax.axvspan(0, 1.0, color="#38a169", alpha=0.14)
    ax.text(0.35, 0.35, "连续区\n（绳子没断）", fontsize=9, color="#276749")
    ax.set_xlabel("r / r0")
    ax.set_ylabel("厚度 h")
    ax.set_title("厚度剖面：中心薄、外围厚\n越过 h_c 就断", fontsize=10)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.25)

    # 2) 临界点密度
    ax = axes[1]
    rs = np.linspace(0, 3, 600)
    d = np.array([pit_density(float(r), 1.0) for r in rs])
    ax.fill_between(rs, 0, d, color="#dd6b20", alpha=0.55)
    ax.plot(rs, d, lw=2.2, color="#c53030")
    ax.axvline(1.0, color="#38a169", ls="--", lw=1.6, label="r_c")
    ax.set_xlabel("r / r_c")
    ax.set_ylabel("临界点密度")
    ax.set_title("泡泡不是均匀分布\n断口集中在 r_c 附近", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    # 3) 平面示意：环带
    ax = axes[2]
    rng = np.random.default_rng(11)
    N = 2600
    r = rng.uniform(0, 3, N)
    keep = rng.uniform(0, 1, N) < np.array([pit_density(float(x), 1.0) for x in r])
    rr = r[keep]
    th = rng.uniform(0, 2 * np.pi, len(rr))
    ax.scatter(rr * np.cos(th), rr * np.sin(th), s=3, color="#dd6b20", alpha=0.75)
    circle = plt.Circle((0, 0), 1.0, fill=False, ls="--", lw=1.6, color="#38a169")
    ax.add_patch(circle)
    ax.text(0, 0, "连续区\n无泡泡", ha="center", va="center", fontsize=9,
            color="#276749")
    ax.set_aspect("equal")
    ax.set_title("空间分布：环带，不是满铺", fontsize=10)
    ax.set_xlabel("x / r_c")
    ax.set_ylabel("y / r_c")
    ax.grid(alpha=0.2)

    fig.suptitle("中心薄且不断、外围厚且断裂 —— 给出的相对结构", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="径向剖面")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    text = build_report()
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text + "\n")
        print("report written: %s" % args.out)
    else:
        print(text)

    if args.plot:
        make_plot(args.plot)


if __name__ == "__main__":
    main()
