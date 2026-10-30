#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
油滴被螺旋搅动：拉成细丝，断成一串
=====================================

作者 2026-10-27：
    **"你觉得油滴被螺旋搅动会怎么破开"**
    **"我倒是想到油膜，一个大球破开了就变成更小的球"**
    **"破裂"**（不是整体拉开）
    **"螺距那个我都不知道啥意思"** → 螺距弃用

物理图像（Rayleigh–Plateau）：
    1. 螺旋搅动 = 剪切 + 旋转 → 油滴被**拉成细丝**（沿螺旋线）
    2. 细丝不是均匀断的 —— 它表面会先出现**周期性的鼓包**
    3. 最不稳定的波长 **λ ≈ 9.02 r**（r = 细丝半径）
    4. 鼓包长大、颈缩、断开 → **一串等间距的小球**，沿原螺旋排列

所以「一次裂成几块」这个问法不对——**细丝断裂天然给出一长串**，
数量由「细丝拉多长 / 波长多长」定。

体积守恒（细丝半径 r、长度 L、原球半径 R）：
    L = k·R          （k = 拉伸比）
    π r² L = (4/3)π R³   →   r = R·√(4/(3k))
    液滴数        n ≈ L / (9.02 r)
    单个液滴半径  r_drop ≈ 1.89 r
    尺度比        r_drop / R ≈ 1.89·√(4/(3k))

关键：**尺度比 ∝ 1/√k**，不是 ∝ 1/k。
所以想把球变小 100 倍，单级要拉 4.8 万倍——**不现实**。
但**级联**很快：每级拉 10 倍、变小 4.6 倍，三级就 97 倍。

用法：
    python droplet_spiral_break.py
    python droplet_spiral_break.py --plot ../figures/36_droplet_break.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np

LAM_FACTOR = 9.02      # 最不稳定波长 / 细丝半径
DROP_FACTOR = 1.89     # 液滴半径 / 细丝半径


def breakup(k):
    """给定拉伸比 k，返回 (细丝半径/R, 液滴数, 液滴半径/R)。"""
    r_over_R = math.sqrt(4.0 / (3.0 * k))
    L = k                      # 以 R 为单位
    n = L / (LAM_FACTOR * r_over_R)
    r_drop = DROP_FACTOR * r_over_R
    return r_over_R, n, r_drop


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("油滴被螺旋搅动：拉成细丝，断成一串")
    add("=" * 78)
    add("")
    add("  Rayleigh–Plateau：最不稳定波长 λ ≈ %.2f r，液滴半径 ≈ %.2f r"
        % (LAM_FACTOR, DROP_FACTOR))
    add("")
    add("  拉伸比 k | 细丝半径/R | 液滴数 | 液滴半径/R | 尺度缩小")
    add("  ---------+------------+--------+------------+---------")
    for k in (2, 5, 10, 20, 50, 100, 1000):
        rr, n, rd = breakup(k)
        add("  %8d | %10.3f | %6.1f | %10.4f | %7.1f 倍"
            % (k, rr, n, rd, 1.0 / rd))

    add("")
    add("-" * 78)
    add("一、单级做不到 100 倍")
    add("-" * 78)
    add("  尺度比 ∝ 1/√k，所以：")
    target = 0.01
    k_need = (DROP_FACTOR ** 2) * 4.0 / (3.0 * target ** 2)
    add("    要单级缩小 100 倍（r_drop/R = 0.01），需要 k ≈ %.0f" % k_need)
    add("    ——把半径 1 的球拉成快 5 万倍长的细丝。不现实。")
    add("")

    add("-" * 78)
    add("二、级联很快就够")
    add("-" * 78)
    add("  每级拉伸比 k  →  每级缩小倍数  →  到 100 倍需要几级：")
    for k in (10, 50, 100, 1000):
        rr, n, rd = breakup(k)
        per = 1.0 / rd
        need = math.log(100.0) / math.log(per)
        add("    k=%4d  →  缩小 %.2f 倍  →  %.1f 级" % (k, per, need))
    add("")
    add("  ★ **每级拉 100 倍，三级就够；每级拉 10 倍，要十二级。**")
    add("    所以「原初 200 Mpc → 现在 2 Mpc」不需要宇宙整体膨胀 100 倍，")
    add("    只需要**同一张膜按某个强度破了若干级**——级数取决于每级拉多狠。")

    add("")
    add("-" * 78)
    add("三、这个图像和「膨胀」的区别")
    add("-" * 78)
    add("  膨胀（拉开）：所有东西等比例变远，**结构本身尺寸不变**。")
    add("               观测：特征尺度应随宇宙学红移按 (1+z) 缩放。")
    add("")
    add("  破裂（本次）：**球本身变小**，数量变多。")
    add("               观测：高红移的\"球\"应该更大、更少；")
    add("                     低红移的球更小、更多，且**沿原来的螺旋排布**。")
    add("")
    add("  ★ 这两条的观测特征不一样，可以分开。")
    add("    「球的数量密度随红移上升、尺寸下降」支持破裂；")
    add("    「只是被拉开、尺寸按 (1+z) 缩放」支持膨胀。")

    add("")
    add("-" * 78)
    add("四、边界")
    add("-" * 78)
    add("  1. Rayleigh–Plateau 是**圆柱液柱**的结果（牛顿流体、无粘修正）。")
    add("     球状液滴被剪切流拉断，细节会不同（tip streaming、end pinching），")
    add("     但\"拉长→断裂→一串\"这个定性图像成立。")
    add("  2. 9.02 和 1.89 这两个系数是**流体**的数，带粘性、密度在里面。")
    add("     按最高原则（不用质量/速度），它们应该被理解为**几何常数**：")
    add("     断裂波长是半径的 9 倍、液滴是半径的 1.9 倍——纯比例关系。")
    add("  3. 级联的\"每级拉 10 倍\"是我取的，不是推出来的。")
    return "\n".join(out)


def make_plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4))

    # 1) 螺旋拉伸示意
    ax = axes[0]
    th = np.linspace(0, 6 * np.pi, 500)
    R0 = 1.0
    # 螺旋半径收缩，模拟拉伸
    rs = R0 * np.exp(-th / (6 * np.pi) * 1.8)
    ax.plot(rs * np.cos(th), rs * np.sin(th), lw=2.4, color="#dd6b20", alpha=0.9)
    ax.plot(0, 0, "o", ms=18, color="#c53030", alpha=0.35)
    ax.set_aspect("equal")
    ax.set_title("螺旋搅动把球拉成细丝\n（螺旋半径逐渐收缩）", fontsize=10)
    ax.grid(alpha=0.2)

    # 2) 细丝断裂示意
    ax = axes[1]
    x = np.linspace(0, 9 * 2 * np.pi, 2000)
    # 半径为 1 的细丝，加周期性扰动（波长 9r）
    r_base = 1.0
    amp = 0.18
    r = r_base + amp * np.sin(x / 9.0 * 2 * np.pi)
    up = r
    dn = -r
    ax.fill_between(x, dn, up, color="#dd6b20", alpha=0.55)
    ax.plot(x, up, lw=1.2, color="#c53030")
    ax.plot(x, dn, lw=1.2, color="#c53030")
    # 标出颈缩处
    for i in range(1, 9):
        xx = i * 9.0
        if xx < x[-1]:
            ax.plot(xx, 0, "|", ms=12, color="black", alpha=0.6)
    ax.set_xlabel("沿细丝的位置（单位：细丝半径 r）")
    ax.set_ylabel("细丝轮廓")
    ax.set_title("Rayleigh–Plateau：波长 ≈ 9r\n颈缩处断开 → 一串小球", fontsize=10)
    ax.set_aspect("equal")
    ax.grid(alpha=0.2)

    # 3) 拉伸比 vs 尺度
    ax = axes[2]
    ks = np.logspace(0.2, 3.5, 200)
    scales = [1.0 / breakup(float(k))[2] for k in ks]
    ax.loglog(ks, scales, lw=2.8, color="#2b6cb0", label="单级缩小倍数")
    ax.axhline(100, color="#c53030", ls="--", lw=1.6, label="目标：100 倍")
    # 级联 10x/级
    cum = 1.0
    rr, n, rd = breakup(10)
    cums, levs = [], []
    for lev in range(1, 5):
        cum *= 1 / rd
        cums.append(cum); levs.append(lev)
    ax.plot([10 * lev for lev in levs], cums, "o-", ms=9, lw=1.8,
            color="#38a169", label="级联（每级拉 10 倍）")
    ax.set_xlabel("拉伸比 k 或 级联累计拉伸")
    ax.set_ylabel("球缩小倍数")
    ax.set_title("单级达不到 100 倍\n级联要看每级拉多狠", fontsize=10)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.25, which="both")

    fig.suptitle("油滴被螺旋搅动：拉成细丝 → 断成一串 → 级联变小", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="油滴螺旋破裂")
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
