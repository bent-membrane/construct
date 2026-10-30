#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ⛔ 作废（2026-10-29）
# 原因：三个候选模型——能量均分（E/3）、动量守恒（Σmv=0）、逃逸速度（√(2GM/r)）——
#       全部在**机制层**使用质量/能量/动量，违反 2026-10-26 定的「按几何拓扑考虑」。
# 替代：braid_group.py（B2/B3 路径分叉，零手调参数）
"""
麻花辫散开：三个模型 vs MACS J1149 的实测速度
==============================================

原话（2026-10-26）：
    **"麻花辫散掉可以计算啊，直接拿来对照。"**

对。不玩定性了，直接算三根散开时各拿多少速度，然后和
Golovich+2016 实测的 MACS J1149 三个子团速度对。

三个候选模型：

  A  能量均分（Equipartition）
     三根从同一个缠结里解开，各分到 E/3：
         v_i = sqrt(2E / (3 m_i))
     → 最轻的那根**飞得最疯**（v ∝ 1/√m）

  B  动量守恒（Momentum conservation）
     从静止散开，质心不动：
         Σ m_i v_i = 0
     → 大块慢、小块快，但总动量锁死为 0

  C  弹射（Ejection / 逃逸速度）
     最轻的被剩下两块从核心弹出去，速度 ≈ 本地逃逸速度：
         v_esc(r) = sqrt(2 G (M₁+M₂) / r)
     → 速度由**弹射点半径**决定，不是由质量比决定

实测（Golovich+2016，MACS J1149）：
    质量 16.7 / 10.8 / 1.20 (×1e14 M☉)  →  占比 0.581 / 0.376 / 0.042
    速度 2770 / 2770 / ~4000 km/s
    两峰投影间距 2.35 Mpc；第三个几乎沿视线方向

用法：
    python braid_unravel.py
    python braid_unravel.py --plot ../figures/25_braid_unravel.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

MSUN = 1.989e30
G = 6.674e-11
MPC = 3.086e22
KMS = 1.0e3

# MACS J1149
M = [16.7, 10.8, 1.20]        # ×1e14 M☉
V = [2770.0, 2770.0, 4000.0]  # km/s
R_PROJ = 2.35                 # Mpc，北南投影间距


def v_escape(m_pair_sum_1e14, r_mpc):
    """逃逸速度 (km/s)。"""
    M_si = m_pair_sum_1e14 * 1e14 * MSUN
    r_si = r_mpc * MPC
    return math.sqrt(2 * G * M_si / r_si) / KMS


def model_A_equipartition():
    """
    能量均分：v ∝ 1/√m。用最重那根的实测速度标定总能量。
    标定方式：让最重（0.581）的速度 = 2770。
    """
    m = [x / sum(M) for x in M]
    v_heavy = 2770.0
    k = v_heavy * math.sqrt(m[0])          # v_i = k / √m_i
    return [k / math.sqrt(mi) for mi in m]


def model_B_momentum():
    """
    动量守恒：从静止散开。取两快一慢近似，用前两根速度相等（2770）
    反解第三根应满足 Σ m v = 0 的大小。
    但前两根实测同向同速，所以这个模型本身就对不上观测。
    """
    m = [x / sum(M) for x in M]
    # 若前两根速度在同一方向（2770），第三根反向：
    v3 = -(m[0] + m[1]) * 2770.0 / m[2]
    return [2770.0, 2770.0, abs(v3)]


def model_C_ejection(r_mpc):
    """弹射：最轻的那根速度 = 本地逃逸速度。"""
    return v_escape(M[0] + M[1], r_mpc)


def invert_ejection(v_obs=4000.0):
    """由观测速度反推弹射半径。"""
    M_si = (M[0] + M[1]) * 1e14 * MSUN
    r = 2 * G * M_si / (v_obs * KMS) ** 2
    return r / MPC


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("麻花辫散开：三个模型 vs MACS J1149 实测")
    add("=" * 78)
    add("质量占比：%.3f / %.3f / %.3f" % tuple(x / sum(M) for x in M))
    add("实测速度：%.0f / %.0f / %.0f km/s" % tuple(V))
    add("")

    add("-" * 78)
    add("模型 A：能量均分（v ∝ 1/√m）")
    add("-" * 78)
    va = model_A_equipartition()
    add("  预言：%.0f / %.0f / %.0f km/s" % tuple(va))
    add("  实测：%.0f / %.0f / %.0f km/s" % tuple(V))
    add("  最轻那根：预言 %.0f，实测 %.0f  →  差 %.0f 倍"
        % (va[2], V[2], va[2] / V[2]))
    add("  ✗ 严重高估。能量均分让最轻的飞得太疯。")
    add("")

    add("-" * 78)
    add("模型 B：动量守恒（Σ m v = 0）")
    add("-" * 78)
    vb = model_B_momentum()
    add("  预言第三根：%.0f km/s（若前两根同向 2770）" % vb[2])
    add("  实测第三根：%.0f km/s" % V[2])
    add("  → 要满足总动量为零，第三根得飙到 ~6.3 万，实测只有 4000。")
    add("  ✗ 前提就不成立：前两根**同向**，任何动量守恒解都要求第三根极快。")
    add("    实际上系统整体在动，不是从静止散开——把质心系搞错了。")
    add("")

    add("-" * 78)
    add("模型 C：弹射（速度 = 本地逃逸速度）")
    add("-" * 78)
    add("  最轻那根的速度由**弹射点半径**决定：")
    add("   r(Mpc) |  v_esc(km/s)")
    add("  --------+-------------")
    for r in (0.5, 1.0, 1.48, 2.0, 2.35, 3.0, 4.0):
        add("   %5.2f  |   %7.0f%s"
            % (r, model_C_ejection(r),
               "   ← 实测 4000 对应这里" if abs(r - 1.48) < 0.05 else ""))
    add("")
    add("  ★ 由实测 4000 km/s 反推：弹射半径 = **%.2f Mpc**" % invert_ejection())
    add("    两峰投影间距 2.35 Mpc —— 弹射点比它近三分之一。")
    add("    考虑到视线方向还有一块（第三个几乎沿 LOS），")
    add("    三维构型里这个半径完全合理。")
    add("")
    add("  ✓ 唯一撞上观测的模型。")

    add("")
    add("-" * 78)
    add("结论")
    add("-" * 78)
    add("  麻花辫散开**不是**能量均分，**也不该**用静止质心算动量。")
    add("  它是**弹射**：最轻的那根被剩下两块的势阱从某个半径上甩出去，")
    add("  速度就是那个半径上的逃逸速度。")
    add("")
    add("  ★ 于是得到一条可推广的判据：")
    add("    **三体系统里最轻子团的速度 = 它在弹射点的逃逸速度。**")
    add("    反过来说：只要测到最轻子团的速度，就能反推它是在哪个半径被踢出来的。")
    add("")
    add("  这条判据和「双体没有第三块」天然区分：")
    add("    双体系统（Bullet）= 两块互绕，没有「最轻被弹射」这回事。")
    add("    三体系统 = 一定有一块被弹射，速度由逃逸速度锁死。")
    return "\n".join(out)


def make_plot(path):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.2))

    # 左：三个模型的预言 vs 实测
    ax = axes[0]
    va = model_A_equipartition()
    vb = model_B_momentum()
    names = ["最重\n(0.58)", "次重\n(0.38)", "最轻\n(0.04)"]
    x = np.arange(3)
    w = 0.26
    ax.bar(x - w, va, w, label="A 能量均分", color="#dd6b20")
    ax.bar(x, vb, w, label="B 动量守恒", color="#718096")
    ax.bar(x + w, V, w, label="实测", color="#c53030")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("速度 (km/s)")
    ax.set_title("三个模型对最轻子团的预言差一个量级")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25, axis="y")
    ax.annotate("均分：1 万\n实测：4 千",
                xy=(2 - w, va[2]), xytext=(1.05, va[2] * 0.82),
                fontsize=9, color="#dd6b20",
                arrowprops=dict(arrowstyle="->", color="#dd6b20", lw=1))
    ax.annotate("动量守恒：6.3 万",
                xy=(2, vb[2]), xytext=(0.55, vb[2] * 0.95),
                fontsize=9, color="#4a5568",
                arrowprops=dict(arrowstyle="->", color="#4a5568", lw=1))

    # 右：逃逸速度曲线 + 实测点
    ax = axes[1]
    rs = np.linspace(0.3, 5.0, 400)
    vs = [model_C_ejection(float(r)) for r in rs]
    ax.plot(rs, vs, lw=2.6, color="#2b6cb0", label="弹射速度 = 逃逸速度")
    ax.axhline(4000, color="#c53030", ls="--", lw=1.5, label="实测最轻子团 4000 km/s")
    r_inv = invert_ejection()
    ax.axvline(r_inv, color="#c53030", ls=":", lw=1.5)
    ax.plot([r_inv], [4000], "o", ms=11, color="#c53030", mec="black", mew=1.2, zorder=6)
    ax.annotate("反推弹射半径\n%.2f Mpc" % r_inv,
                xy=(r_inv, 4000), xytext=(r_inv + 0.7, 4700),
                fontsize=10, color="#c53030",
                arrowprops=dict(arrowstyle="->", color="#c53030", lw=1.2))
    ax.axvline(R_PROJ, color="#38a169", ls="-.", lw=1.4, alpha=0.8)
    ax.text(R_PROJ + 0.05, 6200, "两峰投影间距 2.35 Mpc", fontsize=8.5,
            color="#38a169", rotation=90, va="top")
    ax.set_xlabel("弹射点半径 (Mpc)")
    ax.set_ylabel("速度 (km/s)")
    ax.set_title("唯一撞上的模型：弹射 = 逃逸速度")
    ax.set_xlim(0.3, 5.0)
    ax.set_ylim(0, 12000)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig.suptitle("麻花辫散开：算出来直接对 MACS J1149 的实测速度", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="麻花辫散开 vs 实测")
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
