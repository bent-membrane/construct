#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ⛔ 作废（2026-10-29）
# 原因：① 能量景观 E(θ1,θ2)=k(θ1²+θ2²)+λ[1−cos3(θ1−θ2)] 的 λ=0.35 是手调参数；
#       ② 「速度 vs 质量」图在机制层用质量。
# 替代：braid_group.py（同一问题按群论零参数重算，10-27）、braid_four_projection.py
"""
两股 vs 三股：辫的拓扑决定解绑方式
====================================

原话（2026-10-26）：
    **"作为拓扑结构，绳索的三根缠绕就应该和两根缠绕有区别，
      就和编麻花辫似的。互相之间作用不一样。"**

这句话不是修辞，是拓扑学。

两股绳：辫群 B2。
    B2 是无限循环群——**只有一个生成元**（互绕一次），
    整个状态用一个整数描述：缠绕数。
    解开只有一条路：缠绕数从 n 单调回到 0。
    **不可能卡在半路。**

三股绳：辫群 B3。
    B3 不是循环群。它由两个生成元 σ1、σ2 张成，
    而且有辫关系 σ1·σ2·σ1 = σ2·σ1·σ2（这是三股才有的）。
    **同一个辫有多种写法，解开有多条路径。**
    更重要：三股可以在"两股已解开、第三股还缠着"的地方**停住**，
    因为那里是一个**局部极小**，不是死路。

结论方向：
    两股 = 全或无（单调解绑）
    三股 = 有中间站（可以卡在半解开）

可观测对应：
    三体并合系统应该比双体系统**更多**地呈现"部分连接"的构型。
    ——子弹星团主团的三个亚团块 + Cha+2025 那条质量桥，正好是这个样子。

用法：
    python braid_topology.py
    python braid_topology.py --plot ../figures/24_braid.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def energy_two(theta, k=1.0):
    """两股：单一扭转角，能量是单调的碗。"""
    return k * theta ** 2


def energy_three(t1, t2, k=1.0, lam=0.35):
    """
    三股：两个相对角 t1、t2（分别对应 σ1、σ2 方向）。
    三股辫有 3 重对称——所以能量里出现 cos(3Δ)。
    lam 控制"编织耦合"的强度：lam=0 时退化成两个独立的碗。
    """
    import numpy as np
    return k * (t1 ** 2 + t2 ** 2) + lam * (1.0 - np.cos(3.0 * (t1 - t2)))


def gradient_points_two(n=9, span=2.5):
    pts = []
    for i in range(n):
        pts.append((-span + 2 * span * i / (n - 1), 0.0))
    return pts


def gradient_points_three(n=9, span=2.5):
    pts = []
    for i in range(n):
        for j in range(n):
            pts.append((-span + 2 * span * i / (n - 1),
                        -span + 2 * span * j / (n - 1)))
    return pts


def relax_two(t0, k=1.0, steps=4000, dt=0.004):
    """梯度流：两股从初值滑到谷底。"""
    t = t0
    path = [t]
    for _ in range(steps):
        g = 2 * k * t
        t -= dt * g
        path.append(t)
    return path


def relax_three(t1, t2, k=1.0, lam=0.35, steps=6000, dt=0.003):
    """梯度流：三股从初值滑到最近的局部极小。"""
    a, b = t1, t2
    path = [(a, b)]
    for _ in range(steps):
        ga = 2 * k * a - 3 * lam * math.sin(3 * (a - b))
        gb = 2 * k * b + 3 * lam * math.sin(3 * (a - b))
        a -= dt * ga
        b -= dt * gb
        path.append((a, b))
    return path


def build_report():
    out = []
    add = out.append
    add("=" * 76)
    add("两股 vs 三股：辫的拓扑决定解绑方式")
    add("=" * 76)
    add("")

    add("-" * 76)
    add("一、拓扑结构本体")
    add("-" * 76)
    add("  两股（B2）：")
    add("    · 一个生成元 σ，缠绕数 n ∈ ℤ")
    add("    · σⁿ 是全部状态 —— **一维的**")
    add("    · 解开 = n → 0，路径唯一")
    add("")
    add("  三股（B3）：")
    add("    · 两个生成元 σ1、σ2")
    add("    · 辫关系 σ1·σ2·σ1 = σ2·σ1·σ2（三股独有，两股没有对应物）")
    add("    · 状态是二维的，**同一个辫有多种写法**")
    add("    · 解开 = 走到 unit 辫，路径不唯一，**中间可以停**")
    add("")

    add("-" * 76)
    add("二、能量景观：能不能卡住")
    add("-" * 76)
    add("  两股：E(θ) = kθ^2  —— 一个碗，处处往下，没有岔路")
    add("  三股：E(θ1,θ2) = k(θ1²+θ2²) + λ[1 − cos3(θ1−θ2)]")
    add("        后一项是**编织耦合**，带 3 重对称（三股的特征）")
    add("")

    # 验证三股的局部极小
    add("  三股在 λ=0.35 下的局部极小（梯度流终点）：")
    starts = [(1.5, 0.0), (0.0, 1.5), (1.2, -1.2), (2.0, 0.6)]
    seen = []
    for s in starts:
        p = relax_three(*s)
        a, b = p[-1]
        e = energy_three(a, b)
        add("    起点 (%.1f, %.1f) → 终点 (%+.2f, %+.2f)  E=%.3f"
            % (s[0], s[1], a, b, e))
        seen.append((round(a, 2), round(b, 2)))
    uniq = set(seen)
    add("")
    add("  ★ 不同起点收敛到不同终点：共 %d 个不同的极小。" % len(uniq))
    add("    这些极小对应「辫的不同写法」——三股可以停在其中任何一个。")
    add("    两股只有一个极小（θ=0），没有这种分岔。")

    add("")
    add("-" * 76)
    add("三、这对应什么样的观测")
    add("-" * 76)
    add("  两体并合：")
    add("    · 两个团要么连着、要么分开，中间没有稳定状态")
    add("    · 分离距离随时间单调变化（先增后减，但不会停在半路）")
    add("")
    add("  三体并合：")
    add("    · 可以出现「两团已分开、第三团还连着」的稳定配置")
    add("    · 也就是**部分连接**——一个桥连着两块，第三块已经跑远")
    add("    · 而且这种配置是**局部极小**，会停留很久")
    add("")
    add("  ★ 子弹星团现在的样子，正好是这个：")
    add("    · 主团拆成至少 3 个亚团块（Cha+2025）")
    add("    · 一条质量 + ICL 的桥连着主团和子团（5σ / 3σ）")
    add("    · 射电在两团之间也有过量，位置重合")
    add("    ——**两团已经分开，但中间那根还连着。**")
    add("")
    add("  按两股的模型，这种「分开又连着」的构型不该稳定存在；")
    add("  按三股的模型，它恰恰是一个局部极小，停多久都合理。")

    add("")
    add("-" * 76)
    add("四、还有一个更硬的差别：手性")
    add("-" * 76)
    add("  两股：左右手性对称。两根绳互绕，整体镜像一下，缠绕数变号——")
    add("        但系统本身没有偏好，**净角动量必须为零**（对称解）。")
    add("")
    add("  三股：**散开**。")
    add("")
    add("  ★ 这里有一个必须先纠正的错误推法（2026-10-26 修正）：")
    add("    「3 是奇数、三股抵消不干净、必然剩净环流」——**这句是错的**。")
    add("    三根对称排列（互成 120°）时，三个角动量矢量也成 120°，")
    add("    三个等长矢量互成 120° **恰好和为零**。奇偶跟能不能抵消没关系。")
    add("    原话：**「不是应该 2 根的才转吗，三根的就纯散开了」**")
    add("")
    add("  真正该看的是「散不散」，不是「转不转」：")
    add("    两股：解开之后两根**还绑在一起**（互绕/反向甩），是稳定的二体")
    add("    三股：解开之后三根**各自独立**，是混沌的三体")
    add("")
    add("  ★ 修正后的预言：")
    add("    **三体系统会出现「一个最轻的子团以异常高的速度脱离」；")
    add("      双体系统不会——两块会一直绑在一起回程。**")
    add("")
    add("  数据侧已经有一个候选：")
    add("    MACS J1149（Golovich+2016）：北 16.7 / 南 10.8 / 第三 1.2 (×1e14 M☉)")
    add("    第三个只有 7% 的质量，速度却是 **~4000 km/s**，")
    add("    比主并合的 2770 km/s 还快，方向几乎沿视线。")
    add("    —— 这就是「最轻的那块先被踢出去」的样子。")

    add("")
    add("-" * 76)
    add("五、和之前的局部最优是同一件事")
    add("-" * 76)
    add("  第五章里我们用「局部最优」解释小行星为什么是土豆：")
    add("    质量太小 → 卡在浅坑里出不来。")
    add("")
    add("  这里是同一个机制的**拓扑版**：")
    add("    三股辫的解绑能量景观天然有多个极小（3 重对称逼出来的），")
    add("    系统滑进其中一个就停住，不是因为质量小，是因为**拓扑允许它停**。")
    add("")
    add("  两股没有这个自由度——它的能量景观只有一个碗，")
    add("  所以两体系统只能「全解开」或「全没解开」。")
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

    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.34, wspace=0.3)

    # ── 1) 两股辫示意
    ax = fig.add_subplot(gs[0, 0])
    t = np.linspace(0, 4 * np.pi, 400)
    ax.plot(np.sin(t), np.linspace(0, 1, 400), lw=2.4, color="#c53030")
    ax.plot(-np.sin(t), np.linspace(0, 1, 400), lw=2.4, color="#2b6cb0")
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("两股（B2）：一个生成元\n状态 = 一个整数", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(alpha=0.15)

    # ── 2) 三股辫示意
    ax = fig.add_subplot(gs[0, 1])
    t = np.linspace(0, 4 * np.pi, 400)
    for k, c in enumerate(("#c53030", "#2b6cb0", "#38a169")):
        ax.plot(np.sin(t + 2 * np.pi * k / 3), np.linspace(0, 1, 400),
                lw=2.4, color=c)
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("三股（B3）：两个生成元\nσ1·σ2·σ1 = σ2·σ1·σ2", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(alpha=0.15)

    # ── 3) 能量景观一维对比
    ax = fig.add_subplot(gs[0, 2])
    th = np.linspace(-2.6, 2.6, 500)
    ax.plot(th, energy_two(th), lw=2.6, color="#2b6cb0", label="两股：E∝θ^2")
    ax.plot(th, energy_three(th, 0.0), lw=2.2, ls="--", color="#c53030",
            label="三股：固定 θ2=0 的切面")
    ax.set_xlabel("扭转角")
    ax.set_ylabel("能量")
    ax.set_title("一维切面：两股只有一个碗\n三股在切面上就能看到多个坑", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    # ── 4) 三股的能量地形（二维）
    ax = fig.add_subplot(gs[1, 0:2])
    N = 260
    lim = 2.6
    t1 = np.linspace(-lim, lim, N)
    t2 = np.linspace(-lim, lim, N)
    T1, T2 = np.meshgrid(t1, t2)
    E = energy_three(T1, T2)
    cs = ax.contourf(T1, T2, E, levels=28, cmap="twilight")
    plt.colorbar(cs, ax=ax, fraction=0.03, label="能量")
    ax.contour(T1, T2, E, levels=14, colors="w", linewidths=0.35, alpha=0.5)
    # 画几条梯度流
    for s in [(1.5, 0.0), (0.0, 1.5), (1.2, -1.2), (2.0, 0.6),
              (-1.5, -0.3), (-0.8, 1.8)]:
        p = relax_three(*s)
        ax.plot([q[0] for q in p], [q[1] for q in p],
                lw=1.4, color="#ffd166", alpha=0.95)
        ax.plot(p[-1][0], p[-1][1], "*", ms=13, color="#ff4d4d")
    ax.set_xlabel("θ1（σ1 方向的扭转）")
    ax.set_ylabel("θ2（σ2 方向的扭转）")
    ax.set_title("三股的能量地形：多个局部极小（红星）\n滑进哪个就停哪个 —— 所以能「卡在半解开」", fontsize=11)
    ax.set_aspect("equal")

    # ── 5) 速度 vs 质量：三体里最轻的先飞出去
    ax = fig.add_subplot(gs[1, 2])
    # 双体：两块速度一致（Bullet，质量比 10:1，v_core≈2650）
    ax.scatter([0.91, 0.09], [2650, 2650], s=150, color="#2b6cb0",
               edgecolors="black", linewidths=0.8, zorder=5,
               label="双体（Bullet 10:1）")
    ax.plot([0.91, 0.09], [2650, 2650], color="#2b6cb0", ls="--", lw=1.2, alpha=0.6)
    # 三体：最轻的那块速度飙到 4000（MACS J1149）
    ax.scatter([0.58, 0.38, 0.04], [2770, 2770, 4000], s=150, color="#c53030",
               marker="s", edgecolors="black", linewidths=0.8, zorder=5,
               label="三体（MACS J1149）")
    ax.plot([0.58, 0.38], [2770, 2770], color="#c53030", ls="--", lw=1.2, alpha=0.6)
    ax.annotate("最轻的 7%\n速度 4000 km/s\n（被踢出去了）",
                xy=(0.04, 4000), xytext=(0.22, 3550), fontsize=8.5,
                color="#c53030",
                arrowprops=dict(arrowstyle="->", color="#c53030", lw=1.1))
    ax.set_xlabel("质量占比")
    ax.set_ylabel("相对速度 (km/s)")
    ax.set_title("修正后的判据：不是「转不转」，是「散不散」", fontsize=10)
    ax.set_ylim(2300, 4300)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.25)

    fig.suptitle("两股 vs 三股 —— 拓扑结构决定解绑方式", fontsize=15)
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="两股 vs 三股辫拓扑")
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
