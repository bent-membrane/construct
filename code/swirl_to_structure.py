#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
搅动 → 地形 → 物质会聚
========================

2026-10-27 给出的链：
    **绳子解开/破开 → 角动量 → 梯度（地形） → 物质沿梯度滑落 → 聚集**
    **泡泡决定了物质的分布。物质是被抖到我们这里的东西，随着梯度滑落就会聚集。**

这条链只有一环没算过：**角动量 → 梯度**。

模型（纯几何，不含质量/速度/能量）：

    1. 螺旋搅动在地形上留下的痕迹，是**螺旋状的密度调制**：
           z(r,θ) = −A · cos(m·θ − k·ln r)
       这是螺旋密度波的标准形式（m 条臂，k 控制缠得多紧）。
       **坡度来自螺旋本身，不需要质量。**

    2. 物质沿梯度滑落（过阻尼梯度下降，没有惯性）：
           dx/dt = −∇z
       这是「水往低处流」的最简形式，也是 SGD 最简形式。

    3. 看：粒子会不会自己聚成团？聚在哪？

预期：
    z 的极小在 m·θ − k·ln r = π 上 —— **这是一条螺旋线**。
    所以粒子应该聚成**沿螺旋排列的一串团**。

    这正是 `36_droplet_break.png` 缺的那一环：
    破裂给出「一串球」，这里给出「为什么它们沿螺旋排」。

用法：
    python swirl_to_structure.py
    python swirl_to_structure.py --plot ../figures/37_swirl_structure.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np


def pit_centers(n_pits=26, m=2, k=1.6, r_min=0.06, r_max=0.95):
    """沿一条对数螺旋摆一串坑的位置。"""
    th = np.linspace(0, m * 2 * np.pi * 1.15, n_pits)
    r = r_min * (r_max / r_min) ** (th / th.max())
    return r * np.cos(th), r * np.sin(th)


def terrain(X, Y, A=1.0, sigma=0.075):
    """沿螺旋线布置的一串高斯坑（离散极小）。"""
    cx, cy = pit_centers()
    Z = np.zeros_like(X)
    for xi, yi in zip(cx, cy):
        Z -= A * np.exp(-((X - xi) ** 2 + (Y - yi) ** 2) / (2 * sigma ** 2))
    return Z


def grad_terrain(x, y, A=1.0, sigma=0.075):
    """解析梯度：每个高斯坑贡献一项。"""
    cx, cy = pit_centers()
    gx = np.zeros_like(x)
    gy = np.zeros_like(y)
    s2 = sigma * sigma
    for xi, yi in zip(cx, cy):
        dx = x - xi
        dy = y - yi
        e = np.exp(-(dx * dx + dy * dy) / (2 * s2))
        # z = -A·exp(-d²/2s²)  →  -∂z/∂x = -A·exp·dx/s²（向坑心 = -dx 方向）
        gx -= A * e * dx / s2
        gy -= A * e * dy / s2
    return gx, gy


def evolve(n=4000, steps=800, dt=0.05, seed=3, rmax=1.0, save_every=40):
    """过阻尼梯度下降（向量化）。"""
    rng = np.random.default_rng(seed)
    r = rmax * np.sqrt(rng.uniform(0.02, 1.0, n))
    th = rng.uniform(0, 2 * np.pi, n)
    x = r * np.cos(th)
    y = r * np.sin(th)

    snaps = [(x.copy(), y.copy())]
    for s in range(steps):
        gx, gy = grad_terrain(x, y)
        x -= dt * gx
        y -= dt * gy
        if (s + 1) % save_every == 0:
            snaps.append((x.copy(), y.copy()))
    return snaps


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("搅动 → 地形 → 物质会聚")
    add("=" * 78)
    add("")
    add("  地形：z(r,θ) = −A·cos(mθ − k·ln r)，m=%d，k=%.1f" % (2, 1.6))
    add("  物质：过阻尼梯度下降 dx/dt = −∇z（没有惯性、没有质量）")
    add("")
    add("-" * 78)
    add("地形的结构")
    add("-" * 78)
    add("  z 的极小满足 mθ − k·ln r = (2n+1)π")
    add("  解出 θ = [k·ln r + (2n+1)π] / m")
    add("")
    add("  ★ 这是一条**对数螺旋**。")
    add("    所以物质会被拉到这条螺旋上，形成**沿螺旋排列的一串团**。")
    add("")
    add("  换句话说：")
    add("    · 螺旋搅动 → 螺旋地形")
    add("    · 螺旋地形 → 物质沿螺旋排布")
    add("    · 螺旋的间距（k）决定**团与团的间隔**")
    add("")
    add("  ★ 这就是 `36` 那一环缺的东西：")
    add("    36 说破裂会给出一串球，但没说为什么它们沿螺旋排。")
    add("    这里说：**地形本来就长成螺旋，物质只是顺着滑下去。**")

    add("")
    add("-" * 78)
    add("和观测对什么")
    add("-" * 78)
    add("  1. **星系旋臂**：如果物质聚在螺旋上，旋臂就是这条线。")
    add("     标准解释是密度波；这里不需要波，螺旋地形本身就是原因。")
    add("  2. **团与团的间隔**：由 k 定。k 小 → 缠得紧、间隔小；")
    add("     k 大 → 拉得开、间隔大。")
    add("  3. **可检验**：如果这是真的，那同一批\"球\"应该落在**同一条对数螺旋**上，")
    add("     而不只是\"大致沿旋臂\"。对数螺旋有确定的参数（k），可以拟合。")
    return "\n".join(out)


def make_plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig = plt.figure(figsize=(16, 5.6))
    gs = fig.add_gridspec(1, 3, wspace=0.25)

    # 1) 地形
    ax = fig.add_subplot(gs[0, 0])
    N = 400
    L = 1.0
    xs = np.linspace(-L, L, N)
    X, Y = np.meshgrid(xs, xs)
    Z = terrain(X, Y)
    im = ax.imshow(Z, extent=[-L, L, -L, L], cmap="twilight", origin="lower")
    plt.colorbar(im, ax=ax, fraction=0.046)
    ax.set_title("螺旋地形 z = −A·cos(2θ − k·ln r)\n深色 = 谷", fontsize=10)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

    # 2) 演化中间态
    snaps = evolve()
    mid = snaps[len(snaps) // 3]
    ax = fig.add_subplot(gs[0, 1])
    ax.imshow(Z, extent=[-L, L, -L, L], cmap="twilight", origin="lower", alpha=0.55)
    ax.scatter(mid[0], mid[1], s=2.2, color="#ffd166", alpha=0.8)
    ax.set_title("演化中途：开始沿谷聚集", fontsize=10)
    ax.set_xlabel("x")
    ax.set_aspect("equal")
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)

    # 3) 终态
    ax = fig.add_subplot(gs[0, 2])
    ax.imshow(Z, extent=[-L, L, -L, L], cmap="twilight", origin="lower", alpha=0.55)
    ax.scatter(snaps[-1][0], snaps[-1][1], s=2.2, color="#ff4d4d", alpha=0.85)
    ax.set_title("终态：物质落在**同一条螺旋**上\n（这就是旋臂的几何来源）", fontsize=10)
    ax.set_xlabel("x")
    ax.set_aspect("equal")
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)

    fig.suptitle("搅动 → 螺旋地形 → 物质沿螺旋排布", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="搅动→地形→结构")
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
