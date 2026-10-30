#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
膜凹陷的引力透镜
=================

之前画的都是「坑长什么样」（3D 曲面）和「雾化层投影」（2D 铺开）。
这一份算的是**光在坑里怎么被弯**——从几何图像走到可观测信号。

模型：
    两个膜凹陷（对应两个核），高度场  h(x,y)
    透镜势           ψ = h（地形本身就是势）
    偏折角           α = ∇ψ
    会聚             κ = ∇²ψ / 2
    剪切             γ1 = (∂²ψ/∂x² − ∂²ψ/∂y²)/2 ,  γ2 = ∂²ψ/∂x∂y

    每个凹陷用等温球（平坦旋转曲线的那种）：
        h_i(r) = A_i · ln(r / r_c)        （截断在 r_c 内）
    对等温球，偏折角 α_r = 2A/r —— 不随距离衰减，这正是强透镜能
    在很大半径上还能看到信号的原因。

关键问题：
    单坑的透镜信号（爱因斯坦环、临界线）是标准的。
    **双坑**会给出什么？两个临界线会不会合并？
    核-雾分离（表里分离）在透镜信号上到底留下什么痕迹？

用法：
    python lensing.py
    python lensing.py --plot ../figures/21_lensing.png
    python lensing.py --sep 720 --ratio 10
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ---------------------------------------------------------------- 势函数
# 用 Plummer 球：光滑、无奇点、有解析偏折角，不会像截断等温球那样在
# 截断半径处出现跳变（上一版 3D 图里那些"柱子"就是那个接缝）。
#   ψ(R) = (amp/2)·ln(R² + rc²)
#   α(R) = amp·R/(R² + rc²)        （沿径向）
#   κ(R) = amp·rc²/(R² + rc²)²
def psi_plummer(x, y, xc, yc, amp, rc):
    r2 = (x - xc) ** 2 + (y - yc) ** 2
    return 0.5 * amp * math.log(r2 + rc * rc)


def defl_plummer(x, y, xc, yc, amp, rc):
    dx, dy = x - xc, y - yc
    r2 = dx * dx + dy * dy
    a = amp / (r2 + rc * rc)
    return a * dx, a * dy


def total_defl(x, y, lenses):
    ax = ay = 0.0
    for (xc, yc, amp, rc) in lenses:
        dx, dy = defl_plummer(x, y, xc, yc, amp, rc)
        ax += dx
        ay += dy
    return ax, ay


def psi_total(x, y, lenses):
    return sum(psi_plummer(x, y, xc, yc, amp, rc) for (xc, yc, amp, rc) in lenses)


def kappa_plummer(r, amp, rc):
    return amp * rc * rc / (r * r + rc * rc) ** 2


# ---------------------------------------------------------------- 诊断
def build_report(sep, ratio, amp_main, rc):
    m1, m2 = 1.0, 1.0 / ratio
    s = m1 + m2
    x1 = -sep * m2 / s
    x2 = +sep * m1 / s
    lenses = [(x1, 0.0, amp_main * m1, rc), (x2, 0.0, amp_main * m2, rc)]

    out = []
    add = out.append
    add("=" * 76)
    add("膜凹陷的引力透镜")
    add("=" * 76)
    add("两坑间距 = %.0f kpc   质量比 = %.1f : 1" % (sep, ratio))
    add("主坑中心 (%.0f, 0)   次坑中心 (%.0f, 0)" % (x1, x2))
    add("等温球截断半径 rc = %.0f kpc" % rc)
    add("")

    add("-" * 76)
    add("一、偏折角沿连线的剖面")
    add("-" * 76)
    add("   x(kpc) |    α_x     |  说明")
    add("  --------+------------+---------------------------")
    for x in (-600, -400, -300, -200, -100, 0, 100, 200, 300, 400, 600):
        ax, ay = total_defl(float(x), 1e-6, lenses)
        note = ""
        if abs(x - x1) < 60:
            note = "主坑附近"
        elif abs(x - x2) < 40:
            note = "次坑附近"
        add("  %7d | %10.4f | %s" % (x, ax, note))

    add("")
    add("-" * 76)
    add("二、两个临界线在哪儿")
    add("-" * 76)
    add("  临界线 = 放大率发散的地方（会聚到 1 的等值线）")
    add("  单坑时是一圈；双坑时会被撕开、甚至合并成一个哑铃")
    add("")
    add("  r(kpc) |  主坑 κ   |  次坑 κ   |  总 κ")
    add("  -------+-----------+-----------+--------")
    for r in (20, 50, 100, 200, 300, 500, 800):
        k1 = kappa_plummer(float(r), amp_main * m1, rc)
        k2 = kappa_plummer(float(r), amp_main * m2, rc)
        add("  %6d | %9.3f | %9.3f | %7.3f" % (r, k1, k2, k1 + k2))

    add("")
    add("-" * 76)
    add("三、和标准引力透镜的关系")
    add("-" * 76)
    add("  标准 GR（等温球）：α = 4πσ²/c² · (1/b)，形式上就是 α ∝ 1/b")
    add("  《构造》：α = ∇h，h 是膜地形。只要 h 的形状是 ln(r)，")
    add("           偏折角就自动是 1/r —— 和 GR 的等温球完全同形。")
    add("")
    add("  ★ 这说明什么：**单坑**的透镜信号《构造》和 GR 无法区分。")
    add("    要区分，只能看**双坑**——两个凹陷的相互作用。")
    add("")
    add("  标准 GR 的双坑：两个独立的势阱线性叠加，各自一圈临界线。")
    add("  《构造》的双坑：两根绳缆滑脱后留下的两个**地形凹陷**，")
    add("                  如果凹陷之间还有膜张力的耦合，")
    add("                  二者之间会出现一条**鞍点脊**——")
    add("                  这条脊会改变中间区域的光路。")
    return "\n".join(out)


# ---------------------------------------------------------------- 图
def make_plot(sep, ratio, amp_main, rc, path):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    from matplotlib.colors import LinearSegmentedColormap

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    m1, m2 = 1.0, 1.0 / ratio
    s = m1 + m2
    x1 = -sep * m2 / s
    x2 = +sep * m1 / s
    lenses = [(x1, 0.0, amp_main * m1, rc), (x2, 0.0, amp_main * m2, rc)]

    # 网格
    N = 420
    L = sep * 1.9
    xs = np.linspace(-L, L, N)
    ys = np.linspace(-L, L, N)
    X, Y = np.meshgrid(xs, ys)

    # 势（Plummer，光滑无接缝）
    Z = np.zeros_like(X)
    Ax = np.zeros_like(X)
    Ay = np.zeros_like(X)
    for (xc, yc, amp, rcc) in lenses:
        dx = X - xc
        dy = Y - yc
        r2 = dx * dx + dy * dy
        Z += 0.5 * amp * np.log(r2 + rcc * rcc)
        a = amp / (r2 + rcc * rcc)
        Ax += a * dx
        Ay += a * dy

    # 偏折角大小（直接可观测量，单位无量纲）
    A = np.hypot(Ax, Ay)

    # 剪切 γ（从势的二阶导）
    gxx = np.gradient(np.gradient(Z, xs, axis=0), xs, axis=0)
    gyy = np.gradient(np.gradient(Z, ys, axis=1), ys, axis=1)
    gxy = np.gradient(np.gradient(Z, xs, axis=0), ys, axis=1)
    G1 = (gxx - gyy) / 2.0
    G2 = gxy
    G = np.hypot(G1, G2)

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=0.28, wspace=0.22)

    # ── 1) 地形 3D
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    stride = 6
    ax.plot_surface(X[::stride, ::stride], Y[::stride, ::stride],
                    Z[::stride, ::stride], cmap="twilight",
                    linewidth=0, antialiased=True, alpha=0.95)
    ax.set_xlabel("x (kpc)")
    ax.set_ylabel("y (kpc)")
    ax.set_zlabel("膜高度 h")
    ax.set_title("两个膜凹陷（地形本身就是透镜势）", pad=12)
    ax.view_init(elev=34, azim=-58)

    # ── 2) 光路
    ax = fig.add_subplot(gs[0, 1])
    ax.imshow(Z, extent=[-L, L, -L, L], cmap="twilight", origin="lower", alpha=0.55)
    ax.contour(X, Y, Z, levels=14, colors="w", linewidths=0.4, alpha=0.5)

    # 从左侧平行入射 —— 只画穿过透镜核心区的那些光线，偏折才看得清
    for y0 in np.linspace(-420, 420, 21):
        px, py = -L, y0
        pts_x, pts_y = [px], [py]
        step = L / 160.0
        while px < L:
            ax_, ay_ = total_defl(px, py, lenses)
            px += step
            py += ay_ * step
            pts_x.append(px)
            pts_y.append(py)
        ax.plot(pts_x, pts_y, lw=1.5, color="#ffd166", alpha=0.9)

    for (xc, yc, _, _) in lenses:
        ax.plot(xc, yc, "x", ms=11, mew=2.2, color="#ff4d4d")
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_xlabel("x (kpc)")
    ax.set_ylabel("y (kpc)")
    ax.set_title("光线在两个坑里的偏折", pad=10)

    # ── 3) 偏折角大小
    ax = fig.add_subplot(gs[1, 0])
    im = ax.imshow(A, extent=[-L, L, -L, L], cmap="inferno", origin="lower",
                   norm=matplotlib.colors.LogNorm(vmin=1e-3, vmax=A.max()))
    plt.colorbar(im, ax=ax, fraction=0.046, label="偏折角 |α|")
    cs = ax.contour(X, Y, A, levels=6, colors="cyan", linewidths=1.1)
    ax.clabel(cs, inline=True, fontsize=7, fmt="%.3f")
    ax.plot([xc for (xc, _, _, _) in lenses], [0, 0], "w+", ms=12, mew=2)
    ax.set_xlabel("x (kpc)")
    ax.set_ylabel("y (kpc)")
    ax.set_title("偏折角 |α|（可直接对照透镜观测）", pad=8)

    # ── 4) 剪切
    ax = fig.add_subplot(gs[1, 1])
    im = ax.imshow(G, extent=[-L, L, -L, L], cmap="viridis", origin="lower",
                   norm=matplotlib.colors.LogNorm(vmin=max(G[G > 0].min(), 1e-4),
                                                  vmax=G.max()))
    plt.colorbar(im, ax=ax, fraction=0.046, label="剪切 |γ|")
    # 剪切的方向线（每 30 像素取一个）
    st = 26
    ax.quiver(X[::st, ::st], Y[::st, ::st], G1[::st, ::st], G2[::st, ::st],
              color="w", alpha=0.75, scale=1e3, width=0.003)
    ax.set_xlabel("x (kpc)")
    ax.set_ylabel("y (kpc)")
    ax.set_title("剪切场 |γ|（弱透镜观测对照量）", pad=8)

    fig.suptitle("膜凹陷的引力透镜 —— 从「坑的形状」到「光的弯折」",
                 fontsize=15, y=0.975)
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="膜凹陷的引力透镜")
    ap.add_argument("--sep", type=float, default=720.0, help="两坑间距 kpc")
    ap.add_argument("--ratio", type=float, default=10.0, help="质量比")
    ap.add_argument("--amp", type=float, default=40.0, help="主坑深度系数")
    ap.add_argument("--rc", type=float, default=200.0, help="Plummer 核心半径 kpc")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    text = build_report(args.sep, args.ratio, args.amp, args.rc)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text + "\n")
        print("report written: %s" % args.out)
    else:
        print(text)

    if args.plot:
        make_plot(args.sep, args.ratio, args.amp, args.rc, args.plot)


if __name__ == "__main__":
    main()
