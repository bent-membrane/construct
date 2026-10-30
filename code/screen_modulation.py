#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
屏幕上的概率调制 —— 条纹不是波峰相遇，是落点概率被配对规则改写
==================================================================
对应 9.4「干涉到底是什么」。

模型（标准双缝的分布，按《构造》重新解释）：

    P(x) = 包络(x) × [1 + V·cos(2π d x / (λL))] / 2

    包络(x) = sinc²(π a x / (λL))     单缝自己的落点分布
    V        = 配对规则的可见度（= |⟨cos φ⟩|，见 phase_visibility.py）

三件事分开看：
    ① 包络：一个缝单独开时，粉尘落在哪里（近轴像钟形）
    ② 调制：配对规律在这个包络上刻出的明暗
    ③ 单光子：每个光子按 P(x) 随机落一次；N 个之后直方图才长出条纹

★ 这不是新物理，就是夫琅禾费双缝的强度分布。
   它的作用是让「概率被配对规则改写」这句话有一个能画、能算的落点。

用法：
    python screen_modulation.py
    python screen_modulation.py --plot ../figures/47_screen_modulation.png
"""
import argparse
import os
import numpy as np


def envelope(x, lam, a, L):
    """单缝包络 sinc²(π a x / λL)，归一化到 1。"""
    u = np.pi * a * x / (lam * L)
    out = np.ones_like(x)
    nz = np.abs(u) > 1e-12
    out[nz] = (np.sin(u[nz]) / u[nz]) ** 2
    return out


def prob(x, lam, a, d, L, V):
    """落点概率密度：包络 × 配对调制。"""
    env = envelope(x, lam, a, L)
    mod = 0.5 * (1.0 + V * np.cos(2 * np.pi * d * x / (lam * L)))
    return env * mod


def numbers():
    lam, a, d, L = 0.5e-6, 0.1e-3, 0.5e-3, 1.0
    fringe = lam * L / d          # 条纹间距
    width = lam * L / a           # 包络半宽（第一个零点）
    print("参数：λ=0.5 μm, a=0.1 mm, d=0.5 mm, L=1 m")
    print(f"  条纹间距 = {fringe*1e3:.3f} mm")
    print(f"  包络第一零点 = {width*1e3:.3f} mm")
    print(f"  包络内可见条纹数 ≈ {2*width/fringe:.1f}")
    print("\n可见度与亮暗比（只算调制因子的峰谷）：")
    for V in (1.0, 0.7, 0.5, 0.25, 0.0):
        imax, imin = 0.5 * (1 + V), 0.5 * (1 - V)
        ratio = imax / imin if imin > 0 else float("inf")
        print(f"  V={V:>4}  亮/暗 = {ratio:.1f}" if imin > 0 else f"  V={V:>4}  无非零暗点")


def plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    lam, a, d, L = 0.5e-6, 0.1e-3, 0.5e-3, 1.0
    x = np.linspace(-6e-3, 6e-3, 4000)
    env = envelope(x, lam, a, L)

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9))

    # ① 包络 × 调制
    ax = axes[0][0]
    ax.plot(x * 1e3, env, "--", color="#888888", lw=1.3, label="包络（单缝）")
    ax.plot(x * 1e3, prob(x, lam, a, d, L, 1.0), color="#1f4e79", lw=2.0, label="P(x)（V=1）")
    ax.fill_between(x * 1e3, 0, prob(x, lam, a, d, L, 1.0), color="#1f4e79", alpha=0.15)
    ax.set_title("① 包络被配对规则切开")
    ax.set_xlabel("屏幕位置 x（mm）")
    ax.set_ylabel("相对概率")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    # ② 三个可见度
    ax = axes[0][1]
    for V, c, s in [(1.0, "#1f4e79", "-"), (0.5, "#c07000", "--"), (0.0, "#888888", ":")]:
        ax.plot(x * 1e3, prob(x, lam, a, d, L, V), s, color=c, lw=1.9, label=f"V={V:g}")
    ax.set_title("② 标记越强，调制越浅")
    ax.set_xlabel("屏幕位置 x（mm）")
    ax.set_ylabel("相对概率")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    # ③④ 单光子累积
    rng = np.random.default_rng(7)
    p = prob(x, lam, a, d, L, 1.0)
    p = p / p.sum()
    for ax, N in [(axes[1][0], 200), (axes[1][1], 20000)]:
        hits = rng.choice(x, size=N, p=p)
        ax.hist(hits * 1e3, bins=120, range=(-6, 6), color="#1f4e79", alpha=0.75)
        ax.set_title(f"{'③' if N == 200 else '④'} {N} 个光子落点")
        ax.set_xlabel("屏幕位置 x（mm）")
        ax.set_ylabel("个数")
        ax.grid(alpha=0.2)

    fig.suptitle("条纹 = 配对规则对落点概率的调制（不是波峰相遇）", fontsize=13, y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()
    if args.plot:
        os.makedirs(os.path.dirname(os.path.abspath(args.plot)), exist_ok=True)
        plot(args.plot)
        return
    numbers()


if __name__ == "__main__":
    main()
