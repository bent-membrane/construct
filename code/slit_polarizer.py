#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
双缝的偏振选择度 —— 9.4 的耦合项，先给极限和过渡
====================================================
2026-10-30，要求：先补第九章第 1 条账（双缝耦合项），
「哪怕先给量纲和极限情形」。

问题：9.4 说「双缝本身就是偏振片」。这句话要能算，缺两个量：

    ① 单缝选择度 η —— 这条缝对偏振约束多硬，由 a/λ 决定
    ② 两缝夹角 θ —— 两条缝给的标签对不对得上

标准光学里，狭缝的偏振依赖是现成的：缝宽 a 小到跟波长 λ 可比时，
平行板波导的两个偏振基模截止条件不同（一个无截止、一个有截止 λc=2a），
所以窄缝本身就是一个偏振器；a 远大于 λ 时两模都通，偏振效应消失。

本脚本不推 Maxwell 方程，只给一个标称过渡：

    η(a/λ) = 1 / (1 + (a/λ)²)

    极限是对的：a 远小于 λ → η→1（完全约束）；a 远大于 λ → η→0（不约束）。
    中间的形状是标称的，不是严格解 —— 严格解要解缝的边界条件，
    这一格留给后面。

可见度：两路偏振夹角 θ 的光叠加，条纹可见度按 cos θ 掉。
把它跟 η 接起来：

    V = 1 − η·(1 − cos θ)

    η=0（宽缝、无标签）  → V=1，不管 θ，条纹都在
    η=1（窄缝、完全标签） → V=cos θ；θ=90° 时完全消失
    θ=0（普通双缝）      → V=1，条纹完整

这三个极限正好对上三件事：普通双缝（θ=0）、Shih-Alley 正交偏振片
（θ=90°）、以及「缝太宽就不会偏振」。

★ 必须标明：这是借标准光学的狭缝偏振当对应物，不是从膜张力推出来的。
   书里要写清楚，不能拿外部的尺子冒充自己的推导。

用法：
    python slit_polarizer.py --plot ../figures/45_slit_polarizer.png
"""
import argparse
import os
import numpy as np


def eta(x):
    """单缝偏振选择度。x = a/λ。0 → 1（完全偏振），∞ → 0（无偏振）。"""
    return 1.0 / (1.0 + np.asarray(x, dtype=float) ** 2)


def visibility(k, theta_deg):
    """两路标签夹角 theta 下的条纹可见度。k 是单缝选择度。"""
    th = np.deg2rad(np.asarray(theta_deg, dtype=float))
    return 1.0 - k * (1.0 - np.cos(th))


def plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.6))

    # ---- 左：η 随 a/λ ----
    x = np.logspace(-1.2, 1.2, 400)
    k = eta(x)
    ax1.semilogx(x, k, lw=2.4, color="#1f4e79")
    ax1.axvline(1.0, color="#c00000", ls="--", lw=1.2)
    ax1.text(1.08, 0.92, "a = λ", color="#c00000", fontsize=10)

    for xx, label in [(0.1, "a 远小于 λ\n完全约束"), (10.0, "a 远大于 λ\n不约束")]:
        kk = eta(xx)
        ax1.plot([xx], [kk], "o", color="#c00000", ms=6)
        ax1.annotate(label, (xx, kk), textcoords="offset points",
                     xytext=(8, -34 if xx < 1 else 6), fontsize=9,
                     color="#404040")

    ax1.set_xlabel("缝宽 / 波长  a/λ")
    ax1.set_ylabel("单缝选择度  η")
    ax1.set_ylim(-0.05, 1.08)
    ax1.set_title("① 缝越窄，对转向的约束越硬")
    ax1.grid(alpha=0.25)

    # ---- 右：可见度随两缝夹角 ----
    theta = np.linspace(0, 90, 300)
    for kk, style in [(1.0, "-"), (0.6, "--"), (0.25, ":")]:
        v = visibility(kk, theta)
        ax2.plot(theta, v, style, lw=2.2,
                 label=f"η = {kk:g}")
    ax2.axvline(90, color="#c00000", ls="--", lw=1.0)
    ax2.annotate("Shih-Alley\n两缝正交", (90, 0.55), textcoords="offset points",
                 xytext=(-70, 0), fontsize=9, color="#c00000")
    ax2.annotate("普通双缝\n两缝同向", (0, 1.0), textcoords="offset points",
                 xytext=(10, -32), fontsize=9, color="#404040")

    ax2.set_xlabel("两缝标签夹角  θ（度）")
    ax2.set_ylabel("条纹可见度  V")
    ax2.set_ylim(-0.05, 1.08)
    ax2.set_title("② 两缝对不上，条纹就没了")
    ax2.legend(loc="lower left", fontsize=9, framealpha=0.9)
    ax2.grid(alpha=0.25)

    fig.suptitle("双缝 = 偏振片：两个量、三个极限（标称过渡，非严格解）",
                 fontsize=13, y=1.02)
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

    print("极限核对：")
    for x in (0.1, 0.5, 1.0, 2.0, 10.0):
        print(f"  a/λ={x:>5}  η={eta(x):.3f}")
    print("\n可见度核对（θ=90°，即完全正交）：")
    for k in (0.0, 0.25, 0.6, 1.0):
        print(f"  η={k:>4}  V={visibility(k, 90):.3f}")


if __name__ == "__main__":
    main()
