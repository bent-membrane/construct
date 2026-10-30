#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
相位随机化 → 条纹可见度 —— 把 φ(t)、κ、θ 接成一个统计量
============================================================
2026-10-30，要求：「算一下」。
承接 9.4「双缝是 B2 装置」里那条红字：
「随机相位 φ(t) 是谁加的、加在哪个自由度上，还没有说清；
 它和本节前面那组 κ、θ 是什么关系，也没有接上。」

设定：
    两股绳，一股相位 θ1，另一股 θ2，相位差 Δθ = θ1 − θ2。
    屏幕某点两路的光程差给一个 δ。
    两路偏振夹角 φ（随机变量）。

    强度 I(δ) = 2 + 2·<cos(δ+φ)>
              = 2 + 2·cos δ·<cos φ> − 2·sin δ·<sin φ>

    φ 的分布对称时 <sin φ> = 0，于是

        V = |<cos φ>|

    —— 可见度就是 cos φ 的平均值。三种情形都归这一条：

    ① φ 固定为 θ（完全标记）      V = |cos θ|
    ② φ 均匀分布（完全退相干）    V = 0
    ③ φ 高斯，标准差 σ            V = exp(−σ²/2)
    ④ φ 双峰：κ 概率取 θ，         V = 1 − κ(1 − cos θ)
       1−κ 概率取 0（部分标记）

    ④ 就是 9.4 里那个公式。所以 κ 和 φ 的分布不是两个东西：
    κ = 「被标记的比例」，等价于 φ 分布里那根尖峰的高度。

★ 这不是新物理，是标准退相干的唯象写法。它的作用是让 9.4 那条红字
  有一个能算的落点，不是补上了「φ 是谁加的」。

用法：
    python phase_visibility.py
    python phase_visibility.py --plot ../figures/46_phase_visibility.png
"""
import argparse
import os
import numpy as np


def v_gauss(sigma):
    """高斯相位噪声下的可见度。"""
    return np.exp(-0.5 * np.asarray(sigma, dtype=float) ** 2)


def v_partial(kappa, theta_deg):
    """部分标记：κ 比例被转到 θ，其余不转。"""
    th = np.deg2rad(np.asarray(theta_deg, dtype=float))
    return 1.0 - kappa * (1.0 - np.cos(th))


def plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.4))

    # 左：高斯相位噪声
    sig = np.linspace(0, 3.2, 400)
    ax1.plot(sig, v_gauss(sig), lw=2.4, color="#1f4e79")
    for s in (0.5, 1.0, 2.0):
        v = v_gauss(s)
        ax1.plot([s], [v], "o", color="#c00000", ms=6)
        ax1.annotate(f"σ={s:g}  V={v:.2f}", (s, v), textcoords="offset points",
                     xytext=(8, 8), fontsize=9, color="#404040")
    ax1.set_xlabel("相位噪声标准差 σ（弧度）")
    ax1.set_ylabel("可见度 V")
    ax1.set_ylim(-0.03, 1.05)
    ax1.set_title("① 随机相位把条纹糊掉：V = exp(−σ²/2)")
    ax1.grid(alpha=0.25)

    # 右：部分标记
    k = np.linspace(0, 1, 300)
    for th, style in [(0, "-"), (45, "--"), (90, ":")]:
        ax2.plot(k, v_partial(k, th), style, lw=2.2, label=f"θ = {th}°")
    ax2.annotate("κ=1, θ=90°\n完全消失", (1.0, 0.0), textcoords="offset points",
                 xytext=(-95, 25), fontsize=9, color="#c00000")
    ax2.set_xlabel("被标记的比例 κ")
    ax2.set_ylabel("可见度 V")
    ax2.set_ylim(-0.03, 1.05)
    ax2.set_title("② 部分标记：V = 1 − κ(1 − cos θ)")
    ax2.legend(loc="lower left", fontsize=9, framealpha=0.9)
    ax2.grid(alpha=0.25)

    fig.suptitle("可见度 = cos φ 的平均值：φ(t)、κ、θ 是同一个统计量的三种写法",
                 fontsize=12.5, y=1.02)
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

    print("高斯相位噪声：")
    for s in (0.1, 0.5, 1.0, 2.0, 3.14):
        print(f"  sigma={s:>5}  V={v_gauss(s):.4f}")
    print("\n部分标记（theta=90 度，即完全正交）：")
    for k in (0.0, 0.25, 0.5, 1.0):
        print(f"  kappa={k:>4}  V={v_partial(k, 90):.4f}")
    print("\n部分标记（theta=45 度）：")
    for k in (0.0, 0.5, 1.0):
        print(f"  kappa={k:>4}  V={v_partial(k, 45):.4f}")


if __name__ == "__main__":
    main()
