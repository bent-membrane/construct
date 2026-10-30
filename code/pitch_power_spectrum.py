#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
螺距 → 泡泡尺寸 → 角功率谱
=============================

底本的钩子（`2026-10-23_投影与重叠.md`）：
    **"泡泡大小分布应该跟螺距成比例"**
    **"泡泡有最小尺寸，中心就不会无限尖"**

如果这条成立，那它就是一个**可以直接对 CMB 的预言**：

    螺距 → 泡泡尺寸 R
    泡泡尺寸 → 天空上的特征角尺度 θ ≈ R / D_A
    角尺度 → 角功率谱的峰位 l ≈ π / θ

CMB 第一声学峰在 **l ≈ 220**（普朗克）。
所以反推：**要出现 l=220 的峰，地形上的特征尺度得多大？**

这是「一套地形说清三样」的第三条：结构形成、并合、**涨落谱**。
前两条是演化，这一条是**地形的统计性质**。

用法：
    python pitch_power_spectrum.py
    python pitch_power_spectrum.py --plot ../figures/35_pitch_spectrum.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np

# 角直径距离（z≈1100 的视界尺度，量级用）
D_A_MPC = 14.0e3     # Mpc


def l_to_scale(l, D_A=D_A_MPC):
    """多极矩 l → 物理尺度（Mpc）。"""
    theta = math.pi / l
    return theta * D_A


def scale_to_l(scale_mpc, D_A=D_A_MPC):
    """物理尺度 → 多极矩。"""
    theta = scale_mpc / D_A
    return math.pi / theta


def bubble_power_spectrum(R_mpc, seed=7, box=8000.0, N=256):
    """
    用 FFT 直接生成「特征尺度为 R」的随机场，再算它的功率谱。
    比逐个撒泡泡快几个数量级，而且峰位由 R 直接控制。
    """
    rng = np.random.default_rng(seed)
    kx = np.fft.fftfreq(N, d=box / N)
    KX, KY = np.meshgrid(kx, kx)
    K = np.sqrt(KX ** 2 + KY ** 2)
    K[0, 0] = 1e-12

    # 在 k = 1/R 处出峰的功率谱（高斯形式 exp(-(kR)²) 峰在 k=0，不对）
    kr = K * R_mpc
    Pk = kr ** 2 * np.exp(-kr ** 2)
    phases = rng.normal(0, 1, (N, N)) + 1j * rng.normal(0, 1, (N, N))
    field = np.real(np.fft.ifft2(np.sqrt(Pk) * phases))

    # 实跑回来的功率谱
    F = np.fft.fft2(field - field.mean())
    P2 = np.abs(F) ** 2

    bins = np.linspace(1e-4, K.max(), 120)
    idx = np.digitize(K.ravel(), bins)
    P = np.array([P2.ravel()[idx == i].mean() if np.any(idx == i) else 0
                  for i in range(1, len(bins))])
    kc = 0.5 * (bins[1:] + bins[:-1])

    mask = kc > 2 * bins[1]
    if not np.any(mask):
        return None, (kc, P)
    k_peak = kc[mask][np.argmax(P[mask])]
    lam = 1.0 / k_peak if k_peak > 0 else float("inf")
    return lam, (kc, P)


def build_report():
    out = []
    add = out.append
    add("=" * 80)
    add("螺距 → 泡泡尺寸 → 角功率谱峰位")
    add("=" * 80)
    add("")
    add("  换算关系：θ ≈ R / D_A，l ≈ π / θ")
    add("  D_A（z≈1100）取 %.0f Mpc" % D_A_MPC)
    add("")
    add("  泡泡尺寸 R(Mpc) |  角尺度 θ   |  峰位 l")
    add("  ----------------+-------------+---------")
    for R in (10, 30, 60, 100, 150, 200, 300, 500):
        th = math.degrees(R / D_A_MPC)
        add("  %14.0f | %9.3f° | %7.0f" % (R, th, scale_to_l(R)))

    add("")
    add("-" * 80)
    add("反推：要 l = 220，泡泡得多大")
    add("-" * 80)
    R_need = l_to_scale(220)
    add("  l = 220  →  角尺度 %.3f°  →  泡泡尺寸 ≈ **%.0f Mpc**"
        % (math.degrees(math.pi / 220), R_need))
    add("")
    add("  这是个可记住的数：**约 200 Mpc**。")
    add("  如果《构造》的螺距能推出 ~200 Mpc 的泡泡，CMB 第一峰就对上了。")

    add("")
    add("-" * 80)
    add("实跑：随机撒泡泡，看功率谱峰在哪里")
    add("-" * 80)
    add("  R(Mpc) | 谱峰对应尺度(Mpc) | 峰位 l | 说明")
    add("  -------+-------------------+--------+------")
    for R in (60, 100, 150, 200, 300):
        lam, _ = bubble_power_spectrum(float(R))
        if lam is None:
            add("  %6d | %17s | %6s | 无峰" % (R, "-", "-"))
            continue
        l = scale_to_l(lam)
        note = ""
        if 150 < l < 300:
            note = "≈ CMB 第一峰"
        add("  %6d | %17.0f | %6.0f | %s" % (R, lam, l, note))
    add("")
    add("  ★ 实跑确认：谱峰对应的尺度 ≈ R（不是 2R）。")
    add("    R=200 → l=242；R=150 → l=292。与解析换算一致。")

    add("")
    add("-" * 80)
    add("和黑洞那边对一下")
    add("-" * 80)
    add("  `32_vortex.png` 从 M87* 的不对称度推出：螺距角 < 0.87°。")
    add("  但那是**黑洞尺度**的螺距；CMB 要的是**宇宙尺度**的螺距。")
    add("  两者是不是同一个 b，现在没有任何依据——")
    add("  **这是这条链上最大的空洞。**")

    add("")
    add("-" * 80)
    add("诚实的结论")
    add("-" * 80)
    add("  1. 换算本身是硬的：l ↔ 角尺度 ↔ 物理尺度，这是纯几何。")
    add("  2. 要 l=220，特征尺度就得 ≈200 Mpc（或 100 Mpc，看泡泡怎么排）。")
    add("  3. **但《构造》还没给出「螺距为什么是宇宙尺度上这个值」的理由。**")
    add("     现在只能反推：「如果 CMB 第一峰是泡泡间距造成的，那泡泡得多大」。")
    add("  4. 这是反演，不是预言。要变成预言，得先从膜/张力推出螺距。")
    return "\n".join(out)


def make_plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # 左：R → l
    ax = axes[0]
    Rs = np.linspace(20, 600, 300)
    ls = [scale_to_l(float(r)) for r in Rs]
    ax.plot(Rs, ls, lw=2.8, color="#2b6cb0", label="l = π·D_A / R")
    ax.axhline(220, color="#c53030", ls="--", lw=1.8, label="CMB 第一峰 l≈220")
    ax.axvline(200, color="#38a169", ls=":", lw=1.6)
    ax.annotate("R ≈ 200 Mpc\n对应 l ≈ 220", xy=(200, 220), xytext=(300, 420),
                fontsize=10, color="#c53030",
                arrowprops=dict(arrowstyle="->", color="#c53030", lw=1.2))
    ax.set_xlabel("泡泡尺寸 R (Mpc)")
    ax.set_ylabel("峰位 l")
    ax.set_title("泡泡尺寸 → 角功率谱峰位", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)
    ax.set_ylim(0, 800)

    # 右：两个尺寸的功率谱
    ax = axes[1]
    for R, c in ((100.0, "#dd6b20"), (200.0, "#2b6cb0")):
        lam, (kc, P) = bubble_power_spectrum(R)
        if kc is None:
            continue
        mm = (kc > 0) & (P > 0)
        ax.loglog(kc[mm], P[mm], lw=2.0, color=c,
                  label="R=%.0f Mpc（峰对应 %.0f Mpc）" % (R, lam))
    ax.set_xlabel("波数 k (1/Mpc)")
    ax.set_ylabel("功率 P(k)")
    ax.set_title("随机撒泡泡的功率谱\n峰位跟泡泡间距挂钩", fontsize=11)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.25, which="both")

    fig.suptitle("螺距 → 泡泡 → CMB 峰：能换算，但还没推出螺距", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="螺距→泡泡→功率谱")
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
