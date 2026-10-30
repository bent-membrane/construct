#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
σ 在实验上到底怎么测 —— 用 Hornberger 2003 的 C70 碰撞退相干实验
=====================================================================
文献：Hornberger, Uttenthaler, Brezger, Hackermüller, Arndt, Zeilinger,
      "Collisional Decoherence Observed in Matter Wave Interferometry",
      Phys. Rev. Lett. 90, 160401 (2003); arXiv:quant-ph/0303093

他们测的不是 σ，是**可见度随气体压强的衰减**：

    V(p) = V0 · exp(−p/p0),      p0 = k_B T / (2 L σ_eff)

    V0 = 41%          高真空下的最大可见度（C70，v=117 m/s）
    L  = 0.22 m       干涉仪长度
    σ_eff            有效碰撞截面（由 C6 色散系数算，无可调参数）

结论：可见度确实按指数掉，而且用不同气体（H2 到 Xe，质量跨两个数量级）
测出来的 p0 和理论预测对得上——**没有可调参数**。

★ 本脚本做两件事：
   1. 画文献里那条曲线长什么样（p0 是从 Fig. 2 读图估的，
      正文只给了公式和 41%，没给数据点）
   2. 反推：从可见度 V 到相位噪声 σ 的换算 —— σ = sqrt(−2 ln(V/V0))

★ 关键：σ 是反推出来的，不是直接测的。实验只给你 V；σ 取决于模型
   给 V 配的一个参数。这一步必须有独立扰动源标定，否则就是循环论证。

用法：
    python sigma_measurement.py
    python sigma_measurement.py --plot ../figures/49_sigma_measurement.png
"""
import argparse
import os
import numpy as np

V0 = 0.41                 # 文献给的高真空可见度
P0_MBAR = 5e-7            # 读图估计（Fig.2：V 掉到 1/e 时 p≈0.55e-6）
P_MAX = 2.5e-6            # 文献扫描上限


def visibility(p, p0=P0_MBAR, v0=V0):
    return v0 * np.exp(-p / p0)


def sigma_from_visibility(v_ratio):
    """V/V0 = exp(-σ²/2)  =>  σ = sqrt(-2 ln(V/V0))"""
    r = np.clip(np.asarray(v_ratio, dtype=float), 1e-6, 1.0)
    return np.sqrt(-2.0 * np.log(r))


def numbers():
    print("Hornberger 2003：C70 + CH4，V0 = 41%")
    print(f"  估计 p0 ≈ {P0_MBAR:.1e} mbar（量级估计，非文献精确值）")
    print("\n  压强 (mbar)   可见度 V")
    for p in (0, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 2.5e-6):
        print(f"  {p:>9.1e}   {visibility(p):.3f}")

    print("\n从可见度反推相位噪声 σ（通用换算，不依赖这个实验）：")
    print("  V/V0    σ (弧度)")
    for r in (1.0, 0.9, 0.7, 0.5, 0.3, 0.1, 0.01):
        print(f"  {r:>5.2f}   {sigma_from_visibility(r):.3f}")


def plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4))

    # ① 文献里那条曲线
    p = np.linspace(0, P_MAX, 400)
    ax1.plot(p * 1e6, visibility(p), lw=2.4, color="#1f4e79")
    ax1.axhline(V0, ls=":", color="#888888", lw=1.0)
    ax1.annotate("V0 = 41%（高真空）", (0, V0), textcoords="offset points",
                 xytext=(8, 6), fontsize=9, color="#404040")
    ax1.axvline(P0_MBAR * 1e6, ls="--", color="#c00000", lw=1.1)
    ax1.annotate("p0（可见度掉到 1/e）", (P0_MBAR * 1e6, 0.30),
                 textcoords="offset points", xytext=(6, 0), fontsize=9, color="#c00000")
    ax1.set_xlabel("气体压强 p（10^-6 mbar）")
    ax1.set_ylabel("条纹可见度 V")
    ax1.set_title("① 实验测的是这条曲线：V = V0·exp(−p/p0)")
    ax1.set_ylim(0, 0.46)
    ax1.grid(alpha=0.25)

    # ② 从 V 反推 σ
    r = np.linspace(0.01, 1.0, 400)
    ax2.plot(r, sigma_from_visibility(r), lw=2.4, color="#2e7d32")
    for rr in (0.9, 0.5, 0.1):
        s = sigma_from_visibility(rr)
        ax2.plot([rr], [s], "o", color="#c00000", ms=6)
        ax2.annotate(f"V/V0={rr:g}  σ={s:.2f}", (rr, s), textcoords="offset points",
                     xytext=(8, 4), fontsize=9, color="#404040")
    ax2.set_xlabel("可见度比 V/V0")
    ax2.set_ylabel("相位噪声 σ（弧度）")
    ax2.set_title("② 反推：σ = sqrt(-2 ln(V/V0))")
    ax2.grid(alpha=0.25)

    fig.suptitle("σ 不直接测：测可见度，再换算", fontsize=13, y=1.0)
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
