#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
两条待办账的数值核对
======================
① κ/η 约束度 如何和偏振参数对应
② 驻波频率和 EEG 怎么匹配

① 偏振参数 → 选择度 η
   消光比 ER = T∥ / T⊥
   η = (ER − 1) / (ER + 1)
   正交时可见度 V(90°) = 1 − η = 2 / (ER + 1)

② 驻波频率 → EEG
   若"我"是闭环驻波，频率 f = v / L（传导速度 / 环路长度）
   拿不同大小的脑做对照：人 ~10 cm，鼠 ~1 cm
   若 f ∝ 1/L，人脑频率应比鼠脑低约 10 倍
   —— 实际 EEG 主频跨物种差不到 2 倍
   所以"几何共振"这个版本站不住。

用法：
    python exp_checks.py
    python exp_checks.py --plot ../figures/50_exp_checks.png
"""
import argparse
import os
import numpy as np


def eta_from_ER(ER):
    return (ER - 1.0) / (ER + 1.0)


def numbers():
    print("① 偏振片：消光比 → 选择度 η → 正交可见度")
    print("   ER        η        V(90°)")
    for ER in (10, 100, 300, 1000, 10000):
        print(f"  {ER:>6.0f}   {eta_from_ER(ER):.4f}   {1-eta_from_ER(ER):.4f}")

    print("\n   亚波长缝的 η（η = 1/(1+(a/λ)^2)）相当于多少消光比：")
    for ar in (0.1, 0.3, 0.5, 1.0):
        eta = 1.0 / (1.0 + ar * ar)
        er = (1 + eta) / (1 - eta)
        print(f"   a/λ={ar:>4}  η={eta:.3f}  ≈ ER {er:.0f}")

    print("\n② 驻波频率 f = v/L")
    brains = [("人", 0.10, 10.0), ("猫", 0.03, 10.0), ("鼠", 0.01, 7.0)]
    print("   脑      L(m)    EE 主频(Hz)   v 要多少才能凑出来")
    for name, L, f in brains:
        print(f"   {name}     {L:.2f}     {f:>6.1f}        {f*L:.3f} m/s")
    print("\n   若 v 相同（同类型神经纤维），f 应 ∝ 1/L：")
    v = 1.0
    for name, L, f in brains:
        print(f"   {name}  f_pred={v/L:>7.1f} Hz   实测 {f:>5.1f} Hz   差 {v/L/f:.1f} 倍")


def plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4))

    # ① ER → η
    er = np.logspace(0.5, 4, 400)
    ax1.semilogx(er, eta_from_ER(er), lw=2.4, color="#1f4e79", label="η（选择度）")
    ax1.semilogx(er, 1 - eta_from_ER(er), lw=1.8, ls="--", color="#c00000",
                 label="V(90°)（正交可见度）")
    for e in (100, 1000):
        ax1.plot([e], [eta_from_ER(e)], "o", color="#1f4e79", ms=6)
        ax1.annotate(f"ER={e} η={eta_from_ER(e):.3f}", (e, eta_from_ER(e)),
                     textcoords="offset points", xytext=(6, -14), fontsize=9, color="#404040")
    ax1.set_xlabel("消光比 ER = T∥/T⊥")
    ax1.set_ylabel("数值")
    ax1.set_title("① 偏振参数 → 选择度 η")
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.25)

    # ② 跨物种
    sizes = np.array([0.01, 0.03, 0.10])   # m
    freqs = np.array([7.0, 10.0, 10.0])    # Hz
    ax2.plot(sizes * 100, freqs, "o", ms=9, color="#c00000", label="实测 EEG 主频")
    for s, f, n in zip(sizes * 100, freqs, ["鼠", "猫", "人"]):
        ax2.annotate(n, (s, f), textcoords="offset points", xytext=(8, 6), fontsize=10)
    L = np.logspace(-2, -0.7, 200)
    for v, ls in [(0.1, ":"), (0.3, "--")]:
        ax2.plot(L * 100, v / L, ls, color="#888888", lw=1.6, label=f"f = v/L，v={v} m/s")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("环路尺度 L（cm，对数轴）")
    ax2.set_ylabel("频率（Hz，对数轴）")
    ax2.set_title("② 几何共振要求 f ∝ 1/L——实测不跟")
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.25, which="both")

    fig.suptitle("两条账的数值核对：① 能对上　② 对不上", fontsize=13, y=1.0)
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
