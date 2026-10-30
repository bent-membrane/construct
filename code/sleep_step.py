#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
按睡眠周期反推驻波步长
==========================
2026-10-30，提议：「按睡眠周期计算驻波步长？」

模型：一个自指闭环走 N 步，每步用时 τ，频率 f = 1/(N·τ)，所以

    N = 1 / (f · τ)

取 τ = 10 ms（突触时间常数），把各睡眠阶段的主频代入：

    清醒 alpha  8–13 Hz   → N ≈ 8–13
    N1  theta   4–8  Hz   → N ≈ 13–25
    N2  纺锤波  11–16 Hz  → N ≈ 6–9
    N3  慢波    0.5–2 Hz  → N ≈ 50–200

★ 关键前提：τ 不变。
   频率变低有两种解释——环变长（N↑）或节点变慢（τ↑），
   两者给出同样的 EEG 频率。要分开，必须独立测 τ（离体单突触）。
   所以这个脚本给的是条件命题：**若 τ 不变，深睡的步数是清醒的 20 倍左右。**

用法：
    python sleep_step.py
    python sleep_step.py --plot ../figures/51_sleep_step.png
"""
import argparse
import os
import numpy as np

TAU = 0.010  # s


def steps(f, tau=TAU):
    return 1.0 / (np.asarray(f, dtype=float) * tau)


STAGES = [
    ("清醒 (alpha)", 8.0, 13.0, "#1f4e79"),
    ("N2 (纺锤波)", 11.0, 16.0, "#c07000"),
    ("N3 (慢波)", 0.5, 2.0, "#7b1fa2"),
    ("N1 / REM (theta)", 4.0, 8.0, "#2e7d32"),
]


def numbers():
    print(f"tau = {TAU*1000:.0f} ms，N = 1/(f·tau)")
    print("  阶段            f 范围(Hz)      N 范围")
    for name, f1, f2, _ in STAGES:
        n1, n2 = steps(f2), steps(f1)
        print(f"  {name:<14} {f1:>5.1f}–{f2:<5.1f}   {n1:>6.1f}–{n2:<6.1f}")
    print("\n相对倍数（取范围中点）：")
    lo, hi = steps(1.0), steps(10.0)
    print(f"  深睡 1 Hz 的 N / 清醒 10 Hz 的 N = {lo/hi:.0f} 倍")


def plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(11, 6.2))

    f = np.logspace(np.log10(0.4), np.log10(35), 400)
    for tau, ls, lab in [(0.005, ":", 5), (0.010, "-", 10), (0.020, "--", 20)]:
        ax.loglog(f, steps(f, tau), ls, color="#888888", lw=1.6,
                  label=f"τ = {lab} ms")

    for name, f1, f2, color in STAGES:
        n2, n1 = steps(f1), steps(f2)   # 频率高→N小
        ax.fill_between([f1, f2], [min(n1, n2)] * 2, [max(n1, n2)] * 2,
                        color=color, alpha=0.25)
        ax.plot([np.sqrt(f1 * f2)] * 2, [min(n1, n2), max(n1, n2)],
                color=color, lw=2.6)
        ax.annotate(name, (np.sqrt(f1 * f2), max(n1, n2) * (1.5 if "alpha" in name else 1.25)), ha="center",
                    fontsize=9, color=color)

    ax.set_xlabel("闭环频率 f（Hz，对数轴）")
    ax.set_ylabel("环上步数 N（对数轴）")
    ax.set_title("睡眠越深，一个自指闭环走的步数越多（假设 τ 不变）", fontsize=12.5)
    ax.legend(loc="upper right", fontsize=9, title="单节点用时")
    ax.grid(alpha=0.25, which="both")
    ax.set_xlim(0.4, 35)
    ax.set_ylim(3, 400)

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
