#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
驻波的谐波序列 vs EEG 频段
=============================
2026-10-30，作者：「驻波的边界就是转一圈回到原点吧」

对——闭合条件要求绕一圈的相位变化是 2π 的整数倍，所以驻波频率是量子化的：

    f_n = n · f0

这个脚本拿两个候选基频（2 Hz 和 6 Hz）去套标准 EEG 频段，看哪个覆盖得上。

    delta   0.5–4 Hz
    theta   4–8 Hz
    alpha   8–13 Hz
    beta    13–30 Hz
    gamma   30+ Hz

结论：
    f0 = 2 Hz 时，2/4/6/8/10/12/14… 能把五个频段全部覆盖；
    f0 = 6 Hz 时，6/12/18/24/30… 缺 delta，且 alpha 里只有一个音（12）。

★ 账：真实 EEG 不是纯谐波，频段是人为划的，不是精确倍数。
   这个图只说明"量级上能覆盖"，不说明"就是谐波"。

用法：
    python eeg_harmonics.py
    python eeg_harmonics.py --plot ../figures/52_eeg_harmonics.png
"""
import argparse
import os
import numpy as np

BANDS = [
    ("delta", 0.5, 4, "#7b1fa2"),
    ("theta", 4, 8, "#2e7d32"),
    ("alpha", 8, 13, "#c07000"),
    ("beta", 13, 30, "#1f4e79"),
    ("gamma", 30, 60, "#c00000"),
]


def harmonics(f0, fmax=60):
    return [n * f0 for n in range(1, int(fmax / f0) + 1)]


def numbers():
    for f0 in (2.0, 6.0):
        hs = harmonics(f0)
        print(f"基频 f0 = {f0:g} Hz，谐波：{', '.join(f'{h:g}' for h in hs[:16])}")
        print("  覆盖：")
        for name, lo, hi, _ in BANDS:
            hit = [h for h in hs if lo <= h < hi]
            print(f"    {name:<6} {lo:>4}–{hi:<3} Hz  " +
                  (f"命中 {len(hit)} 个：{', '.join(f'{h:g}' for h in hit)}" if hit else "无"))
        print()


def plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, opts = plt.subplots(2, 1, figsize=(12, 6.4), sharex=True,
                             gridspec_kw={"height_ratios": [1, 2.2]})

    # 上：频段色带
    ax = opts[0]
    for name, lo, hi, color in BANDS:
        ax.axvspan(lo, hi, color=color, alpha=0.25)
        ax.text(np.sqrt(lo * hi), 0.5, name, ha="center", va="center",
                fontsize=10, color=color)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_title("标准 EEG 频段（人为划分，不是精确倍数）", fontsize=11)

    # 下：谐波序列（上下分区，避免重叠）
    ax = opts[1]
    for name, lo, hi, color in BANDS:
        ax.axvspan(lo, hi, color=color, alpha=0.08)

    rows = [(2.0, "#1f4e79", "f0 = 2 Hz", 0.5, 1.0),
            (6.0, "#c07000", "f0 = 6 Hz", 0.0, 0.5)]
    for f0, color, label, y0, y1 in rows:
        hs = [h for h in harmonics(f0) if h <= 40]
        for h in hs:
            ax.plot([h, h], [y0, y1], color=color, lw=1.5, alpha=0.7)
        ax.axhline(y1, color=color, lw=0.8, alpha=0.3)
        ax.text(hs[-1] + 1.0, (y0 + y1) / 2, label, color=color,
                fontsize=10, va="center")

    ax.set_xscale("log")
    ax.set_xlim(0.5, 48)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("频率（Hz，对数轴）")
    ax.set_title("驻波只能站在 n·f0 上：上排 f0=2，下排 f0=6", fontsize=11)
    ax.grid(axis="x", alpha=0.2, which="both")

    fig.suptitle("转一圈回到原点 → 频率量子化 → 谐波序列", fontsize=13, y=0.98)
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
