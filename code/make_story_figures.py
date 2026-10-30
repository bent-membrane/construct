#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
为《当星系团撞在一起》配两张解释图
=====================================
2026-10-27。不是数据图，是把概念摆出来给人看。

图 43：两根 vs 三根 —— 解绑路线从哪一刻开始分叉
图 44：四根的粗粒化 —— 三根合成一根粗绳，远看就是两根

用法：
    python make_story_figures.py
    python make_story_figures.py --outdir ../figures
"""

import argparse
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


COLORS = ["#c53030", "#2b6cb0", "#38a169", "#b7791f"]


def fig_two_vs_three(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig = plt.figure(figsize=(15.5, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], wspace=0.22)

    # ── A：两根
    ax = fig.add_subplot(gs[0, 0])
    t = np.linspace(0, 6 * np.pi, 700)
    y = np.linspace(0, 1, 700)
    ax.plot(np.sin(t), y, lw=3.0, color=COLORS[0])
    ax.plot(-np.sin(t), y, lw=3.0, color=COLORS[1])
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-0.15, 1.18)
    ax.set_title("两根：只有一条路", fontsize=12)
    ax.annotate("", xy=(0, 1.16), xytext=(0, 1.04),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#4a5568"))
    ax.text(0, 1.14, "解开", ha="center", va="bottom", fontsize=9, color="#4a5568")
    ax.text(0, -0.12, "全解 或 全不解", ha="center", fontsize=10.5, color="#2d3748")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    # ── B：三根
    ax = fig.add_subplot(gs[0, 1])
    for k, c in enumerate(COLORS[:3]):
        ax.plot(np.sin(t + 2 * np.pi * k / 3), y, lw=3.0, color=c)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-0.15, 1.18)
    ax.set_title("三根：解到一半会分岔", fontsize=12)

    # 分叉示意
    ax.annotate("", xy=(-0.55, 1.10), xytext=(0, 1.03),
                arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#4a5568"))
    ax.annotate("", xy=(0.55, 1.10), xytext=(0, 1.03),
                arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#4a5568"))
    ax.text(-0.62, 1.11, "先拆这边", fontsize=8.5, ha="right", color="#4a5568")
    ax.text(0.62, 1.11, "先拆那边", fontsize=8.5, ha="left", color="#4a5568")
    ax.text(0, -0.12, "可以停在「两根已分开、第三根还别着」",
            ha="center", fontsize=10.5, color="#2d3748")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    # ── C：路线数
    ax = fig.add_subplot(gs[0, 2])
    d = np.arange(0, 6)
    b2 = np.ones_like(d, dtype=float)
    b3 = np.array([1, 1, 1, 2, 3, 4], dtype=float)
    ax.plot(d, b2, "-o", lw=2.4, ms=8, color="#718096", label="两根")
    ax.plot(d, b3, "-s", lw=2.4, ms=8, color=COLORS[0], label="三根")
    ax.axvline(2.5, color="#cbd5e0", ls="--", lw=1.2)
    ax.text(2.55, 3.6, "从这一步起\n出现岔路", fontsize=8.5, color="#4a5568")
    ax.set_xlabel("缠的长度")
    ax.set_ylabel("同一状态有几条解绑路线")
    ax.set_title("分叉从「长度 3」开始", fontsize=12)
    ax.set_ylim(0.6, 4.4)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.25)

    fig.suptitle("两根绳子和三根绳子：解绑的自由度不一样", fontsize=14.5, y=1.02)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


def fig_coarse_grain(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig = plt.figure(figsize=(15.5, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.2)

    t = np.linspace(0, 6 * np.pi, 700)
    y = np.linspace(0, 1, 700)

    # ── A：四根平等
    ax = fig.add_subplot(gs[0, 0])
    for k, c in enumerate(COLORS):
        ax.plot(np.sin(t + 2 * np.pi * k / 4) * 0.95, y, lw=2.8, color=c)
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-0.15, 1.18)
    ax.set_title("四根：595 种不同的缠法", fontsize=12)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    # ── B：一根散开，三根成束
    ax = fig.add_subplot(gs[0, 1])
    # 三根紧贴成束（细小的相位差，看得出是三股）
    for k, c in enumerate(COLORS[:3]):
        ax.plot(0.07 * np.sin(t + 1.6 * k), y, lw=7.0, color=c, alpha=0.85)
    # 第四根在外侧绕
    ax.plot(0.95 * np.sin(t), y, lw=2.8, color=COLORS[3])
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-0.15, 1.18)
    ax.set_title("一根散开，三根编成一束", fontsize=12)
    ax.text(0.0, -0.14, "细股 + 粗束", ha="center", fontsize=10.5, color="#2d3748")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    # ── C：远看
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(np.sin(t) * 0.14, y, lw=18.0, color="#4a5568", alpha=0.8, label="粗束（三根）")
    ax.plot(0.95 * np.sin(t), y, lw=2.8, color=COLORS[3], label="细的那根")
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-0.15, 1.18)
    ax.set_title("远看：两根", fontsize=12)
    ax.text(0.0, -0.14, "595 种 → 只能分出 9 种读数",
            ha="center", fontsize=10.5, color="#c53030")
    ax.legend(fontsize=8.5, loc="upper left", framealpha=0.9)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    fig.suptitle("粗粒化：三根当成一根粗绳，四根就伪装成两根",
                 fontsize=14.5, y=1.02)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser(description="故事配图")
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "..", "figures"))
    args = ap.parse_args()
    outdir = os.path.abspath(args.outdir)

    fig_two_vs_three(os.path.join(outdir, "43_story_two_vs_three.png"))
    fig_coarse_grain(os.path.join(outdir, "44_story_coarse_grain.png"))


if __name__ == "__main__":
    main()
