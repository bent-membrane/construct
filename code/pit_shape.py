#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
坑的形状决定坡度：同深同宽，坡度差一倍
=======================================

原话（2026-10-26）：
    **"凹陷深不深和它坑周围坡度有关吗？
      锥子一样和碗一样深，坡度根本不同。"**

对。深度 A 和半径 R 相同，坡度可以完全不同：

    锥形    z = -A(1 - r/R)        坡度 = A/R，**处处相同**
    碗形    z = -A(1 - (r/R)²)     坡度 = 2Ar/R²，边缘最陡
    Plummer z = -A/√(1+(r/rc)²)    坡度先增后减
    高斯    z = -A·exp(-(r/R)²)    中心平、边缘最陡、外面又平

所以：
    · **不能用"凹陷有多深"这一个数描述坑。**
    · 可观测的是**坡度**（它直接决定力的面内/向上之比）
    · 坡度由**形状**定，深度只是形状的一个参数

对《构造》的意义：
    面内/向上 = |∇z| —— 这个比值决定「引力 vs 张力」的占比。
    同样的深度，换个形状，这个比值能差一倍以上。
    所以「反推凹陷深度」是错的提法，应该反推**形状 + 尺度**。

用法：
    python pit_shape.py
    python pit_shape.py --plot ../figures/30_pit_shape.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np

A = 1.0      # 深度
R = 3.0      # 边缘半径

SHAPES = [
    ("锥形", "#c53030"),
    ("碗形(抛物)", "#dd6b20"),
    ("Plummer", "#2b6cb0"),
    ("高斯", "#38a169"),
]


def profile(shape, r, rc=None):
    """返回 z(r)。rc 只对 Plummer 用。"""
    r = np.asarray(r, dtype=float)
    if shape == "锥形":
        return -A * np.clip(1 - r / R, 0, None)
    if shape == "碗形(抛物)":
        return -A * np.clip(1 - (r / R) ** 2, 0, None)
    if shape == "Plummer":
        rc = rc if rc is not None else R / 3.0
        return -A / np.sqrt(1 + (r / rc) ** 2)
    if shape == "高斯":
        return -A * np.exp(-(r / R) ** 2)
    raise ValueError(shape)


def slope(shape, r):
    """|∇z| 的解析式。"""
    r = np.asarray(r, dtype=float)
    if shape == "锥形":
        return np.where(r < R, A / R, 0.0)
    if shape == "碗形(抛物)":
        return np.where(r < R, 2 * A * r / R ** 2, 0.0)
    if shape == "Plummer":
        rc = R / 3.0
        return A * r / (rc ** 2 * (1 + r ** 2 / rc ** 2) ** 1.5)
    if shape == "高斯":
        return 2 * A * r / R ** 2 * np.exp(-(r / R) ** 2)
    raise ValueError(shape)


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("坑的形状决定坡度（同深 A=%.1f，同边缘 R=%.1f）" % (A, R))
    add("=" * 78)
    add("")
    add("  %-12s %12s %12s %10s" %
        ("形状", "最大坡度", "最大倾角", "最陡位置"))
    add("  " + "-" * 52)

    rr = np.linspace(0.001, R * 1.6, 4000)
    maxima = {}
    for name, _ in SHAPES:
        s = slope(name, rr)
        i = int(np.argmax(s))
        maxima[name] = (float(s[i]), rr[i])
        add("  %-12s %12.3f %11.1f° %10.2f" %
            (name, s[i], math.degrees(math.atan(s[i])), rr[i]))

    add("")
    add("  ★ 最大坡度相差 %.1f 倍（%s %.2f vs %s %.2f）"
        % (max(v[0] for v in maxima.values()) / min(v[0] for v in maxima.values()),
           "碗形", maxima["碗形(抛物)"][0],
           "高斯" if maxima["高斯"][0] < maxima["碗形(抛物)"][0] else "锥形",
           min(v[0] for v in maxima.values())))
    add("")

    add("-" * 78)
    add("坡度分布的形状差异（这才是可观测的）")
    add("-" * 78)
    add("  r/R  |  锥形  |  碗形  | Plummer |  高斯")
    add("  -----+--------+--------+---------+--------")
    for r in (0.1, 0.3, 0.5, 0.7, 1.0, 1.3, 1.6):
        add("  %4.1f | %6.3f | %6.3f |  %6.3f | %6.3f"
            % (r,
               float(slope("锥形", r)),
               float(slope("碗形(抛物)", r)),
               float(slope("Plummer", r)),
               float(slope("高斯", r))))
    add("")
    add("  · 锥形：坡度**处处相同** → 面内/向上 在整条斜面上是常数")
    add("  · 碗形：坡度随 r **线性增加** → 越靠边，面内占比越大")
    add("  · Plummer：先增后减 → 有一个「最陡圈」")
    add("  · 高斯：和 Plummer 类似，但外面衰减更快")
    add("")

    add("-" * 78)
    add("这对「反推」意味着什么")
    add("-" * 78)
    add("  可观测的是**坡度**（=面内/向上之比 = 引力与张力的占比）。")
    add("  但同一个坡度，可以由：")
    add("    · 一个很深的锥子，")
    add("    · 或一个很浅但很陡的碗，")
    add("  给出来——**两者不可区分。**")
    add("")
    add("  所以只能反推「坡度分布」，反推不出「深度」。")
    add("  深度是形状的参数，不是观测量。")
    add("")
    add("  ★ 底本里有一条可能是线索：**「暗物质 = 褶皱」**。")
    add("    褶皱是**尖的、不规则**的——形状上更接近锥形/折线，")
    add("    不是光滑的碗或 Plummer 球。")
    add("    如果按「褶皱」来定形状，最陡处会落在折痕上，")
    add("    坡度可以远高于光滑模型的 0.385。")
    return "\n".join(out)


def make_plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6))

    rr = np.linspace(0, R * 1.6, 500)

    ax = axes[0]
    for name, color in SHAPES:
        ax.plot(rr, profile(name, rr), lw=2.4, color=color, label=name)
    ax.set_xlabel("r")
    ax.set_ylabel("z")
    ax.set_title("四个坑：深度一样、边缘一样\n但形状完全不同", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)
    ax.set_ylim(-1.15, 0.1)

    ax = axes[1]
    for name, color in SHAPES:
        ax.plot(rr, slope(name, rr), lw=2.4, color=color, label=name)
    ax.set_xlabel("r")
    ax.set_ylabel("坡度（单位深度/单位半径）")
    ax.set_title("坡度分布：同一个深度，峰值差 %.1f 倍" %
                 (max(float(np.max(slope(n, rr))) for n, _ in SHAPES) /
                  min(float(np.max(slope(n, rr))) for n, _ in SHAPES)),
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig.suptitle("坑的形状决定坡度 —— 深度不是关键参数", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="坑的形状与坡度")
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
