#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
面条化：不是物质被拉长，是成像被拉长
=======================================

原话（2026-10-26）：
    **"还有什么面条化，我觉得都是观察方向的光子因为所在方向力小了，
      所以流速变慢了。"**

标准物理：潮汐力把物体真的拉长（径向拉力差）。
这里：**物体本身没动，是它各部分的 photons 速度不同，
      到达时间被拉开，成像就长了。**

这个说法在《构造》自己的设定下能算，因为底本早就写了：

    · **光速是粉尘速度，不是宇宙上限**
    · **光在横截面螺旋走 → 光速是旋转带来的**
    · **黑 = 那儿没有投影可看**

模型：
    坡度       |∇z|(r)           （坑底最小，最陡圈最大）
    局部光速   v(r) ∝ |∇z|(r)    （力小 → 光慢）
    成像放大   dL_app / dL_true = v_max / v(r)

    ★ 关键推论：
      · 坑底 |∇z| → 0  →  v → 0  →  放大率 → ∞（无限拉长）
      · 但坑底 v=0 意味着**那里没有光出来** → 黑洞的"黑"
      · 所以能看到的是**坑底外围那一圈**，那里放大率极大 → 面条化

    这跟标准面条化的区别：
      标准：物质密度被拉伸（真的变形）
      这里：物质不动，是像被拉伸（纯成像效应）

用法：
    python noodle_imaging.py
    python noodle_imaging.py --plot ../figures/33_noodle.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np

R_C = 1.0     # 核心半径（归一化）


def slope(r):
    """Plummer 坑的坡度，r 以 rc 为单位，取 A = rc。"""
    return np.where(r > 0, r / (1.0 + r * r) ** 1.5, 0.0)


def light_speed(r):
    """局部光速 ∝ 坡度。峰值归一化为 1。"""
    s = slope(r)
    # 最陡处 r = 1/sqrt(2) ≈ 0.707
    s_max = slope(np.array([1 / math.sqrt(2)]))[0]
    return s / s_max


def magnification(r):
    """成像放大率 = v_max / v(r)。"""
    v = light_speed(r)
    with np.errstate(divide="ignore"):
        return np.where(v > 1e-12, 1.0 / v, np.inf)


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("面条化：成像被拉长，不是物质被拉长")
    add("=" * 78)
    add("")
    add("  r/rc  |  坡度  |  局部光速  |  成像放大率")
    add("  ------+--------+------------+-------------")
    for r in (0.02, 0.05, 0.1, 0.2, 0.4, 0.707, 1.0, 1.5, 2.0, 3.0):
        s = float(slope(np.array([r]))[0])
        v = float(light_speed(np.array([r]))[0])
        mag = 1.0 / v if v > 1e-12 else float("inf")
        note = ""
        if r < 0.1:
            note = "  ← 坑底附近，几乎拉成无限长"
        elif abs(r - 0.707) < 0.01:
            note = "  ← 最陡圈，放大率最小"
        add("  %5.3f | %6.3f | %10.3f | %11.2f%s" % (r, s, v, mag, note))

    add("")
    add("-" * 78)
    add("三条推论")
    add("-" * 78)
    add("  1. **拉长发生在「漏出去之前」**：光还在膜面上时，")
    add("     坡度小的地方光速慢，成像被拉长。")
    add("     （作者 10-26 原话：「没漏出去之前就是会拉长」）")
    add("")
    add("  ⚠️ 但**「黑」不是这一段的事**：底本里黑洞是**针眼通道**，")
    add("     光被带出膜面、投影消失——那发生在**应力破维**之后，")
    add("     由应力临界决定，不由坡度决定。本图不覆盖那一段。")
    add("")
    add("  2. **坑底外围被拉长**：放大率在靠近坑底时发散。")
    add("     一个球形物体，靠近坑底的那一面放大率远大于远端，")
    add("     **成像被拉成面条**。")
    add("")
    add("  3. **最陡圈反而不拉长**：r≈0.7rc 处放大率最小（=1）。")
    add("     所以面条化有一个**方向性**：朝向坑底的方向被拉长，")
    add("     沿着等高线的方向不被拉长。")
    add("")

    add("-" * 78)
    add("和标准面条化怎么区分")
    add("-" * 78)
    add("  标准（潮汐力）：")
    add("    · 物质真的被拉长，密度下降，结构被破坏")
    add("    · 拉伸方向 = 径向（朝黑洞）")
    add("    · 拉伸是**物理过程**，需要时间")
    add("")
    add("  这里（成像）：")
    add("    · 物体本身不动，只是不同部分的 photons 到达时间不同")
    add("    · 拉伸方向 = **坡度最缓的方向**（靠近坑底的一侧）")
    add("    · 是**瞬时几何效应**，不需要时间演化")
    add("")
    add("  ★ 可分辨的地方：")
    add("    潮汐拉伸会改变物体的**光谱**（密度变了、温度变了）；")
    add("    成像拉伸只改变**几何**，光谱不变（只是红移被拉长）。")
    add("    不过这在黑洞附近几乎测不到——两边都淹没在强引力红移里。")

    add("")
    add("-" * 78)
    add("这个模型的问题")
    add("-" * 78)
    add("  1. **v ∝ 坡度是我假设的**，底本里只说了「光速是粉尘速度」，")
    add("     没说它和坡度成什么关系。我取线性是为了算得出来。")
    add("  2. 如果 v 和坡度的关系不是线性，放大率曲线会变，")
    add("     但「坑底 v→0 → 无限拉长」这个定性结论不变——")
    add("     只要 v 随坡度单调增。")
    add("  3. 真实的黑洞成像里，还有引力透镜、多普勒增亮、红移，")
    add("     这些都没算进去。")
    return "\n".join(out)


def make_plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))

    rs = np.linspace(0.005, 3.0, 600)

    # 1) 坡度
    ax = axes[0]
    ax.plot(rs, slope(rs), lw=2.6, color="#2b6cb0")
    ax.set_xlabel("r / rc")
    ax.set_ylabel("坡度")
    ax.set_title("坑的坡度\n坑底为 0，最陡圈在 0.71", fontsize=10)
    ax.grid(alpha=0.25)
    ax.axvline(0.707, color="gray", ls=":", lw=1.2)

    # 2) 局部光速
    ax = axes[1]
    ax.plot(rs, light_speed(rs), lw=2.6, color="#dd6b20")
    ax.set_xlabel("r / rc")
    ax.set_ylabel("局部光速（归一化）")
    ax.set_title("光速 ∝ 坡度\n坑底光速为 0 = 黑", fontsize=10)
    ax.grid(alpha=0.25)
    ax.axvspan(0, 0.1, color="black", alpha=0.25)
    ax.text(0.15, 0.5, "这里没有光\n（黑）", fontsize=9, color="black")

    # 3) 成像放大率
    ax = axes[2]
    mag = magnification(rs)
    mag = np.clip(mag, 0, 30)
    ax.plot(rs, mag, lw=2.6, color="#c53030")
    ax.set_xlabel("r / rc")
    ax.set_ylabel("成像放大率")
    ax.set_title("放大率在坑底发散\n= 面条化", fontsize=10)
    ax.grid(alpha=0.25)
    ax.set_ylim(0, 12)
    ax.axvline(0.707, color="gray", ls=":", lw=1.2)
    ax.text(0.75, 10, "最陡圈\n放大率最小", fontsize=8.5, color="gray")

    fig.suptitle("面条化 = 光速变慢造成的成像拉伸，不是物质被拉长", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="面条化的成像解释")
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
