#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
斜坡上的力怎么分解：向上 vs 向左右
=====================================

原话（2026-10-26）：
    **"就是在斜坡上有向上和向左或者右的分量。"**

所以底本那句「凹陷/凸起处有张力，且有向观测方向的力矩」
说的**不是转动矩**，是**力的方向分解**：

    平的地方：张力纯向上（第四维）→ 只有引力，膜面内没有东西
    斜坡上　：法向倾斜 → 张力一分为二
                 向上  = 第四维分量  = 引力
                 向左右 = 膜面内分量 = 指向坑底的拖拽

关键问题：
    这个膜面内分量，和标准引力的径向分量是不是同一件事？
    如果是，等价；如果是不同比例，那就是可分辨的差异。

模型：
    Plummer 凹陷   z(r) = -A / sqrt(1 + (r/rc)²)
    坡度           |∇z| = A·r / (rc²·(1+r²/rc²)^{3/2})
    法向倾角       θ = arctan(|∇z|)
    分量比         面内 / 向上 = |∇z|

对照标准引力：
    Plummer 势     Φ(r) = -GM / sqrt(r²+rc²)
    加速度         g(r) = GM·r / (r²+rc²)^{3/2}   （指向中心）
    这是一个**纯径向**的力，没有"向上"的分量。

    ★ 于是出现一个可以直接对的东西：
      《构造》说 面内分量 / 向上分量 = |∇z|
      标准引力说 全部是面内（径向），向上分量为 0

    如果《构造》的"向上分量"就是引力，那两边描述的是同一个力，
    只是分解方式不同——但**比例**可以被检验。

用法：
    python membrane_tilt.py
    python membrane_tilt.py --plot ../figures/29_tilt.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np

# 归一化单位：rc = 1
A_DEPTH = 1.0      # 凹陷深度（单位 rc）


def dzdr(r):
    """Plummer 凹陷的径向坡度 |∇z|，r 以 rc 为单位。"""
    return A_DEPTH * r / (1.0 + r * r) ** 1.5


def tilt_deg(r):
    return math.degrees(math.atan(dzdr(r)))


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("斜坡上的力分解：向上 vs 向左右")
    add("=" * 78)
    add("凹陷：Plummer z = -A/sqrt(1+(r/rc)²)，取 A = rc")
    add("")
    add("  r/rc |  坡度 |∇z| |  倾角 θ |  面内/向上 |  说明")
    add("  -----+------------+----------+------------+------------------")
    for r in (0.0, 0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0):
        s = dzdr(r)
        note = ""
        if abs(r - 0.0) < 1e-9:
            note = "坑底，完全垂直"
        elif abs(r - 0.707) < 0.4:
            note = "最陡处附近"
        elif r >= 3.0:
            note = "远处，几乎又平了"
        add("  %4.1f | %10.3f | %7.1f° | %10.3f | %s"
            % (r, s, tilt_deg(r), s, note))

    add("")
    # 找最陡点
    rs = np.linspace(0.01, 5, 2000)
    slopes = np.array([dzdr(float(r)) for r in rs])
    imax = int(np.argmax(slopes))
    add("  ★ 最陡点在 r = %.2f rc，倾角 %.1f°，面内/向上 = %.3f"
        % (rs[imax], math.degrees(math.atan(slopes[imax])), slopes[imax]))
    add("")

    add("-" * 78)
    add("和标准引力对一下")
    add("-" * 78)
    add("  标准引力（Plummer 势）：")
    add("    g(r) = GM·r/(r²+rc²)^{3/2}   —— 纯径向，指向中心")
    add("    **没有「向上」的分量。全部都在膜面内。**")
    add("")
    add("  《构造》：")
    add("    向上分量 / 面内分量 = 1 / |∇z|")
    add("    在坑底：向上占 100%（面内为 0）")
    add("    在最陡处：面内占 %.0f%%，向上占 %.0f%%"
        % (100 * slopes[imax] / math.hypot(1, slopes[imax]),
           100 * 1 / math.hypot(1, slopes[imax])))
    add("    在远处：又回到几乎纯向上")
    add("")
    add("  ★ 这里有一处**根本性的不同**：")
    add("    标准引力在任何地方都是**纯面内**（径向），")
    add("    《构造》里**面内分量只是斜坡上的副产品**，坑底和远处都消失。")
    add("")
    add("    换句话说：")
    add("      标准引力：质量在哪里，力就从哪里发出，径向。")
    add("      《构造》  ：力来自膜的法向倾斜，坑底（最陡为零）反而没有面内分量。")
    add("")
    add("    ★ 可检验的差别：**在凹陷最深处，膜面内的力应该趋于零。**")
    add("      标准引力在那里是**最强**的（r=0 处 g=0，但这个零是因为")
    add("      对称性，不是因为没有质量）。")

    add("")
    add("-" * 78)
    add("但这需要小心")
    add("-" * 78)
    add("  1. Plummer 势在 r=0 处引力也是 0（球对称，内部无质量贡献的抵消），")
    add("     所以「坑底无力」这一点两边**都能对上**，不构成区分。")
    add("  2. 真正的区别在**比例**：构造给的是 面内/向上 = |∇z|，")
    add("     标准给的是 面内/向上 = ∞（没有向上分量）。")
    add("     但「向上分量」在标准物理里根本不存在——没法比。")
    add("  3. 所以这个分解本身**不能直接对照观测**，除非能找到")
    add("     「向上分量」的独立观测效应（比如它是否影响光的传播）。")
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

    # 左：剖面 + 力分解
    ax = axes[0]
    r = np.linspace(-3, 3, 400)
    z = -A_DEPTH / np.sqrt(1 + r ** 2)
    ax.plot(r, z, lw=2.6, color="#2b6cb0", label="膜凹陷 z(r)")

    for rr in (0.0, 0.7, 1.5, 2.5):
        s = dzdr(abs(rr)) * np.sign(rr if rr != 0 else 1)
        zz = -A_DEPTH / math.sqrt(1 + rr * rr)
        # 法向
        nx = -s
        nz = 1.0
        norm = math.hypot(nx, nz)
        nx, nz = nx / norm, nz / norm
        ax.annotate("", xy=(rr + nx * 0.8, zz + nz * 0.8),
                    xytext=(rr, zz),
                    arrowprops=dict(arrowstyle="->", color="#c53030", lw=2))
        # 分解
        ax.plot([rr, rr], [zz, zz + nz * 0.55], color="#dd6b20", lw=2.2, ls="--")
        ax.plot([rr, rr + nx * 0.55], [zz, zz], color="#38a169", lw=2.2, ls="--")

    ax.set_xlabel("r / rc")
    ax.set_ylabel("膜高度 z")
    ax.set_title("法向在斜坡上倾斜\n红=张力方向  橙=向上(引力)  绿=向左右(面内)", fontsize=10)
    ax.grid(alpha=0.25)
    ax.set_ylim(-1.6, 0.35)

    # 右：倾角 vs r
    ax = axes[1]
    rr = np.linspace(0.01, 5, 500)
    th = [tilt_deg(float(x)) for x in rr]
    ax.plot(rr, th, lw=2.6, color="#c53030")
    ax.axhline(45, color="gray", ls="--", lw=1.2, label="45° 线")
    imax = int(np.argmax(th))
    ax.plot(rr[imax], th[imax], "o", ms=10, color="#dd6b20", mec="black")
    ax.annotate("最陡：%.1f°\n(r=%.2f rc)" % (th[imax], rr[imax]),
                xy=(rr[imax], th[imax]), xytext=(rr[imax] + 0.8, th[imax] - 8),
                fontsize=10,
                arrowprops=dict(arrowstyle="->", color="#dd6b20", lw=1.2))
    ax.set_xlabel("r / rc")
    ax.set_ylabel("法向倾角 (°)")
    ax.set_title("斜坡越陡，面内分量越大\n坑底和远处都回到纯垂直")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig.suptitle("斜坡上的力怎么分：向上（引力）vs 向左右（面内）", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="斜坡上的力分解")
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
