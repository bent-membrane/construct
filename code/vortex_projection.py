#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
漩涡的投影：为什么黑洞看起来像圆
===================================

底本（`2026-10-11_宇宙模型续篇.md` §1.1，10-11 深夜）：
    **黑洞 = 螺旋盘面 + 一扇门**
    **"面上坠落" = 树叶在漩涡里螺旋（面把"穿门"稀释成漫游，每圈近一点）**
    **漩涡 = 旋转 + 凹陷的儿子**

底本（`2026-10-23_投影与重叠.md`）：
    **泡泡大小分布应该跟螺距成比例**
    **泡泡有最小尺寸，中心就不会无限尖**

所以黑洞周围的「环」不是静态圆盘被压扁——
它是**螺旋下降通道的一圈**。决定它长什么样的，不是倾角，是**螺距**。

模型：
    对数螺旋      r(θ) = r0 · exp(b·θ)
    螺距角         α = arctan(b)        （b 小 = 缠得紧）
    沿轴下降       z = -h·θ
    投影（倾角 i） x' = x , y' = y·cos i − z·sin i

关键量：
    **一圈内的半径变化**  Δr/r = exp(2πb) − 1

    · b 小（缠紧）→ 一圈转完半径几乎不变 → **看上去是个圆环**
    · b 大（松）  → 一圈内半径涨一大截 → **看到明显的螺旋臂**

    ★ 所以「黑洞看起来像圆」这个事实，直接给 b 一个上限：
      一圈内的半径涨幅必须小于观测能分辨的不对称度。

用法：
    python vortex_projection.py
    python vortex_projection.py --plot ../figures/32_vortex.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np


def spiral(b, r0=1.0, turns=1.6, n=800):
    theta = np.linspace(0, 2 * np.pi * turns, n)
    r = r0 * np.exp(b * theta)
    return r * np.cos(theta), r * np.sin(theta), theta


def project(x, y, z, inc_deg):
    i = math.radians(inc_deg)
    return x, y * math.cos(i) - z * math.sin(i)


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("漩涡的投影：为什么黑洞看起来像圆")
    add("=" * 78)
    add("")
    add("  对数螺旋 r = r0·exp(b·θ)，螺距角 α = arctan(b)")
    add("")
    add("  b      |  螺距角 α |  一圈半径涨  | 投影看起来")
    add("  -------+-----------+-------------+------------------")
    for b in (0.005, 0.01, 0.02, 0.05, 0.10, 0.20):
        growth = math.exp(2 * math.pi * b) - 1
        if growth < 0.03:
            look = "几乎正圆"
        elif growth < 0.10:
            look = "略扁/略开"
        elif growth < 0.30:
            look = "能看出螺旋"
        else:
            look = "明显螺旋臂"
        add("  %5.3f  |  %5.1f°  |  %7.1f%%   | %s"
            % (b, math.degrees(math.atan(b)), 100 * growth, look))

    add("")
    add("-" * 78)
    add("用 EHT 的观测卡一下")
    add("-" * 78)
    add("  M87*：环的不对称度 < 10%（2019 年 Paper VI）")
    add("        → 一圈内半径涨幅 < 10%")
    add("        → exp(2πb) < 1.10")
    add("        → **b < %.4f，螺距角 < %.2f°**"
        % (math.log(1.10) / (2 * math.pi),
           1000 * math.degrees(math.atan(math.log(1.10) / (2 * math.pi))) / 1000))
    add("")
    add("  这个数意味着：漩涡**缠得非常紧**。")
    add("  一圈转下来，半径几乎不变——所以投影上就是一圈细环。")
    add("")

    add("-" * 78)
    add("这和「缠紧的绳结」是同一个图像")
    add("-" * 78)
    add("  第三章说：四维是绳子被极限压缩、绞死的超螺旋块。")
    add("  这里算出来：黑洞周围的漩涡螺距角 < 0.87°，缠得极紧。")
    add("  **两处对上了**——越靠近黑洞，缠得越紧，螺距越小。")
    add("")
    add("  而且底本自己写过一个推论：")
    add("    「泡泡有最小尺寸，中心就不会无限尖」")
    add("  最小尺寸由螺距定：螺距小 → 泡泡小 → 但**不为零**。")
    add("  → **核心是平的（core），不是尖的（cusp）。**")
    add("")
    add("  ★ 这是《构造》在这个问题上明确的站位：")
    add("    站 **core** 那边。")
    add("    而观测上确实有一批团的核心是平的（cusp-core 问题里 core 一方）。")

    add("")
    add("-" * 78)
    add("但要说清楚还差什么")
    add("-" * 78)
    add("  1. 「不对称度 < 10%」是 EHT 2019 的上限，不是精确测量；")
    add("     后来 Tiede & Broderick 2024 给椭率 0.09±0.07，仍不是硬数字。")
    add("  2. 这个 b 上限是**把整个不对称度都算给螺距**得到的。")
    add("     实际的不对称还来自多普勒增亮（环一边亮一边暗），")
    add("     所以真正的 b 上限比算出来的更松。")
    add("  3. 倾角在这里完全没用上——因为漩涡投影的形状主要由螺距定，")
    add("     不是由倾角定。这跟「椭圆环」的图像正好相反。")
    return "\n".join(out)


def make_plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.4))

    for ax, b, title in zip(
            axes,
            (0.01, 0.06, 0.18),
            ("螺距紧 b=0.01（约 3.6%）\n投影：几乎正圆",
             "中等 b=0.06（约 46%）\n投影：看得出不对称",
             "松 b=0.18（约 210%）\n投影：明显螺旋臂")):
        x, y, th = spiral(b, turns=1.5)
        z = -0.15 * th
        xp, yp = project(x, y, z, 20.0)
        ax.plot(xp, yp, lw=2.2, color="#c53030", alpha=0.9)
        # 标出起点终点
        ax.plot(xp[0], yp[0], "o", ms=8, color="#2b6cb0")
        ax.plot(xp[-1], yp[-1], "s", ms=8, color="#38a169")
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.2)
        ax.set_xlabel("投影 x")
        ax.set_ylabel("投影 y")

    fig.suptitle("漩涡的投影：螺距决定看到什么，不是倾角", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="漩涡的投影")
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
