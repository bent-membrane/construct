#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
四维矢量场拓扑：先分类临界点，再看谁分隔
==========================================

按 Hofmann/Rieck/Sadlo 2018 (CGF 37:301) 的框架重跑。
**不再编地形**——先找临界点，用 Jacobian 特征值分类，判断分隔与否。

《构造》的对应：
    临界点 = 泡泡
    不变流形 = 绳子
    鞍点的余维一流形 = **表里的边界**
    source 类（全排斥/全吸引）= 没有表里边界

分类规则（论文）：
    按实部符号 + 复共轭对数 g（旋转平面数）
    4D 共 14 种情况，归 8 类：
        source / 1:3 saddle / 2:2 saddle
        1-spiral source / 1:spiral-3 saddle / 2:spiral-2 saddle
        2-spiral source / 2-spiral saddle

    分隔的：**所有 saddle 类**（1:3、2:2、1:s3、2:s2）
    不分隔的：**所有 source/sink 类**

本脚本用四个具体场演示：一个纯螺旋、一个螺旋鞍点、一个 1:3 鞍点、一个纯源。

用法：
    python fourd_topology.py
    python fourd_topology.py --plot ../figures/38_fourd_topology.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np


def classify(evals):
    """按实部符号与复共轭对数分类（4D）。"""
    re = np.sort(np.real(evals))
    n_pos = int((re > 1e-9).sum())
    n_neg = int((re < -1e-9).sum())
    g = int((np.abs(np.imag(evals)) > 1e-9).sum() // 2)
    if n_pos > 0 and n_neg == 0:
        base = "source"
    elif n_neg > 0 and n_pos == 0:
        base = "sink"
    else:
        base = "%d:%d saddle" % (n_pos, n_neg)
    if g == 1:
        base = "1-spiral " + base
    elif g == 2:
        base = "2-spiral " + base
    separating = ("saddle" in base)
    return base, g, separating


def make_fields():
    """四个 4x4 矩阵，各代表一类。"""
    fields = {}

    # 1) 纯 2-spiral source：两个平面都在旋转，全部排斥
    a, b = 0.35, 1.0
    fields["2-spiral source"] = np.array([
        [a, -b, 0, 0],
        [b, a, 0, 0],
        [0, 0, a, -b],
        [0, 0, b, a],
    ], dtype=float)

    # 2) 2-spiral saddle：一个平面排斥、一个吸引，都在转
    c = -0.35
    fields["2-spiral saddle"] = np.array([
        [a, -b, 0, 0],
        [b, a, 0, 0],
        [0, 0, c, -b],
        [0, 0, b, c],
    ], dtype=float)

    # 3) 1:3 saddle：一个排斥、三个吸引（无旋转）
    fields["1:3 saddle"] = np.diag([0.6, -0.5, -0.7, -0.9])

    # 4) 1-spiral source：一个平面转，全部排斥
    fields["1-spiral source"] = np.array([
        [0.4, -0.8, 0, 0],
        [0.8, 0.4, 0, 0],
        [0, 0, 0.5, 0],
        [0, 0, 0, 0.6],
    ], dtype=float)

    return fields


def integrate(A, x0, steps=400, dt=0.02):
    """简单欧拉积分（线性场稳）。"""
    xs = np.zeros((steps + 1, 4))
    x = np.array(x0, dtype=float)
    xs[0] = x
    for i in range(steps):
        x = x + dt * (A @ x)
        xs[i + 1] = x
    return xs


def build_report():
    out = []
    add = out.append
    add("=" * 80)
    add("四维矢量场拓扑：临界点分类与分隔判定")
    add("=" * 80)
    add("")
    add("  %-20s %30s %8s" % ("场", "特征值", "分类"))
    add("  " + "-" * 62)

    rows = []
    for name, A in make_fields().items():
        ev = np.linalg.eigvals(A)
        cls, g, sep = classify(ev)
        evs = ", ".join("%+.2f%+.2fi" % (r, i) if abs(i) > 1e-9 else "%+.2f" % r
                        for r, i in zip(np.real(ev), np.imag(ev)))
        add("  %-20s %30s %8s" % (name, evs, cls))
        rows.append((name, cls, g, sep))

    add("")
    add("-" * 80)
    add("分隔判定（《构造》里的「表里边界」）")
    add("-" * 80)
    for name, cls, g, sep in rows:
        mark = "★ **有表里边界**" if sep else "无"
        add("  %-20s %s" % (name, mark))
    add("")
    add("  ★ 规律：**所有 saddle 类分隔，所有 source/sink 类不分隔。**")
    add("    所以「一个地方有没有表里边界」= 「那里的临界点是不是鞍点」。")
    add("    这是可判定的，不需要看质量、速度、光度。")

    add("")
    add("-" * 80)
    add("对应到《构造》")
    add("-" * 80)
    add("  · 纯 source（全排斥）：只是一个「源」，物质被甩出去，没有表里")
    add("  · 2-spiral source：边转边甩出去——**这就是「螺旋搅动」的拓扑身份**")
    add("  · 2-spiral saddle：一个平面转着排斥、另一个平面转着吸引")
    add("                       → **表里边界就在这两个平面之间**")
    add("  · 1:3 saddle：一个方向出、三个方向进，无旋转")
    add("")
    add("  ★ 「螺旋搅动」在 4D 拓扑里不是源，是**鞍点**——")
    add("    因为它同时有出有进，而那正是「表里」的来源。")
    return "\n".join(out)


def make_plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fields = make_fields()
    names = list(fields.keys())

    fig = plt.figure(figsize=(16, 8))
    for k, name in enumerate(names):
        A = fields[name]
        ev = np.linalg.eigvals(A)
        cls, g, sep = classify(ev)

        ax = fig.add_subplot(2, 2, k + 1, projection="3d")
        # 从几个初值积分，投影到 (x,y,z)，丢掉 w
        rng = np.random.default_rng(5)
        starts = rng.normal(0, 0.6, (14, 4))
        for s in starts:
            xs = integrate(A, s, steps=260, dt=0.02)
            # 投影：x,y,z（丢掉 w）
            ax.plot(xs[:, 0], xs[:, 1], xs[:, 2], lw=1.0,
                    color="#c53030" if sep else "#2b6cb0", alpha=0.75)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        tag = "★ 分隔" if sep else "不分隔"
        ax.set_title("%s\n%s   [%s]" % (name, cls, tag), fontsize=10)
        ax.tick_params(labelsize=7)

    fig.suptitle("四维矢量场：临界点分类与分隔判定（投影丢掉 w 维）\n"
                 "红色 = 分隔（有表里边界）　蓝色 = 不分隔", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="4D 矢量场拓扑")
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
