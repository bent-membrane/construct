#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
回程 = 动量：用优化算法重跑「冲过谷底再荡回来」
================================================

2026-10-26 定的认知框架（最高优先级）：
    **"就记住怎么优化算法，梯度下降就行了。"**

对照表：
    质量        → 无此变量
    引力        → 梯度
    引力势      → 损失函数地形
    物体下落    → 参数沿负梯度更新
    **回程      → 动量冲过谷底再荡回来**
    束缚        → 收敛到某个局部极小
    暗物质      → 为解释「梯度为什么这么算」而虚构的项

所以「回程」根本不需要质量、不需要能量收支——它是**动量项**的必然结果：

    v ← β·v − η·∇L(θ)
    θ ← θ + v

    β = 0      无动量，单调滑向谷底，不回头
    β > 0      有动量，冲过谷底，荡回来
    β 越大     冲得越远

这正是底本里「解螺旋之后会旋回来」的优化版表述。

用法：
    python momentum_return.py
    python momentum_return.py --plot ../figures/34_momentum.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np


def loss(theta):
    """损失函数：一个碗（对称的坑）。"""
    return 0.5 * theta ** 2


def grad(theta):
    return theta


def run(beta, eta=0.02, theta0=3.0, steps=600):
    """动量法轨迹。返回 theta 序列。"""
    theta = theta0
    v = 0.0
    traj = [theta]
    for _ in range(steps):
        v = beta * v - eta * grad(theta)
        theta = theta + v
        traj.append(theta)
    return np.array(traj)


def overshoot(traj):
    """过冲量：越过 0 之后到的最远负值（绝对值）。"""
    neg = traj[traj < 0]
    return abs(neg.min()) if len(neg) else 0.0


def settle_time(traj, tol=0.05):
    """收敛时间：最后一次 |theta| > tol 的位置。"""
    idx = np.where(np.abs(traj) > tol)[0]
    return int(idx[-1]) if len(idx) else 0


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("回程 = 动量（用优化算法重跑）")
    add("=" * 78)
    add("")
    add("  地形：L(θ) = ½θ²，起点 θ₀ = 3.0，学习率 η = 0.02")
    add("")
    add("   β    | 过冲量 | 收敛步数 | 行为")
    add("  ------+--------+----------+------------------")
    for beta in (0.0, 0.3, 0.6, 0.8, 0.9, 0.95):
        tr = run(beta)
        ov = overshoot(tr)
        st = settle_time(tr)
        if beta == 0:
            beh = "单调滑落，从不回头"
        elif ov < 0.05:
            beh = "几乎不过冲"
        elif ov < 1.0:
            beh = "过冲，荡一次就稳"
        else:
            beh = "大幅过冲，来回荡"
        add("  %4.2f  | %6.2f | %8d | %s" % (beta, ov, st, beh))

    add("")
    add("-" * 78)
    add("这就是「回程」")
    add("-" * 78)
    add("  参数点滑到谷底时**并没有停**——它带着动量冲了过去，")
    add("  然后被梯度拽回来，再冲过去一点，最后才收敛。")
    add("")
    add("  **整个过程没有一个「力」，没有「质量」，没有「能量收支」。**")
    add("  只有：地形、梯度、动量。")
    add("")
    add("  ★ 对照今天白算的那些：")
    add("    「回来的能量从哪来？」——优化里没有这个问题。")
    add("    动量是**状态变量**，它自己带着速度，不需要谁来供能。")
    add("")
    add("  ★ 底本的「解螺旋之后会旋回来」，在这个语言里就是：")
    add("    **β ≠ 0 的系统一定会过冲，过冲就是回程。**")

    add("")
    add("-" * 78)
    add("顺带：暗物质在这里是什么")
    add("-" * 78)
    add("  假设有个观察者，**看不到地形，只能看到 θ 的轨迹**。")
    add("  它会看到参数的运动会加速、会过冲、会来回荡。")
    add("")
    add("  为了解释这些，它有两种做法：")
    add("    (a) 假设有一股看不见的「拉力」——于是需要一个看不见的「质量」；")
    add("    (b) 承认有地形和动量——不需要任何看不见的东西。")
    add("")
    add("  ★ **暗物质 = (a) 的产物。** 它不是被证伪，是没必要。")
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

    # 左：过冲量 vs β
    ax = axes[0]
    bs = np.linspace(0.0, 0.98, 120)
    ovs = [overshoot(run(float(b), steps=600)) for b in bs]
    ax.plot(bs, ovs, lw=2.8, color="#c53030")
    ax.axhline(0.0, color="black", ls="--", lw=1.0, alpha=0.5)
    for b, c in ((0.0, "#38a169"), (0.6, "#dd6b20"), (0.9, "#c53030")):
        v = overshoot(run(b, steps=600))
        ax.plot([b], [v], "o", ms=10, color=c, mec="black", zorder=6)
    ax.annotate("β<0.75：几乎不过冲\n单调滑落，不回头",
                xy=(0.35, 0.02), xytext=(0.30, 0.9), fontsize=9,
                arrowprops=dict(arrowstyle="->", color="gray", lw=1))
    ax.annotate("β>0.75：开始过冲\n= 回程出现",
                xy=(0.9, 0.89), xytext=(0.45, 1.5), fontsize=9,
                arrowprops=dict(arrowstyle="->", color="#c53030", lw=1.1))
    ax.set_xlabel("动量系数 β")
    ax.set_ylabel("过冲量（越过谷底多远）")
    ax.set_title("过冲量 vs 动量系数\n动量足够大，才「荡回来」", fontsize=10)
    ax.grid(alpha=0.25)

    # 右：θ(t)
    ax = axes[1]
    for beta, color, lbl in ((0.0, "#38a169", "β=0"),
                             (0.6, "#dd6b20", "β=0.6"),
                             (0.9, "#c53030", "β=0.9")):
        tr = run(beta, steps=300)
        ax.plot(tr, lw=1.8, color=color, label=lbl)
    ax.axhline(0, color="black", ls="--", lw=1.2)
    ax.set_xlabel("迭代步")
    ax.set_ylabel("θ")
    ax.set_title("θ 随迭代的变化\n过冲 = 穿过 0 再回来 = 回程", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig.suptitle("回程 = 动量，不需要质量也不需要能量收支", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="回程=动量")
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
