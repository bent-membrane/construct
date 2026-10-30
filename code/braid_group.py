#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
按群论重算 B2 vs B3 —— 零手调参数
====================================
2026-10-27，要求：「按群论再算一遍」

昨天（10-26）图 24 用的能量景观是手调的：
    E(θ1,θ2) = k(θ1²+θ2²) + λ[1 − cos 3(θ1−θ2)],  λ = 0.35
结论（B2 单极小、B3 多极小）方向对，但那个 λ 是手调的，不是群算出来的。

这次不用能量函数。只问群本身：

    B2 = ⟨σ⟩  ≅ ℤ
        一个生成元、零关系。
        每个元素唯一写成 σ^n。
        从任何元素走到单位元，只有一条约化路径。

    B3 = ⟨σ1, σ2 | σ1σ2σ1 = σ2σ1σ2⟩
        两个生成元、一个关系。
        同一个辫有多种写法（等价类）。
        从某个元素走到单位元，可以有多条等长路径。
        —— 路径分叉的地方，就是「中间站」。

算的量全部由群结构决定，没有可调参数：
    1. 每个长度上有多少个元素（等价类）
    2. 每个元素到单位元的最短路径数
    3. 有多条最短路径的元素占比
    4. 多条路径经过的不同中间态

B2 的答案：路径数恒为 1 —— 全或无。
B3 的答案：从长度 3 开始出现分叉 —— 可以卡在半解开。

注意（2026-10-26 修正过一次）：
    不要再说「3 是奇数 → 抵消不干净 → 剩净环流」。
    三股对称时三个角动量矢量互成 120°，和恰好为零。
    这里算的是「路径分不分叉」，不是「转不转」。

用法：
    python braid_group.py
    python braid_group.py --plot ../figures/40_braid_group.png
"""

import argparse
import sys
from collections import defaultdict, deque
from functools import lru_cache

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ═══════════════════════════════════════════════════════════════
#  词与等价类
# ═══════════════════════════════════════════════════════════════

B2_GENS = ((1, 1), (1, -1))
B3_GENS = ((1, 1), (1, -1), (2, 1), (2, -1))

# 生成元写成 (编号, 符号)：σ1 → (1, +1)，σ1⁻¹ → (1, -1)，σ2 → (2, +1) …


def free_reduce(word):
    """自由群约化：相邻的 σ σ⁻¹ 相消。"""
    out = []
    for g in word:
        if out and out[-1][0] == g[0] and out[-1][1] == -g[1]:
            out.pop()
        else:
            out.append(g)
    return tuple(out)


def braid_step(word):
    """所有能由一次辫关系 σ1σ2σ1 ↔ σ2σ1σ2 得到的词（正逆生成元各自成立）。"""
    res = set()
    for i in range(len(word) - 2):
        a, b, c = word[i], word[i + 1], word[i + 2]
        if a == (1, 1) and b == (2, 1) and c == (1, 1):
            res.add(word[:i] + ((2, 1), (1, 1), (2, 1)) + word[i + 3:])
        elif a == (2, 1) and b == (1, 1) and c == (2, 1):
            res.add(word[:i] + ((1, 1), (2, 1), (1, 1)) + word[i + 3:])
        elif a == (1, -1) and b == (2, -1) and c == (1, -1):
            res.add(word[:i] + ((2, -1), (1, -1), (2, -1)) + word[i + 3:])
        elif a == (2, -1) and b == (1, -1) and c == (2, -1):
            res.add(word[:i] + ((1, -1), (2, -1), (1, -1)) + word[i + 3:])
    return res


@lru_cache(maxsize=None)
def normalize_b2(word):
    """B2 没有关系，只需自由约化。"""
    return free_reduce(word)


@lru_cache(maxsize=None)
def normalize_b3(word):
    """B3 的等价类代表：自由约化 → 辫关系闭包 → 取字典序最小的词。"""
    word = free_reduce(word)
    seen = {word}
    stack = [word]
    while stack:
        cur = stack.pop()
        for nxt in braid_step(cur):
            nxt = free_reduce(nxt)
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    return min(seen)


# ═══════════════════════════════════════════════════════════════
#  由群结构直接建图：节点 = 等价类，边 = 乘一个生成元
# ═══════════════════════════════════════════════════════════════

def build_graph(normalize, gens, max_depth):
    start = ()
    dist = {start: 0}
    adj = defaultdict(set)
    q = deque([start])

    while q:
        v = q.popleft()
        if dist[v] >= max_depth:
            continue
        for g in gens:
            u = normalize(v + (g,))
            adj[v].add(u)
            adj[u].add(v)
            if u not in dist:
                dist[u] = dist[v] + 1
                q.append(u)

    # 最短路径计数：按距离从小到大累加
    paths = defaultdict(int)
    paths[start] = 1
    for v in sorted(dist, key=lambda x: dist[x]):
        for u in adj[v]:
            if dist.get(u) == dist[v] + 1:
                paths[u] += paths[v]

    return dist, paths, adj


def summarize(dist, paths, max_depth):
    """按深度统计：元素数、路径总数、多路径元素数、最大路径数。"""
    rows = []
    for d in range(max_depth + 1):
        nodes = [v for v in dist if dist[v] == d]
        if not nodes:
            continue
        multi = [v for v in nodes if paths[v] > 1]
        rows.append({
            "depth": d,
            "nodes": len(nodes),
            "multi": len(multi),
            "max_paths": max(paths[v] for v in nodes),
            "total_paths": sum(paths[v] for v in nodes),
        })
    return rows


def word_str(word):
    if not word:
        return "1"
    parts = []
    for gen, sign in word:
        parts.append(f"σ{gen}" + ("" if sign > 0 else "^-1"))
    return "·".join(parts)


def gen_str(g):
    return f"σ{g[0]}" + ("" if g[1] > 0 else "^-1")


def path_str(nodes, gens):
    """把一条路径写成 1 --σ2--> σ2 --σ1--> … 的样子。"""
    parts = [word_str(nodes[0])]
    for g, node in zip(gens, nodes[1:]):
        parts.append(f"--{gen_str(g)}--> {word_str(node)}")
    return " ".join(parts)


def shortest_paths_to(normalize, gens, target, max_depth):
    """
    列出从单位元到 target 的所有最短路径。
    返回 [(节点列表, 生成元列表)]，两者长度差 1，gi 把 nodes[i-1] 送到 nodes[i]。
    """
    dist, _, adj = build_graph(normalize, gens, max_depth)
    if target not in dist:
        return []
    out = []
    stack = [(target, [target], [])]
    while stack:
        v, nodes, gs = stack.pop()
        if v == ():
            # nodes 已经是 () → … → target；这里翻成「解绑」方向 target → 1，
            # 每步乘的生成元取逆（v = u·g  ⇒  u = v·g⁻¹）。
            out.append((list(reversed(nodes)),
                        [(g[0], -g[1]) for g in reversed(gs)]))
            continue
        for u in adj[v]:
            if dist.get(u) == dist[v] - 1:
                for g in gens:
                    if normalize(u + (g,)) == v:
                        stack.append((u, [u] + nodes, [g] + gs))
                        break
    return out


# ═══════════════════════════════════════════════════════════════
#  报告
# ═══════════════════════════════════════════════════════════════

MAX_DEPTH = 5


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("按群论重算 B2 vs B3 —— 零手调参数")
    add("=" * 78)
    add("")

    add("-" * 78)
    add("一、两个群的严格写法")
    add("-" * 78)
    add("  B2 = ⟨σ⟩            （零关系，无限循环群）")
    add("       元素 ↔ 整数 n      唯一约化词 = σⁿ")
    add("       线性序，没有分支。")
    add("")
    add("  B3 = ⟨σ1, σ2 | σ1σ2σ1 = σ2σ1σ2⟩")
    add("       两个生成元，一个关系。")
    add("       等价类代表词 = 自由约化 + 辫关系闭包后字典序最小的那个。")
    add("       非交换：σ1σ2 ≠ σ2σ1。")
    add("")

    add("-" * 78)
    add("二、核心量：到单位元的约化路径数")
    add("-" * 78)
    add("  在『节点=等价类、边=乘一个生成元』的 Cayley 图上：")
    add("  一个元素的路径数 = 它有多少种不同的最短约化写法。")
    add("  路径分叉 = 解绑时可以走不同的中间站。")
    add("")

    dist2, paths2, _ = build_graph(normalize_b2, B2_GENS, MAX_DEPTH)
    rows2 = summarize(dist2, paths2, MAX_DEPTH)
    add("  B2：")
    add("    长度  元素数  多路径元素  最大路径数")
    for r in rows2:
        add("    %4d  %6d  %10d  %10d"
            % (r["depth"], r["nodes"], r["multi"], r["max_paths"]))
    add("    → 每一层都只有 1 个元素、路径数恒为 1：全或无。")
    add("")

    dist3, paths3, adj3 = build_graph(normalize_b3, B3_GENS, MAX_DEPTH)
    rows3 = summarize(dist3, paths3, MAX_DEPTH)
    add("  B3：")
    add("    长度  元素数  多路径元素  最大路径数")
    for r in rows3:
        add("    %4d  %6d  %10d  %10d"
            % (r["depth"], r["nodes"], r["multi"], r["max_paths"]))
    add("")

    # Δ 的例子
    delta = normalize_b3(((1, 1), (2, 1), (1, 1)))
    paths_delta = shortest_paths_to(normalize_b3, B3_GENS, delta, MAX_DEPTH)
    add("  ★ 最小分叉点：Δ = σ1σ2σ1 = σ2σ1σ2")
    add("    到单位元的最短路径数 = %d" % len(paths_delta))
    for i, (nodes, gs) in enumerate(paths_delta, 1):
        add("      路径 %d：%s" % (i, path_str(nodes, gs)))
    add("")
    mids = {nodes[1] for nodes, _ in paths_delta if len(nodes) > 2}
    add("    两条路径的第一个中间站不同：%s"
        % "、".join(sorted(word_str(m) for m in mids)))
    add("    这两个中间态彼此不等价（B3 非交换）——")
    add("    所以三股解到一半，可以停在两种不同的『半解开』上。")
    add("")

    add("-" * 78)
    add("三、和昨天图 24 的关系")
    add("-" * 78)
    add("  图 24 的结论（B2 单极小 / B3 多极小）是对的，")
    add("  但那个多极小是从手调能量 E = kΣθ² + λ[1−cos3Δ] 里看出来的，λ=0.35。")
    add("  这次不设能量：多极小 = 路径分叉，由辫关系自己决定。")
    add("")
    add("  用词说清楚：")
    add("    B2 的 Cayley 图是一条直线 → 每个态只有一条出路 → 只能全解或全不解。")
    add("    B3 的 Cayley 图有环/分叉 → 同一个态有多条等长出路 → 可以停在中间。")
    add("")

    add("-" * 78)
    add("四、边界")
    add("-" * 78)
    add("  · 这里算的是纯辫群（拓扑），不是物理系统。")
    add("    它回答的是『拓扑允不允许卡住』，不回答『实际卡多久』。")
    add("  · 深度截到 %d 层，长辫的等价类更大，但分叉从长度 3 就已经出现，结论不依赖深度。" % MAX_DEPTH)
    add("  · 手性那段不再重算：2026-10-26 已修正——")
    add("    三股对称时角动量抵消，奇偶无关；判据是『散不散』。")
    add("")

    return "\n".join(out)


# ═══════════════════════════════════════════════════════════════
#  图
# ═══════════════════════════════════════════════════════════════

def make_plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    dist3, paths3, adj3 = build_graph(normalize_b3, B3_GENS, MAX_DEPTH)
    rows2 = summarize(*build_graph(normalize_b2, B2_GENS, MAX_DEPTH)[:2], MAX_DEPTH)
    rows3 = summarize(dist3, paths3, MAX_DEPTH)

    fig = plt.figure(figsize=(15, 9.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.36, wspace=0.32)

    # ── 1) B2：一条直线
    ax = fig.add_subplot(gs[0, 0])
    xs = list(range(-3, 4))
    ax.plot(xs, [0] * len(xs), "-o", color="#2b6cb0", lw=2.4, ms=9)
    for x in xs:
        ax.annotate(f"σ^{x}" if x != 1 else "σ", (x, 0), textcoords="offset points",
                    xytext=(0, 12), ha="center", fontsize=8.5, color="#2b6cb0")
    ax.annotate("0", (0, 0), textcoords="offset points", xytext=(0, -18),
                ha="center", fontsize=9, color="#c53030", fontweight="bold")
    ax.set_xlim(-3.6, 3.6)
    ax.set_ylim(-0.6, 0.6)
    ax.set_title("B2 = <σ>  (无限循环群)：一条直线\n每个态只有一条出路", fontsize=10)
    ax.axis("off")

    # ── 2) B3：分层 Cayley 图
    ax = fig.add_subplot(gs[0, 1:])
    by_depth = defaultdict(list)
    for v, d in dist3.items():
        by_depth[d].append(v)
    pos = {}
    for d in sorted(by_depth):
        nodes = sorted(by_depth[d])
        for i, v in enumerate(nodes):
            x = (i - (len(nodes) - 1) / 2) * (1.0 if len(nodes) < 8 else 6.0 / len(nodes))
            pos[v] = (x, -d)
    for v in pos:
        for u in adj3[v]:
            if dist3.get(u) == dist3[v] + 1:
                ax.plot([pos[v][0], pos[u][0]], [pos[v][1], pos[u][1]],
                        color="#cbd5e0", lw=0.7, zorder=1)
    for d in sorted(by_depth):
        nodes = by_depth[d]
        multi = [v for v in nodes if paths3[v] > 1]
        single = [v for v in nodes if paths3[v] == 1]
        if single:
            ax.scatter([pos[v][0] for v in single], [pos[v][1] for v in single],
                       s=46, color="#2b6cb0", zorder=3, edgecolors="white", linewidths=0.5)
        if multi:
            ax.scatter([pos[v][0] for v in multi], [pos[v][1] for v in multi],
                       s=64, color="#e53e3e", zorder=4, edgecolors="white", linewidths=0.5)
        ax.text(-3.9, -d, f"长度 {d}", fontsize=7.5, va="center", color="#4a5568")
    ax.set_title("B3 = <σ1,σ2 | σ1σ2σ1=σ2σ1σ2>：分层 Cayley 图\n"
                 "红点 = 有多条最短路径的态（可以卡在中间）", fontsize=10)
    ax.axis("off")

    # ── 3) 路径数分布
    ax = fig.add_subplot(gs[1, 0])
    d2 = [r["depth"] for r in rows2]
    m2 = [r["max_paths"] for r in rows2]
    d3 = [r["depth"] for r in rows3]
    m3 = [r["max_paths"] for r in rows3]
    ax.plot(d2, m2, "-o", color="#2b6cb0", lw=2.2, ms=7, label="B2")
    ax.plot(d3, m3, "-s", color="#e53e3e", lw=2.2, ms=7, label="B3")
    ax.set_xlabel("词长")
    ax.set_ylabel("最大最短路径数")
    ax.set_title("同一节点有几条等长解绑路线", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    # ── 4) 多路径元素占比
    ax = fig.add_subplot(gs[1, 1])
    share2 = [0.0 if r["nodes"] == 0 else r["multi"] / r["nodes"] for r in rows2]
    share3 = [0.0 if r["nodes"] == 0 else r["multi"] / r["nodes"] for r in rows3]
    w = 0.38
    ax.bar([d - w / 2 for d in d2], share2, width=w, color="#2b6cb0", label="B2")
    ax.bar([d + w / 2 for d in d3], share3, width=w, color="#e53e3e", label="B3")
    ax.set_xlabel("词长")
    ax.set_ylabel("多路径元素占比")
    ax.set_ylim(0, 1.05)
    ax.set_title("分叉出现的位置", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25, axis="y")

    # ── 5) Δ 的两条路
    ax = fig.add_subplot(gs[1, 2])
    delta = normalize_b3(((1, 1), (2, 1), (1, 1)))
    paths_delta = shortest_paths_to(normalize_b3, B3_GENS, delta, MAX_DEPTH)
    for i, (nodes, gs) in enumerate(paths_delta):
        ys = list(range(len(nodes)))
        ax.plot([i] * len(nodes), ys, "-o", lw=2.0, ms=8,
                color="#e53e3e" if i == 0 else "#38a169")
        for j, wd in enumerate(nodes):
            label = word_str(wd) if j == 0 else f"{gen_str(gs[j - 1])}· {word_str(wd)}"
            ax.annotate(label, (i, j), textcoords="offset points",
                        xytext=(10, 0), fontsize=6.8, va="center")
    ax.set_xlim(-0.6, 1.9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["路径 1", "路径 2"])
    ax.set_yticks([])
    ax.set_title("Δ = σ1σ2σ1 = σ2σ1σ2 的两条最短路径\n第一个中间站不同 → 两种『半解开』",
                 fontsize=10)

    fig.suptitle("按群论重算：B2 只有一条路，B3 从长度 3 开始分叉", fontsize=15)
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="B2 vs B3 群论重算")
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
