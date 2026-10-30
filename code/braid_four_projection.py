#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
四股辫的粗粒化：一根 + 三根看起来像两根
==========================================
2026-10-27，观察：

    「四股辫如果一根散开了，也可能和另外三根编在一起的成为双螺旋，
      也就是说一个单体和三体也能当成两个独立个体看」

群论上的操作叫粗粒化（投影）：

    π: B4 → B2 ,  σ1 ↦ σ, σ2 ↦ σ, σ3 ↦ σ

投影是合法同态（B4 的辫关系投过去仍然成立）。
投影之后只剩「总共绕了几圈」一个数 —— 就是两股的样子。

被丢掉的是投影的核：三根之间怎么编的那些细节。
本脚本量化一件事：

    从两体的观测反推四体，到底丢了多少信息？

具体算：
    1. B4 每一层有多少个元素（等价类）
    2. 这些元素投影到 B2 后，每个 B2 元素底下压着多少个 B4 元素
       （原像越大 = 两体观测越看不见四体的差别）
    3. 核（投影到 1 的那一堆）有多大
    4. 不含 σ3 的子群 ≅ B3 —— 就是「三根编成一根」的构型

B4 = ⟨σ1,σ2,σ3 | σ1σ2σ1=σ2σ1σ2,  σ2σ3σ2=σ3σ2σ3,  σ1σ3=σ3σ1⟩

注意 σ1 与 σ3 对易（不相交的两对），这是"两头可以独立"的硬依据；
但 σ2 把两边连起来，所以 1+3 的独立只是部分的，不是直积。

用法：
    python braid_four_projection.py
    python braid_four_projection.py --plot ../figures/41_four_projection.png
    python braid_four_projection.py --depth 5
"""

import argparse
import sys
from collections import Counter, defaultdict, deque
from functools import lru_cache

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


B4_GENS = ((1, 1), (1, -1), (2, 1), (2, -1), (3, 1), (3, -1))


def free_reduce(word):
    out = []
    for g in word:
        if out and out[-1][0] == g[0] and out[-1][1] == -g[1]:
            out.pop()
        else:
            out.append(g)
    return tuple(out)


def braid_step_b4(word):
    """B4 的三种关系：σ1σ2σ1=σ2σ1σ2、σ2σ3σ2=σ3σ2σ3、σ1σ3=σ3σ1。"""
    res = set()
    n = len(word)

    # 相邻对易：σ1 与 σ3
    for i in range(n - 1):
        a, b = word[i], word[i + 1]
        if {a[0], b[0]} == {1, 3}:
            res.add(word[:i] + (b, a) + word[i + 2:])

    # 三段辫关系
    for i in range(n - 2):
        a, b, c = word[i], word[i + 1], word[i + 2]
        for x, y, z in (((1, 1), (2, 1), (1, 1)), ((2, 1), (1, 1), (2, 1)),
                        ((2, 1), (3, 1), (2, 1)), ((3, 1), (2, 1), (3, 1)),
                        ((1, -1), (2, -1), (1, -1)), ((2, -1), (1, -1), (2, -1)),
                        ((2, -1), (3, -1), (2, -1)), ((3, -1), (2, -1), (3, -1))):
            if (a, b, c) == (x, y, z):
                swap = ((2, 1), (1, 1), (2, 1)) if x == (1, 1) else \
                       ((1, 1), (2, 1), (1, 1)) if x == (2, 1) and y == (1, 1) else \
                       ((3, 1), (2, 1), (3, 1)) if x == (2, 1) else \
                       ((2, 1), (3, 1), (2, 1)) if x == (3, 1) else \
                       ((2, -1), (1, -1), (2, -1)) if x == (1, -1) else \
                       ((1, -1), (2, -1), (1, -1)) if x == (2, -1) and y == (1, -1) else \
                       ((3, -1), (2, -1), (3, -1)) if x == (2, -1) else \
                       ((2, -1), (3, -1), (2, -1))
                res.add(word[:i] + swap + word[i + 3:])
    return res


@lru_cache(maxsize=None)
def normalize_b4(word):
    word = free_reduce(word)
    seen = {word}
    stack = [word]
    while stack:
        cur = stack.pop()
        for nxt in braid_step_b4(cur):
            nxt = free_reduce(nxt)
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    return min(seen)


def projection(word):
    """π: σi ↦ σ，所以投影就是指数和（B4 的 Abel 化）。"""
    return sum(sign for _, sign in word)


def build_b4(max_depth):
    start = ()
    dist = {start: 0}
    q = deque([start])
    while q:
        v = q.popleft()
        if dist[v] >= max_depth:
            continue
        for g in B4_GENS:
            u = normalize_b4(v + (g,))
            if u not in dist:
                dist[u] = dist[v] + 1
                q.append(u)
    return dist


def analyze(max_depth):
    dist = build_b4(max_depth)
    by_depth = defaultdict(list)
    for v, d in dist.items():
        by_depth[d].append(v)

    rows = []
    for d in sorted(by_depth):
        nodes = by_depth[d]
        proj = Counter(projection(v) for v in nodes)
        rows.append({
            "depth": d,
            "nodes": len(nodes),
            "proj": proj,
            "with3": sum(1 for v in nodes if any(g[0] == 3 for g in v)),
            "no3": sum(1 for v in nodes if not any(g[0] == 3 for g in v)),
        })

    # 全深度累计的原像大小
    all_proj = Counter(projection(v) for v in dist)
    # 核：投影到 0 的元素（纯辫那块看不见的部分）
    kernel = all_proj.get(0, 0)

    return dist, rows, all_proj, kernel


def build_report(max_depth):
    out = []
    add = out.append
    add("=" * 78)
    add("四股辫的粗粒化：一根 + 三根看起来像两根")
    add("=" * 78)
    add("")
    add("  π: B4 → B2 ,  σ1↦σ, σ2↦σ, σ3↦σ")
    add("  投影 = 指数和（总共绕了几圈）—— 这就是两体观测能看到的那一个数。")
    add("")

    dist, rows, all_proj, kernel = analyze(max_depth)

    add("-" * 78)
    add("一、每一层的情况")
    add("-" * 78)
    add("  词长  B4元素数  含σ3  不含σ3(≅B3子群)  投影分布(指数和→个数)")
    for r in rows:
        pj = " ".join(f"{k:+d}:{v}" for k, v in sorted(r["proj"].items()))
        add("  %4d %9d %6d %14d      %s"
            % (r["depth"], r["nodes"], r["with3"], r["no3"], pj))
    add("")

    add("-" * 78)
    add("二、粗粒化到底丢了什么")
    add("-" * 78)
    add("  累计（到深度 %d）：" % max_depth)
    add("    B4 元素总数        = %d" % len(dist))
    add("    投影到 B2 的像     = %d 个不同的指数和" % len(all_proj))
    add("    ★ 平均每个 B2 点数 = %.1f 个 B4 元素" % (len(dist) / len(all_proj)))
    add("    核（指数和 = 0）   = %d 个元素完全看不见" % kernel)
    add("")
    add("  逐个指数和的原像大小（这就是『看起来一样』的族）：")
    for n in sorted(all_proj):
        bar = "#" * min(60, all_proj[n])
        add("    %+3d : %4d  %s" % (n, all_proj[n], bar))
    add("")

    add("-" * 78)
    add("三、『一根散开 + 三根编成一根』的构型")
    add("-" * 78)
    total_no3 = sum(r["no3"] for r in rows)
    add("  不含 σ3 的元素（股 4 完全不参与）= %d / %d = %.1f%%"
        % (total_no3, len(dist), 100.0 * total_no3 / len(dist)))
    add("  这些元素全部落在 B3 子群 ⟨σ1,σ2⟩ 里 —— 三根内部怎么编都在，")
    add("  但股 4 一参与，它们就和别的元素混进同一个投影类。")
    add("")
    add("  ★ 这个情况就是：股 4 只做『整体绕』时，")
    add("    投影把『细股 + 三股束』压成一个数，看起来和真两股一样。")
    add("    区别全在核里 —— 观测窗里拿不到。")
    add("")

    add("-" * 78)
    add("四、边界")
    add("-" * 78)
    add("  · 这里只做 Abel 化投影（σi↦σ）。更细的粗粒化（比如只看置换部分）")
    add("    会保留更多信息，但那就不是『两体』这个观测了。")
    add("  · 深度截到 %d，指数增长很快；分叉和压缩比在低层就稳定。" % max_depth)
    add("  · 结论是拓扑层面的『不可区分』，不是动力学上的『一定认错』。")
    add("")

    return "\n".join(out)


def make_plot(path, max_depth):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    dist, rows, all_proj, kernel = analyze(max_depth)

    fig = plt.figure(figsize=(15, 9.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.36, wspace=0.32)

    # ── 1) B4 分层增长
    ax = fig.add_subplot(gs[0, 0])
    ds = [r["depth"] for r in rows]
    ns = [r["nodes"] for r in rows]
    ax.plot(ds, ns, "-o", color="#2b6cb0", lw=2.2, ms=7, label="B4 元素数")
    ax.plot(ds, [2 * d + 1 for d in ds], "--s", color="#718096", lw=1.8, ms=6,
            label="B2 元素数 (2d+1)")
    ax.set_yscale("log")
    ax.set_xlabel("词长")
    ax.set_ylabel("元素数（对数）")
    ax.set_title("B4 指数增长 vs B2 线性增长\n差距就是被压缩掉的信息", fontsize=10)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.25)

    # ── 2) 原像大小
    ax = fig.add_subplot(gs[0, 1:])
    ks = sorted(all_proj)
    vals = [all_proj[k] for k in ks]
    bars = ax.bar(ks, vals, color="#e53e3e", width=0.75)
    ax.axvline(0, color="#2b6cb0", ls="--", lw=1.4)
    ax.set_xlabel("投影到 B2 的指数和（总共绕了几圈）")
    ax.set_ylabel("压着多少个 B4 元素")
    ax.set_title("粗粒化的原像：同一个『两体』观测底下压着多少四体构型", fontsize=10)
    ax.grid(alpha=0.25, axis="y")

    # ── 3) 含 σ3 / 不含 σ3
    ax = fig.add_subplot(gs[1, 0])
    w = 0.38
    ax.bar([d - w / 2 for d in ds], [r["no3"] for r in rows], width=w,
           color="#2b6cb0", label="不含 σ3（就是 B3 子群）")
    ax.bar([d + w / 2 for d in ds], [r["with3"] for r in rows], width=w,
           color="#e53e3e", label="含 σ3（股 4 参与）")
    ax.set_xlabel("词长")
    ax.set_ylabel("元素数")
    ax.set_title("「三根自己编」与「第四根参与」\n蓝柱 1/4/12/34/96 正是图 40 的 B3 序列", fontsize=10)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.25, axis="y")

    # ── 4) 投影分布（每层）
    ax = fig.add_subplot(gs[1, 1])
    for r in rows:
        if r["depth"] < 2:
            continue
        ks2 = sorted(r["proj"])
        ax.plot(ks2, [r["proj"][k] for k in ks2], "-o", ms=4,
                lw=1.4, label=f"深度 {r['depth']}")
    ax.set_xlabel("指数和")
    ax.set_ylabel("该层元素数")
    ax.set_title("每一层投影后都向中间聚拢\n两体观测的『分辨率』随深度下降", fontsize=10)
    ax.legend(fontsize=7.5, ncol=2)
    ax.grid(alpha=0.25)

    # ── 5) 结论框
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    ax.text(0.02, 0.95,
            "粗粒化账本\n"
            "────────────\n"
            f"B4 元素（深度 {max_depth}）：{len(dist)}\n"
            f"B2 像的大小：{len(all_proj)}\n"
            f"平均原像：{len(dist)/len(all_proj):.1f} 个 B4 / 1 个 B2\n"
            f"核（指数和 = 0）：{kernel} 个\n\n"
            "读法：\n"
            "· 投影 = 只剩『绕了几圈』\n"
            "· 核 = 三根内部编织，全部\n"
            "  落在两体观测的盲区\n"
            "· 所以『一根 + 三根』可以\n"
            "  伪装成『两根』",
            va="top", ha="left", fontsize=9.5, family="sans-serif",
            bbox=dict(boxstyle="round,pad=0.6", fc="#f7fafc", ec="#cbd5e0"))

    fig.suptitle("四股辫的粗粒化：三根合成一根粗绳，看起来就是双螺旋", fontsize=15)
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="B4 → B2 粗粒化")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--depth", type=int, default=4)
    args = ap.parse_args()

    text = build_report(args.depth)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text + "\n")
        print("report written: %s" % args.out)
    else:
        print(text)

    if args.plot:
        make_plot(args.plot, args.depth)


if __name__ == "__main__":
    main()
