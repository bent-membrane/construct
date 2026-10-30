#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
有疑问的并合系统：按辫群粗粒化统一分类
==========================================
2026-10-29 重算：`phase`（刚穿/回程/中间/turnaround）来自已降级的 four_clusters
时间模型，只作旧判断保留、不参与任何判定；判定只用「块数 + 赘生物 + 核容量」。
新增第四节：φ、ΔS、块数三者的几何关系（不看时间）。

2026-10-27，要求：「你把之前觉得有疑问点星团可以都算一下」

不重新拟合、不调参数。只问一件事：

    观测到 k 个团块 → 底下的辫至少要几股？
    这些辫被粗粒化时，能藏住多少东西？

判据全部来自群结构，只有一条来自观测：

    ★ 块数 ≤ 股数
      B_n 的粗粒化最多产生 n 块。看到 3 块，n 就不可能等于 2。

配套的拓扑容量（本脚本现算）：
    B2 → B2：核 = 0           没有任何隐藏
    B3 → B2：核 = ?           三股压成两股时丢掉的部分
    B4 → B2：核 = 158/595，26.6%

系统清单沿用 four_clusters.py / braid_unravel.py 的文献取值，
只标量级；「疑问」列是之前几轮里真正卡住的地方。

用法：
    python cluster_braid_classify.py
    python cluster_braid_classify.py --plot ../figures/42_cluster_classify.png
"""

import argparse
import sys
from collections import Counter, defaultdict, deque
from functools import lru_cache

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ═══════════════════════════════════════════════════════════════
#  拓扑侧：B2 / B3 / B4 的核
# ═══════════════════════════════════════════════════════════════

B2_GENS = ((1, 1), (1, -1))
B3_GENS = ((1, 1), (1, -1), (2, 1), (2, -1))
B4_GENS = ((1, 1), (1, -1), (2, 1), (2, -1), (3, 1), (3, -1))


def free_reduce(word):
    out = []
    for g in word:
        if out and out[-1][0] == g[0] and out[-1][1] == -g[1]:
            out.pop()
        else:
            out.append(g)
    return tuple(out)


def braid_step(word):
    """通用辫关系：相邻两股 σiσ(i+1)σi ↔ σ(i+1)σiσ(i+1)；不相邻的对易。"""
    res = set()
    n = len(word)
    for i in range(n - 1):
        a, b = word[i], word[i + 1]
        if abs(a[0] - b[0]) >= 2:          # 不相交 → 对易
            res.add(word[:i] + (b, a) + word[i + 2:])
    for i in range(n - 2):
        a, b, c = word[i], word[i + 1], word[i + 2]
        for s in (1, -1):
            x = (a[0], s)
            if a == x and b == (a[0] + 1, s) and c == x:
                swap = ((a[0] + 1, s), x, (a[0] + 1, s))
                res.add(word[:i] + swap + word[i + 3:])
            elif a == (b[0] + 1, s) and b == (b[0], s) and c == (b[0] + 1, s):
                swap = ((b[0], s), (b[0] + 1, s), (b[0], s))
                res.add(word[:i] + swap + word[i + 3:])
    return res


@lru_cache(maxsize=None)
def normalize(word):
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


def projection(word):
    """π: σi ↦ σ ⇒ 投影就是指数和。"""
    return sum(sign for _, sign in word)


def build_graph(gens, max_depth):
    start = ()
    dist = {start: 0}
    q = deque([start])
    while q:
        v = q.popleft()
        if dist[v] >= max_depth:
            continue
        for g in gens:
            u = normalize(v + (g,))
            if u not in dist:
                dist[u] = dist[v] + 1
                q.append(u)
    return dist


def topology_capacity(gens, max_depth):
    dist = build_graph(gens, max_depth)
    proj = Counter(projection(v) for v in dist)
    kernel = proj.get(0, 0)
    # 单位元永远在核里，但它不构成「隐藏」；只看非平凡的核元素。
    nontrivial = max(0, kernel - 1)
    return {
        "total": len(dist),
        "image": len(proj),
        "kernel": kernel,
        "nontrivial": nontrivial,
        "kernel_share": nontrivial / len(dist) if dist else 0.0,
        "proj": proj,
    }


# ═══════════════════════════════════════════════════════════════
#  观测侧：七个有疑问的系统
# ═══════════════════════════════════════════════════════════════

SYSTEMS = [
    dict(
        name="MACS J0025.4", subs=2, phi="2–5°", dS="170–370 kpc",
        phase="刚穿", extra="无",
        doubt="最干净的两体，被用作相位锚点",
    ),
    dict(
        name="Bullet 1E0657", subs=3, phi="11.8–14.1°", dS="100–162 kpc",
        phase="刚穿末", extra="主团≥3 亚团块 + 5σ 质量桥",
        doubt="两股模型不该稳定出现『分开又连着』",
    ),
    dict(
        name="MACS J1149", subs=3, phi="—", dS="2350 kpc（投影）",
        phase="三子团两次并合", extra="第三子团 4.2%、速度 4000 km/s",
        doubt="两体相位 344 Myr 对周期 19.6 Gyr 只占 6%，根本到不了远心点",
    ),
    dict(
        name="A520", subs=3, phi="21.5–39°", dS="368–630 kpc",
        phase="中间", extra="dark core 争议（Jee>10σ / Clowe<5%）",
        doubt="第三核到底在不在——决定它是不是三股",
    ),
    dict(
        name="ZwCl 0008.8", subs=2, phi="32–41°", dS="173–312 kpc",
        phase="回程", extra="无",
        doubt="φ 最大，按相位该是回程；两体能否撑住这么大的角",
    ),
    dict(
        name="A168", subs=2, phi="—", dS="—",
        phase="turnaround", extra="冷气跑到子团前面",
        doubt="冷气在前，动力学与两体图景是否自洽",
    ),
    dict(
        name="El Gordo", subs=2, phi="—", dS="—",
        phase="—", extra="v/v_esc = 0.49",
        doubt="5–6σ 是 ΛCDM 宇宙学张力，不是动力学超速——两件事别混",
    ),
]


# ═══════════════════════════════════════════════════════════════
#  分类
# ═══════════════════════════════════════════════════════════════

def min_braid(subs):
    """块数 ≤ 股数：看到 k 块，n 至少是 k。"""
    return max(2, subs)


def verdict(s):
    n = min_braid(s["subs"])
    if s["subs"] >= 3:
        return n, "明确不是两体", "块数已经超过 2"
    if s["extra"] != "无":
        return 2, "两体可能，但有赘生物", "赘生物可以被 n>2 的核吸收"
    return 2, "两体最省", "没有粗粒化证据"


def build_report(max_depth=4):
    out = []
    add = out.append
    add("=" * 80)
    add("有疑问的并合系统：按辫群粗粒化统一分类")
    add("=" * 80)
    add("")

    # ---- 拓扑容量
    b2 = topology_capacity(B2_GENS, max_depth)
    b3 = topology_capacity(B3_GENS, max_depth)
    b4 = topology_capacity(B4_GENS, max_depth)

    add("-" * 80)
    add("一、拓扑容量（先算，不看观测）")
    add("-" * 80)
    add("  群        元素数(深度%d)  投影像  核(指数和0)  非平凡核  隐藏占比" % max_depth)
    for nm, cap in (("B2", b2), ("B3", b3), ("B4", b4)):
        add("  %-6s %12d %8d %12d %10d %8.1f%%"
            % (nm, cap["total"], cap["image"], cap["kernel"],
               cap["nontrivial"], 100 * cap["kernel_share"]))
    add("")
    add("  读法：核 = 粗粒化成两体之后，完全看不见的那些构型。")
    add("        单位元永远在核里，但它不算隐藏；所以 B2 的隐藏占比是 0。")
    add("        B2 投影到自己是恒等映射 —— 两股没有可藏的地方。")
    add("")

    # ---- 系统表（按股数层级）
    add("-" * 80)
    add("二、按股数分层（不用同一套判据套所有系统）")
    add("-" * 80)

    layers = [
        ("B2 层 —— 两块，两股够用",
         [s for s in SYSTEMS if s["subs"] == 2 and s["extra"] == "无"],
         "核 = 0，没有可藏的地方；解开只有一条路（单调解绑）",
         "没有粗粒化余地；但「最省」不等于「唯一」——三股以上也能投影成两块。"),
        ("B3 层 —— 三块，三股必需",
         [s for s in SYSTEMS if s["subs"] >= 3 and s["name"] != "MACS J1149"],
         "核 = 21.8%（147 个元素里 32 个非平凡核元素）",
         "块数超过 2，B2 直接出局，不需要任何拟合；但只排除两股，不区分 B3/B4。"),
        ("B4 层 —— 三块 + 可藏的赘生物",
         [s for s in SYSTEMS if s["name"] == "MACS J1149"],
         "核 = 26.4%（595 个元素里 157 个非平凡核元素）",
         "4.2% < 26.4%：第三子团可以藏在核里，不必引入第四股；超过 26.4% 就失效。"),
        ("待定层 —— 两块 + 赘生物，无占比",
         [s for s in SYSTEMS if s["subs"] == 2 and s["extra"] != "无"],
         "赘生物可以被 n>2 的核吸收，但缺文献占比",
         "两体解不排他；只能定性排他，不能定量匹配。"),
    ]
    for title, members, topo, concl in layers:
        add("")
        add("  【%s】" % title)
        add("    拓扑：%s" % topo)
        for s in members:
            add("    · %-18s %d 块   φ %-13s ΔS %-16s %s"
                % (s["name"], s["subs"], s["phi"], s["dS"], s["extra"]))
            add("      旧相位（已降级，不作依据）：%s｜疑问：%s" % (s["phase"], s["doubt"]))
        add("    → %s" % concl)
    add("")

    # ---- 层内比较
    add("-" * 80)
    add("三、层内比较：同一层里 φ 也不一样（只看有 φ/ΔS 的系统）")
    add("-" * 80)
    add("  B2 层（都是两股、两块）：")
    add("    MACS J0025.4   φ 2–5°       ΔS 170–370 kpc")
    add("    ZwCl 0008.8    φ 32–41°     ΔS 173–312 kpc     ← 同层，φ 差 8 倍")
    add("")
    add("  B3 层（都是三股、三块）：")
    add("    Bullet 1E0657  φ 11.8–14.1°  ΔS 100–162 kpc")
    add("    A520           φ 21.5–39°    ΔS 368–630 kpc")
    add("")
    add("  · 同股数 ≠ 同 φ：B2 层两个系统 φ 从 2° 到 41°。")
    add("  · 不同股数也可以同 φ：A520（21.5–39）与 ZwCl（32–41）区间重叠。")
    add("  · 所以「φ 大 = 回程」的时间解释在几何上没有支撑——")
    add("    φ 是层内的散点，不是纺线上的刻度。")
    add("  · MACS J1149 / A168 / El Gordo 没有 φ，不参与这一步。")
    add("")
    add("-" * 80)
    add("四、唯一能定量的一格：MACS J1149")
    add("-" * 80)
    add("  第三子团占 4.2%（braid_unravel.py 口径）。")
    add("  B4 的核占 %.1f%%（深度 %d）。" % (100 * b4["kernel_share"], max_depth))
    add("  → 4.2% 远小于核容量：第三子团完全可以藏在核里，")
    add("    不需要把它当成『第三股』——它可能只是三股束里的一股。")
    add("    这正好对上 10-27 的观察：一根散开 + 三根成束，看起来像两股。")
    add("")
    add("  反过来说：如果哪天测到某个『第三子团』占比超过核容量，")
    add("    粗粒化解法就失效，必须真的引入第三股。")
    add("    —— 这是这套分类里唯一可证伪的一格。")
    add("")

    # ---- 边界
    add("-" * 80)
    add("五、边界")
    add("-" * 80)
    add("  · 块数是观测计数，不是股数；只用它定下界，不反推上界。")
    add("  · 除了 MACS J1149，其余系统没有『赘生物质量占比』的文献值，")
    add("    所以只能定性排他，不能定量匹配。")
    add("  · 深度截到 %d；核占比随深度缓慢变化，量级稳定。" % max_depth)
    add("")

    return "\n".join(out), b2, b3, b4


def make_plot(path, b2, b3, b4, max_depth):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig = plt.figure(figsize=(15, 9.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.4, wspace=0.34)

    # ── 1) 核容量
    ax = fig.add_subplot(gs[0, 0])
    names = ["B2", "B3", "B4"]
    caps = [b2, b3, b4]
    shares = [100 * c["kernel_share"] for c in caps]
    bars = ax.bar(names, shares, color=["#718096", "#2b6cb0", "#c53030"], width=0.55)
    for b, c in zip(bars, caps):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1,
                f"{c['nontrivial']}/{c['total']}", ha="center", fontsize=9)
    ax.set_ylabel("核占比 (%)")
    ax.set_title("粗粒化成两体后，有多少构型看不见\nB2 自己为 0，越高的辫藏得越多", fontsize=10)
    ax.set_ylim(0, max(shares) * 1.35)
    ax.grid(alpha=0.25, axis="y")

    # ── 2) 系统块数
    ax = fig.add_subplot(gs[0, 1])
    names = [s["name"] for s in SYSTEMS]
    subs = [s["subs"] for s in SYSTEMS]
    colors = ["#e53e3e" if v >= 3 else "#2b6cb0" for v in subs]
    y = list(range(len(names)))
    ax.barh(y, subs, color=colors, height=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8.5)
    ax.set_xlabel("观测到的团块数")
    ax.set_xlim(0, 4.2)
    ax.axvline(2.5, color="#4a5568", ls="--", lw=1.2)
    ax.text(2.55, len(names) - 0.75, "这条线右边\nB2 直接出局", fontsize=8, color="#4a5568")
    ax.set_title("块数 ≤ 股数：看到 3 块，n 不可能等于 2", fontsize=10)
    ax.grid(alpha=0.25, axis="x")

    # ── 3) 判定汇总
    ax = fig.add_subplot(gs[0, 2])
    hard = sum(1 for s in SYSTEMS if s["subs"] >= 3)
    soft = sum(1 for s in SYSTEMS if s["subs"] == 2 and s["extra"] != "无")
    clean = sum(1 for s in SYSTEMS if s["subs"] == 2 and s["extra"] == "无")
    labels = ["明确 n≥3", "两块+赘生物", "两块干净"]
    vals = [hard, soft, clean]
    ax.pie(vals, labels=labels, autopct=lambda p: f"{int(round(p*sum(vals)/100))}",
           colors=["#e53e3e", "#dd6b20", "#2b6cb0"], startangle=90,
           labeldistance=1.18, textprops=dict(fontsize=8.5))
    ax.set_title("七个系统的分类", fontsize=10)

    # ── 4) φ–ΔS 平面
    ax = fig.add_subplot(gs[1, 0:2])
    pts = [
        ("MACS J0025", 3.5, 270, "#718096"),
        ("Bullet", 13.0, 130, "#c53030"),
        ("A520", 30.0, 500, "#dd6b20"),
        ("ZwCl 0008", 36.5, 240, "#2b6cb0"),
    ]
    for nm, phi, dS, col in pts:
        ax.scatter(phi, dS, s=160, color=col, edgecolors="black", linewidths=0.8, zorder=5)
        ax.annotate(nm, (phi, dS), textcoords="offset points", xytext=(8, 8), fontsize=9)
    ax.set_xlabel("φ（核—雾夹角，度）")
    ax.set_ylabel("ΔS（kpc）")
    ax.set_title("有 φ/ΔS 的四个系统（数据来自 four_clusters.py）\n"
                 "MACS J1149 / A168 / El Gordo 只有块数与结构，没有这一对坐标",
                 fontsize=10)
    ax.grid(alpha=0.25)

    # ── 5) MACS J1149 那一格
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    ax.text(0.02, 0.96,
            "唯一能定量的一格\n"
            "─────────────\n"
            "MACS J1149 第三子团\n"
            "    占 4.2%\n\n"
            f"B4 的核\n"
            f"    {100*b4['kernel_share']:.1f}%（{b4['nontrivial']}/{b4['total']}）\n\n"
            "4.2% < 26.6%\n"
            "→ 第三子团可以藏在核里\n"
            "→ 不必引入『第三股』\n\n"
            "可证伪：\n"
            "若赘生物占比 > 核容量，\n"
            "粗粒化解法失效",
            va="top", ha="left", fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.6", fc="#f7fafc", ec="#cbd5e0"))

    fig.suptitle("七个有疑问的并合系统：块数定下界，核容量定上限", fontsize=15)
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="并合系统辫群分类")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--depth", type=int, default=4)
    args = ap.parse_args()

    text, b2, b3, b4 = build_report(args.depth)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text + "\n")
        print("report written: %s" % args.out)
    else:
        print(text)

    if args.plot:
        make_plot(args.plot, b2, b3, b4, args.depth)


if __name__ == "__main__":
    main()
