# -*- coding: utf-8 -*-
"""
图 12：通讯过程——卷曲穿层与闭环
左：螺距对上 → 信号一路穿到第 4 层 → 回传 → 闭环闭合
右：差半个螺距 → 信号卡在第 3 层 → 闭环失败

每一层画成一条卷曲（正弦侧视）。信号沿层前进，到层间尝试穿透：
对齐（螺距对上）才通，错开（差半个螺距）就卡住——全或无。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "figures"))

LAYERS = 4
DY = 2.8
LAM = 12.0
AMP = 0.34
X_MAX = 100.0
CHECKPOINTS = [(30.0, 1), (55.0, 2), (80.0, 3)]   # (位置, 从第几层往上穿)

C_LAYER = ["#3fa9f5", "#39d0c8", "#b07cff", "#ff9f43"]
C_SIG   = "#ffd166"
C_OK    = "#27ae60"
C_BAD   = "#c0392b"
C_DIM   = "#666666"


def y_of(k):
    return (k - 1) * DY


def state(p, phase_fail):
    """p: 0→1 去程, 1→1.6 回程, 1.6→2 停留。返回 (x, y, 层, 状态)"""
    if phase_fail and p >= 0.80:
        return 80.0 + 1.5 * np.sin(26.0 * p), y_of(3), 3, "stuck"
    if p <= 1.0:
        x = p * X_MAX
        layer = 1
        for cx, lf in CHECKPOINTS:
            if x >= cx:
                layer = lf + 1
        return x, y_of(layer), layer, "go"
    if p <= 1.6:
        s = (p - 1.0) / 0.6
        x = X_MAX * (1 - s)
        y = y_of(4) + (y_of(1) - y_of(4)) * s + 4.2 * np.sin(np.pi * s)
        return x, y, 4, "back"
    return 0.0, y_of(1), 1, "done"


def build(ax, phase_fail, title):
    ax.set_facecolor("#05060a")
    ax.set_xlim(-9, 109)
    ax.set_ylim(-1.8, 11.8)
    ax.axis("off")
    ax.set_title(title, color="white", fontsize=13.5, pad=8)

    xs = np.linspace(0, X_MAX, 900)
    for k in range(1, LAYERS + 1):
        ph = np.pi if (phase_fail and k == 4) else 0.0
        ax.plot(xs, y_of(k) + AMP * np.sin(2*np.pi*xs/LAM + ph),
                color=C_LAYER[k-1], lw=1.7, alpha=0.85, zorder=2)
        ax.text(-4.5, y_of(k), f"第{k}层", color=C_LAYER[k-1], fontsize=10,
                va="center", ha="right")

    # 穿层位置
    for cx, lf in CHECKPOINTS:
        ok = not (lf == 3 and phase_fail)
        ax.plot([cx, cx], [y_of(lf)+0.5, y_of(lf+1)-0.5],
                color=C_OK if ok else C_BAD, lw=1.1, ls=":", alpha=0.6, zorder=1)
        ax.text(cx, y_of(lf)+0.5+0.12, "通" if ok else "断",
                color=C_OK if ok else C_BAD, fontsize=9.5, ha="center", va="bottom")

    trail, = ax.plot([], [], color=C_SIG, lw=2.0, alpha=0.75, zorder=4)
    halo = ax.scatter([], [], s=340, c=C_SIG, alpha=0.15, linewidths=0, zorder=5)
    dot  = ax.scatter([], [], s=75,  c=C_SIG, alpha=1.0,  linewidths=0, zorder=6)
    txt  = ax.text(0.02, 0.97, "", transform=ax.transAxes, color="white",
                   fontsize=12, va="top")

    ring = ax.scatter([0.92], [0.80], s=1700, transform=ax.transAxes,
                      c="#242424", edgecolors="#444444", linewidths=1.2, zorder=3)
    ring_txt = ax.text(0.92, 0.80, "", transform=ax.transAxes,
                       ha="center", va="center", color="#777777",
                       fontsize=10.5, zorder=4)
    return trail, halo, dot, txt, ring, ring_txt


def make_update(objs, phase_fail):
    trail, halo, dot, txt, ring, ring_txt = objs

    def update(p):
        x, y, layer, status = state(p, phase_fail)

        ps = np.linspace(0, p, 70)
        xs, ys = [], []
        for q in ps:
            qx, qy, _, _ = state(q, phase_fail)
            xs.append(qx); ys.append(qy)
        trail.set_data(xs, ys)
        halo.set_offsets([[x, y]])
        dot.set_offsets([[x, y]])

        if status == "go":
            txt.set_text(f"去程 · 第 {layer} 层")
            txt.set_color("white")
        elif status == "stuck":
            txt.set_text("卡住：差半个螺距 → 无（不通）")
            txt.set_color(C_BAD)
            ring.set_facecolors(C_BAD)
            ring.set_edgecolors(C_BAD)
            ring_txt.set_text("断")
            ring_txt.set_color("white")
        elif status == "back":
            txt.set_text("回程 · 出门和进门，同一个口")
            txt.set_color("#9fd3ff")
            ring.set_facecolors(C_OK)
            ring.set_edgecolors(C_OK)
            ring_txt.set_text("闭环")
            ring_txt.set_color("white")
        elif status == "done":
            txt.set_text("闭环建立 · 硬开关合上，软信息允许有损")
            txt.set_color(C_OK)
            ring.set_facecolors(C_OK)
            ring.set_edgecolors(C_OK)
            ring_txt.set_text("闭环")
            ring_txt.set_color("white")
        return trail, halo, dot, txt, ring, ring_txt

    return update


fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.3), facecolor="#05060a")
fig.subplots_adjust(left=0.04, right=0.98, top=0.855, bottom=0.04, wspace=0.10)
fig.suptitle("通讯过程：卷曲 → 穿层（全或无）→ 闭环", color="white", fontsize=16, y=0.975)

objs_ok  = build(axes[0], False, "螺距对上：全（通）")
objs_bad = build(axes[1], True,  "差半个螺距：无（不通）")

FRAMES = np.linspace(0, 2.0, 120)


def update_all(p):
    objs = []
    objs += list(make_update(objs_ok,  False)(p))
    objs += list(make_update(objs_bad, True)(p))
    return objs


ani = FuncAnimation(fig, update_all, frames=FRAMES, interval=80, blit=False)
p_out = os.path.join(OUT, "12_communication_process.gif")
ani.save(p_out, writer=PillowWriter(fps=14), dpi=105)
plt.close()
print("saved", p_out, f"{os.path.getsize(p_out)/1e6:.2f} MB")

# ---------------------------------------------------------------- 关键时刻数值
print("\n关键时刻：")
for q in [0.0, 0.30, 0.55, 0.80, 0.95, 1.0, 1.3, 1.6, 2.0]:
    a = state(q, False)
    b = state(q, True)
    print(f"  p={q:4.2f}  通: x={a[0]:6.1f} 层{a[2]} {a[3]:6s} | "
          f"断: x={b[0]:6.1f} 层{b[2]} {b[3]}")
