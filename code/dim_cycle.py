# -*- coding: utf-8 -*-
"""
图 15：维度循环（按手绘图重排）
四维体（内含卷曲）→ 三维体 → 二维面 → 一维线 → 回到四维体
三次投影 + 一次反推，原点每层都在、位置不变。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Polygon
from matplotlib.animation import FuncAnimation, PillowWriter

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "figures"))

BG = "#05060a"
INK = "#e8eef5"
RED = "#e74c3c"
GRAY = "#8899aa"
GREEN = "#27ae60"


def cube_pts(cx, cy, w=1.7, h=1.7, d=0.8):
    x0, y0 = cx - w/2, cy - h/2
    P0 = (x0, y0); P1 = (x0+w, y0); P2 = (x0+w, y0+h); P3 = (x0, y0+h)
    ox, oy = d*0.78, d*0.55
    P4 = (P3[0]+ox, P3[1]+oy); P5 = (P2[0]+ox, P2[1]+oy)
    P6 = (P1[0]+ox, P1[1]+oy); P7 = (P0[0]+ox, P0[1]+oy)
    return dict(P0=P0, P1=P1, P2=P2, P3=P3, P4=P4, P5=P5, P6=P6, P7=P7)


def draw_cube(ax, cx, cy, textured=False, lw=2.0, alpha=1.0):
    p = cube_pts(cx, cy)
    P0,P1,P2,P3,P4,P5,P6,P7 = (p[k] for k in ["P0","P1","P2","P3","P4","P5","P6","P7"])
    ax.plot(*zip(P0,P1,P2,P3,P0), color=INK, lw=lw, alpha=alpha, solid_capstyle="round")
    ax.plot(*zip(P3,P4,P5,P2), color=INK, lw=lw, alpha=alpha)
    ax.plot(*zip(P1,P6,P5), color=INK, lw=lw, alpha=alpha)
    ax.plot(*zip(P0,P7), color=INK, lw=lw*0.55, ls=":", alpha=alpha*0.55)
    ax.plot(*zip(P7,P6), color=INK, lw=lw*0.55, ls=":", alpha=alpha*0.55)
    ax.plot(*zip(P7,P4), color=INK, lw=lw*0.55, ls=":", alpha=alpha*0.55)

    if textured:
        # 顶面：一排泡泡
        for u in np.linspace(0.18, 0.82, 3):
            for v in np.linspace(0.18, 0.82, 3):
                px = P3[0] + u*(P2[0]-P3[0]) + v*(P4[0]-P3[0])
                py = P3[1] + u*(P2[1]-P3[1]) + v*(P4[1]-P3[1])
                ax.add_patch(Circle((px, py), 0.085, fill=False, color=INK, lw=0.9, alpha=0.75))
        # 前面：斜纹
        for t in np.linspace(0.12, 0.88, 5):
            x1 = P0[0] + t*(P1[0]-P0[0]); y1 = P0[1] + t*(P1[1]-P0[1])
            ax.plot([x1, x1-0.28], [y1, y1+(P3[1]-P0[1])], color=INK,
                    lw=0.9, alpha=0.55)
    return p


def draw_plane(ax, cx, cy, w=1.05, h=2.3, shear=0.45, lw=2.0):
    """一张斜着的薄片"""
    P0 = (cx - w/2, cy - h/2)
    P1 = (cx + w/2, cy - h/2)
    P2 = (cx + w/2 + shear, cy + h/2)
    P3 = (cx - w/2 + shear, cy + h/2)
    ax.add_patch(Polygon([P0, P1, P2, P3], closed=True, facecolor=INK,
                         alpha=0.12, edgecolor=INK, lw=lw))
    # 薄边
    t = 0.10
    ax.plot([P1[0], P1[0]+t], [P1[1], P1[1]-t*0.35], color=INK, lw=lw*0.8)
    ax.plot([P2[0], P2[0]+t], [P2[1], P2[1]-t*0.35], color=INK, lw=lw*0.8)
    ax.plot([P1[0]+t, P2[0]+t], [P1[1]-t*0.35, P2[1]-t*0.35], color=INK, lw=lw*0.8)
    return dict(P0=P0, P1=P1, P2=P2, P3=P3)


fig, ax = plt.subplots(figsize=(12.5, 9), facecolor=BG)
ax.set_facecolor(BG)
ax.set_xlim(0, 16); ax.set_ylim(0, 11)
ax.set_aspect("equal"); ax.axis("off")

# 四个节点
c4 = draw_cube(ax, 3.4, 6.6, textured=True)     # 四维体
c3 = draw_cube(ax, 8.0, 8.6, textured=False)    # 三维体
pl = draw_plane(ax, 12.6, 5.6)                  # 二维面
ax.plot([8.0, 8.0], [1.8, 3.6], color=INK, lw=2.6, solid_capstyle="round")  # 一维线

# 原点红点（每层都在）
dots = [
    ((3.4+0.85+0.35), 6.6+0.15),
    ((8.0+0.85+0.35), 8.6+0.15),
    (12.6+0.15, 5.6+0.15),
    (8.0, 2.7),
]
for (dx, dy) in dots:
    ax.add_patch(Circle((dx, dy), 0.10, color=RED, zorder=6))

# 循环箭头
arrows = [
    ((4.7, 7.4), (6.7, 8.4), "投影", -0.28),
    ((9.4, 8.3), (11.7, 6.7), "投影", -0.28),
    ((11.9, 4.6), (9.1, 3.4), "投影", -0.28),
    ((7.0, 3.2), (4.4, 5.6), "反推", -0.28),
]
for (a, b, lab, rad) in arrows:
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=18,
                                 color=RED, lw=1.8,
                                 connectionstyle=f"arc3,rad={rad}"))
    mx, my = (a[0]+b[0])/2, (a[1]+b[1])/2
    ax.text(mx, my+0.35, lab, color=RED, fontsize=10.5, ha="center")

# 标签
ax.text(3.4, 4.9, "四维体\n（内含卷曲）", color=INK, ha="center", fontsize=11)
ax.text(8.0, 7.3, "三维体", color=INK, ha="center", fontsize=11)
ax.text(12.9, 3.95, "二维面", color=INK, ha="center", fontsize=11)
ax.text(8.0, 1.35, "一维线", color=INK, ha="center", fontsize=11)

# 环心
ax.text(8.0, 5.55, "闭环", color=GREEN, ha="center", fontsize=15)
ax.text(8.0, 5.05, "原点每层都在，位置不变", color=GRAY, ha="center", fontsize=10)

ax.set_title("维度循环：三次投影，一次反推，回到四维",
             color=INK, fontsize=16, pad=14)

p_out = os.path.join(OUT, "15_维度循环.png")
plt.savefig(p_out, dpi=170, bbox_inches="tight", facecolor=BG)
plt.close()
print("saved", p_out)
