# -*- coding: utf-8 -*-
"""
图 11：通讯是怎么走的
A 投影链：球 → 圆 → 简谐，原点跟着走
B 卷曲与平行：各卷各的，螺距对上才通
C 闭环：硬开关 + 软信息
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "figures"))

BLUE   = "#2980b9"
ORANGE = "#e67e22"
RED    = "#c0392b"
GREEN  = "#27ae60"
GRAY   = "#7f8c8d"
INK    = "#2c3e50"

# 三个面板统一坐标框，保证标题对齐
XL = (-3.6, 3.6)
YL = (-2.0, 2.0)

fig = plt.figure(figsize=(16.5, 5.8))
fig.suptitle("信息怎么走：投影 · 卷曲 · 闭环", fontsize=17, color=INK, y=0.98)

# ============================ A ============================
ax = fig.add_subplot(1, 3, 1)
ax.set_xlim(*XL); ax.set_ylim(*YL)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("A. 投影链：球 → 圆 → 简谐", fontsize=13, color=INK, pad=10)

# 球（左）
cx = -2.4
ax.add_patch(Circle((cx, 0), 0.85, fill=False, lw=2.2, color=BLUE))
ax.add_patch(Ellipse((cx, 0), 1.70, 0.44, fill=False, lw=1.0, color=BLUE, alpha=0.40))
ax.add_patch(Ellipse((cx, 0), 0.44, 1.70, fill=False, lw=1.0, color=BLUE, alpha=0.40))
t = np.linspace(0, 2*np.pi, 200)
ax.plot(cx + 0.85*np.cos(t), 0.85*np.sin(t), color=ORANGE, lw=1.1, ls="--", alpha=0.75)
ax.plot(cx + 0.85, 0, "o", color=ORANGE, ms=7, zorder=5)
ax.annotate("", xy=(cx+0.15, 0.86), xytext=(cx+0.15, 1.12),
            arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.3))
ax.text(cx+0.24, 1.14, "自转", fontsize=9, color=BLUE)
ax.text(cx, -1.14, "三维：球", ha="center", fontsize=10.5, color=BLUE)

# 圆（中）
cx2 = 0.0
ax.add_patch(Circle((cx2, 0), 0.85, fill=False, lw=2.2, color=BLUE))
ax.plot(cx2 + 0.85*np.cos(t), 0.85*np.sin(t), color=ORANGE, lw=1.1, ls="--", alpha=0.75)
ax.plot(cx2 - 0.85, 0, "o", color=ORANGE, ms=7, zorder=5)
ax.text(cx2, -1.14, "二维：圆", ha="center", fontsize=10.5, color=BLUE)

# 一维：线段 + 简谐波（右）
cx3 = 2.4
ax.plot([cx3-1.05, cx3+1.05], [0.55, 0.55], color=BLUE, lw=2.4, solid_capstyle="round")
ax.plot(cx3-0.62, 0.55, "o", color=ORANGE, ms=7, zorder=5)
ax.annotate("", xy=(cx3-0.30, 0.55), xytext=(cx3-0.52, 0.55),
            arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.3))
xs = np.linspace(cx3-1.05, cx3+1.05, 300)
ax.plot(xs, -0.30 + 0.24*np.sin(2*np.pi*(xs-cx3)/0.70), color=ORANGE, lw=1.6)
ax.text(cx3, -1.14, "一维：简谐", ha="center", fontsize=10.5, color=BLUE)

# 原点：三层都在，红点连虚线
ax.plot(cx, 0, "o", color=RED, ms=10, zorder=6)
ax.plot(cx2, 0, "o", color=RED, ms=10, zorder=6)
ax.plot(cx3, 0.55, "o", color=RED, ms=10, zorder=6)
ax.plot([cx, cx2, cx3], [0, 0, 0.55], ls=":", color=RED, lw=1.3, alpha=0.75)
ax.text(cx2, 1.55, "原点跟着投影走（三层都在）", ha="center", fontsize=10.5, color=RED)

# 投影箭头
for x0, x1 in [(cx+0.95, cx2-0.95), (cx2+0.95, cx3-1.14)]:
    ax.annotate("", xy=(x1, 0), xytext=(x0, 0),
                arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.6))
ax.text((cx+cx2)/2, 0.16, "投影", ha="center", fontsize=9, color=GRAY)
ax.text((cx2+cx3)/2, 0.16, "投影", ha="center", fontsize=9, color=GRAY)
ax.text(-0.6, -1.80, "丢一维，不丢相对位置", ha="center", fontsize=10.5, color=INK)

# ============================ B ============================
ax = fig.add_subplot(1, 3, 2)
ax.set_xlim(*XL); ax.set_ylim(*YL)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("B. 卷曲与平行：各卷各的", fontsize=13, color=INK, pad=10)

x = np.linspace(-2.6, 2.6, 800)
T = 1.5
# 上半：螺距对上（同相）
ax.plot(x, 0.95 + 0.24*np.sin(2*np.pi*x/T), color=BLUE, lw=2.2)
ax.plot(x, 0.42 + 0.24*np.sin(2*np.pi*x/T), color=BLUE, lw=2.2, alpha=0.50)
ax.text(0.0, 1.52, "螺距对上：全（通）", ha="center", fontsize=11, color=GREEN)

# 分隔线
ax.plot([-2.6, 2.6], [0.02, 0.02], color=GRAY, lw=0.8, ls=":")

# 下半：差半个螺距（反相）
ax.plot(x, -0.62 + 0.24*np.sin(2*np.pi*x/T), color=BLUE, lw=2.2)
ax.plot(x, -1.15 + 0.24*np.sin(2*np.pi*x/T + np.pi), color=RED, lw=2.2, alpha=0.85)
ax.text(0.0, -1.62, "差半个螺距：无（不通）", ha="center", fontsize=11, color=RED)

# 穿层箭头
ax.annotate("", xy=(-2.78, -1.15), xytext=(-2.78, 0.95),
            arrowprops=dict(arrowstyle="<|-|>", color=GRAY, lw=1.3))
ax.text(-2.66, 0.10, "穿层", fontsize=9, color=GRAY, rotation=90)
ax.text(-2.60, 1.55, "一束=一维", fontsize=9, color=GRAY, ha="left")

# ============================ C ============================
ax = fig.add_subplot(1, 3, 3)
ax.set_xlim(*XL); ax.set_ylim(*YL)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("C. 闭环：硬开关 + 软信息", fontsize=13, color=INK, pad=10)

R, r = 1.95, 1.10
ax.add_patch(Ellipse((0, 0), 2*R, 2*r, fill=False, lw=2.6, color=BLUE))

# 去 / 回
ax.annotate("", xy=(0.85, r-0.05), xytext=(0.15, r-0.05),
            arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.9))
ax.text(0.95, r+0.06, "去", ha="left", fontsize=10, color=BLUE)
ax.annotate("", xy=(-0.85, -r+0.05), xytext=(-0.15, -r+0.05),
            arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.9))
ax.text(-0.95, -r-0.12, "回", ha="right", fontsize=10, color=BLUE)

# 硬开关：左顶点画一个闸刀（断开态）
sx = -R
ax.plot([sx], [0.22], "o", color=GREEN, ms=8, zorder=6)
ax.plot([sx], [-0.22], "o", color=GREEN, ms=8, zorder=6)
ax.plot([sx, sx+0.42], [-0.22, 0.34], color=GREEN, lw=3.0, zorder=6)
ax.text(sx-0.08, -0.66, "全或无\n（断=不通）", ha="center", va="top", fontsize=9.5, color=GREEN, linespacing=1.5)

# 层内软信息：噪声波
seg = np.linspace(0.0, 1.0, 160)
ax.plot(-0.95 + 1.90*seg, -0.55 + 0.08*np.sin(34*np.pi*seg), color=ORANGE, lw=1.5)
ax.text(0.0, -0.86, "软信息（允许有损）", ha="center", fontsize=9.5, color=ORANGE)

# 环心
ax.text(0, 0.26, "闭环才是活的", ha="center", fontsize=13.5, color=INK)
ax.text(0, -0.12, "出门和进门，同一个口", ha="center", fontsize=10, color=GRAY)

plt.tight_layout(rect=[0, 0, 1, 0.945])
p = os.path.join(OUT, "11_communication_schematic.png")
plt.savefig(p, dpi=170, bbox_inches="tight")
plt.close()
print("saved", p)
