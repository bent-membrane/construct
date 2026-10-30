# -*- coding: utf-8 -*-
"""
图 14：线 → 面 → 体 → 绳索（垂直关系），绳索上跑闭环
左：点、线、面、体，每一级沿"垂直于上一级"的方向展开；反推升维、投影降维，都有损
右：四维绳索（环面螺旋）上的信息回路，回到起点才算通
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.animation import FuncAnimation, PillowWriter

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "figures"))

C_PT   = "#ffffff"
C_LINE = "#e74c3c"
C_FACE = "#3fa9f5"
C_BODY = "#27ae60"
C_ROPE = "#ffd166"
C_OK   = "#27ae60"
BG     = "#05060a"
GRAY   = "#8899aa"

fig = plt.figure(figsize=(16, 7), facecolor=BG)

# ================= 左：线面体（2D 示意，垂直关系） =================
ax = fig.add_axes([0.02, 0.12, 0.64, 0.74])
ax.set_xlim(0, 11); ax.set_ylim(0, 5)
ax.axis("off")

# 点
ax.scatter([0.9], [2.6], s=150, c=C_PT, zorder=6)
ax.text(0.9, 1.75, "点\n0 维", color=C_PT, ha="center", fontsize=11.5)

# 线
ax.plot([2.8, 2.8], [2.0, 3.7], color=C_LINE, lw=3.4, solid_capstyle="round", zorder=5)
ax.scatter([2.8], [2.0], s=45, c=C_LINE, zorder=6)
ax.plot([2.8, 3.02], [2.0, 2.0], color=C_LINE, lw=1.3)
ax.plot([2.8, 2.8], [2.0, 2.22], color=C_LINE, lw=1.3)
ax.text(2.8, 1.75, "线\n1 维", color=C_LINE, ha="center", fontsize=11.5)

# 面
ax.add_patch(Rectangle((4.5, 2.0), 1.5, 1.6, facecolor=C_FACE, alpha=0.25,
                       edgecolor=C_FACE, lw=1.6))
ax.scatter([4.5], [2.0], s=45, c=C_FACE, zorder=6)
ax.plot([4.5, 4.72], [2.0, 2.0], color=C_FACE, lw=1.2)
ax.plot([4.5, 4.5], [2.0, 2.22], color=C_FACE, lw=1.2)
ax.text(5.25, 1.75, "面\n2 维", color=C_FACE, ha="center", fontsize=11.5)

# 体（等轴测立方体）
x0, y0, w, h, d = 7.3, 2.0, 1.25, 1.45, 0.55
ax.add_patch(Rectangle((x0, y0), w, h, facecolor=C_BODY, alpha=0.10,
                       edgecolor=C_BODY, lw=1.4))
ax.add_patch(Polygon([(x0, y0+h), (x0+w, y0+h), (x0+w+d, y0+h+d), (x0+d, y0+h+d)],
                     closed=True, facecolor=C_BODY, alpha=0.20, edgecolor=C_BODY, lw=1.4))
ax.add_patch(Polygon([(x0+w, y0), (x0+w, y0+h), (x0+w+d, y0+h+d), (x0+w+d, y0+d)],
                     closed=True, facecolor=C_BODY, alpha=0.07, edgecolor=C_BODY, lw=1.4))
ax.scatter([x0], [y0], s=45, c=C_BODY, zorder=6)
ax.plot([x0, x0+0.22], [y0, y0], color=C_BODY, lw=1.2)
ax.plot([x0, x0], [y0, y0+0.22], color=C_BODY, lw=1.2)
ax.text(x0+w/2, 1.75, "体\n3 维", color=C_BODY, ha="center", fontsize=11.5)

# 升维箭头（⊥）
for x1, x2, yy in [(1.20, 2.50, 4.35), (3.15, 4.25, 4.35), (6.20, 7.05, 4.35)]:
    ax.annotate("", xy=(x2, yy), xytext=(x1, yy),
                arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.6))
    ax.text((x1+x2)/2, yy+0.14, "⊥", color=GRAY, ha="center", fontsize=14)

# 降维箭头（投影）
for x1, x2, yy in [(2.50, 1.20, 1.15), (4.25, 3.15, 1.15), (7.05, 6.20, 1.15)]:
    ax.annotate("", xy=(x2, yy), xytext=(x1, yy),
                arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.3,
                                linestyle=(0, (4, 3))))
ax.text(0.9, 4.85, "→ 反推（有损）", color=GRAY, fontsize=12)
ax.text(4.6, 0.65, "← 投影（有损）", color="#c0392b", fontsize=12)

# ================= 右：绳索（3D） =================
ax3 = fig.add_axes([0.68, 0.14, 0.30, 0.70], projection="3d")
ax3.set_facecolor(BG)
ax3.set_box_aspect([1, 1, 0.78])
ax3.view_init(elev=16, azim=-70)
ax3.set_xlim(-1.25, 1.25); ax3.set_ylim(-1.25, 1.25); ax3.set_zlim(-0.15, 1.75)
ax3.set_axis_off()
ax3.set_title("四维绳索：绕四圈，回原点", color="white", fontsize=12, pad=-4)

RP, HGT, TURNS = 0.80, 1.55, 4


def rope_pt(s):
    """s ∈ [0,1]：沿曲面（圆柱）绕四圈上升到顶；s ∈ [1,1.35]：外侧回程"""
    if s <= 1.0:
        phi = TURNS * 2*np.pi * s
        z = HGT * s
        return np.array([RP*np.cos(phi), RP*np.sin(phi), z])
    t = (s - 1.0) / 0.35
    phi = 0.45 * np.pi * np.sin(np.pi * t)
    rad = RP + 0.12 * np.sin(np.pi * t)
    return np.array([rad*np.cos(phi), rad*np.sin(phi), HGT*(1-t)])


# 曲面：圆柱的母线与上下圈
for ang in np.linspace(0, 2*np.pi, 12, endpoint=False):
    ax3.plot([RP*np.cos(ang)]*2, [RP*np.sin(ang)]*2, [0, HGT],
             color="#6b7c93", lw=0.5, alpha=0.22)
zz = np.linspace(0, 2*np.pi, 120)
for z0 in (0.0, HGT):
    ax3.plot(RP*np.cos(zz), RP*np.sin(zz), np.full_like(zz, z0),
             color="#6b7c93", lw=0.8, alpha=0.35)

th = np.linspace(0, 1.0, 800)
allp = np.array([rope_pt(t) for t in th])
ax3.plot(allp[:, 0], allp[:, 1], allp[:, 2], color=C_ROPE, lw=1.2, alpha=0.28)
bp = np.linspace(1.0, 1.35, 100)
backp = np.array([rope_pt(t) for t in bp])
ax3.plot(backp[:, 0], backp[:, 1], backp[:, 2], color="#8899aa", lw=1.0,
         ls="--", alpha=0.45)

p0 = rope_pt(0.0)
ax3.scatter(*p0, s=90, c="white", marker="*", zorder=9)
ax3.text(p0[0]*1.25, p0[1]*1.25, p0[2]-0.15, "起点", color="white", fontsize=8.5)

trail, = ax3.plot([], [], [], color=C_ROPE, lw=3.0, zorder=10)
dot  = ax3.scatter([], [], [], s=90, c=C_ROPE, zorder=12)
halo = ax3.scatter([], [], [], s=480, c=C_ROPE, alpha=0.15, zorder=11)
ring = ax3.scatter([], [], [], s=650, c=C_OK, alpha=0.85, zorder=13)

fig.text(0.5, 0.035,
         "每一级都垂直于上一级；两条路都有损，能闭环就行",
         color="#9fb3c8", fontsize=12, ha="center")
fig.suptitle("线 → 面 → 体 → 绳索：垂直展开，绳索上跑回路",
             color="white", fontsize=16, y=0.955)

FRAMES = np.linspace(0, 1.35, 150)


def update(s_cur):
    ss = np.linspace(0, min(s_cur, 1.0), 160)
    pts = np.array([rope_pt(q) for q in ss])
    if s_cur > 1.0:
        bp = np.linspace(1.0, s_cur, 40)
        pts = np.vstack([pts, np.array([rope_pt(q) for q in bp])])
    trail.set_data(pts[:, 0], pts[:, 1])
    trail.set_3d_properties(pts[:, 2])

    p = rope_pt(s_cur)
    dot._offsets3d = ([p[0]], [p[1]], [p[2]])
    halo._offsets3d = ([p[0]], [p[1]], [p[2]])

    if s_cur >= 1.35 - 1e-6:
        ring._offsets3d = ([p[0]], [p[1]], [p[2]])
    else:
        ring._offsets3d = ([], [], [])
    return trail, dot, halo, ring


ani = FuncAnimation(fig, update, frames=FRAMES, interval=80, blit=False)
p_out = os.path.join(OUT, "14_rope_loop.gif")
ani.save(p_out, writer=PillowWriter(fps=16), dpi=100)
plt.close()
print("saved", p_out, f"{os.path.getsize(p_out)/1e6:.2f} MB")
