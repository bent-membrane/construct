# -*- coding: utf-8 -*-
"""
图 13：三维卷曲闭环（环面版 · (1,4) 螺旋）
信号沿环面走一圈大圆（θ: 0→2π），同时绕小圆四圈（φ = 4θ）。
四圈 = 四次卷曲；走完回到起点 = 第五条回到第一条。

左：四圈螺距均匀 → 终点与起点重合（闭环）
右：第四圈只转半圈 → 终点错开（闭不上）
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "figures"))

R, r = 3.0, 1.0
COLORS = ["#3fa9f5", "#39d0c8", "#b07cff", "#ff9f43"]
C_SIG  = "#ffd166"
C_OK   = "#27ae60"
C_BAD  = "#c0392b"
C_GHOST = "#6b7c93"
BG     = "#05060a"

TWO_PI = 2.0 * np.pi


def P(th, fail=False):
    """环面上的 (1,4) 螺旋。th ∈ [0, 2π]。
    fail=True 时第四段（最后 90°）小圆只转半圈，终点错开。"""
    if fail and th > 1.5 * np.pi:
        phi = 6.0 * np.pi + (th - 1.5 * np.pi)
    else:
        phi = 4.0 * th
    return np.array([
        (R + r * np.cos(phi)) * np.cos(th),
        (R + r * np.cos(phi)) * np.sin(th),
        r * np.sin(phi),
    ])


def setup_ax(ax, fail, title):
    ax.set_facecolor(BG)
    ax.set_box_aspect([1, 1, 0.52])
    ax.view_init(elev=30, azim=-62)
    ax.set_xlim(-4.4, 4.4)
    ax.set_ylim(-4.4, 4.4)
    ax.set_zlim(-1.6, 1.6)
    ax.set_axis_off()
    ax.set_title(title, color="white", fontsize=13.5, pad=-2)

    # 环面参考网（经线：固定 φ；纬线：固定 θ）
    th = np.linspace(0, TWO_PI, 240)
    for phi0 in np.linspace(0, TWO_PI, 12, endpoint=False):
        ax.plot((R + r*np.cos(phi0))*np.cos(th),
                (R + r*np.cos(phi0))*np.sin(th),
                np.full_like(th, r*np.sin(phi0)),
                color=C_GHOST, lw=0.6, alpha=0.22, zorder=1)
    ph = np.linspace(0, TWO_PI, 100)
    for t0 in np.linspace(0, TWO_PI, 16, endpoint=False):
        ax.plot((R + r*np.cos(ph))*np.cos(t0),
                (R + r*np.cos(ph))*np.sin(t0),
                r*np.sin(ph),
                color=C_GHOST, lw=0.6, alpha=0.14, zorder=1)

    # 完整路径（淡）
    ss = np.linspace(0, TWO_PI, 1200)
    pts = np.array([P(q, fail) for q in ss])
    for k in range(4):
        lo, hi = k*np.pi/2, (k+1)*np.pi/2
        m = (ss >= lo) & (ss <= hi)
        ax.plot(pts[m, 0], pts[m, 1], pts[m, 2], color=COLORS[k],
                lw=1.1, alpha=0.20, zorder=2)

    p0 = P(0.0, False)
    ax.scatter(*p0, s=120, c="white", marker="*", zorder=9)
    ax.text(p0[0]*1.14, p0[1]*1.14, p0[2]+0.32, "起点", color="white", fontsize=9)

    trails = []
    for k in range(4):
        ln, = ax.plot([], [], [], color=COLORS[k], lw=2.6, alpha=1.0, zorder=5)
        trails.append(ln)

    dot = ax.scatter([], [], [], s=95, c=C_SIG, zorder=10)
    halo = ax.scatter([], [], [], s=520, c=C_SIG, alpha=0.15, zorder=8)
    ring = ax.scatter([], [], [], s=760, c=C_OK, alpha=0.9, zorder=11)
    txt = ax.text2D(0.03, 0.95, "", transform=ax.transAxes, color="white", fontsize=11.5)

    return trails, dot, halo, ring, txt


def update_factory(objs, fail):
    trails, dot, halo, ring, txt = objs

    def update(th_cur):
        for k in range(4):
            lo = k*np.pi/2
            hi = min((k+1)*np.pi/2, th_cur)
            if hi <= lo + 1e-9:
                trails[k].set_data([], [])
                trails[k].set_3d_properties([])
            else:
                ss = np.linspace(lo, hi, 70)
                pts = np.array([P(q, fail) for q in ss])
                trails[k].set_data(pts[:, 0], pts[:, 1])
                trails[k].set_3d_properties(pts[:, 2])

        p = P(th_cur, fail)
        dot._offsets3d = ([p[0]], [p[1]], [p[2]])
        halo._offsets3d = ([p[0]], [p[1]], [p[2]])

        if th_cur >= TWO_PI:
            if fail:
                txt.set_text("错开 · 第四圈只转了半圈，闭不上")
                txt.set_color(C_BAD)
                ring.set_facecolors(C_BAD)
            else:
                txt.set_text("闭环 · 第五条 = 第一条，回到起点")
                txt.set_color(C_OK)
                ring.set_facecolors(C_OK)
            ring._offsets3d = ([p[0]], [p[1]], [p[2]])
        else:
            seg = min(int(th_cur / (np.pi/2)), 3)
            txt.set_text(f"第 {seg+1} 圈 · 第 {seg+1} 层")
            txt.set_color("white")
            ring._offsets3d = ([], [], [])
        return trails + [dot, halo, ring, txt]

    return update


fig = plt.figure(figsize=(13.5, 7.4), facecolor=BG)
ax1 = fig.add_subplot(1, 2, 1, projection="3d")
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
fig.subplots_adjust(left=0.0, right=1.0, top=0.88, bottom=0.0, wspace=0.02)
fig.suptitle("三维卷曲：绕小圆四圈，回到起点", color="white", fontsize=16, y=0.965)

o1 = setup_ax(ax1, False, "四圈螺距均匀：全（通）")
o2 = setup_ax(ax2, True,  "第四圈差半圈：无（不通）")

FRAMES = np.linspace(0, TWO_PI, 150)


def update_all(th_cur):
    return update_factory(o1, False)(th_cur) + update_factory(o2, True)(th_cur)


ani = FuncAnimation(fig, update_all, frames=FRAMES, interval=80, blit=False)
p_out = os.path.join(OUT, "13_3d_curl_loop.gif")
ani.save(p_out, writer=PillowWriter(fps=16), dpi=100)
plt.close()
print("saved", p_out, f"{os.path.getsize(p_out)/1e6:.2f} MB")
