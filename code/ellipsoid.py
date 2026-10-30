# -*- coding: utf-8 -*-
"""
图 10：三维椭球密度场 + 倾斜投影
把"雾"从粒子云换成连续密度：轴对称长椭球，轴向沿喷流。
任意倾角下投影的轴比与面密度轮廓都有解析解。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- 椭球参数
A = 1.0        # 横向半轴
B = 1.0
C = 3.0        # 轴向半轴（沿喷流拉长）
Q = C / A      # 长宽比


def proj_axes(alpha):
    """视线与喷流轴成 alpha；返回投影椭圆两个半轴。"""
    # 垂直于 n 和 z 的方向：半轴 A
    # n-z 平面内的垂直方向：半轴 sqrt(A^2 cos^2 + C^2 sin^2)
    s1 = A
    s2 = np.sqrt(A**2 * np.cos(alpha)**2 + C**2 * np.sin(alpha)**2)
    return s1, s2


def axis_ratio(alpha):
    s1, s2 = proj_axes(alpha)
    return min(s1, s2) / max(s1, s2)


# ---------------------------------------------------------------- 图
fig = plt.figure(figsize=(17, 5.4))

# A: 三维椭球等密度面
ax = fig.add_subplot(1, 3, 1, projection="3d")
u = np.linspace(0, 2 * np.pi, 60)
v = np.linspace(0, np.pi, 40)
x = A * np.outer(np.cos(u), np.sin(v))
y = B * np.outer(np.sin(u), np.sin(v))
z = C * np.outer(np.ones_like(u), np.cos(v))
ax.plot_surface(x, y, z, color="#3fa9f5", alpha=0.35, linewidth=0)
ax.plot([0, 0], [0, 0], [-C * 1.2, C * 1.2], color="#c0392b", lw=2, ls="--")
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_zlim(-3, 3)
ax.set_xlabel("面内 x")
ax.set_ylabel("面内 y")
ax.set_zlabel("喷流 z")
ax.set_title("A. 雾 = 沿喷流拉长的椭球\n（虚线 = 喷流轴）", fontsize=12)
ax.view_init(elev=18, azim=-60)

# B: 解析轴比 vs 倾角
alphas = np.linspace(0, np.pi / 2, 61)
ratios = np.array([axis_ratio(a) for a in alphas])

ax2 = fig.add_subplot(1, 3, 2)
ax2.plot(np.rad2deg(alphas), ratios, color="#c0392b", lw=2.2,
         label=f"椭球 c/a = {Q}")
# 对比：厚饼近似 cos(alpha)
ax2.plot(np.rad2deg(alphas), np.abs(np.cos(alphas)), color="#7f8c8d",
         ls=":", lw=1.5, label="薄饼 cos α")
ax2.axhline(1.0, color="#2980b9", ls="--", lw=1.2, label="球对称")
ax2.set_xlabel("视线与喷流轴夹角 α（度）")
ax2.set_ylabel("投影轴比 b/a")
ax2.set_title("B. 椭球模型的轴比曲线\n比粒子模型平缓得多", fontsize=12)
ax2.legend(frameon=False, fontsize=9)
ax2.grid(alpha=0.25)
ax2.set_ylim(0, 1.1)

# C: 不同倾角的投影面密度（解析高斯椭球）
# 三维高斯：exp(-(x^2/A^2 + y^2/B^2 + z^2/C^2)/2)
# 投影后仍是二维高斯，半轴由 proj_axes 给出
ax3 = fig.add_subplot(1, 3, 3)
tt = np.linspace(0, 2 * np.pi, 400)
for adeg, color, ls in [(0, "#c0392b", "-"), (30, "#e67e22", "--"), (60, "#2980b9", ":")]:
    a = np.deg2rad(adeg)
    s1, s2 = proj_axes(a)
    ax3.plot(s2 * np.cos(tt), s1 * np.sin(tt), color=color, ls=ls, lw=2,
             label=f"α={adeg}°  b/a={axis_ratio(a):.2f}")
ax3.set_aspect("equal")
ax3.set_xlim(-3.4, 3.4)
ax3.set_ylim(-3.4, 3.4)
ax3.axhline(0, color="gray", lw=0.6, alpha=0.5)
ax3.axvline(0, color="gray", lw=0.6, alpha=0.5)
ax3.set_xlabel("投影平面的长轴方向")
ax3.set_ylabel("投影平面的短轴方向")
ax3.set_title("C. 同一个椭球，不同角度的投影轮廓", fontsize=12)
ax3.legend(frameon=False, fontsize=9)
ax3.grid(alpha=0.25)

plt.tight_layout()
p = os.path.join(OUT, "10_ellipsoid.png")
plt.savefig(p, dpi=170, bbox_inches="tight")
plt.close()
print("saved", p)

# ---------------------------------------------------------------- 数
print(f"\n椭球 c/a = {Q}")
print("倾角 -> 轴比：")
for adeg in [0, 15, 30, 45, 60, 75, 90]:
    print(f"  α={adeg:>2}°  b/a={axis_ratio(np.deg2rad(adeg)):.3f}")

print("\n和粒子模型对比（粒子模型数值来自沙盘 07）：")
particle = {0: 0.997, 15: 0.600, 30: 0.362, 45: 0.265, 60: 0.219, 75: 0.197, 90: 0.190}
for adeg in [0, 15, 30, 45, 60, 75, 90]:
    e = axis_ratio(np.deg2rad(adeg))
    pv = particle[adeg]
    print(f"  α={adeg:>2}°  椭球={e:.3f}  粒子={pv:.3f}  差={abs(e-pv):.3f}")

# 反推：要匹配观测常见的 b/a=0.7，需要多大倾角
target = 0.7
best = min(alphas, key=lambda a: abs(axis_ratio(a) - target))
print(f"\n椭球模型：要得到 b/a={target}，需要倾角 α ≈ {np.rad2deg(best):.1f}°")
