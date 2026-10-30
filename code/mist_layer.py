# -*- coding: utf-8 -*-
"""
图 02：雾化层 —— "挤出去的吧，和雾化似的"
内层准直射流 + 外层雾化，看投影到观测面上是什么轮廓。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(7)

# ---------------------------------------------------------------- 粒子
N = 300_000
frac_core = 0.18
n_core = int(N * frac_core)
n_mist = N - n_core

theta_core = rng.uniform(0, np.deg2rad(4), n_core)
theta_mist = rng.exponential(np.deg2rad(11), n_mist)
keep = theta_mist < np.deg2rad(70)
theta_mist = theta_mist[keep]
n_mist = len(theta_mist)

# 内层快、外层慢——"挤出去"的速度差
v_core = np.clip(rng.normal(5.0, 0.2, n_core), 2.0, None)
v_mist = np.clip(rng.normal(2.0, 0.8, n_mist), 0.2, None)

phi_core = rng.uniform(0, 2 * np.pi, n_core)
phi_mist = rng.uniform(0, 2 * np.pi, n_mist)


# 持续喷流：每个粒子已飞了不同时长 tau，形成一整根而不是"一团"
tau_core = rng.uniform(0, 1, n_core)
tau_mist = rng.uniform(0, 1, n_mist)


def to_xyz(theta, phi, v, tau):
    r = v * tau * np.sin(theta)
    return r * np.cos(phi), r * np.sin(phi), v * tau * np.cos(theta)


xc, yc, zc = to_xyz(theta_core, phi_core, v_core, tau_core)
xm, ym, zm = to_xyz(theta_mist, phi_mist, v_mist, tau_mist)

x = np.concatenate([xc, xm])
y = np.concatenate([yc, ym])
z = np.concatenate([zc, zm])
is_mist = np.concatenate([np.zeros(n_core, bool), np.ones(n_mist, bool)])


def project(x, y, z, alpha):
    """观测方向与喷流轴夹角 alpha；返回天球平面坐标。"""
    u = y
    v = x * np.cos(alpha) - z * np.sin(alpha)
    return u, v


# ---------------------------------------------------------------- 图
fig = plt.figure(figsize=(17, 5.4))

# --- A：侧视，内层 / 外层分色
ax = fig.add_subplot(1, 3, 1)
ax.scatter(xc, zc, s=0.4, c="#c0392b", alpha=0.25, linewidths=0, label="内层：准直射流")
ax.scatter(xm, zm, s=0.4, c="#2980b9", alpha=0.08, linewidths=0, label="外层：雾化")
ax.set_xlabel("面内 x")
ax.set_ylabel("离面高度 z")
ax.set_title("A. 侧看：一根芯，外面一团雾", fontsize=12)
ax.legend(frameon=False, markerscale=12, fontsize=9, loc="upper right")
ax.grid(alpha=0.2)
ax.set_xlim(-3, 3)
ax.set_ylim(-1.0, 6.2)

# --- B：投影面密度（沿轴俯视）
ax2 = fig.add_subplot(1, 3, 2)
u, v = project(x, y, z, 0.0)
H, xe, ye = np.histogram2d(u, v, bins=260, range=[[-3, 3], [-3, 3]])
ax2.imshow(
    H.T, origin="lower", extent=[xe[0], xe[-1], ye[0], ye[-1]],
    cmap="inferno", norm=LogNorm(vmin=1, vmax=H.max()), aspect="equal",
)
ax2.set_title("B. 投影到膜上：中心一个点，外面一圈雾", fontsize=12)
ax2.set_xlabel("投影 x")
ax2.set_ylabel("投影 y")

# --- C：径向轮廓，不同倾角
ax3 = fig.add_subplot(1, 3, 3)
bins = np.linspace(0, 3, 120)
centers = 0.5 * (bins[1:] + bins[:-1])
for alpha_deg, color in [(0, "#c0392b"), (30, "#e67e22"), (60, "#2980b9")]:
    uu, vv = project(x, y, z, np.deg2rad(alpha_deg))
    R = np.hypot(uu, vv)
    cnt, _ = np.histogram(R, bins=bins)
    area = np.pi * (bins[1:] ** 2 - bins[:-1] ** 2)
    ax3.plot(centers, cnt / area, color=color, lw=2, label=f"倾角 {alpha_deg}°")

# 只芯 / 只雾
for mask, color, ls, lab in [
    (~is_mist, "#8e44ad", "--", "只有内层芯"),
    (is_mist, "#16a085", ":", "只有外层雾"),
]:
    uu, vv = project(x[mask], y[mask], z[mask], 0.0)
    R = np.hypot(uu, vv)
    cnt, _ = np.histogram(R, bins=bins)
    area = np.pi * (bins[1:] ** 2 - bins[:-1] ** 2)
    ax3.plot(centers, cnt / area, color=color, ls=ls, lw=1.6, label=lab)

ax3.set_yscale("log")
ax3.set_xlabel("投影半径")
ax3.set_ylabel("投影面密度（相对）")
ax3.set_title("C. 径向轮廓：芯给尖峰，雾给外沿", fontsize=12)
ax3.legend(frameon=False, fontsize=8.5)
ax3.grid(alpha=0.2)

plt.tight_layout()
p = os.path.join(OUT, "02_mist_layer.png")
plt.savefig(p, dpi=170, bbox_inches="tight")
plt.close()
print("saved", p)

# ---------------------------------------------------------------- 数
for label, mask in [("全部", np.ones_like(is_mist)), ("只有内层", ~is_mist), ("只有外层", is_mist)]:
    uu, vv = project(x[mask], y[mask], z[mask], 0.0)
    R = np.hypot(uu, vv)
    print(f"{label:6s}  中位半径={np.median(R):.3f}  90%分位={np.quantile(R,0.9):.3f}  最远={R.max():.2f}")

u0, v0 = project(x, y, z, 0.0)
R0 = np.hypot(u0, v0)
for r in [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]:
    print(f"  R={r:.2f}  面密度 ~ {np.sum((R0>=r*0.9)&(R0<r*1.1))/(np.pi*((r*1.1)**2-(r*0.9)**2)):.1f}")
