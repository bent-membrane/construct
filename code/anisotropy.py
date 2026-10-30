# -*- coding: utf-8 -*-
"""
图 07：雾化投影的各向异性
喷流沿 z 轴喷出，从不同倾角看，投影晕会被压扁/拉长多少？
能不能和球对称的 CDM 晕区分开？
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(11)

# ---------------------------------------------------------------- 粒子（同沙盘 02）
N = 300_000
frac_core = 0.18
n_core = int(N * frac_core)
n_mist = N - n_core

tc = rng.uniform(0, np.deg2rad(4), n_core)
tm = rng.exponential(np.deg2rad(11), n_mist)
tm = tm[tm < np.deg2rad(70)]

vc = np.clip(rng.normal(5.0, 0.2, n_core), 2.0, None)
vm = np.clip(rng.normal(2.0, 0.8, len(tm)), 0.2, None)
tauc = rng.uniform(0, 1, n_core)
taum = rng.uniform(0, 1, len(tm))

pc = rng.uniform(0, 2 * np.pi, n_core)
pm = rng.uniform(0, 2 * np.pi, len(tm))


def xyz(theta, phi, v, tau):
    r = v * tau * np.sin(theta)
    return r * np.cos(phi), r * np.sin(phi), v * tau * np.cos(theta)


xc, yc, zc = xyz(tc, pc, vc, tauc)
xm, ym, zm = xyz(tm, pm, vm, taum)

x = np.concatenate([xc, xm])
y = np.concatenate([yc, ym])
z = np.concatenate([zc, zm])


def project(alpha):
    u = y
    v = x * np.cos(alpha) - z * np.sin(alpha)
    return u, v


def axis_ratio(u, v):
    cov = np.cov(np.vstack([u, v]))
    eig = np.linalg.eigvalsh(cov)
    return np.sqrt(max(eig[0], 1e-12) / eig[1]), eig


alphas = np.linspace(0, np.pi / 2, 61)
ratios = np.array([axis_ratio(*project(a))[0] for a in alphas])

# ---------------------------------------------------------------- 图
from matplotlib.colors import LogNorm

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5.2))

ax1.plot(np.rad2deg(alphas), ratios, color="#c0392b", lw=2.2, label="雾化投影")
ax1.plot(np.rad2deg(alphas), np.cos(alphas), color="#7f8c8d", ls=":", lw=1.5, label="薄饼：cos α")
ax1.axhline(1.0, color="#2980b9", ls="--", lw=1.5, label="球对称晕：1.00")
ax1.set_xlabel("观测方向与喷流轴的夹角 α（度）")
ax1.set_ylabel("投影轴比 b/a")
ax1.set_title("A. 倾斜看，晕会被压扁", fontsize=12)
ax1.legend(frameon=False, fontsize=9)
ax1.grid(alpha=0.25)
ax1.set_ylim(0, 1.15)

# B：α=30° 的投影密度
a = np.deg2rad(30)
u, v = project(a)
H, xe, ye = np.histogram2d(u, v, bins=180, range=[[-2.5, 2.5], [-2.5, 2.5]])
ax2.imshow(H.T, origin="lower", extent=[xe[0], xe[-1], ye[0], ye[-1]],
           cmap="inferno", norm=LogNorm(vmin=1, vmax=H.max()), aspect="equal")
ax2.set_title("B. α = 30° 时的投影", fontsize=12)
ax2.set_xlabel("投影 x")
ax2.set_ylabel("投影 y")

# C：等效椭圆
for adeg, c in [(0, "#c0392b"), (30, "#e67e22"), (60, "#2980b9")]:
    a = np.deg2rad(adeg)
    u, v = project(a)
    cov = np.cov(np.vstack([u, v]))
    eig, vec = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(vec[1, 1], vec[0, 1]))
    ratio = np.sqrt(eig[0] / eig[1])
    t = np.linspace(0, 2 * np.pi, 300)
    aa, bb = np.sqrt(eig[1]), np.sqrt(eig[0])
    px = aa * np.cos(t) * np.cos(np.radians(ang)) - bb * np.sin(t) * np.sin(np.radians(ang))
    py = aa * np.cos(t) * np.sin(np.radians(ang)) + bb * np.sin(t) * np.cos(np.radians(ang))
    ax3.plot(px, py, color=c, lw=2, label=f"α={adeg}°  b/a={ratio:.2f}")
ax3.set_aspect("equal")
ax3.set_xlim(-1.2, 1.2)
ax3.set_ylim(-1.2, 1.2)
ax3.set_title("C. 等效椭圆", fontsize=12)
ax3.legend(frameon=False, fontsize=9)
ax3.grid(alpha=0.25)

plt.tight_layout()
p = os.path.join(OUT, "07_anisotropy.png")
plt.savefig(p, dpi=170, bbox_inches="tight")
plt.close()
print("saved", p)

print("\n倾角 -> 轴比：")
for adeg in [0, 15, 30, 45, 60, 75, 90]:
    a = np.deg2rad(adeg)
    r, _ = axis_ratio(*project(a))
    print(f"  α={adeg:>2}°  b/a={r:.3f}")
