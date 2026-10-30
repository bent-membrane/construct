# -*- coding: utf-8 -*-
"""
图 01：膜凹陷 → 喷流 → 投影晕
把"没有暗物质，是挤出去的，和雾化似的"变成一张能看的几何沙盘。
输出：三张图，不写公式。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import splu

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- 1. 膜
def build_operator(N=257, L=10.0):
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]
    X, Y = np.meshgrid(x, x)
    R = np.hypot(X, Y)

    main = -2.0 * np.ones(N)   # 一维二阶差分对角；二维由 kron 相加得 -4
    off = np.ones(N - 1)
    D2 = sparse.diags([off, main, off], [-1, 0, 1], shape=(N, N), format="lil") / dx**2
    I = sparse.eye(N)
    Lap = (sparse.kron(D2, I) + sparse.kron(I, D2)).tolil()

    idx = np.arange(N * N).reshape(N, N)
    bd = np.unique(np.concatenate([idx[0, :], idx[-1, :], idx[:, 0], idx[:, -1]]))
    Lap[bd, :] = 0
    Lap[bd, bd] = 1

    return x, X, Y, R, Lap.tocsr(), bd


def membrane_shape(lu, bd, R, M=1.0, T=1.0, sigma=0.8, N=257):
    rho = M / (2 * np.pi * sigma**2) * np.exp(-(R**2) / (2 * sigma**2))
    f = (rho / T).ravel().copy()
    f[bd] = 0
    return lu.solve(f).reshape(N, N)


x, X, Y, R, Lap, bd = build_operator()
lu = splu(Lap)

h = membrane_shape(lu, bd, R)


# ---------------------------------------------------------------- 2. 画图
fig = plt.figure(figsize=(17, 5.6))

# --- 图 A：凹陷本体（3D）
ax = fig.add_subplot(1, 3, 1, projection="3d")
sk = slice(None, None, 6)
ax.plot_surface(
    X[sk, sk], Y[sk, sk], h[sk, sk],
    cmap="viridis", linewidth=0, antialiased=True, rstride=1, cstride=1,
)
ax.set_title("A. 质量压出的凹陷\n（引力 = 膜的垂直化）", fontsize=12, pad=8)
ax.set_zlim(h.min() * 1.05, 0.3)
ax.view_init(elev=32, azim=-58)
ax.set_xlabel("面内 x")
ax.set_ylabel("面内 y")


# --- 图 B：不同张力下的剖面（可拨的旋钮）
ax2 = fig.add_subplot(1, 3, 2)
colors = ["#c0392b", "#e67e22", "#2980b9", "#27ae60"]
for T, c in zip([0.5, 1.0, 2.0, 4.0], colors):
    hh = membrane_shape(lu, bd, R, T=T)
    ax2.plot(x, hh[257 // 2, :], color=c, lw=2, label=f"张力 T = {T}")
ax2.axhline(0, color="gray", lw=0.8, ls=":")
ax2.set_title("B. 同一个质量，不同的膜张力\n张力越大，坑越浅越窄", fontsize=12)
ax2.set_xlabel("离中心的距离")
ax2.set_ylabel("下沉深度")
ax2.legend(frameon=False, fontsize=9)
ax2.grid(alpha=0.25)

# --- 图 C：喷流锥 + 投影（"漏斗状的当然看着圈大"）
ax3 = fig.add_subplot(1, 3, 3, projection="3d")

theta = np.deg2rad(18)      # 喷流半张角
H = 6.0                     # 喷流高度
phi = np.linspace(0, 2 * np.pi, 60)
z = np.linspace(0, H, 30)
PHI, Z = np.meshgrid(phi, z)
Rj = Z * np.tan(theta)
Xj = Rj * np.cos(PHI)
Yj = Rj * np.sin(PHI)

ax3.plot_surface(Xj, Yj, Z, color="#7fb3d5", alpha=0.55, linewidth=0)
ax3.plot(np.cos(phi) * H * np.tan(theta), np.sin(phi) * H * np.tan(theta),
         np.full_like(phi, H), color="#2980b9", lw=1.5)

# 投影到膜平面（z=0）
ax3.plot(np.cos(phi) * H * np.tan(theta), np.sin(phi) * H * np.tan(theta),
         np.zeros_like(phi), color="#c0392b", lw=2, ls="--")

# 膜面
sk2 = slice(None, None, 10)
ax3.plot_surface(X[sk2, sk2], Y[sk2, sk2], np.zeros_like(X[sk2, sk2]),
                 color="gray", alpha=0.12, linewidth=0)

ax3.set_title("C. 喷流锥与它在膜上的投影\n外圈是投影，不是物质真铺那么广", fontsize=12, pad=8)
ax3.set_xlim(-6, 6)
ax3.set_ylim(-6, 6)
ax3.set_zlim(-0.2, H * 1.1)
ax3.view_init(elev=22, azim=-60)
ax3.set_xlabel("面内 x")
ax3.set_ylabel("面内 y")

plt.tight_layout()
p = os.path.join(OUT, "01_membrane_jet.png")
plt.savefig(p, dpi=170, bbox_inches="tight")
plt.close()
print("saved", p)


# ---------------------------------------------------------------- 3. 数值小结（不进图）
h_all = {T: membrane_shape(lu, bd, R, T=T) for T in [0.5, 1.0, 2.0, 4.0]}
for T, hh in h_all.items():
    depth = -hh.min()
    # 45° 等坡度点（视界候选）
    prof = hh[257 // 2, 128:]
    rr = x[128:]
    slope = np.abs(np.gradient(prof, rr))
    idx = np.where(slope > 1.0)[0]
    r_h = rr[idx[-1]] if len(idx) else float("nan")
    print(f"T={T:>4}  最深={depth:.3f}   |dh/dr|>1 的最外半径={r_h:.2f}")
