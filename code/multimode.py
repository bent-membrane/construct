# -*- coding: utf-8 -*-
"""
图 09：多模 —— 连续谱 + 离散模
模式竞争：耦合太强只剩一个（纯音，没信息）；太弱全是噪声。
中间地带，几个模式共存 = 既有结构又有变化。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(33)

N = 240
DT = 0.01
STEPS = 6000


def evolve(beta, seed=1):
    r = np.random.default_rng(seed)
    g = r.normal(1.5, 0.45, N)      # 各模式增益（固定结构里的"强度"）
    A = r.uniform(0, 0.01, N)
    hist = np.zeros((STEPS // 20, N))
    for i in range(STEPS):
        total = np.sum(A ** 2)
        dA = (g - 1.0 - A ** 2 - beta * (total - A ** 2)) * A
        A = np.clip(A + DT * dA, 0, None)
        if i % 20 == 0:
            hist[i // 20] = A
    return g, A, hist


betas = [0.0, 0.2, 0.5, 0.8, 1.0]
alive = []
for b in betas:
    g, A, _ = evolve(b)
    alive.append(np.sum(A > 0.05))

fig = plt.figure(figsize=(16.5, 5.2))

# A: 演化轨迹
ax1 = fig.add_subplot(1, 3, 1)
for b, c in [(0.0, "#27ae60"), (0.5, "#e67e22"), (1.0, "#c0392b")]:
    g, A, h = evolve(b)
    idx = np.argsort(g)[-8:]
    for j, k in enumerate(idx):
        ax1.plot(np.arange(h.shape[0]) * 20 * DT, h[:, k],
                 color=c, alpha=0.35 if j else 0.9, lw=1.2 if j else 2)
    ax1.plot([], [], color=c, lw=2, label=f"β = {b}")
ax1.set_xlabel("时间")
ax1.set_ylabel("模式振幅")
ax1.set_title("A. 竞争：β 越大，共存越少", fontsize=12)
ax1.legend(frameon=False, fontsize=9)
ax1.set_yscale("log")
ax1.grid(alpha=0.25)

# B: 存活模式数 vs 交叉饱和
ax2 = fig.add_subplot(1, 3, 2)
bb = np.linspace(0, 1, 21)
counts = []
for b in bb:
    g, A, _ = evolve(b, seed=2)
    counts.append(np.sum(A > 0.05))
ax2.plot(bb, counts, "-o", ms=4, color="#2980b9")
ax2.set_xlabel("交叉饱和 β（1 = 完全竞争）")
ax2.set_ylabel("存活模式数")
ax2.set_title("B. 从'全都要'到'只剩一个'", fontsize=12)
ax2.grid(alpha=0.25)

# C: 最终频谱
ax3 = fig.add_subplot(1, 3, 3)
for b, c in [(0.0, "#27ae60"), (0.5, "#e67e22"), (1.0, "#c0392b")]:
    g, A, _ = evolve(b)
    ax3.bar(np.arange(N), A, color=c, alpha=0.55, width=1.0,
            label=f"β = {b}（存活 {np.sum(A>0.05)}）")
ax3.set_xlabel("模式编号（相当于频率）")
ax3.set_ylabel("最终振幅")
ax3.set_title("C. 最终频谱：多峰 vs 单峰", fontsize=12)
ax3.legend(frameon=False, fontsize=9)
ax3.set_xlim(0, N)

plt.tight_layout()
p = os.path.join(OUT, "09_multimode.png")
plt.savefig(p, dpi=170, bbox_inches="tight")
plt.close()
print("saved", p)

print("\nβ -> 存活模式数：")
for b, n in zip(bb, counts):
    print(f"  β={b:.2f}  存活 {n}")
