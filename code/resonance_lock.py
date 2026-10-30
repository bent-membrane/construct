# -*- coding: utf-8 -*-
"""
图 08：共振对齐
膜在震动/旋转 -> 每个模式的固有频率有一个宽度 sigma。
耦合强度 K 要多大，才能把这一群振子锁成一个"闭环"？
（Kuramoto 模型）
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(21)

N = 600
DT = 0.02
T_TOTAL = 120.0
STEPS = int(T_TOTAL / DT)


def run(K, sigma, omega0=0.0, seed=1):
    r = np.random.default_rng(seed)
    omega = r.normal(omega0, sigma, N)
    theta = r.uniform(0, 2 * np.pi, N)
    rs = np.empty(STEPS)
    for i in range(STEPS):
        # 平均场形式，O(N)
        z = np.mean(np.exp(1j * theta))
        theta += DT * (omega + K * np.abs(z) * np.sin(np.angle(z) - theta))
        rs[i] = np.abs(np.mean(np.exp(1j * theta)))
    return rs


sigmas = [0.3, 1.0, 2.0]
Ks = np.linspace(0, 6, 31)

fig = plt.figure(figsize=(16.5, 5.2))

# A: 序参量 vs 耦合强度
ax1 = fig.add_subplot(1, 3, 1)
for sig, c in zip(sigmas, ["#2980b9", "#e67e22", "#c0392b"]):
    rs_final = []
    for K in Ks:
        rr = run(K, sig, seed=3)
        rs_final.append(np.mean(rr[-int(10 / DT):]))
    ax1.plot(Ks, rs_final, "-o", ms=3, color=c, label=f"膜抖动 σ = {sig}")
ax1.set_xlabel("耦合强度 K")
ax1.set_ylabel("序参量 r（0=散，1=锁死）")
ax1.set_title("A. 抖动越大，要锁住需要的耦合越强", fontsize=12)
ax1.legend(frameon=False, fontsize=9)
ax1.grid(alpha=0.25)

# B: 锁定过程的时间曲线
ax2 = fig.add_subplot(1, 3, 2)
for K, c in zip([0.8, 2.0, 4.0], ["#c0392b", "#e67e22", "#27ae60"]):
    rr = run(K, 1.0, seed=5)
    ax2.plot(np.arange(STEPS) * DT, rr, color=c, lw=1.8,
             label=f"K = {K}")
ax2.set_xlabel("时间")
ax2.set_ylabel("序参量 r")
ax2.set_title("B. K 太小锁不住，K 过临界才起飞", fontsize=12)
ax2.legend(frameon=False, fontsize=9)
ax2.grid(alpha=0.25)

# C: 相位分布（极坐标）
ax3 = fig.add_subplot(1, 3, 3, projection="polar")
for K, c, lab in [(0.5, "#c0392b", "K=0.5 散"), (4.0, "#27ae60", "K=4.0 锁")]:
    rr = run(K, 1.0, seed=7)
    r = np.random.default_rng(7)
    omega = r.normal(0, 1.0, N)
    theta = r.uniform(0, 2 * np.pi, N)
    for i in range(STEPS):
        z = np.mean(np.exp(1j * theta))
        theta += DT * (omega + K * np.abs(z) * np.sin(np.angle(z) - theta))
    ax3.hist(theta, bins=36, color=c, alpha=0.45, label=lab)
ax3.set_title("C. 相位分布：散 vs 聚成一束", fontsize=12, pad=16)
ax3.legend(frameon=False, fontsize=8, loc="lower left")

plt.tight_layout()
p = os.path.join(OUT, "08_resonance_lock.png")
plt.savefig(p, dpi=170, bbox_inches="tight")
plt.close()
print("saved", p)

# ---------------------------------------------------------------- 数
print("\n序参量（终态平均）：")
for sig in sigmas:
    vals = []
    for K in [0, 0.5, 1, 1.5, 2, 3, 4]:
        rr = run(K, sig, seed=3)
        vals.append(np.mean(rr[-int(10 / DT):]))
    print(f"  σ={sig}: " + "  ".join(f"K={k}:{v:.2f}" for k, v in zip([0,0.5,1,1.5,2,3,4], vals)))

# 临界耦合（近似：r 首次超过 0.5）
for sig in sigmas:
    for K in Ks:
        rr = run(K, sig, seed=3)
        if np.mean(rr[-int(10 / DT):]) > 0.5:
            print(f"  σ={sig}: r>0.5 的临界 K ≈ {K:.1f}")
            break
