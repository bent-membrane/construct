# -*- coding: utf-8 -*-
"""
图 03：反着来 —— 从暗物质圈的轮廓反推喷射速度
正问题：给速度 → 算轮廓
逆问题：给轮廓 → 反推速度（能反推到哪一步？）
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = os.path.dirname(os.path.abspath(__file__))

THETA_CORE_MAX = np.deg2rad(4.0)
THETA_MIST_MEAN = np.deg2rad(11.0)
FRAC_CORE = 0.18
BINS = np.linspace(0, 3.0, 90)
CENTERS = 0.5 * (BINS[1:] + BINS[:-1])
AREA = np.pi * (BINS[1:] ** 2 - BINS[:-1] ** 2)


def simulate_R(v_core, v_mist, n=60_000, seed=1):
    """正模型：给定速度，返回投影半径样本（tau 归一化为 1）。"""
    rng = np.random.default_rng(seed)
    n_core = int(n * FRAC_CORE)
    n_mist = n - n_core

    tc = rng.uniform(0, THETA_CORE_MAX, n_core)
    tm = rng.exponential(THETA_MIST_MEAN, n_mist)
    tm = tm[tm < np.deg2rad(70)]

    vc = np.clip(rng.normal(v_core, 0.04 * v_core, n_core), 0.2 * v_core, None)
    vm = np.clip(rng.normal(v_mist, 0.40 * v_mist, len(tm)), 0.05 * v_mist, None)

    tauc = rng.uniform(0, 1, n_core)
    taum = rng.uniform(0, 1, len(tm))

    Rc = vc * tauc * np.sin(tc)
    Rm = vm * taum * np.sin(tm)
    return np.concatenate([Rc, Rm])


def profile(R):
    h, _ = np.histogram(R, bins=BINS)
    return h / AREA / len(R)


# ---------------------------------------------------------------- 观测
TRUE_VC, TRUE_VM = 5.0, 2.0
R_obs = simulate_R(TRUE_VC, TRUE_VM, n=400_000, seed=42)
obs = profile(R_obs)
rng_noise = np.random.default_rng(3)
obs_noisy = obs * (1 + 0.05 * rng_noise.standard_normal(len(obs)))
obs_noisy = np.clip(obs_noisy, 1e-12, None)

# ---------------------------------------------------------------- 反推
vcs = np.linspace(2.0, 8.0, 13)
vms = np.linspace(0.8, 3.2, 13)
resid = np.zeros((len(vcs), len(vms)))
best = (np.inf, None, None)

for i, vc in enumerate(vcs):
    for j, vm in enumerate(vms):
        Rm = simulate_R(vc, vm, n=40_000, seed=11)
        mod = profile(Rm)
        m = mod > 0
        err = np.mean((np.log10(mod[m]) - np.log10(obs_noisy[m])) ** 2)
        resid[i, j] = err
        if err < best[0]:
            best = (err, vc, vm)

err, vc_hat, vm_hat = best
print(f"真值       v_core={TRUE_VC:.2f}  v_mist={TRUE_VM:.2f}  比值={TRUE_VC/TRUE_VM:.2f}")
print(f"反推       v_core={vc_hat:.2f}  v_mist={vm_hat:.2f}  比值={vc_hat/vm_hat:.2f}")

# ---------------------------------------------------------------- 图
fig = plt.figure(figsize=(17, 5.4))

ax = fig.add_subplot(1, 3, 1)
ax.plot(CENTERS, obs_noisy, ".", ms=3, color="#7f8c8d", label="观测到的轮廓（加了 5% 噪声）")
ax.plot(CENTERS, profile(simulate_R(TRUE_VC, TRUE_VM, 200_000, 7)), color="#c0392b", lw=2, label="真相")
ax.plot(CENTERS, profile(simulate_R(vc_hat, vm_hat, 200_000, 7)), color="#2980b9", lw=2, ls="--", label="反推出来的")
ax.set_yscale("log")
ax.set_xlabel("投影半径")
ax.set_ylabel("投影面密度（相对）")
ax.set_title("A. 只给这条曲线，能把速度反推到哪一步？", fontsize=12)
ax.legend(frameon=False, fontsize=9)
ax.grid(alpha=0.2)

ax2 = fig.add_subplot(1, 3, 2)
cs = ax2.contourf(vcs, vms, resid.T, levels=25, cmap="viridis")
ax2.plot(TRUE_VC, TRUE_VM, "*", ms=18, color="red", label="真相")
ax2.plot(vc_hat, vm_hat, "o", ms=9, mfc="none", mec="white", mew=2, label="反推")
ax2.set_xlabel("内层速度（以飞行时间为单位）")
ax2.set_ylabel("外层速度")
ax2.set_title("B. 速度是能反推出来的\n（谷底就在真值附近）", fontsize=12)
ax2.legend(frameon=False, fontsize=9)
fig.colorbar(cs, ax=ax2, label="轮廓误差（越小越像）")

ax3 = fig.add_subplot(1, 3, 3)
taus = np.logspace(-1, 1, 200)
for v_hat, name, c in [(vc_hat, "内层", "#c0392b"), (vm_hat, "外层", "#2980b9")]:
    ax3.plot(taus, v_hat / taus, color=c, lw=2, label=f"{name}：速度 × 年龄 = 常数")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlabel("喷流年龄（以沙盘的“1”为单位）")
ax3.set_ylabel("喷射速度")
ax3.set_title("C. 观测只给一条双曲线\n要落成一个点，得有一个时钟", fontsize=12)
ax3.legend(frameon=False, fontsize=9)
ax3.grid(alpha=0.2, which="both")
for tau in [0.5, 1.0, 2.0]:
    ax3.axvline(tau, color="gray", ls=":", lw=0.8)

plt.tight_layout()
p = os.path.join(OUT, "03_inverse_speed.png")
plt.savefig(p, dpi=170, bbox_inches="tight")
plt.close()
print("saved", p)

# 换算表
print("\n如果知道年龄，速度就是：")
for tau in [0.5, 1.0, 2.0, 5.0]:
    print(f"  年龄={tau:>4} → 内层 {vc_hat/tau:6.2f}   外层 {vm_hat/tau:6.2f}")
