"""
四维投影 → 反推折叠模数 k
对应：2026-10-23_构造.md 第 9 节「刻度：模」

设定
----
四维点同时绕两个正交平面旋转（U154 的"公转加自转"）：
    (x, y) 平面转 θ = k·φ
    (z, w) 平面转 φ
转数比 k 就是折叠模数。星星散布在管半径 TUBE 内。

投影：丢掉 w，得到三维点云。再叠加观测噪声（按结构尺度 R1 的百分比）。

反推
----
给定候选 k，用 z 估计相位 φ，xy 的角度应当是 k·φ。
两支相位（φ 和 2π−φ）都试，取贴合的那支。
所有点的相位残差聚合起来，对齐度接近 1 说明 k 对。

对照
----
1. 打乱：只重排 z，(x,y) 与 z 的对应断掉
2. 聚类：手写 k-means + 轮廓系数，看它能不能自己找到 k

运行：python fourd_k_inverse.py
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
rng = np.random.default_rng(20261023)

# ---------------- 参数 ----------------
N = 1_000_000
K_TRUE = 5
R1 = R2 = 1.0
TUBE = 0.06
KS = np.arange(2, 31)
NOISE_LEVELS = (0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20)


def make(k, noise_frac, n=N):
    """四维双旋转螺旋 → 投影到三维 → 加观测噪声"""
    t = rng.uniform(0, 1, n)
    phi = 2 * np.pi * t
    theta = k * phi
    x = R1 * np.cos(theta) + rng.normal(0, TUBE, n)
    y = R1 * np.sin(theta) + rng.normal(0, TUBE, n)
    z = R2 * np.cos(phi) + rng.normal(0, TUBE, n)
    _w = R2 * np.sin(phi) + rng.normal(0, TUBE, n)      # 投影丢掉
    s = noise_frac * R1
    x += rng.normal(0, s, n)
    y += rng.normal(0, s, n)
    z += rng.normal(0, s, n)
    return x, y, z


def score(x, y, z, k):
    """候选 k 的对齐度：1 = 完全对齐，0 = 完全散开"""
    p1 = np.arccos(np.clip(z, -1, 1))
    p2 = 2 * np.pi - p1
    th = np.arctan2(y, x)
    d1 = np.angle(np.exp(1j * (th - k * p1)))
    d2 = np.angle(np.exp(1j * (th - k * p2)))
    d = np.where(np.abs(d1) < np.abs(d2), d1, d2)
    return np.abs(np.exp(1j * d).mean())


def scan_noise():
    print(f"真 k = {K_TRUE}，{N:,} 颗星，四维双旋转 → 三维投影\n")
    print(f"{'噪声':>5} | {'反推':>4} | {'真数据峰':>8} | {'打乱峰':>7} | "
          f"{'比值':>5} | {'信号':>6} | 判定")
    print("-" * 70)
    for nf in NOISE_LEVELS:
        x, y, z = make(K_TRUE, nf)
        s = np.array([score(x, y, z, k) for k in KS])
        k_hat = int(KS[np.argmax(s)])

        zr = rng.permutation(z)                      # 真打乱：只重排 z
        sr = np.array([score(x, y, zr, k) for k in KS])

        ratio = s.max() / sr.max()
        signal = s.max() - sr.max()                  # 扣掉地板之后剩下的
        ok = "✅" if (k_hat == K_TRUE and ratio > 1.5) else "❌"
        print(f"{nf*100:>4.0f}% | {k_hat:>4} | {s.max():>8.3f} | {sr.max():>7.3f} | "
              f"{ratio:>5.2f} | {signal:>6.3f} | {ok}")

    # 15% 下的完整排名
    x, y, z = make(K_TRUE, 0.15)
    s = np.array([score(x, y, z, k) for k in KS])
    print("\n15% 噪声下前五名：")
    for i in np.argsort(s)[::-1][:5]:
        print(f"  k={KS[i]:>2}  对齐度 {s[i]:.3f}")


# ---------------- 聚类对照 ----------------
def kmeans(P, k, iters=12):
    cent = P[rng.choice(len(P), k, replace=False)].copy()
    for _ in range(iters):
        d = ((P[:, None, :] - cent[None, :, :]) ** 2).sum(-1)
        lab = d.argmin(1)
        for j in range(k):
            m = lab == j
            if m.any():
                cent[j] = P[m].mean(0)
    return lab, ((P - cent[lab]) ** 2).sum()


def silhouette(P, lab, sample=1500):
    ii = rng.choice(len(P), min(sample, len(P)), replace=False)
    Ps, Ls = P[ii], lab[ii]
    D = np.sqrt(((Ps[:, None, :] - Ps[None, :, :]) ** 2).sum(-1))
    out = np.zeros(len(ii))
    for j in range(len(ii)):
        same = Ls == Ls[j]
        same[j] = False
        if not same.any():
            continue
        a = D[j, same].mean()
        b = min(D[j, Ls == c].mean() for c in np.unique(Ls) if c != Ls[j])
        out[j] = (b - a) / max(a, b)
    return out.mean()


def cluster_control():
    print("\n聚类对照（手写 k-means，抽 2 万点，噪声 2%，真 k=5）")
    x, y, z = make(K_TRUE, 0.02)
    idx = rng.choice(N, 20000, replace=False)
    P = np.column_stack([x[idx], y[idx], z[idx]])

    print(f"{'簇数':>4} | {'惯性':>10} | {'轮廓系数':>8}")
    best_k, best_s = None, -2
    for kk in range(2, 11):
        lab, inert = kmeans(P, kk)
        sc = silhouette(P, lab)
        print(f"{kk:>4} | {inert:>10.1f} | {sc:>8.3f}")
        if sc > best_s:
            best_s, best_k = sc, kk
    print(f"\n轮廓系数选出的簇数 = {best_k}（真 k = {K_TRUE}）")
    print("→ 聚类找团，不找周期；它没有能力自己知道该分几堆")


if __name__ == "__main__":
    scan_noise()
    cluster_control()
