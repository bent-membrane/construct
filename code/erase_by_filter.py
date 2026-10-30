#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
擦除在图上的样子 —— 不改任何一个落点，只分堆
================================================
对应 9.4 与 9.7：擦除不是倒流，是分类。

做法：
    1. 造 20000 个光子。每个光子带一个隐藏标签（被标记 / 没被标记）。
    2. 被标记的那部分（比例 κ）按 V=0 落点（无条纹）；
       没被标记的那部分（1-κ）按 V=1 落点（有条纹）。
    3. 混在一起看：可见度是 V = 1-κ(1-cos θ)，落在中间。
    4. 按隐藏标签分成两堆重画：一堆条纹回来，一堆是平缓曲线。

关键：两堆用的落点，跟混在一起时**一模一样**。没有改过去，没有倒流，
只是把已有的数据按标签分了一次组。

这直接演示：延迟擦除 = 事后分类。

用法：
    python erase_by_filter.py
    python erase_by_filter.py --plot ../figures/48_erase_by_filter.png
"""
import argparse
import os
import numpy as np

LAM, A, D, L = 0.5e-6, 0.1e-3, 0.5e-3, 1.0
KAPPA = 0.5          # 被标记的比例
N = 20000


def envelope(x):
    u = np.pi * A * x / (LAM * L)
    out = np.ones_like(x)
    nz = np.abs(u) > 1e-12
    out[nz] = (np.sin(u[nz]) / u[nz]) ** 2
    return out


def prob(x, V):
    return envelope(x) * 0.5 * (1.0 + V * np.cos(2 * np.pi * D * x / (LAM * L)))


def sample(n, V, rng, grid):
    p = prob(grid, V)
    p = p / p.sum()
    return rng.choice(grid, size=n, p=p)


def visibility_of(xs, half=3.5e-3):
    """按包络加权后取条纹频率分量：V = 2|Σ w·e^{-ikx}| / Σw，w=1/包络。"""
    m = np.abs(xs) < half
    if m.sum() < 100:
        return float("nan")
    x2 = xs[m]
    w = 1.0 / np.maximum(envelope(x2), 1e-3)
    k = 2 * np.pi * D / (LAM * L)
    return float(2 * abs(np.sum(w * np.exp(-1j * k * x2))) / np.sum(w))


def numbers():
    rng = np.random.default_rng(11)
    grid = np.linspace(-6e-3, 6e-3, 6000)

    tagged = sample(int(N * KAPPA), 0.0, rng, grid)      # 被标记：无条纹
    untagged = sample(N - int(N * KAPPA), 1.0, rng, grid)  # 没被标记：条纹
    mixed = np.concatenate([tagged, untagged])

    print(f"N={N}, κ={KAPPA}")
    print(f"  混在一起   V≈{visibility_of(mixed):.3f}")
    print(f"  被标记堆   V≈{visibility_of(tagged):.3f}")
    print(f"  未标记堆   V≈{visibility_of(untagged):.3f}")
    print(f"  两堆计数   {len(tagged)} / {len(untagged)}")


def plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    rng = np.random.default_rng(11)
    grid = np.linspace(-6e-3, 6e-3, 6000)
    n_tag = int(N * KAPPA)
    tagged = sample(n_tag, 0.0, rng, grid)
    untagged = sample(N - n_tag, 1.0, rng, grid)
    mixed = np.concatenate([tagged, untagged])

    fig, axes = plt.subplots(3, 1, figsize=(11, 10.5), sharex=True)
    bins = np.linspace(-6, 6, 121)

    panels = [
        (mixed, "① 全部 20000 个光子混着看：可见度居中", "#1f4e79"),
        (tagged, f"② 按标签挑出「被标记」的 {n_tag} 个：条纹没了", "#c07000"),
        (untagged, f"③ 挑出「没被标记」的 {N-n_tag} 个：条纹回来", "#2e7d32"),
    ]
    for ax, (data, title, color) in zip(axes, panels):
        ax.hist(data * 1e3, bins=bins, color=color, alpha=0.78)
        ax.set_title(title, fontsize=11.5)
        ax.set_ylabel("个数")
        ax.grid(alpha=0.2)
    axes[-1].set_xlabel("屏幕位置 x（mm）")

    fig.suptitle("擦除 = 分堆。三张图用的是同一批落点，一个都没改", fontsize=13, y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()
    if args.plot:
        os.makedirs(os.path.dirname(os.path.abspath(args.plot)), exist_ok=True)
        plot(args.plot)
        return
    numbers()


if __name__ == "__main__":
    main()
