# Global crystallization of dilute alternating Gaussian chains

[English paper (PDF)](paper-en.pdf) · [English proof (Markdown)](paper-en.md) · [中文论文 (PDF)](paper-zh.pdf) · [中文全文](paper-zh.md)

For an alternating chain with like-species attraction and opposite-species repulsion of Gaussian magnitude $e^{-ar^2}$, we prove that **equal spacing uniquely minimizes the unrestricted periodic energy whenever $ah^2\ge4$**, where $h$ is the mean spacing. The result covers every even period and allows gaps arbitrarily close to zero.

**Rigidity.** At fixed density in this range, if the energy per particle tends to its minimum as the period grows, the average absolute deviation of the gaps from $h$ tends to zero. The proof gives a lower bound on the energy excess in terms of the number of particles removed in closest-pair deletions and the deviations of the remaining gaps.

**Where the boundary lies.** Combined with lower-energy competitors in the [companion manuscript](https://github.com/gamehhh7/alternating-chain-stability/tree/3a7012e285a7482a9caadc5130cc21bbe8057172), this brackets the onset of eventual global crystallization between approximately 2.474167 and 4. The precise transition is open. The question originates in [Bétermin, Knüpfer and Nolte, Conjecture 2.7](https://doi.org/10.1007/s10955-020-02603-2).

## Reproduce the checks

Run `python verify.py` with Python 3.9 or later. It checks the proof's elementary constants with exact fractions. A separate floating-point calculation compares two independent sums for 607 cells, tests 1,120 closest-pair deletions including cells crossing the origin, and checks 240 gap Hessians. The proofs in the manuscripts cover all periods and configurations.

## 中文摘要

对于同类相吸、异类相斥的交替高斯链，只要 $ah^2\ge4$，等距排列就是任意偶数周期下唯一的最低能排列。即使两个粒子可以任意接近，结论仍然成立。在这一密度范围内，若每粒子能量趋于最低值，全部间距与平均值的绝对偏差的平均也趋于零。

AI-generated research manuscript, 26 November 2026. Independent peer review and literature-priority confirmation are outstanding.
