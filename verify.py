"""Reproduce the rational certificates and exercise the cell estimates (stdlib only)."""

from fractions import Fraction as Q
from itertools import accumulate
from math import ceil, exp, factorial, fsum, log, sqrt
from random import Random


D = sqrt(log(3))
GAMMA = 1 / 3 - 18 * exp(-4)


def certify_constants():
    # Lower partial sum of log(3) = 2 atanh(1/2).
    lower_log = 2 * (Q(1, 2) + Q(1, 2) ** 3 / 3)
    assert lower_log == Q(13, 12)
    second = 4 * (8 * lower_log - 1) / ((2 * lower_log - 1) * 27)
    ratio = Q(5, 4) ** 4 / 3**9
    tail = Q(4 * 4**4, 3**15)
    assert second == Q(184, 189)
    assert ratio < Q(1, 2) and second + tail < Q(49, 50)
    assert sum(Q(4**m, factorial(m)) for m in range(10)) > 54
    assert 1 + Q(3, 4) + Q(3, 4) ** 2 / 2 > 2


def kernel(r):
    return exp(-r * r)


def curvature(r):
    return (4 * r * r - 2) * kernel(r)


def image_energy(points, length):
    # Independent signed sum over neighboring image cells, with a Gaussian cutoff.
    if not points:
        return 0.0
    reach = ceil(10 / length) + 1
    total = []
    for i, x in enumerate(points):
        for j, y in enumerate(points):
            sign = 1 if (i - j) % 2 == 0 else -1
            total.append(sign * fsum(kernel(x - y + k * length)
                                     for k in range(-reach, reach + 1)))
    return len(points) - fsum(total)


def gaps_for(points, length):
    if not points:
        return []
    return [y - x for x, y in zip(points, points[1:] + [points[0] + length])]


def windows(gaps, start, limit):
    dist = 0.0
    outer = 1
    while dist < limit:
        dist += gaps[(start + outer - 1) % len(gaps)]
        yield outer, dist
        outer += 1


def gap_energy(gaps):
    return 2 * fsum((1 if m % 2 else -1) * kernel(r)
                    for i in range(len(gaps)) for m, r in windows(gaps, i, 10))


def uniform(h):
    e = 2 * fsum((-1)**(m + 1) * kernel(m * h) for m in range(1, 12))
    slope = -4 * h * fsum((-1)**(m + 1) * m**2 * kernel(m * h)
                          for m in range(1, 12))
    return e, slope


def test_cell(points, length):
    n = len(points)
    h = length / n
    original = image_energy(points, length)
    assert abs(original - gap_energy(gaps_for(points, length))) < 2e-11 * n

    sites = points[:]
    energy, steps = original, 0
    while sites:
        gaps = gaps_for(sites, length)
        idx = min(range(len(gaps)), key=gaps.__getitem__)
        if gaps[idx] >= D:
            break
        remove = {idx, (idx + 1) % len(sites)}
        sites = [x for j, x in enumerate(sites) if j not in remove]
        after = image_energy(sites, length)
        lower = 1 if gaps[idx] < 0.5 else 2 * kernel(gaps[idx])
        assert energy - after >= lower - 1e-10
        energy = after
        steps += 1

    m = len(sites)
    e, slope = uniform(h)
    bound = (n - m) * (1 / 3 - e + h * slope)
    if m:
        s = length / m
        tangents = fsum(kernel(t) - kernel(s) + 2 * s * kernel(s) * (t - s)
                        for t in gaps_for(sites, length))
        bound += tangents / 25
    assert original - n * e >= bound - 2e-10 * n
    return steps


def check_hessian(gaps, direction):
    weighted = fsum(curvature(t) * u**2 for t, u in zip(gaps, direction))
    terms = []
    for i in range(len(gaps)):
        change = 0.0
        for m, r in windows(gaps, i, 12):
            change += direction[(i + m - 1) % len(gaps)]
            terms.append((1 if m % 2 else -1) * curvature(r) * change**2)
    assert fsum(terms) >= weighted / 50 - 1e-12


def main():
    certify_constants()
    rng = Random(170803)
    special = [
        ([0.0, 1e-9], 4.0),
        ([0.0, 3.999999], 4.0),  # wraparound deletion
        ([0.0, 0.5], 4.0),
        ([0.0, D - 1e-10], 4.0),
        ([0.0, D + 1e-10], 4.0),
        ([0.0, 1.48, 1.52, 3.0, 6.0, 9.0], 12.0),
        ([0.0, 1e-8, 2e-8, 3e-8, 4e-8, 5e-8], 12.0),
    ]
    count, steps = 0, 0
    for points, length in special:
        steps += test_cell(points, length)
        count += 1
    for n in (2, 4, 6, 10, 16):
        for h in (2.0, 2.03, 2.5, 4.0):
            for shape in (0.12, 1.0, 4.0):
                for _ in range(10):
                    weights = [max(rng.gammavariate(shape, 1), 1e-8) for _ in range(n)]
                    gaps = [n * h * t / sum(weights) for t in weights]
                    points = [0.0, *accumulate(gaps[:-1])]
                    steps += test_cell(points, n * h)
                    count += 1
    hessians = 0
    for n in (2, 4, 10, 18):
        for scale in (0, 0.01, 0.5, 2.0):
            for _ in range(15):
                gaps = [D + scale * rng.random() for _ in range(n)]
                direction = [rng.uniform(-1, 1) for _ in gaps]
                check_hessian(gaps, direction)
                hessians += 1
    print("Exact rational inequalities passed.")
    print(f"Floating-point cross-checks: {count} cells, {steps} pair deletions,"
          f" {hessians} gap Hessians.")
    print(f"Uniform energy at h=2: {uniform(2)[0]:.12f}")
    print(f"Uniform deletion coefficient lower bound: {GAMMA:.12f}")


if __name__ == "__main__":
    main()
