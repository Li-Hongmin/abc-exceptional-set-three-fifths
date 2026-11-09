"""BLT v2 MILP (arXiv:2410.12234v2, Prop. 6.1) + the 24 two-letter (2,3)/(3,4) collision constraints (star).

Run with:
    uv run --no-project --with 'pulp==3.3.0' python units/abc/computations/blt-lp/lp_flagship.py 1
"""
import itertools
import sys
import pulp

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import lp as base


def build_tight(lam=1.0, thue=True, det=True, extra_geo=()):
    """Build the BLT exponent MILP with tight Big-M constants (M=1, 2)."""
    P = pulp.LpProblem("blt_flagship", pulp.LpMaximize)
    V = {}
    for s in "abc":
        for i in base.R:
            V[s, i] = pulp.LpVariable(f"{s}{i}", 0)
        V[s] = pulp.LpVariable(s, 0)
    a, b, c = V["a"], V["b"], V["c"]
    D = pulp.LpVariable("D", 0, 1)
    L = a + b + c
    P += D
    P += L <= lam
    for s in "abc":
        P += pulp.lpSum(V[s, i] for i in base.R) <= V[s]
        P += 7 * V[s] - 1 <= pulp.lpSum((7 - i) * V[s, i] for i in base.R)
    P += D <= a + b
    P += D <= a + c
    P += D <= b + c
    zc = [0]

    def geo(I, J, K):
        nonlocal P
        s1 = (
            pulp.lpSum(V["a", i] for i in I)
            + pulp.lpSum(V["b", i] for i in J)
            + pulp.lpSum(V["c", i] for i in K)
        )
        s2 = (
            pulp.lpSum((i - 1) * V["a", i] for i in I)
            + pulp.lpSum((i - 1) * V["b", i] for i in J)
            + pulp.lpSum((i - 1) * V["c", i] for i in K)
        )
        zc[0] += 1
        z = pulp.LpVariable(f"z{zc[0]}", cat="Binary")
        P += s1 <= L - D + 1.0 * z
        P += 1 - s2 <= L - D + 2.0 * (1 - z)

    base_cuts = [((1, 2, 4),) * 3]
    S12, S123, S13 = (1, 2), (1, 2, 3), (1, 3)
    for t in set(itertools.permutations([S12, S12, S123])):
        base_cuts.append(t)
    for t in set(itertools.permutations([S12, S123, S123])):
        base_cuts.append(t)
    for t in set(itertools.permutations([S12, S12, S13])):
        base_cuts.append(t)
    for t in set(itertools.permutations([S12, S13, S13])):
        base_cuts.append(t)
    for t in base_cuts + list(extra_geo):
        geo(*t)

    F = {}
    for s in "abc":
        F[s] = pulp.LpVariable("F" + s, 0)
        for j in (2, 3, 5):
            P += pulp.lpSum(V[s, i] for i in base.R if i % j == 0) <= F[s]
    P += F["a"] + F["b"] + F["c"] <= 4 * L - 6 * D

    if thue:
        for x, y in [("a", "b"), ("a", "c"), ("b", "c")]:
            for p in range(2, 7):
                P += D <= L - pulp.lpSum(
                    V[x, i] + V[y, i] for i in base.R if i % p == 0
                )
    if det:
        for x, y in itertools.permutations("abc", 2):
            for p in base.R:
                for q in base.R:
                    P += D <= L - V[x, p] - V[y, q] + V[x, p] * (1.0 / q)
                    P += D <= L - V[x, p] - V[y, q] + V[y, q] * (1.0 / p)
    return P, V, D


def add_one_two_letter(P, V, D, A, pA, B, pB):
    """Add the two-letter collision constraint for ordered letters (A, B) freeing pA, pB in {(2,3),(3,4)}."""
    C = ({"a", "b", "c"} - {A, B}).pop()
    tag = f"{A}{pA[0]}{pA[1]}_{B}{pB[0]}{pB[1]}"
    p, q = V[A, pA[0]], V[A, pA[1]]
    r, s = V[B, pB[0]], V[B, pB[1]]
    outer = V[A] - p - q + V[B] - r - s
    w = pulp.lpSum(i * V[B, i] for i in base.R if i not in pB) + 7 * (
        V[B] - pulp.lpSum(V[B, i] for i in base.R)
    )
    zc = pulp.LpVariable(f"fl_cauchy_{tag}", cat="Binary")
    positives = []
    for idx, coordinate in enumerate((p, q)):
        h = pulp.LpVariable(f"fl_pos_{tag}_{idx}", 0)
        z = pulp.LpVariable(f"fl_sign_{tag}_{idx}", cat="Binary")
        P += h >= 2 * coordinate - w
        P += h <= 2 * coordinate - w + 1.0 * (1 - z)
        P += h <= 1.0 * z
        positives.append(h)
    kappa = p + q + 0.5 * pulp.lpSum(positives)
    P += D <= outer + V[C] + 1.0 * zc
    P += D <= outer + 0.5 * (V[C] + 2 * r + kappa) + 1.0 * (1 - zc)
    P += D <= outer + 0.5 * (V[C] + 2 * s + kappa) + 1.0 * (1 - zc)
    return P


def add_flagship_two_letter(P, V, D, pairs_exp=((2, 3), (3, 4))):
    """Add all ordered two-letter cuts with pA, pB in pairs_exp."""
    for A, B in itertools.permutations("abc", 2):
        for pA in pairs_exp:
            for pB in pairs_exp:
                add_one_two_letter(P, V, D, A, pA, B, pB)
    return P


def separate_two_letter(x, pairs_exp=((2, 3), (3, 4))):
    """Find the most violated (2,3)/(3,4) two-letter collision cut at x."""
    best = None
    for A, B in itertools.permutations("abc", 2):
        C = ({"a", "b", "c"} - {A, B}).pop()
        for pA in pairs_exp:
            for pB in pairs_exp:
                p, q = x[A, pA[0]], x[A, pA[1]]
                r, s = x[B, pB[0]], x[B, pB[1]]
                outer = x[A] - p - q + x[B] - r - s
                w = sum(i * x[B, i] for i in base.R if i not in pB) + 7 * (
                    x[B] - sum(x[B, i] for i in base.R)
                )
                kappa = p + q + 0.5 * (max(0.0, 2 * p - w) + max(0.0, 2 * q - w))
                cap = min(r + s, 2 * s, 2 * r)
                v = outer + max(x[C], 0.5 * (x[C] + cap + kappa))
                if best is None or v < best[0]:
                    best = (v, (A, pA, B, pB))
    return best


INITIAL_GEO_CUTS = [
    ((1,), (4,), (1,)),
    ((1,), (5,), (1, 4)),
    ((1,), (1,), (4,)),
    ((1,), (1, 4), (1,)),
    ((1,), (1,), (1, 4)),
    ((1, 5), (2, 4), (1, 3, 5)),
    ((1, 2, 4), (1, 5), (3,)),
    ((1, 5), (1, 2, 3, 4), ()),
    ((), (1, 5), (1, 2, 3, 4)),
    ((), (1, 2, 3, 4), (1, 5)),
    ((1, 5), (), (1, 2, 3, 4)),
    ((1, 5), (1, 2, 4), (1, 3)),
    ((1, 2, 4), (1,), (1, 3, 5)),
    ((1,), (1, 3, 5), (1, 2, 4)),
    ((1, 4), (1, 2), (1,)),
    ((1,), (1, 4), (1, 2)),
    ((1, 2), (1,), (1, 4)),
    ((1, 2, 4), (1,), (1, 3)),
    ((1,), (1,), (1, 2, 3, 4)),
    ((1,), (1, 2, 4), (1, 3)),
    ((1, 2, 3), (1, 2, 3, 5), (1, 4)),
]

INITIAL_TL_CUTS = [
    ("b", (3, 4), "c", (2, 3)),
    ("c", (2, 3), "b", (3, 4)),
    ("c", (3, 4), "b", (2, 3)),
    ("a", (2, 3), "b", (3, 4)),
    ("b", (2, 3), "a", (2, 3)),
    ("c", (2, 3), "a", (2, 3)),
    ("a", (3, 4), "b", (2, 3)),
    ("b", (2, 3), "c", (3, 4)),
    ("a", (2, 3), "c", (3, 4)),
    ("b", (2, 3), "a", (3, 4)),
    ("a", (2, 3), "b", (2, 3)),
    ("c", (2, 3), "a", (3, 4)),
    ("a", (3, 4), "c", (2, 3)),
    ("b", (3, 4), "a", (2, 3)),
    ("c", (3, 4), "a", (2, 3)),
]


def run(lam=1.0, th=True, de=True, maxit=4, verbose=False):
    """Run lazy cutting-plane separation over geometric and (2,3)/(3,4) two-letter cuts."""
    extra_geo = list(INITIAL_GEO_CUTS)
    extra_tl = list(INITIAL_TL_CUTS)
    for it in range(maxit):
        P, V, D = build_tight(lam, thue=th, det=de, extra_geo=extra_geo)
        for tl_cut in extra_tl:
            add_one_two_letter(P, V, D, *tl_cut)
        P.solve(pulp.PULP_CBC_CMD(msg=False, threads=1))
        if pulp.LpStatus[P.status] != "Optimal":
            raise RuntimeError(f"MILP failed: {pulp.LpStatus[P.status]}")
        x = {k: V[k].value() for k in V}
        Dv = D.value()
        vg, ttg = base.separate(x, Dv, lam)
        vtl, tttl = separate_two_letter(x)
        if verbose:
            print(
                f"  it={it:2d} D={Dv:.9f} vg={vg:.6f} vtl={vtl:.6f} geo={len(extra_geo)} tl={len(extra_tl)}",
                flush=True,
            )
        if vg < Dv - 1e-7 or vtl < Dv - 1e-7:
            if vg <= vtl:
                if ttg not in extra_geo:
                    extra_geo.append(ttg)
            else:
                if tttl not in extra_tl:
                    extra_tl.append(tttl)
        else:
            break
    return pulp.LpStatus[P.status], Dv, x, extra_geo, extra_tl, P


if __name__ == "__main__":
    lams = [float(a) for a in sys.argv[1:] if not a.startswith("-")] or [1.0]
    verbose = "--verbose" in sys.argv
    for lam in lams:
        st, Dv, x, ex_geo, ex_tl, P = run(lam, verbose=verbose)
        print(
            f"lam={lam} flagship_two_letter: {st} D={Dv:.9f} (< 0.6: {Dv < 0.6 - 1e-6}) geo_cuts={len(ex_geo)} tl_cuts={len(ex_tl)}",
            flush=True,
        )
        print(
            "   ",
            {
                f"{k[0]}{k[1]}" if isinstance(k, tuple) else k: round(v, 5)
                for k, v in x.items()
                if v and v > 1e-6
            },
        )
