import pulp, itertools, sys
R = range(1, 7)
def build(lam=1.0, thue=False, det=False, z13fix=True, extra_geo=()):
    P = pulp.LpProblem("x", pulp.LpMaximize)
    V = {}
    for s in 'abc':
        for i in R: V[s, i] = pulp.LpVariable(f"{s}{i}", 0)
        V[s] = pulp.LpVariable(s, 0)
    a, b, c = V['a'], V['b'], V['c']
    D = pulp.LpVariable("D"); L = a + b + c
    P += D
    P += L <= lam
    for s in 'abc':
        P += pulp.lpSum(V[s, i] for i in R) <= V[s]
        P += 7 * V[s] - 1 <= pulp.lpSum((7 - i) * V[s, i] for i in R)
    P += D <= a + b; P += D <= a + c; P += D <= b + c
    M = 100; zc = [0]
    def geo(I, J, K):
        nonlocal P
        s1 = pulp.lpSum(V['a', i] for i in I) + pulp.lpSum(V['b', i] for i in J) + pulp.lpSum(V['c', i] for i in K)
        s2 = pulp.lpSum((i - 1) * V['a', i] for i in I) + pulp.lpSum((i - 1) * V['b', i] for i in J) + pulp.lpSum((i - 1) * V['c', i] for i in K)
        zc[0] += 1; z = pulp.LpVariable(f"z{zc[0]}", cat="Binary")
        P += s1 <= L - D + M * z          # term1: L - s1 >= D fails ...
        P += 1 - s2 <= L - D + M * (1 - z)
    base = [((1,2,4),)*3]
    S12, S123, S13 = (1,2), (1,2,3), (1,3)
    for t in set(itertools.permutations([S12, S12, S123])): base.append(t)
    for t in set(itertools.permutations([S12, S123, S123])): base.append(t)
    for t in set(itertools.permutations([S12, S12, S13])): base.append(t)
    for t in set(itertools.permutations([S12, S13, S13])): base.append(t)
    for t in base + list(extra_geo): geo(*t)
    if not z13fix:  # paper's appendix line 116-117 (typo'd) instead of the correct 12,13,13 (a)
        pass
    F = {}
    for s in 'abc':
        F[s] = pulp.LpVariable("F" + s, 0)
        for j in (2, 3, 5): P += pulp.lpSum(V[s, i] for i in R if i % j == 0) <= F[s]
    P += F['a'] + F['b'] + F['c'] <= 4 * L - 6 * D
    if thue:  # v1 Prop 3.6 for each pair of letters, p>=2: D <= L - sum_{p|i}(x_i+y_i)
        for x, y in [('a','b'), ('a','c'), ('b','c')]:
            for p in range(2, 7):
                P += D <= L - pulp.lpSum(V[x, i] + V[y, i] for i in R if i % p == 0)
    if det:  # v1 Prop 3.5: D <= L - x_p - y_q + min(x_p/q, y_q/p), both orders of letters
        for x, y in itertools.permutations('abc', 2):
            for p in R:
                for q in R:
                    P += D <= L - V[x, p] - V[y, q] + V[x, p] * (1.0 / q)
                    P += D <= L - V[x, p] - V[y, q] + V[y, q] * (1.0 / p)
    return P, V, D

def sol(P, V, D):
    P.solve(pulp.PULP_CBC_CMD(msg=False))
    return pulp.LpStatus[P.status], D.value(), {k: V[k].value() for k in V}

def separate(x, Dv, lam):
    """find subsets I,J,K of {1..6} with both geometry terms < D (violated)."""
    L = x['a'] + x['b'] + x['c']
    best = None
    opts = []
    for s in 'abc':
        o = []
        for m in range(64):
            I = tuple(i for i in R if m >> (i - 1) & 1)
            o.append((I, sum(x[s, i] for i in I), sum((i - 1) * x[s, i] for i in I)))
        opts.append(o)
    for A in opts[0]:
        for B in opts[1]:
            for C in opts[2]:
                s1 = A[1] + B[1] + C[1]; s2 = A[2] + B[2] + C[2]
                v = max(L - s1, L - 1 + s2)
                if best is None or v < best[0]: best = (v, (A[0], B[0], C[0]))
    return best

if __name__ == "__main__":
    lam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    for th, de in [(0,0),(1,0),(0,1),(1,1)]:
        extra = []
        for it in range(40):
            st, Dv, x = sol(*build(lam, th, de, extra_geo=extra))
            v, t = separate(x, Dv, lam)
            if v < Dv - 1e-7: extra.append(t)
            else: break
        print(f"lam={lam} thue={th} det={de}: {st} D={Dv:.6f}  geo_cuts_added={len(extra)} best_geo_slack={v-Dv:+.2e}")
        if th == 0 and de == 0 and not extra: pass
        print("   ", {f"{k[0]}{k[1]}" if isinstance(k, tuple) else k: round(val, 5) for k, val in x.items() if val and val > 1e-6})
