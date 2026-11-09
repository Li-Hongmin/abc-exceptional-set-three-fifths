"""Independent verification, in rational arithmetic, of the branch-and-bound certificate D <= 59969/100033.

Standard-library verification:
    python3 verify_certificate.py
Optional comparison with both actual PuLP builders (no solver invocation):
    python3 verify_certificate.py --compare-models   (requires pulp)

This checks algebra, coverage, dual bounds and explicit witnesses, not analytic
number theory. No code is imported from the model builders except in --compare-models.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product
import json
from math import lcm
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
NAMES = [v for s in 'abc' for v in [*(f'{s}{i}' for i in range(1, 7)), s]] + ['Fa', 'Fb', 'Fc', 'D']
POS = {v: i for i, v in enumerate(NAMES)}
LIMIT = F(59969, 100033)
ZERO = (F(0),) * len(NAMES)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def var(name):
    return tuple(F(i == POS[name]) for i in range(len(NAMES)))


def add(*vectors):
    return tuple(map(sum, zip(*vectors))) if vectors else ZERO


def scale(n, vector):
    return tuple(F(n) * x for x in vector)


def subtract(a, b):
    return add(a, scale(-1, b))


def row(lhs, rhs=0):
    return lhs, F(rhs)


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


L = add(*(var(s) for s in 'abc'))
D = var('D')


def base_rows():
    # Order is the certificate's public row-index convention. All x >= 0.
    rows = [row(L, 1)]
    for s in 'abc':
        coords = add(*(var(f'{s}{i}') for i in range(1, 7)))
        tail = subtract(var(s), coords)
        height = add(*(scale(i, var(f'{s}{i}')) for i in range(1, 7)), scale(7, tail))
        rows += [row(subtract(coords, var(s))), row(height, 1)]
    rows += [row(subtract(D, add(var(a), var(b)))) for a, b in [('a', 'b'), ('a', 'c'), ('b', 'c')]]
    for s in 'abc':
        for divisor in (2, 3, 5):
            rows.append(row(subtract(add(*(var(f'{s}{i}') for i in range(1, 7) if i % divisor == 0)), var(f'F{s}'))))
    rows.append(row(add(*(var(f'F{s}') for s in 'abc'), scale(-4, L), scale(6, D))))
    # Every orbit of the full symmetric constraint family has this ordering.
    for a, b in [('a', 'b'), ('b', 'c')]:
        rows.append(row(subtract(add(var(f'{b}4'), scale(2, var(f'{b}5'))), add(var(f'{a}4'), scale(2, var(f'{a}5'))))))
    return rows


def subsets(node):
    sets = []
    for key in ('I', 'J', 'K'):
        values = node[key]
        require(isinstance(values, list) and all(type(i) is int and 1 <= i <= 6 for i in values), f'Invalid subset {key}')
        require(len(set(values)) == len(values), 'Repeated subset index')
        sets.append(values)
    return sets


def geometric_branches(sets):
    selected = [(s, i) for s, indices in zip('abc', sets) for i in indices]
    first_bound = subtract(L, add(*(var(f'{s}{i}') for s, i in selected)))
    # Second affine bound is L + sum((i-1)*s_i) - 1.
    second_linear = add(L, *(scale(i - 1, var(f'{s}{i}')) for s, i in selected))
    return [[row(subtract(D, first_bound))], [row(subtract(D, second_linear), -1)]]


def collision_expressions(A, pairA, B, pairB):
    require(A in tuple('abc') and B in tuple('abc') and A != B, 'Invalid letters')
    require(tuple(pairA) in ((2, 3), (3, 4)) and tuple(pairB) in ((2, 3), (3, 4)), 'Invalid exponent pair')
    C = next(s for s in 'abc' if s not in (A, B))
    p, q = (var(f'{A}{i}') for i in pairA)
    r, s = (var(f'{B}{i}') for i in pairB)
    frozen = subtract(add(var(A), var(B)), add(p, q, r, s))
    tail = subtract(var(B), add(*(var(f'{B}{i}') for i in range(1, 7))))
    modulus_lower = add(scale(7, tail), *(scale(i, var(f'{B}{i}')) for i in range(1, 7) if i not in pairB))
    positive_arguments = [subtract(scale(2, p), modulus_lower), subtract(scale(2, q), modulus_lower)]
    return C, p, q, r, s, frozen, positive_arguments


def collision_branches(A, pairA, B, pairB):
    C, p, q, r, s, frozen, arguments = collision_expressions(A, pairA, B, pairB)
    branches = [[row(subtract(D, add(frozen, var(C))))]]
    # max(c, h) permits the diagonal alternative OR the energy alternative.
    # Four CLOSED sign quadrants exhaust the latter, including all boundaries.
    for signs in product((0, 1), repeat=2):
        signs_rows = [row(scale(1 if sign == 0 else -1, arg)) for sign, arg in zip(signs, arguments)]
        kappa = add(p, q, *(scale(F(sign, 2), arg) for sign, arg in zip(signs, arguments)))
        # min(2r,2s) <= r+s; the omitted cap is redundant for all r,s.
        energy_rows = [row(subtract(D, add(frozen, scale(F(1, 2), add(var(C), kappa, scale(2, cap)))))) for cap in (r, s)]
        branches.append(signs_rows + energy_rows)
    return branches


def read_certificate(path):
    def unique_object(pairs):
        obj = {}
        for key, value in pairs:
            require(key not in obj, f'Duplicate JSON key {key}')
            obj[key] = value
        return obj
    return json.loads(path.read_text(), object_pairs_hook=unique_object)


def verify_tree(tree):
    counts = Counter()
    bounds = []
    max_depth = 0
    root = base_rows()

    def visit(node, rows, path):
        nonlocal max_depth
        kind = node['type']
        counts[kind] += 1
        max_depth = max(max_depth, len(path))
        if kind == 'leaf':
            require(set(node) == {'type', 'y', 'bound'}, f'Unknown leaf fields at {path}')
            combined = list(ZERO)
            bound = F(0)
            require(isinstance(node['y'], dict) and node['y'], f'Empty dual at {path}')
            for key, value in node['y'].items():
                i = int(key)
                require(str(i) == key and 0 <= i < len(rows), f'Invalid dual index at {path}')
                y = F(value)
                require(y > 0, f'Nonpositive dual at {path}')
                coefficients, rhs = rows[i]
                for j, coefficient in enumerate(coefficients):
                    combined[j] += y * coefficient
                bound += y * rhs
            require(combined[POS['D']] == 1, f'Wrong objective coefficient at {path}')
            require(all(value >= 0 for j, value in enumerate(combined) if j != POS['D']), f'Negative remainder at {path}')
            require(bound == F(node['bound']), f'Incorrect declared dual objective at {path}')
            require(bound <= LIMIT < F(3, 5), f'Insufficient bound at {path}')
            bounds.append(bound)
            return
        if kind == 'geo':
            require(set(node) == {'type', 'I', 'J', 'K', 'children'}, 'Unknown geometry fields')
            choices = geometric_branches(subsets(node))
        elif kind == 'tl':
            require(set(node) == {'type', 'A', 'pA', 'B', 'pB', 'children'}, 'Unknown collision fields')
            choices = collision_branches(node['A'], node['pA'], node['B'], node['pB'])
        else:
            raise ValueError(f'Unknown node kind {kind}')
        require(isinstance(node['children'], list) and len(node['children']) == len(choices), f'Incomplete branching at {path}')
        for branch, (extra_rows, child) in enumerate(zip(choices, node['children'])):
            visit(child, rows + extra_rows, path + (branch,))

    visit(tree, root, ())
    require(counts == {'leaf': 777, 'geo': 520, 'tl': 64}, f'Unexpected tree {counts}')
    require(max(bounds) == LIMIT, 'Incorrect maximum bound')
    print(f'PASS: {counts["leaf"]} exact duals; {counts["geo"]} exhaustive geometric and {counts["tl"]} exhaustive collision splits; depth {max_depth}.')
    print(f'PASS: conditional D <= {LIMIT} < 1199/2000 < 3/5.')
    require(LIMIT < F(1199, 2000), 'Incorrect rounded exponent')


def point(letters, exponent):
    values = {}
    for s, coordinates in zip('abc', letters):
        require(len(coordinates) == 6, 'Six coordinates required')
        values.update({f'{s}{i}': F(x) for i, x in enumerate(coordinates, 1)})
        values[s] = sum(coordinates)
        values[f'F{s}'] = max(sum(coordinates[i - 1] for i in range(1, 7) if i % j == 0) for j in (2, 3, 5))
    values['D'] = F(exponent)
    return tuple(F(values[name]) for name in NAMES)


def full_old_constraints(x, include_v1=False):
    # Scale the finite subset scan to integers; no floating arithmetic.
    require(all(v >= 0 for v in x), 'Negative witness coordinate')
    require(all(dot(a, x) <= b for a, b in base_rows()[:20]), 'Witness fails base constraint')
    den = lcm(*(v.denominator for v in x))
    v = {name: int(value * den) for name, value in zip(NAMES, x)}
    options = []
    for s in 'abc':
        options.append([(sum(v[f'{s}{i}'] for i in range(1, 7) if mask >> (i - 1) & 1), sum((i - 1) * v[f'{s}{i}'] for i in range(1, 7) if mask >> (i - 1) & 1)) for mask in range(64)])
    total = sum(v[s] for s in 'abc')
    for terms in product(*options):
        require(v['D'] <= max(total - sum(t[0] for t in terms), total - den + sum(t[1] for t in terms)), 'Witness fails geometric subset')
    if include_v1:
        for a, b in permutations('abc', 2):
            for divisor in range(2, 7):
                require(v['D'] <= total - sum(v[f'{a}{i}'] + v[f'{b}{i}'] for i in range(1, 7) if i % divisor == 0), 'Witness fails Thue bound')
            for i, j in product(range(1, 7), repeat=2):
                require(v['D'] <= total - v[f'{a}{i}'] - v[f'{b}{j}'] + min(F(v[f'{a}{i}'], j), F(v[f'{b}{j}'], i)), 'Witness fails determinant bound')


def collision_value(x, A, pa, B, pb):
    C, p, q, r, s, frozen, args = collision_expressions(A, pa, B, pb)
    kappa = dot(add(p, q), x) + sum(max(F(0), dot(arg, x)) for arg in args) / 2
    return dot(frozen, x) + max(dot(var(C), x), (dot(var(C), x) + min(2 * dot(r, x), 2 * dot(s, x)) + kappa) / 2)


def verify_witnesses():
    v1 = [(0, F(1, 10), F(1, 10), 0, F(1, 10), 0)] * 2 + [(0, F(1, 5), F(1, 5), 0, 0, 0)]
    full_old_constraints(point(v1, F(3, 5)), include_v1=True)
    lower = (F(1, 66), F(4, 33), F(4, 33), 0, F(5, 66), 0)
    x = point([lower] * 3, F(13, 22))
    full_old_constraints(x, include_v1=True)
    for A, B in permutations('abc', 2):
        for pa, pb in product(((2, 3), (3, 4)), repeat=2):
            require(x[POS['D']] <= collision_value(x, A, pa, B, pb), 'Lower witness fails collision bound')
    print('PASS: a D = 3/5 point satisfies all BLT v2 constraints (so 3/5 is attained without (star)); 13/22 witness satisfies all old and 24 new bounds.')


def compare_models():
    # Builders are ONLY comparison targets. They do not construct proof rows.
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(HERE))
    import lp
    import lp_flagship

    def normalized(constraint):
        require(constraint.sense in (-1, 1), 'Unexpected inequality direction')
        sign = -constraint.sense
        return tuple(sorted((v.name, sign * F(str(c))) for v, c in constraint.items() if c)), sign * F(str(-constraint.constant))

    def sparse(vector, rhs=0):
        return tuple(sorted((name, c) for name, c in zip(NAMES, vector) if c)), F(rhs)

    paper_sets = [((1, 2, 4),) * 3]
    for patterns in [((1, 2), (1, 2), (1, 2, 3)), ((1, 2), (1, 2, 3), (1, 2, 3)), ((1, 2), (1, 2), (1, 3)), ((1, 2), (1, 3), (1, 3))]:
        paper_sets.extend(sorted(set(permutations(patterns))))
    require(len(paper_sets) == 13, 'Wrong Proposition 6.1 subset count')

    for builder, M0, M1 in [(lp.build, 100, 100), (lp_flagship.build_tight, 1, 2)]:
        problem, variables, objective = builder(lam=1, thue=False, det=False)
        actual = [normalized(c) for c in problem.constraints.values()]
        # 10 initial linear rows, 26 geometric rows, then 10 Fourier rows.
        expected_linear = [sparse(a, b) for a, b in base_rows()[:20]]
        require(Counter(actual[:10] + actual[36:]) == Counter(expected_linear), 'Builder differs from independently reconstructed base rows')
        seen = []
        for n in range(13):
            first, second = actual[10 + 2 * n:12 + 2 * n]
            first_terms, rhs0 = dict(first[0]), first[1]
            zname = f'z{n + 1}'
            require(first_terms.pop(zname) == -M0 and rhs0 == 0, 'Wrong geometric first Big-M')
            sets = tuple(tuple(i for i in range(1, 7) if first_terms.get(f'{s}{i}') == 1) for s in 'abc')
            branches = geometric_branches(sets)
            require((tuple(sorted(first_terms.items())), rhs0) == sparse(*branches[0][0]), 'Wrong first geometric branch')
            expected_second = dict(sparse(*branches[1][0])[0])
            expected_second[zname] = F(M1)
            require(second == (tuple(sorted(expected_second.items())), F(M1 - 1)), 'Wrong second geometric branch')
            seen.append(sets)
        require(Counter(seen) == Counter(paper_sets), 'Subsets differ from Proposition 6.1 (6.5)-(6.9)')
        require(all(v.lowBound == 0 for k, v in variables.items()), 'Missing coordinate nonnegativity')
        print(f'PASS: {builder.__module__}.{builder.__name__} matches all Proposition 6.1 base constraints, with the corrected (6.9) subsets.')

    for A, B in permutations('abc', 2):
        for pa, pb in product(((2, 3), (3, 4)), repeat=2):
            problem, variables, objective = lp_flagship.build_tight(lam=1, thue=False, det=False)
            before = len(problem.constraints)
            lp_flagship.add_one_two_letter(problem, variables, objective, A, pa, B, pb)
            actual = [normalized(c) for c in list(problem.constraints.values())[before:]]
            C, p, q, r, s, frozen, args = collision_expressions(A, pa, B, pb)
            tag = f'{A}{pa[0]}{pa[1]}_{B}{pb[0]}{pb[1]}'
            zc = f'fl_cauchy_{tag}'
            expected = []
            hnames = []
            def extend(vector, extras, rhs=0):
                terms = dict(sparse(vector)[0])
                terms.update({k: F(v) for k, v in extras.items()})
                return tuple(sorted((k, v) for k, v in terms.items() if v)), F(rhs)
            for i, arg in enumerate(args):
                h, z = f'fl_pos_{tag}_{i}', f'fl_sign_{tag}_{i}'
                hnames.append(h)
                expected += [extend(arg, {h: -1}), extend(scale(-1, arg), {h: 1, z: 1}, 1), extend(ZERO, {h: 1, z: -1})]
            expected.append(extend(subtract(D, add(frozen, var(C))), {zc: -1}))
            for cap in (r, s):
                lhs = subtract(D, add(frozen, scale(F(1, 2), add(var(C), scale(2, cap), p, q))))
                expected.append(extend(lhs, {zc: 1, **{h: F(-1, 4) for h in hnames}}, 1))
            require(Counter(actual) == Counter(expected), f'Collision encoding mismatch {tag}')
            for v in problem.variables():
                if v.name.startswith(('z', 'fl_sign_', 'fl_cauchy_')):
                    require(v.cat == 'Integer' and v.lowBound == 0 and v.upBound == 1, 'Branch variable not binary')
                if v.name.startswith('fl_pos_'):
                    require(v.lowBound == 0, 'Positive part unbounded below')
    print('PASS: all 24 PuLP collision encodings match the independently derived positive parts and max/min branches.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path, default=HERE / 'cert_flagship.json')
    parser.add_argument('--compare-models', action='store_true')
    args = parser.parse_args()
    verify_tree(read_certificate(args.certificate))
    verify_witnesses()
    if args.compare_models:
        compare_models()
    print('Scope: exact conditional optimization proof; analytic inputs are proved in the manuscript.')


if __name__ == '__main__':
    main()
