Supplement to "Breaking the three-fifths barrier for the exceptional set in the abc conjecture"

Files
  verify_certificate.py  Independent verifier (Python >= 3.9, standard library only).
                         Rebuilds every constraint row from the formulas of Sections 2 and 5,
                         checks tree coverage (2 children per geometry node, 5 per collision node)
                         and every leaf's Farkas multiplier in exact Fraction arithmetic,
                         and checks the D = 13/22 feasible point of Section 7.
  cert_flagship.json     Branch-and-bound tree, 777 leaves, rational multipliers.
  lp.py                  Mixed-integer model of BLT v2 Proposition 6.1 (PuLP), used for discovery only.
  lp_flagship.py         The same model with the 24 constraints (star) added (PuLP), discovery only.

Reproduce
  python3 verify_certificate.py
  python3 verify_certificate.py --compare-models     # optional, needs pulp (tested 3.3.0); no solver call

The proof uses only the first command. The PuLP models are not needed to verify the bound.
