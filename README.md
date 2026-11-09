# Breaking the three-fifths barrier for the exceptional set in the abc conjecture

[![verify](https://github.com/Li-Hongmin/abc-exceptional-set-three-fifths/actions/workflows/verify.yml/badge.svg)](https://github.com/Li-Hongmin/abc-exceptional-set-three-fifths/actions/workflows/verify.yml) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24054868-blue)](https://doi.org/10.5281/zenodo.24054868)

**Hongmin Li** (Institute of Science Tokyo; The University of Tokyo) · ORCID [0000-0003-0228-0600](https://orcid.org/0000-0003-0228-0600)

Main result of the manuscript: the number of abc triples a + b = c ≤ X with rad(abc) < c^λ is ≪ X^{59969/100033+ε} for all 0 < λ ≤ 1, improving the exponent 3/5 of Bernert–Browning–Lichtman–Teräväinen; the linear-programming bound is certified by an exact rational Farkas certificate.

- Manuscript on rStation: <https://rstation.ai/problems/math/papers/6/v3> (version 3)
- Archived on Zenodo: concept DOI [10.5281/zenodo.24054868](https://doi.org/10.5281/zenodo.24054868) (all versions); this repository mirrors version 3, DOI [10.5281/zenodo.24054892](https://doi.org/10.5281/zenodo.24054892)

## Contents

- `paper/paper.md`, `paper/abc-exceptional-set-three-fifths.pdf` — the manuscript (version 3).
- `verification/` — the supplement of that version, unpacked unchanged; its own README describes each file.
- `.github/workflows/verify.yml` — reruns every check below on each push.

## What the programs check

| Check | Command (run inside `verification/`) |
|---|---|
| Branch-and-bound certificate (777 leaves, exact Farkas multipliers) and the D = 13/22 point | `python verify_certificate.py` |
| Model comparison (no solver call) | `python verify_certificate.py --compare-models` |

The programs check the finite and exact computations used in the manuscript; the mathematical arguments are those written in the paper.

## Reproducing

```sh
cd verification
pip install pulp==3.3.0
python verify_certificate.py
python verify_certificate.py --compare-models
```

Requirements: Python 3.9+. Run the commands from a copy if you want to keep the recorded outputs untouched.

## Cite

```bibtex
@misc{li2026abcexceptionalsetthreefi,
  author    = {Li, Hongmin},
  title     = {Breaking the three-fifths barrier for the exceptional set in the abc conjecture},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24054868},
  url       = {https://doi.org/10.5281/zenodo.24054868}
}
```

## License

Manuscript and data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: MIT (see `LICENSE`).
