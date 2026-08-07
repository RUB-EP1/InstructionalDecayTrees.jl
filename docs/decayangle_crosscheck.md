# DecayAngle Cross-Check Reproducibility

This document records exactly how `test/fixtures/decayangle_crosscheck.json` was produced.

## Source Provenance

- Python package: `decayangle`
- Repository: `https://github.com/KaiHabermann/decayangle`
- Checked-out commit:
  - `6fbeb1ec090373a13434506e2c555a40c0b849d5`
  - commit date: `2025-10-24 09:02:35 +0200`
  - message: `Better more controlled testing method`
- Declared package version in `pyproject.toml`: `1.1.3`

## Python Environment Used

- Python: `3.14.0`
- Virtual environment path: `/tmp/venvb`
- Installed packages used for generation:
  - `jax==0.9.0.1`
  - `jaxlib==0.9.0.1`
  - `ml_dtypes==0.5.4`
  - `networkx==3.6.1`
  - `numpy==2.4.2`
  - `opt_einsum==3.4.0`
  - `scipy==1.17.1`
  - `tqdm==4.67.3`

## Generation Procedure

1. Clone `decayangle`:

```bash
cd /tmp
git clone https://github.com/KaiHabermann/decayangle.git
cd /tmp/decayangle
git checkout --detach 6fbeb1ec090373a13434506e2c555a40c0b849d5
git rev-parse --verify HEAD
```

2. Create and populate virtual environment:

```bash
python3.14 -m venv /tmp/venvb
/tmp/venvb/bin/python -m pip install \
    jax==0.9.0.1 \
    jaxlib==0.9.0.1 \
    ml_dtypes==0.5.4 \
    networkx==3.6.1 \
    numpy==2.4.2 \
    opt_einsum==3.4.0 \
    scipy==1.17.1 \
    tqdm==4.67.3
```

3. Generate fixture from this repository:

```bash
cd /path/to/your/InstructionalDecayTrees.jl/checkout
DECAYANGLE_SRC=/tmp/decayangle/src /tmp/venvb/bin/python test/generate_decayangle_fixture.py
```

This writes:

- `test/fixtures/decayangle_crosscheck.json`

Verify that the regenerated fixture is byte-for-byte identical:

```bash
shasum -a 256 test/fixtures/decayangle_crosscheck.json
```

Expected SHA-256:

```text
b5b41f02fb812e7a682071e00c3288a7eee41c12bb2460f934c05613b25522c6
```

## Fixture Content Notes

- Convention: `helicity`
- Four-vectors are stored as `(px, py, pz, E)` arrays.
- Relative transform follows:
  - `relative = other @ inverse(reference)` on the Python side.
- Cases currently included:
  - `four_body` (seed `12345`)
  - `five_body` (seed `67890`)

## Julia Validation

The generated JSON is validated by:

- `test/crosscheck_json.jl`

Run with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

This test reconstructs Julia instruction paths from the JSON and checks:
- relative matrix agreement,
- per-path Lorentz decode agreement,
- relative Wigner angle agreement.
