# QCD Kinetic Energy Attractor

The adjacent `qcd_kinetic_attractor_lambda10_Cinf0p87.tsv` is extracted from
the public data release for:

- G. Giacalone, A. Mazeliauskas, and S. Schlichting,
  arXiv:1908.02866.
- Data DOI: `10.4119/unibi/2939684`.
- License: CC-BY 4.0.

Pablos and Takacs cite this work when defining the hydrodynamic-attractor
extrapolation in arXiv:2509.19430v2, Eqs. (2)-(4). The data release's
`data/figure1.py`
constructs the `QCD kinetics`, `C_infinity=0.87` curve from the three-flavor
QCD EKT gluon and fermion files at `lambda=10`. The extraction uses the same
effective degrees of freedom, fixed `eta/s`, late-time normalization fit, and
column definitions as that script.

Checksums:

```text
e8130755e14a9dee30739a6d6d5b0612c02cc425262bfa63873cccd43d93c655  published-data.zip
1bea7289d3dc8ed95819eaa86cf4c489442a054c14aae47eff010cf45155eba0  qcd_kinetic_attractor_lambda10_Cinf0p87.tsv
```

The table has 200 points over
`0.034830472 <= omega <= 6.2900796`. The default OO Plan B interval lies within
this range.

Regenerate it from an extracted copy of the DOI archive with:

```bash
python3 analysis/extract_qcd10_attractor.py \
  --data-root /path/to/published-data/data/QCD \
  --output reference_data/qcd_kinetic_attractor_lambda10_Cinf0p87.tsv
```

This table, together with the implicit temperature solve, the viscous anchor
at `tau_hyd`, `eta/s=0.12`, linear transverse pre-flow, Bjorken longitudinal
flow, and the `tau_min` onset, reproduces every pre-equilibrium-medium detail
specified publicly in arXiv:2509.19430v2. The paper does not publish its
runtime attractor object or explicitly state that the `lambda=10` Figure 1
data file was loaded in its production. Author confirmation or a reference
output remains necessary before claiming bitwise or author-code identity.
