# OO Plan B: arXiv:2509.19430v2

This directory contains the reproducible pre-equilibrium-medium construction
used for the paired O16+O16 HYBRID campaign.

## Published Prescription

The implementation follows arXiv:2509.19430v2, Eqs. (2)-(4):

- solve the implicit effective-temperature equation using the QCD kinetic
  energy attractor;
- anchor to the first hydro slice with the published viscous correction and
  `eta/s=0.12`;
- use the three-flavor conformal equation of state;
- scale transverse flow linearly with `tau/tau_hyd` and use Bjorken
  longitudinal flow;
- apply no quenching before `tau_min`.

The required attractor table is
`reference_data/qcd_kinetic_attractor_lambda10_Cinf0p87.tsv`. Its SHA256 is:

```text
1bea7289d3dc8ed95819eaa86cf4c489442a054c14aae47eff010cf45155eba0
```

It is reconstructed from the public QCD EKT data associated with
arXiv:1908.02866, DOI `10.4119/unibi/2939684`. The production runner rejects a
different checksum, a missing table, `eta/s != 0.12`, or a disabled viscous
anchor.

## Validation

From this directory, with a staged hydro tree available at `staged_hydro/`:

```bash
python3 analysis/test_2509_prehydro_public.py -v
```

The test covers the attractor checksum and known values, implicit-equation
residual, interpolation accuracy, onset, flow, conformal EOS, and pinned
production options.

## Reproducibility Boundary

This reproduces every pre-equilibrium-medium operation stated publicly in the
paper. It does not reproduce the authors' complete semi-analytic energy-loss
calculation, and the paper does not publish its exact runtime attractor object.
An author reference output is required before claiming author-code or bitwise
identity.
