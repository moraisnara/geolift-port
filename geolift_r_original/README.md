# geolift_r_original — vendored GeoLift (trimmed)

This is an **unmodified, trimmed copy** of Meta's
[GeoLift](https://github.com/facebookincubator/GeoLift) R package, **v2.7.5**,
included here as the *original* baseline engine for the port study in
[`../REPORT.md`](../REPORT.md).

It is the upstream source code verbatim — no functions were changed. Only
non-runtime material was removed to keep the repo small. **Kept:** `R/`,
`DESCRIPTION`, `NAMESPACE`, `data/`, `LICENSE.md`. **Removed:** `website/`,
`vignettes/`, `Whitepapers/`, `doc/`, `man/`, `.github/`, the upstream `.git`
history, and community-health files.

The `geolift_r_modified/` folder at the repo root layers parameter presets on
top of this unchanged package (it does **not** fork the code).

## Install

```r
# augsynth must be pinned to the pre-PR#88 commit (see ../REPORT.md):
remotes::install_github("ebenmichael/augsynth", ref = "06415b4")
# then install this trimmed package from the local folder:
remotes::install_local("geolift_r_original")
```

Licensed under the MIT License — see [`LICENSE.md`](LICENSE.md).
Copyright (c) Meta Platforms, Inc. and its affiliates.
