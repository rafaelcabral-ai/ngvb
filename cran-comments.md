## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes

* The package builds on 'INLA', which is not on CRAN. It is listed in `Suggests`
  with the standard `Additional_repositories: https://inla.r-inla-download.org/R/stable`,
  every entry point checks for it with `requireNamespace("INLA")` and fails with an
  actionable install message otherwise, and all examples, tests and vignette code are
  conditional on INLA being available. The package therefore installs and checks
  cleanly without INLA.

## Test environments

* local Windows 11, R 4.4.2
* GitHub Actions: ubuntu-latest (release, devel), macOS-latest (release),
  windows-latest (release)
