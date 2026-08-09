## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes

* The package builds on 'INLA', which is not on CRAN. It is listed in `Suggests`
  with the standard `Additional_repositories: https://inla.r-inla-download.org/R/stable`,
  and every entry point checks for it with `requireNamespace("INLA")` and fails with an
  actionable install message otherwise.

* Nothing in the submitted tarball runs INLA, so a missing or broken INLA
  installation on the check machine cannot affect the results:

  - The shipped vignette (`ngvb_pkg.Rmd`) contains no evaluated code chunks. It is an
    overview that links to the worked examples, which are hosted on the package website
    <https://rafaelcabral-ai.github.io/ngvb/>. The sources for those live in
    `vignettes/articles/` in the GitHub repository and are excluded from the build via
    `.Rbuildignore`.
  - The examples that call INLA are wrapped in `\dontrun{}` rather than `\donttest{}`,
    since they cannot be executed without a package that is not on CRAN.
  - Every test that needs INLA calls a `skip_if_no_inla()` helper, which calls
    `skip_on_cran()` before it even looks for the package. INLA is never attached by the
    test suite.

## Test environments

* local Windows 11, R 4.4.2
* GitHub Actions: ubuntu-latest (release, devel), macOS-latest (release),
  windows-latest (release)
