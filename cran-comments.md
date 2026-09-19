Resubmission

This is a resubmission addressing both points raised by Leonore Hochhauser
on 2026-08-27.

1. Missing \value tag.

   fitted.ngvb.Rd now documents its return value: the class (a data.frame),
   what a row is (one element of the linear predictor), what the columns
   hold (the posterior summaries INLA computed), and the NULL case. It was
   the only exported object missing \value; the rest of the reference index
   was rechecked.

2. Modifying the global environment.

   Two separate places were changed.

   - R/detect.R used <<- in two internal formula walkers. Both now keep
     their state in an explicit local environment.

   - R/api.R was the substantive one: ngvb() assigned its engine objects
     into the user's workspace with assign(..., envir = globalenv()), and
     removed them again on exit. INLA resolves a formula's model = argument
     in its .parent.frame argument, and the argument list recovered from a
     fitted INLA object carries the .parent.frame captured at the user's
     original inla() call, which is .GlobalEnv. Overriding that to a private
     environment whose parent is the global environment makes the lookup
     resolve there instead, so the objects never enter the workspace while
     everything else the formula references still resolves. The two on.exit()
     cleanups this required are gone.

   Nothing in the package writes to .GlobalEnv now.

   Verified against a real INLA installation: ngvb(), ngvb_sample() and
   bayes.factor() all run, and ls(globalenv()) is unchanged across a full
   ngvb() call. The test suite passes with INLA present and the INLA-backed
   tests actually running: 108 passed, 0 failed, 0 skipped.

R CMD check results

0 errors | 0 warnings | 1 note

- This is a new release.

Notes

- The package builds on INLA, which is not on CRAN. It is listed in Suggests
  with Additional_repositories: https://inla.r-inla-download.org/R/stable,
  and every entry point checks for it with requireNamespace("INLA") and fails
  with an actionable install message otherwise.

- Nothing in the submitted tarball runs INLA, so a missing or broken INLA
  installation on the check machine cannot affect the results:

  - The shipped vignette (ngvb_pkg.Rmd) contains no evaluated code chunks.
    It is an overview that links to the worked examples hosted on the package
    website (https://rafaelcabral-ai.github.io/ngvb/). The sources live in
    vignettes/articles/ in the GitHub repository and are excluded from the
    build via .Rbuildignore.
  - The examples that call INLA are wrapped in \dontrun{} rather than
    \donttest{}, since they cannot be executed without a package that is
    not on CRAN.
  - Every test that needs INLA calls a skip_if_no_inla() helper, which calls
    skip_on_cran() before it even looks for the package. INLA is never
    attached by the test suite.

Test environments

- local Windows 11, R 4.4.2: 0 errors, 0 warnings, 1 note.
- local Windows 11, R 4.4.2, with INLA present and the INLA-backed tests
  actually running (NOT_CRAN=true): 108 passed, 0 failed, 0 skipped.
- win-builder R-devel. INLA is not installed there, so this also exercises
  the "INLA absent" path: examples, tests and vignette all pass without it.
- GitHub Actions: ubuntu-latest (release, devel), macOS-latest (release),
  windows-latest (release)
