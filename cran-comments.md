## Resubmission
This is a resubmission.  In this version I have:

* Removed the calls to graphics::par() from the functions plot.choose_uk(), plot.choose_b() and plot.confint_spm().

* Wrapped the examples in the documentation for choose_b() with \donttest{}, rather than \dontrun{}.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- osx (on travis-ci), R-oldrel, R-release            
- ubuntu 12.04 + GCC (on travis-ci), R-release, R-devel
- ubuntu 12.04 + clang (on travis-ci), R-release, R-devel
- win-builder (R-devel and R-release)

## Downstream dependencies

None. This is a new submission
