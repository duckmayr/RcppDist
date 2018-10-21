## Resubmission
This is a resubmission. In this version I have:
* Added a function, bayeslm(), that demonstrates the use of RcppDist functions
  in C++ code, and has R-runnable examples, as well as displays the C++ code
  used for the function in the help file.
* In the package documentation file, push users toward the vignette more,
  and explain what kind of information they'll find there. As the package
  provides C++ headers with functions for users extending R with C++ code,
  it is perhaps more helpful to provide a more detailed explanation of the
  package features in the vignette.

## Test environments
* local Ubuntu 18.04 install, R 3.5.1
* Ubuntu 14.04 (on Travis CI), R 3.4.4, 3.5.1, and devel
* macOS High Sierra 10.13.3 (on Travis CI), R 3.4.4, 3.5.0, and devel
* Windows Server 2012 R2 (on AppVeyor) R 3.5.1 Patched (2018-10-20 r75474)
* win-builder (R-release and R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* The only note was the "New submission" note.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
