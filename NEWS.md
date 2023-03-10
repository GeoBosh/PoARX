# Changes in Version 0.3.4

- first update since 2018

- fixed the arXiv reference in DESCRIPTION.

- changed remaining "http://" to "https://".

- commented out `\SweaveOpts{concordance=TRUE}` since it now leads to errors
  from latex.

# Changes in Version 0.3.3.9000

- .dlambdajdthetaj() changed to avoid mixing positive and negative subscripts
- .predictLambda() process now correct (0 initialisation not in first row 
     previously)

# Changes in Version 0.3.3

- fitPoARX output changed as rhoFit can't be found in independent case


# Changes in Version 0.3.0.9000

- Renamed NEWS to NEWS.md and edited as markdown file.


# Changes in Version 0.2.9

- data set changed from cranlogs countData to building data.


# Changes in Version 0.2.8

- numerical hessian involving `fullPoARXLogLikelihood()` removed, now the
      variance matrix is made non-singular by increasing the `zeroFix` variable


# Changes in Version 0.2.7

- `NULL` exogenous now fixed
- Frank's copula case variance updated with previous changes


# Changes in Version 0.2.6

- init moved to make more clear that an input is required
- different formulae and different lags now accommodated for


# Changes in Version 0.2.5

- vignette removed due to changes (and speeds up build)
- Length of 'init' warning in `fitPoARX()` changed to include lengths of
  variables
- `warnings` in `PoARXPredict()` caused by multiple logical statements fixed using
  `any()`
- `df.residual.PoARX()` edited for `sameParams = TRUE`
- `summary.PoARX()` edited for `sameParams = TRUE` (same as `df.residual`)
- `fu()` can now deal with `NULL` yinit/muinit
