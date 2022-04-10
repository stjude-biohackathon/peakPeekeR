# peakPeekeR

## Installation

```R
devtools::install_github("stjude-biohackathon/peakPeakeR")
```

## Test

This will install the peak caller environments the first time run and print their versions (SICER2 doesn't have a version command). It will take several minutes to run the first time, but will be very quick thereafter, as the environments are retained.

```R
library("peakPeekeR")
test()
```