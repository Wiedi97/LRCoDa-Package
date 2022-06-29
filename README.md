# {LRCoDa-Package}

Linear Regression to combine the use of compositional and non-compositional data in the same model (Planned as extension of the lmCoDaX function from the robComposition package)

```
using LRCoDa
```

data(gemas) 

X <- dplyr::select(gemas, c(MeanTemp, soilclass, Al:Zr))

lr <- LRCoDa(y = gemas$sand, X, external = c('MeanTemp', 'soilclass'))

## What is newly released regarding to the lmCoDaX function from the robComposition Package
- Use linear regression for compositional and non-compositional in the same model
- Include factor variable in the linear regression with compositional data
- Use lmrob for robust regression
- Can specify the norm to calculate the pivot coordinates
- Prefix Notation of the variables available in the output

## Goals
- Extend the possibilities of using the linear regression for compositional data.

## Getting Started

### Dependencies

The package has dependencies on

```
R (>= 2.10), dplyr, robCompositions, robustbase, sjmisc
```

### Installation

Installion of `LRCoDa` is easy when the R-tools are installed. Just use

```
library(devtools)
install_github("Wiedi97/LRCoDa-Package")
```
