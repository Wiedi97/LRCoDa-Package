# {LRCoDa-Package}

Linear Regression to combine the use of compositional and non-compositional data in the same model (Planned as extension of the lmCoDaX function from the robComposition package)

> using LRCoDa

data(gemas) 

X <- dplyr::select(gemas, c(MeanTemp, soilclass, Al:Zr))

lr <- LRCoDa(y = gemas$sand, X, external = c('MeanTemp'), factor_column = 'soilclass')

## What is added to the lmCoDaX function from the robComposition Package
- Use linear regression for compositional and non-compositional in the same model
- Include factor variable in the linear regression with compositional data
- Use lmrob for robust regression
- Can specify the method to calculate the pivot coordinates
- Evt allow formula notation

## Getting Started

### Dependencies

> R (>= 2.10)

### Installation
