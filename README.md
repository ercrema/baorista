# baorista
_baorista_ is an R package that provides a Bayesian inferential tool for analysing time-frequencies of archaeological events associated with time spans typically obtained from relative chronological sequences (e.g. periods and phases) and often analyses using aoristic sums. At its core _baorista_ is a frontend for fitting Bayesian models via the [NIMBLE probabilistic programming language](https://r-nimble.org/). The package is currently in beta development and can be installed using `devtools`:

```
library(devtools)
install_github('ercrema/baorista')
```

## Sample Script
### Data Preparation
_baorista_ can read two types of data: a) tables containing two columns defining the start and end date of time-spans of existence recorded in BP; or b) matrices containing probability mass for each event (row) at each time-block (column). Sample datasets are provided.

```
data(sample_df)
data(sample_mat)
```

The function `createProbMat` creates an object of class `probMat` which standardise the data format and includes additional information such as chronological range and resolution:

```
df.x <- createProbMat(x=sample_df,timeRange=c(5000,3001),resolution=100)
df.mat <- createProbMat(pmat=sample_mat,timeRange=c(5000,3001),resolution=50)
```

The standard `plot()` function displays `probMat` class object using aoristic sums:

```
plot(df.x)
plot(df.mat)
```
### Bayesian Inference

##### Estimating Exponential Growth Rate

##### Non-Parametric Modelling via Random Walk ICAR Model
