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
data(sampledf)
data(samplemat)
```

The function `createProbMat` creates an object of class `probMat` which standardise the data format and includes additional information such as chronological range and resolution:

```
x.df <- createProbMat(x=sampledf,timeRange=c(6500,4001),resolution=50)
x.mat <- createProbMat(pmat=samplemat,timeRange=c(5000,3001),resolution=100)
```

The standard `plot()` function displays `probMat` class object using aoristic sums:

```
plot(x.df)
plot(x.mat)
```
### Bayesian Inference

##### Estimating Exponential Growth Rate

```
exponential.fit <- expfit(x.mat)
```

```
plot(exponential.fit)
```

```
coda::HPDinterval(mcmc(exponential.fit$posterior.r))
```

##### Non-Parametric Modelling via Random Walk ICAR Model

```
non.param.fit <- icarfit(x.df,niter=2000000,nburnin=1000000,thin=100,nchains=4,parallel=TRUE)
```

```
plot(non.param.fit)
```
