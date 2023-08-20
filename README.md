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
par(mfrow=c(1,2))
plot(x.df,main='x.df')
plot(x.mat,main='x.mat')
```
![image](https://github.com/ercrema/baorista/blob/main/readme_figs/fig1.png)

### Bayesian Inference

##### Estimating Exponential Growth Rate

Growth rates can be estimated using the function `expfit()`. The function fits an exponential growth model via NIMBLE, and allows user defined MCMC settings (e.g. number of iterations, number of chains, choice of sampler) and choice of the prior. Convergence is automatically checked via Gelman-Rubin statistics (using `coda::gelman.diag()` function) and a warning is issued if the Rhat is higher than 1.01. Default settings is generally sufficient, and the model is fitted in a couple of minutes.

```
exponential.fit <- expfit(x.mat)
```
Posterior parameters of the growth rate can be extracted from the output:

```
#Compute 95% HPDi
coda::HPDinterval(mcmc(exponential.fit$posterior.r))
```
A dedicated plot function can also be used to visualise the fitted exponential model:

```
plot(exponential.fit)
```
![image](https://github.com/ercrema/baorista/blob/main/readme_figs/fig2.png)

##### Non-Parametric Modelling via Random Walk ICAR Model

_baorista_ offers also a non-parametric model based Random Walk ICAR. The model estimates the probability mass for each time-block accounting for the information from adjacent blocks (and hence temporal autocorrelation). The model requires longer chains and users are reccommended to run the model in parallel chains, and be prepared to adjust MCMC settings as well the prior for the parameter sigma (see details in the help documentation). The following example took approximately 3 hours to be completed, but all parameters reached convergence.

```
non.param.fit <- icarfit(x.df,niter=2000000,nburnin=1000000,thin=100,nchains=4,parallel=TRUE)
```

A dedicated function can display the HPDI of the probability mass of each time-block:  

```
plot(non.param.fit)
```
![image](https://github.com/ercrema/baorista/blob/main/readme_figs/fig3.png)

