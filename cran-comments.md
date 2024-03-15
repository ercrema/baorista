## R CMD check results

0 errors | 0 warnings | 3 note

## NOTE 1
* checking CRAN incoming feasibility ... [4s/11s] NOTE
Maintainer: ‘Enrico Crema <enrico.crema@gmail.com>’

#### RESPONSE
This is a new submission

## NOTE 2

* checking R code for possible problems ... NOTE
expfit : runfun : <anonymous>: no visible global function definition
  for ‘returnType’
expfit : runfun: no visible binding for global variable ‘rAExp’
icarfit : runfun : <anonymous>: no visible global function definition
  for ‘returnType’
icarfit : runfun: no visible binding for global variable ‘rAOG’
logisticfit : runfun : <anonymous>: no visible global function
  definition for ‘returnType’
logisticfit : runfun: no visible binding for global variable ‘rALog’
rAoristicExponentialGrowth_vector: no visible global function
  definition for ‘nimStop’
rAoristicExponentialGrowth_vector: no visible global function
  definition for ‘nimMatrix’
rAoristicGeneral_vector: no visible global function definition for
  ‘nimStop’
rAoristicGeneral_vector: no visible global function definition for
  ‘nimMatrix’
rAoristicLogisticGrowth_vector: no visible global function definition
  for ‘nimStop’
rAoristicLogisticGrowth_vector: no visible global function definition
  for ‘nimMatrix’
Undefined global functions or variables:
  nimMatrix nimStop rAExp rALog rAOG returnType

#### RESPONSES
All instances of 'no visible global function' are associated with the functions `expfit()`, `logisticfit()`, and `icarfit()`. 
Each of these functions consists of generating and fitting (via MCMC) a custom Bayesian model with the nimble package. 
The submitted package effectively creates several of these custom models as statistical distributions (e.g. `dAoristicExponentialGrowth_vector`); the `nimbleModel()` function generates automatically relevant random generator functions (e.g. in this case `rAoristicExponentialGrowth_vector`), which includes a number of sub-functions (e.g. `nimStop`,`nimMatrix`). 
Given the computational cost for the MCMC model fitting, I have also written a subroutine involving the `parLapply()` function. This required the definition of a function (`runfun()`) inside `expfit()`, `logisticfit()`, and `icarfit()`.




