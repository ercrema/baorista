# NOTES

> The Description field is intended to be a (one paragraph) description of what the package does and why it may be useful. Please add more details about the package functionality and implemented methods in your Description text.


Thanks, I have updated the DESCRIPTION file of the package with additional details about its functionality. 

> Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar) -> Missing Rd-tags:
>     plot.fittedExp.Rd: \value
>     plot.fittedICAR.Rd: \value
>     plot.fittedLogistic.Rd: \value
>     plot.probMat.Rd: \value

Thanks, given these are plot functions with no return values I have added the following: \value{No return value (plot function)}.


> \dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user. Does not seem necessary. Please replace \dontrun with \donttest.
> Please unwrap the examples if they are executable in < 5 sec, or replace dontrun{} with \donttest{}.

These have now been removed.

> You write information messages to the console that cannot be easily suppressed.
> It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object. Instead of print()/cat() rather use message()/warning() or if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console. (except for print, summary, interactive functions) -> R/expfit.R, R/icarfit.R, R/logisticfit.R

I have replaced `print()` with `message()`

> Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited.
> e.g.:
> ...
> oldpar <- par(no.readonly = TRUE) # code line i
> on.exit(par(oldpar)) # code line i + 1
> ...
> par(mfrow=c(2,2)) # somewhere after
> ...
> If you're not familiar with the function, please check ?on.exit. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.
> ->  R/plotfittedExp.R, R/plotfittedICAR.R , R/plotfittedLogistic.R

Thank you, I have followed the suggested workaround for ensuring that the par settings are brought restored on exit.

> Please do not modifiy the .GlobalEnv. This is not allowed by the CRAN policies. -> R/expfit.R, R/icarfit.R, R/logisticfit.R 

Several objects are created in the .GlobalEnv by functions derived by the Nimble package and cannot be avoided. To restore the .GlobalEnv to the condition prior to the execution of the function the following lines were added:

```
#Handle cleaning of GlobalEnv on exit
envobj <- ls(envir=.GlobalEnv)
on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))
```

this ensures that the .GlobalEnv is not modified by these functions.

## Test environments
* local Debian 12
* rhub check_for_cran

## R CMD check results

0 errors | 0 warnings | 1 notes

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Enrico Crema  <enrico.crema@gmail..com>'
