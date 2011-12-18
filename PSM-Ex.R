pkgname <- "PSM"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('PSM')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("PSM.estimate")
### * PSM.estimate

flush(stderr()); flush(stdout())

### Name: PSM.estimate
### Title: Estimate population parameters
### Aliases: PSM.estimate
### Keywords: htest models multivariate ts

### ** Examples

cat("\nExamples are included in the package vignette.\n")



cleanEx()
nameEx("PSM.plot")
### * PSM.plot

flush(stderr()); flush(stdout())

### Name: PSM.plot
### Title: Basic plots of data and output
### Aliases: PSM.plot
### Keywords: htest models multivariate ts

### ** Examples

cat("\nExamples are included in the package vignette.\n")



cleanEx()
nameEx("PSM.simulate")
### * PSM.simulate

flush(stderr()); flush(stdout())

### Name: PSM.simulate
### Title: Create simulation data for multiple individuals
### Aliases: PSM.simulate
### Keywords: htest models multivariate ts

### ** Examples

cat("\nExamples are included in the package vignette.\n")



cleanEx()
nameEx("PSM.smooth")
### * PSM.smooth

flush(stderr()); flush(stdout())

### Name: PSM.smooth
### Title: Smoothing of model states based on estimated population
###   parameters.
### Aliases: PSM.smooth
### Keywords: htest models multivariate ts

### ** Examples

cat("\nExamples are included in the package vignette.\n")



cleanEx()
nameEx("PSM.template")
### * PSM.template

flush(stderr()); flush(stdout())

### Name: PSM.template
### Title: Creates a template for a model in PSM
### Aliases: PSM.template
### Keywords: htest models multivariate ts

### ** Examples

# Linear model with input, random effects and dose
PSM.template(Linear=TRUE,dimX=1,dimY=2,dimU=3,dimEta=4)

# Non-linear model without input, random effects and dose
PSM.template(Linear=FALSE,dimX=1,dimY=2,dimU=0,dimEta=0)



cleanEx()
nameEx("matexp")
### * matexp

flush(stderr()); flush(stdout())

### Name: matexp
### Title: Matrix exponential
### Aliases: matexp
### Keywords: math

### ** Examples

##
## The test cases have been taken directly from David Firths MEXP package.
##
##
## ----------------------------
## Test case 1 from Ward (1977)
## ----------------------------
test1 <- t(matrix(c(
    4, 2, 0,
    1, 4, 1,
    1, 1, 4), 3, 3))
matexp(test1)
## Results on Power Mac G3 under Mac OS 10.2.8
##                    [,1]               [,2]               [,3]
## [1,] 147.86662244637000 183.76513864636857  71.79703239999643
## [2,] 127.78108552318250 183.76513864636877  91.88256932318409
## [3,] 127.78108552318204 163.67960172318047 111.96810624637124
## -- these agree with ward (1977, p608)
##
## A naive alternative to mexp, using spectral decomposition:
mexp2 <- function(matrix){
             z <- eigen(matrix,sym=FALSE)
             Re(z$vectors %*% diag(exp(z$values)) %*%
                    solve(z$vectors))
}
try(
mexp2(test1)
) ## now gives an error from solve !
##
## older result was
##                   [,1]                [,2]               [,3]
##[1,] 147.86662244637003  88.500223574029647 103.39983337000028
##[2,] 127.78108552318220 117.345806155250600  90.70416537273444
##[3,] 127.78108552318226  90.384173332156763 117.66579819582827
## -- hopelessly inaccurate in all but the first column.
##
##
## ----------------------------
## Test case 2 from Ward (1977)
## ----------------------------
test2 <- t(matrix(c(
    29.87942128909879, .7815750847907159, -2.289519314033932,
    .7815750847907159, 25.72656945571064,  8.680737820540137,
   -2.289519314033932, 8.680737820540137,  34.39400925519054),
           3, 3))
matexp(test2)
##                   [,1]               [,2]               [,3]
##[1,]   5496313853692357 -18231880972009844 -30475770808580828
##[2,] -18231880972009852  60605228702227024 101291842930256144
##[3,] -30475770808580840 101291842930256144 169294411240859072
## -- which agrees with Ward (1977) to 13 significant figures
mexp2(test2)
##                   [,1]               [,2]               [,3]
##[1,]   5496313853692405 -18231880972009100 -30475770808580196
##[2,] -18231880972009160  60605228702221760 101291842930249376
##[3,] -30475770808580244 101291842930249200 169294411240850880
## -- in this case a very similar degree of accuracy.
##
## ----------------------------
## Test case 3 from Ward (1977)
## ----------------------------
test3 <- t(matrix(c(
    -131, 19, 18,
    -390, 56, 54,
    -387, 57, 52), 3, 3))
matexp(test3)
##                    [,1]                [,2]                [,3]
##[1,] -1.5096441587713636 0.36787943910439874 0.13533528117301735
##[2,] -5.6325707997970271 1.47151775847745725 0.40600584351567010
##[3,] -4.9349383260294299 1.10363831731417195 0.54134112675653534
## -- agrees to 10dp with Ward (1977), p608.
mexp2(test3)
##                   [,1]               [,2]                [,3]
##[1,] -1.509644158796182 0.3678794391103086 0.13533528117547022
##[2,] -5.632570799902948 1.4715177585023838 0.40600584352641989
##[3,] -4.934938326098410 1.1036383173309319 0.54134112676302582
## -- in this case, a similar level of agreement with Ward (1977).
##




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
