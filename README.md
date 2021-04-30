coxmeg 1.0.13
=============

Overview
--------

Time-to-event is one of the most important phenotypes in genetic
epidemiology. The R-package, “coxmeg”, provides a set of utilities to
fit a Cox mixed-effects model and to efficiently perform genome-wide
association analysis of time-to-event phenotypes using a Cox
mixed-effects model. More details can be found in (He and Kulminski
2020).

Installation
------------

### Most-recent version

    install.packages("devtools")
    library(devtools)
    install_github("lhe17/coxmeg")

Functions
---------

The current version provides five functions.

-   `coxmeg`: Fit a Cox mixed-effects model.
-   `coxmeg_m`: Perform a GWAS using a genotype matrix.
-   `coxmeg_plink`: Perform a GWAS using plink files.
-   `coxmeg_gds`: Perform a GWAS using a GDS file. Read more details
    [here](../doc/coxmeg_gds_example.html).
-   `fit_ppl`: Estimate hazard ratios (HRs) given a variance component.

Fit a Cox mixed-effects model with a sparse relatedness matrix
--------------------------------------------------------------

We illustrate how to use coxmeg to fit a Cox mixed-effects model with a
sparse relatedness matrix. We first simulate a block-diagonal
relatedness matrix for a cohort consisting of 200 families, each of
which has five members.

    library(coxmeg)

    ## Loading required package: Rcpp

    ## Warning: package 'Rcpp' was built under R version 3.6.2

    library(MASS)

    ## Warning: package 'MASS' was built under R version 3.6.3

    library(Matrix)

    ## Warning: package 'Matrix' was built under R version 3.6.3

    n_f <- 200
    mat_list <- list()
    size <- rep(5,n_f)
    offd <- 0.5
    for(i in 1:n_f)
    {
      mat_list[[i]] <- matrix(offd,size[i],size[i])
      diag(mat_list[[i]]) <- 1
    }
    sigma <- as.matrix(bdiag(mat_list))
    sigma = as(sigma,'dgCMatrix')

We use ‘dgCMatrix’ to save memory. Next, we simulate random effects and
time-to-event outcomes assuming a constant baseline hazard function. We
assume that the variance component is 0.2. We also simulate a risk
factor with log(HR)=0.1.

    n = nrow(sigma)
    tau_var <- 0.2
    x <- mvrnorm(1, rep(0,n), tau_var*sigma)
    pred = rnorm(n,0,1)
    myrates <- exp(x+0.1*pred-1)
    y <- rexp(n, rate = myrates)
    cen <- rexp(n, rate = 0.02 )
    ycen <- pmin(y, cen)
    outcome <- cbind(ycen,as.numeric(y <= cen))
    head(outcome)

    ##           ycen  
    ## [1,] 3.9563596 1
    ## [2,] 7.1771349 1
    ## [3,] 7.5607081 1
    ## [4,] 0.5005451 1
    ## [5,] 4.0746380 1
    ## [6,] 0.5177455 1

    sigma[1:5,1:5]

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##                         
    ## [1,] 1.0 0.5 0.5 0.5 0.5
    ## [2,] 0.5 1.0 0.5 0.5 0.5
    ## [3,] 0.5 0.5 1.0 0.5 0.5
    ## [4,] 0.5 0.5 0.5 1.0 0.5
    ## [5,] 0.5 0.5 0.5 0.5 1.0

We fit a Cox mixed-effects model using the `coxmeg` function. Here, we
set `type='bd'` because the relatedness matrix is a block-diagonal
matrix. Note that `type='bd'` should be used only for a block-diagonal
matrix or a sparse matrix of which its inverse matrix is also highly
sparse. A sparse kinship matrix can be converted to a block-diagonal
matrix using [kingToMatrix](../doc/coxmeg_gds_example.html). For a
general sparse relatedness matrix of which its inverse is not sparse, it
is recommended that `type='sparse'` be used. When `type='sparse'` is
specified, the relatedness matrix will not be inverted during the
estimation procedure. The function will automatically treat the
relatedness matrix as dense if there are more than 50% non-zero elements
in the matrix. The argument `X` is a design matrix of the predictors,
which can be produced by e.g., the function `model.matrix`. The design
matrix for the Cox model does not include the intercept term. The
columns in `X` should be linear independent; otherwise the function will
stop with an error indicating sigularity.

    re = coxmeg(outcome,sigma,type='bd',X=pred,order=1,detap='diagonal')

    ## Remove 0 subjects censored before the first failure.

    ## There is/are 1 predictors. The sample size included is 1000.

    ## The relatedness matrix is treated as sparse.

    ## The relatedness matrix is inverted.

    ## The method for computing the determinant is 'diagonal'.

    ## Solver: Cholesky decomposition (RcppEigen=TRUE).

    re

    ## $beta
    ## [1] 0.02482134
    ## 
    ## $HR
    ## [1] 1.025132
    ## 
    ## $sd_beta
    ## [1] 0.03666048
    ## 
    ## $p
    ## [1] 0.4983679
    ## 
    ## $tau
    ## [1] 0.1731624
    ## 
    ## $iter
    ## [1] 16
    ## 
    ## $rank
    ## [1] 1000
    ## 
    ## $nsam
    ## [1] 1000
    ## 
    ## $int_ll
    ## [1] 11574.44

In the above result, `tau` is the estimated variance component, and
`int_ll` is -2\*log(lik) of the integrated/marginal likelihood of tau.

We give more details about specifying `order` and `detap`. We set
`order=1` (also by default) to use the first-order approximation of the
inverse Hessian matrix in the optimization. By `detap='diagonal'`, we
tell `coxmeg` to use a diagonal approximation to compute the
determinant, which is much faster under this setting, when estimating
the variance component. By default (`detap='NULL'`), `coxmeg` will
automatically select a method for computing the determinant based on
`type`, the sample size, and whether the relatedness matrix is symmetric
positive definite (SPD).

It should be noted that when the relatedness matrix is SPD, `coxmeg`
will make use of the sparsity by setting `type='sparse'` or `type='bd'`
regardless of whether the relatedness matrix or its inverse is sparse.
However, when the relatedness matrix is symmetric positive semidefinite
(SPSD), `coxmeg` can make use of the sparsity only when its inverse is
sparse. When the relatedness matrix is SPSD and its inverse is dense,
setting `type='sparse'` may result in worse performance. In such a case,
it would be better to use `type='dense'` or to convert the relatedness
matrix to SPD or block-diagonal if possible.

We compare the results with coxme, which are slightly different due to
different approximation of the log-determinant used in the estimation of
the variance component. Also, the integrated log-likelihoods cannot be
compared directly because different approximation of log-determinant is
used.

    library(coxme)

    ## Warning: package 'coxme' was built under R version 3.6.2

    ## Loading required package: survival

    ## Warning: package 'survival' was built under R version 3.6.3

    ## Loading required package: bdsmatrix

    ## Warning: package 'bdsmatrix' was built under R version 3.6.2

    ## 
    ## Attaching package: 'bdsmatrix'

    ## The following object is masked from 'package:base':
    ## 
    ##     backsolve

    bls <- c(1)
    for(i in (size[1]-1):1)
    {bls <- c(bls, c(rep(offd,i),1))}
    tmat <- bdsmatrix(blocksize=size, blocks=rep(bls,n_f),dimnames=list(as.character(1:n),as.character(1:n)))
    re_coxme = coxme(Surv(outcome[,1],outcome[,2])~as.matrix(pred)+(1|as.character(1:n)), varlist=list(tmat),ties='breslow')
    re_coxme

    ## Cox mixed-effects model fit by maximum likelihood
    ## 
    ##   events, n = 949, 1000
    ##   Iterations= 22 94 
    ##                     NULL Integrated   Fitted
    ## Log-likelihood -5624.431  -5619.731 -5491.11
    ## 
    ##                    Chisq     df          p   AIC     BIC
    ## Integrated loglik   9.40   2.00 9.0903e-03  5.40   -4.31
    ##  Penalized loglik 266.64 119.16 2.6634e-13 28.32 -550.25
    ## 
    ## Model:  Surv(outcome[, 1], outcome[, 2]) ~ as.matrix(pred) + (1 | as.character(1:n)) 
    ## Fixed coefficients
    ##                      coef exp(coef)   se(coef)    z   p
    ## as.matrix(pred) 0.0248857  1.025198 0.03671398 0.68 0.5
    ## 
    ## Random effects
    ##  Group             Variable Std Dev   Variance 
    ##  as.character.1.n. Vmat.1   0.4183460 0.1750134

In GWAS, we may split the procedure into two separate steps, (1)
estimate the variance component under the null model, and (2) estimate
the coefficients for the predictors using the estimated variance
component. This can be carried out in the following way.

    re = coxmeg(outcome,sigma,type='bd',order=1,detap='diagonal')

    ## Remove 0 subjects censored before the first failure.

    ## There is/are 0 predictors. The sample size included is 1000.

    ## The relatedness matrix is treated as sparse.

    ## The relatedness matrix is inverted.

    ## The method for computing the determinant is 'diagonal'.

    ## Solver: Cholesky decomposition (RcppEigen=TRUE).

    tau = re$tau
    print(tau)

    ## [1] 0.1704599

    re2 = fit_ppl(pred,outcome,sigma,type='bd',tau=tau,order=1)

    ## Remove 0 subjects censored before the first failure.

    ## The sample size included is 1000.

    ## The relatedness matrix is treated as sparse.

    ## The relatedness matrix is inverted.

    ## Solver: Cholesky decomposition (RcppEigen=TRUE).

    re2

    ## $beta
    ## [1] 0.02472851
    ## 
    ## $HR
    ## [1] 1.025037
    ## 
    ## $sd_beta
    ## [1] 0.03661347
    ## 
    ## $p
    ## [1] 0.4994256
    ## 
    ## $iter
    ## [1] 4
    ## 
    ## $ppl
    ##           [,1]
    ## [1,] -5552.299

Perform GWAS of an age-at-onset phenotype with a sparse relatedness matrix
--------------------------------------------------------------------------

We illustrate how to perform a GWAS using the `coxmeg_plink` function.
This function supports plink bed files. We provide example files in the
package. The example plink files include 20 SNPs and 3000 subjects from
600 families. The following code performs a GWAS for all SNPs in the
example bed files. The `coxmeg_plink` function will write a temporary
.gds file for the SNPs in the folder specified by `tmp_dir`. The user
needs to specify a `tmp_dir` to store the temporary file when `bed` is
provided. The temporary file is removed after the analysis is done.

    library(coxmeg)
    bed = system.file("extdata", "example_null.bed", package = "coxmeg")
    bed = substr(bed,1,nchar(bed)-4)
    pheno = system.file("extdata", "ex_pheno.txt", package = "coxmeg")
    cov = system.file("extdata", "ex_cov.txt", package = "coxmeg")

    ## building a relatedness matrix
    n_f <- 600
    mat_list <- list()
    size <- rep(5,n_f)
    offd <- 0.5
    for(i in 1:n_f)
    {
      mat_list[[i]] <- matrix(offd,size[i],size[i])
      diag(mat_list[[i]]) <- 1
    }
    sigma <- as.matrix(bdiag(mat_list))

    re = coxmeg_plink(pheno,sigma,type='bd',bed=bed,tmp_dir=tempdir(),cov_file=cov,verbose=FALSE)

    ## Excluding 0 SNP on non-autosomes
    ## Excluding 0 SNP (monomorphic: TRUE, MAF: 0.05, missing rate: 0)

    re

    ## $summary
    ##     snp.id chromosome position allele      afreq  afreq_inc   index
    ## 1   null_0          1        1    d/D 0.30983333 0.30983333  null_0
    ## 2   null_1          1        2    d/D 0.23466667 0.23466667  null_1
    ## 3   null_2          1        3    D/d 0.14033333 0.14033333  null_2
    ## 4   null_3          1        4    D/d 0.16183333 0.16183333  null_3
    ## 5   null_4          1        5    d/D 0.19933333 0.19933333  null_4
    ## 6   null_5          1        6    D/d 0.11800000 0.11800000  null_5
    ## 7   null_6          1        7    d/D 0.09483333 0.09483333  null_6
    ## 8   null_7          1        8    D/d 0.49683333 0.49683333  null_7
    ## 9   null_8          1        9    d/D 0.31366667 0.31366667  null_8
    ## 10  null_9          1       10    D/d 0.49183333 0.49183333  null_9
    ## 11 null_10          1       11    d/D 0.34833333 0.34833333 null_10
    ## 12 null_11          1       12    D/d 0.25100000 0.25100000 null_11
    ## 13 null_12          1       13    d/D 0.17500000 0.17500000 null_12
    ## 14 null_13          1       14    D/d 0.06333333 0.06333333 null_13
    ## 15 null_14          1       15    D/d 0.20833333 0.20833333 null_14
    ## 16 null_15          1       16    d/D 0.17050000 0.17050000 null_15
    ## 17 null_16          1       17    D/d 0.33550000 0.33550000 null_16
    ## 18 null_17          1       18    d/D 0.26633333 0.26633333 null_17
    ## 19 null_18          1       19    D/d 0.09433333 0.09433333 null_18
    ## 20 null_19          1       20    d/D 0.11650000 0.11650000 null_19
    ##            beta        HR    sd_beta           p
    ## 1   0.015672101 1.0157956 0.02938524 0.593803537
    ## 2   0.019439150 1.0196293 0.03222054 0.546298835
    ## 3  -0.049845757 0.9513762 0.03860368 0.196628160
    ## 4   0.044130767 1.0451190 0.03701019 0.233106387
    ## 5   0.028473176 1.0288824 0.03432500 0.406811816
    ## 6  -0.114319159 0.8919732 0.04234095 0.006934636
    ## 7  -0.017981231 0.9821795 0.04655562 0.699325464
    ## 8  -0.004207897 0.9958009 0.02717805 0.876957699
    ## 9  -0.063741849 0.9382472 0.02958441 0.031195036
    ## 10 -0.008409562 0.9916257 0.02730686 0.758108827
    ## 11 -0.013581479 0.9865103 0.02859980 0.634872392
    ## 12  0.037508301 1.0382206 0.03113254 0.228282858
    ## 13 -0.017215848 0.9829315 0.03628637 0.635183349
    ## 14 -0.068207724 0.9340664 0.05698849 0.231357835
    ## 15 -0.013965386 0.9861317 0.03431600 0.684034201
    ## 16  0.002172773 1.0021751 0.03685682 0.952990554
    ## 17  0.004762350 1.0047737 0.02859957 0.867749134
    ## 18  0.001786995 1.0017886 0.03098518 0.954009439
    ## 19 -0.016052310 0.9840758 0.04731969 0.734435643
    ## 20 -0.022398126 0.9778508 0.04231689 0.596600710
    ## 
    ## $tau
    ## [1] 0.04028041
    ## 
    ## $rank
    ## [1] 3000
    ## 
    ## $nsam
    ## [1] 3000

The above code first retrieves the full path of the files. If the full
path is not given, `coxmeg_plink` will search the current working
directory. The file name of the bed file should not include the suffix
(.bed). The phenotype and covariate files have the same format as used
in plink, and the IDs must be consistent with the bed files.
Specifically, the phenotype file should include four columns including
family ID, individual ID, time, and status. The covariate file always
starts with two columns, family ID and individual ID. Missing values in
the phenotype and covariate files are denoted by -9 and NA,
respectively. In the current version, the `coxmeg_plink` function does
not impute genotypes itself, and only SNPs without missing values will
be analyzed, so it will be better to use imputed genotype data.

The `coxmeg_plink` function fist estimates the variance component with
only the covariates, and then uses it to analyze each SNP after
filtering. These two steps can be done separately as follows. The first
command without `bed` only esitmates the variance component tau, and the
second command uses the estimated tau to analyze the SNPs.

    re = coxmeg_plink(pheno,sigma,type='bd',cov_file=cov,verbose=FALSE)
    re = coxmeg_plink(pheno,sigma,type='bd',bed=bed,tmp_dir=tempdir(),tau=re$tau,cov_file=cov,verbose=FALSE)

When the genotypes of a group of SNPs are stored in a matrix, the
function `coxmeg_m` can be used to perform GWAS for each of the SNPs.
Similarly, `coxmeg_m` first estimates the variance component without the
SNPs. In the following example, we simulate 10 independent SNPs, and use
`coxmeg_m` to perform association analysis.

    geno = matrix(rbinom(nrow(sigma)*10,2,runif(nrow(sigma)*10,0.05,0.5)),nrow(sigma),10)
    pheno_m = read.table(pheno)
    re = coxmeg_m(geno,pheno_m[,3:4],sigma,type='bd',verbose=FALSE)
    re

    ## $summary
    ##            beta        HR    sd_beta          p
    ## 1  -0.013836758 0.9862585 0.02995662 0.64415798
    ## 2  -0.004709216 0.9953019 0.03023657 0.87623347
    ## 3   0.029468921 1.0299074 0.02963562 0.32004028
    ## 4  -0.017503736 0.9826486 0.02945657 0.55236378
    ## 5   0.048644826 1.0498474 0.02961778 0.10050308
    ## 6  -0.053620077 0.9477921 0.02957838 0.06986063
    ## 7   0.074566608 1.0774171 0.03021080 0.01357904
    ## 8   0.023469183 1.0237468 0.02952778 0.42671989
    ## 9   0.023355308 1.0236302 0.02998626 0.43605809
    ## 10 -0.025082039 0.9752299 0.02937368 0.39316324
    ## 
    ## $tau
    ## [1] 0.04052206
    ## 
    ## $rank
    ## [1] 3000
    ## 
    ## $nsam
    ## [1] 3000

By default, `coxmeg_m` and `coxmeg_plink` will choose an optimal `order`
between 1 and 10 for analyzing the SNPs when `order` is not specified.

Perform GWAS of an age-at-onset phenotype with a dense relatedness matrix
-------------------------------------------------------------------------

When the relatedness matrix is dense and large (&gt;5000),
`type='dense'` should be used. In thise case, it will be more efficient
to use preconditioned conjugate gradiant (PCG) (e.g., by explicitly
specifying `solver=2`) and stochastic lanczos quadrature (SLQ)
`detap='slq'` in the optimization. These can be specified as follows.

    re = coxmeg_plink(pheno,sigma,type='dense',bed=bed,tmp_dir=tempdir(),cov_file=cov,detap='slq',verbose=FALSE,solver=2)

If `solver` is not specified, `coxmeg_plink` will by default choose PCG
as a solver when `type='dense'`. If `detap` is not specified,
`coxmeg_plink` will by default use `detap='slq'` for a dense matrix when
the sample size exceeds 3000. The number of Monte Carlo samples in the
SLQ can be specified by `mc` (by default `mc=100`).

The above command estimates HRs and reports p-values. Instead, a score
test, which is computationally much more efficient, can be used by
specifying `score=TRUE`.

    re = coxmeg_plink(pheno,sigma,type='dense',bed=bed,tmp_dir=tempdir(),tau=re$tau,cov_file=cov,detap='slq',verbose=FALSE,solver=2,score=TRUE)

    ## Excluding 0 SNP on non-autosomes
    ## Excluding 0 SNP (monomorphic: TRUE, MAF: 0.05, missing rate: 0)

    re

    ## $summary
    ##     snp.id chromosome position allele      afreq  afreq_inc   index
    ## 1   null_0          1        1    d/D 0.30983333 0.30983333  null_0
    ## 2   null_1          1        2    d/D 0.23466667 0.23466667  null_1
    ## 3   null_2          1        3    D/d 0.14033333 0.14033333  null_2
    ## 4   null_3          1        4    D/d 0.16183333 0.16183333  null_3
    ## 5   null_4          1        5    d/D 0.19933333 0.19933333  null_4
    ## 6   null_5          1        6    D/d 0.11800000 0.11800000  null_5
    ## 7   null_6          1        7    d/D 0.09483333 0.09483333  null_6
    ## 8   null_7          1        8    D/d 0.49683333 0.49683333  null_7
    ## 9   null_8          1        9    d/D 0.31366667 0.31366667  null_8
    ## 10  null_9          1       10    D/d 0.49183333 0.49183333  null_9
    ## 11 null_10          1       11    d/D 0.34833333 0.34833333 null_10
    ## 12 null_11          1       12    D/d 0.25100000 0.25100000 null_11
    ## 13 null_12          1       13    d/D 0.17500000 0.17500000 null_12
    ## 14 null_13          1       14    D/d 0.06333333 0.06333333 null_13
    ## 15 null_14          1       15    D/d 0.20833333 0.20833333 null_14
    ## 16 null_15          1       16    d/D 0.17050000 0.17050000 null_15
    ## 17 null_16          1       17    D/d 0.33550000 0.33550000 null_16
    ## 18 null_17          1       18    d/D 0.26633333 0.26633333 null_17
    ## 19 null_18          1       19    D/d 0.09433333 0.09433333 null_18
    ## 20 null_19          1       20    d/D 0.11650000 0.11650000 null_19
    ##           score  score_test          p
    ## 1   0.015731431 0.284915628 0.59349729
    ## 2   0.019537166 0.364147765 0.54621165
    ## 3  -0.049037261 1.669161029 0.19637095
    ## 4   0.044750829 1.421986351 0.23307674
    ## 5   0.028694042 0.688246877 0.40676132
    ## 6  -0.110065139 7.299184776 0.00689859
    ## 7  -0.017846384 0.148832197 0.69965386
    ## 8  -0.004209841 0.023985738 0.87692121
    ## 9  -0.063092716 4.642530182 0.03118899
    ## 10 -0.008409498 0.094820647 0.75813589
    ## 11 -0.013558821 0.225502351 0.63487901
    ## 12  0.037870134 1.452836527 0.22807335
    ## 13 -0.017137550 0.225089257 0.63518922
    ## 14 -0.066389447 1.432826764 0.23130366
    ## 15 -0.013905384 0.165261969 0.68435744
    ## 16  0.002166564 0.003450204 0.95316044
    ## 17  0.004762768 0.027683654 0.86785472
    ## 18  0.001804748 0.003388839 0.95357838
    ## 19 -0.015963997 0.115030323 0.73448829
    ## 20 -0.022237828 0.280375767 0.59645503
    ## 
    ## $tau
    ## [1] 0.04052206
    ## 
    ## $rank
    ## [1] 3000
    ## 
    ## $nsam
    ## [1] 3000

In this result, the column `score_test` is the score test statistic,
which follows a *χ*<sup>2</sup> distribution with 1 d.f. The column
`score` is the score function divided by its variance. Note that `score`
is not the estimate of log(HR) under the full model. It is actually a
one-step update of the Newton-Raphson algoritm starting from the values
estimated under the null model. Comparing `score` with log(HR) estimated
in the previous section, we see that they are close to each other in
this example. However, the difference can be large if the genotype is
highly correlated with the covariates or random effects.

Handle positive semidefinite relatedness matrices
-------------------------------------------------

We now assume that the first two subjects in the sample are monozygotic
twins. In this case, the relatedness matrix becomes positive
semidefinite. Specifying `spd=FALSE` will let `coxmeg_plink` handle a
positive semidefinite relatedness matrix.

    sigma[2,1] = sigma[1,2] = 1
    re = coxmeg_plink(pheno,sigma,type='bd',cov_file=cov,verbose=FALSE,spd=FALSE)

    ## Warning in chol.default(x, pivot = TRUE): the matrix is either rank-deficient or
    ## indefinite

    re

    ## $tau
    ## [1] 0.04024134
    ## 
    ## $iter
    ## [1] 15
    ## 
    ## $rank
    ## [1] 3000
    ## 
    ## $nsam
    ## [1] 3000

The warning indicates that the relatedness matrix is not full rank.
Because there is a twin pair in the sample, the rank of the relatedness
matrix is less than the sample size. If the user is not sure whether the
relatedness matrix is positive definite or positive semidefinite, it is
better to use `spd=FALSE` although this is slower because coxmeg will
perform an eigenvalue decomposition under this setting. In the current
version, instead of using the previously proposed GPPL in (He and
Kulminski 2020), coxmeg uses a modified PPL by turning all zero
eigenvalues of the relatedness matrix to a small value (1e-6). This
modification makes coxmeg suitable for twin cohorts.

References
----------

He, Liang, and Alexander M. Kulminski. 2020. “Fast Algorithms for
Conducting Large-Scale GWAS of Age-at-Onset Traits Using Cox
Mixed-Effects Models.” *Genetics*, March, genetics.302940.2019.
<https://doi.org/10.1534/genetics.119.302940>.
