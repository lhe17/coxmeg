## Overview

In addition to the functions described in the main coxmeg vignette, the
`coxmeg_gds` function allows running coxmeg directly on Genomic Data
Structure (GDS) file objects. The basic structure of GDS files is
described in the gdsfmt package.

There are two implementations of GDS files. The original format was
developed for SNP arrays, and an interface to files in this format is
defined in the SNPRelate package. A newer format developed for
sequencing was designed to import all data stored in VCF files. An
interface to the sequencing format is defined in the SeqArray package.
`coxmeg_gds` supports both file types.

## SNPRelate example

In the first example, we use an original GDS file and open it with the
SNPRelate package.

    library(coxmeg)
    library(gdsfmt)
    library(SNPRelate)
    snpfile <- snpgdsExampleFileName()
    snp <- snpgdsOpen(snpfile)

We use the `snpgdsGRM` function in SNPRelate to create a genetic
relationship matrix (GRM) with the GCTA method.

    grm <- snpgdsGRM(snp, method="GCTA", verbose=FALSE)
    sigma <- grm$grm
    dimnames(sigma) <- list(grm$sample.id, grm$sample.id)
    sigma[1:5, 1:5]

    ##            NA19152    NA19139    NA18912    NA19160    NA07034
    ## NA19152  1.3946382  0.2396834  0.2516490  0.2578472 -0.1005009
    ## NA19139  0.2396834  1.3284340  0.2188917  0.2378597 -0.1076276
    ## NA18912  0.2516490  0.2188917  1.3785072  0.2632486 -0.1027841
    ## NA19160  0.2578472  0.2378597  0.2632486  1.4305388 -0.1342964
    ## NA07034 -0.1005009 -0.1076276 -0.1027841 -0.1342964  0.9723706

We create a data.frame of simulated time-to-event outcomes. The first
two columns of the data.frame are expected to be family id and sample
(individual) id.

    sample.id <- read.gdsn(index.gdsn(snp, "sample.id"))
    family.id <- read.gdsn(index.gdsn(snp, "sample.annot/family.id"))
    n <- length(sample.id)
    set.seed(5)
    time <- rnorm(n, mean=100, sd=10)
    set.seed(6)
    status <- rbinom(n, 1, 0.4)
    pheno <- data.frame(family.id, sample.id, time, status,
                        stringsAsFactors=FALSE)
    head(pheno)

    ##   family.id sample.id      time status
    ## 1        72   NA19152  91.59145      1
    ## 2        43   NA19139 113.84359      1
    ## 3        28   NA18912  87.44508      0
    ## 4        56   NA19160 100.70143      0
    ## 5      1341   NA07034 117.11441      1
    ## 6      1341   NA07055  93.97092      1

We will adjust for sex and population group as covariates in the model.
As in the outcome data.frame, family id and sample id are the first two
columns. Categorical variables need to be converted to dummy variables,
so we utilize the `model.matrix` function to prepare the covariates.

    sex <- read.gdsn(index.gdsn(snp, "sample.annot/sex"))
    pop <- read.gdsn(index.gdsn(snp, "sample.annot/pop.group"))
    cov <- data.frame(sex, pop)
    head(cov)

    ##   sex pop
    ## 1   F YRI
    ## 2   M YRI
    ## 3   F YRI
    ## 4   M YRI
    ## 5   M CEU
    ## 6   F CEU

    mm <- model.matrix(~ sex + pop, cov)
    cov <- cbind(pheno[,1:2], mm[,-1])
    head(cov)

    ##   family.id sample.id sexM popHCB popJPT popYRI
    ## 1        72   NA19152    0      0      0      1
    ## 2        43   NA19139    1      0      0      1
    ## 3        28   NA18912    0      0      0      1
    ## 4        56   NA19160    1      0      0      1
    ## 5      1341   NA07034    1      0      0      0
    ## 6      1341   NA07055    0      0      0      0

The coxmeg\_gds function fits a Cox mixed-effects model to variants
stored in a GDS file. The `snp.id` argument allows selecting a subset of
variants.

The GRM is a dense matrix and not postive definite, so we set
`type='dense'` and `spd=FALSE`.

    snp.id <- read.gdsn(index.gdsn(snp, "snp.id"))
    re <- coxmeg_gds(snp, pheno, sigma, type='dense', cov=cov, snp.id=snp.id[1:100], spd=FALSE)

    ## There are 279 subjects who have genotype data and have no missing phenotype or covariates.

    ## Remove 0 subjects censored before the first failure.

    ## There is/are 4 covariates. The sample size included is 279.

    ## Warning in coxmeg_gds(snp, pheno, sigma, type = "dense", cov = cov, snp.id =
    ## snp.id[1:100], : The relatedness matrix has negative eigenvalues. Please use a
    ## positive (semi)definite matrix.

    ## The correlation matrix is treated as dense.

    ## The relatedness matrix is inverted.

    ## The method for computing the determinant is 'exact'.

    ## Solver: PCG (RcppEigen:dense).

    ## Warning in check_tau(tau_e, min_tau, max_tau, ncm): The estimated variance
    ## component equals the lower bound (1e-04), probably suggesting no random effects.

    ## The variance component is estimated. Start analyzing SNPs...

    ## Excluding 28 SNPs (monomorphic: TRUE, MAF: 0.05, missing rate: 0)

    head(re$summary)

    ##   snp.id chromosome position allele     afreq afreq_inc index        beta
    ## 1      1          1  1355433    G/T 0.6989247 0.6989247     1 -0.23752776
    ## 2      2          1  3114015    C/T 0.2043011 0.2043011     2 -0.03767414
    ## 5      5          1  4271809    C/T 0.8817204 0.8817204     5 -0.19487844
    ## 6      6          1  4358909    C/T 0.7956989 0.7956989     6  0.20594275
    ## 8      8          1  4515962    A/G 0.3369176 0.3369176     8 -0.01275018
    ## 9      9          1  4639285    C/G 0.2007168 0.2007168     9 -0.14976807
    ##          HR   sd_beta         p
    ## 1 0.7885750 0.2321316 0.3061916
    ## 2 0.9630267 0.2036705 0.8532479
    ## 5 0.8229347 0.2174205 0.3700817
    ## 6 1.2286829 0.1798058 0.2520591
    ## 8 0.9873308 0.1555024 0.9346519
    ## 9 0.8609076 0.1750789 0.3923122

    snpgdsClose(snp)

## SeqArray example

In this example, we use a sequencing GDS file and open it with the
SeqArray package.

    library(SeqArray)
    seqfile <- seqExampleFileName()
    seq <- seqOpen(seqfile)

As an alternate method of creating a random effects matrix, we use the
KING algorithm to estimate pairwise relatedness, and use the GENESIS
package to transform the results into a sparse block-diagonal matrix.
The threshold chosen corresponds to setting values for pairs less
closely related than second-degree relatives to 0. We multiply the
matrix by 2 so diagonal elements are 1 rather than 0.5 (the latter is
the kinship coefficient for identical genomes).

    library(GENESIS)
    if(requireNamespace('GENESIS', quietly = TRUE))
    {
      king <- snpgdsIBDKING(seq, verbose=FALSE)
      sigma <- GENESIS::kingToMatrix(king, thresh=0.177) * 2
      sigma[1:5,1:5]
    }

    ## Using 90 samples provided

    ## Identifying clusters of relatives...

    ##     16 relatives in 2 clusters; largest cluster = 13

    ## Creating block matrices for clusters...

    ## 74 samples with no relatives included

    ## Putting all samples together into one block diagonal matrix

    ## 5 x 5 sparse Matrix of class "dsCMatrix"
    ##           NA06984   NA06989   NA12156  NA11832  NA12249
    ## NA06984 1.0000000 0.3580247 0.2153846 .        .       
    ## NA06989 0.3580247 1.0000000 0.3870968 .        .       
    ## NA12156 0.2153846 0.3870968 1.0000000 .        .       
    ## NA11832 .         .         .         1.000000 0.238255
    ## NA12249 .         .         .         0.238255 1.000000

We create a data.frame of simulated time-to-event outcomes. The first
two columns of the data.frame are expected to be family id and sample
(individual) id.

    sample.id <- seqGetData(seq, "sample.id")
    family.id <- seqGetData(seq, "sample.annotation/family")
    family.id[family.id == ""] <- sample.id[family.id == ""]
    n <- length(sample.id)
    set.seed(35)
    time <- rnorm(n, mean=100, sd=10)
    set.seed(36)
    status <- rbinom(n, 1, 0.4)
    pheno <- data.frame(family.id, sample.id, time, status,
                        stringsAsFactors=FALSE)
    head(pheno)

    ##   family.id sample.id      time status
    ## 1      1328   NA06984 110.65125      1
    ## 2   NA06985   NA06985 101.32881      1
    ## 3     13291   NA06986  99.65956      1
    ## 4      1328   NA06989  99.55024      0
    ## 5      1340   NA06994 133.37838      1
    ## 6      1340   NA07000  96.07082      0

The SeqArray package allows for pre-selecting variants using the
`seqSetFilter` function. We set `type='bd'` since `sigma` is a block
diagonal matrix.

    seqSetFilter(seq, variant.sel=1:100)

    ## # of selected variants: 100

    re <- coxmeg_gds(seq, pheno, sigma, type='bd')

    ## There are 90 subjects who have genotype data and have no missing phenotype or covariates.

    ## Remove 3 subjects censored before the first failure.

    ## There is/are 0 covariates. The sample size included is 87.

    ## as(<dsCMatrix>, "dgCMatrix") is deprecated since Matrix 1.5-0; do as(., "generalMatrix") instead

    ## The correlation matrix is treated as sparse/block diagonal.

    ## The relatedness matrix is inverted.

    ## The method for computing the determinant is 'diagonal'.

    ## Solver: Cholesky decomposition (RcppEigen=TRUE).

    ## The variance component is estimated. Start analyzing SNPs...

    ## # of selected samples: 87
    ## [..................................................]  0%, ETC: ---    [==================================================] 100%, completed, 0s
    ## # of selected variants: 3

    ## The order is set to be 2.

    head(re$summary)

    ##    snp.id chromosome  position allele      afreq  afreq_inc index       beta
    ## 9       9          1  34435767    G/A 0.85555556 0.85057471     9 -0.4738578
    ## 56     56          1  89618740    C/T 0.90000000 0.89655172    56  0.6003386
    ## 93     93          1 114749804    A/G 0.09444444 0.09770115    93  0.9905076
    ##           HR   sd_beta           p
    ## 9  0.6225958 0.3507818 0.176739691
    ## 56 1.8227359 0.4426871 0.175059724
    ## 93 2.6926010 0.3786899 0.008906643

    seqClose(seq)

In this vignette we presented two examples, each using a different GDS
file type and method of creating a random effects matrix. We note that
the choice of random effects matrix in each case was arbitrary; and
either method of generating the matrix would work with either GDS file
type.

The code in `coxmeg_gds` is very similar to the code in `coxmeg_plink`,
as the latter function converts the plink file to GDS before reading
genotypes. `coxmeg_gds` takes R objects as arguments for phenotypes,
covariates, and a connection to the GDS file; while `coxmeg_plink` takes
paths to files on disk. If the user intends to utilize any other
functions requiring a GDS object, or to run multiple analyses on the
same dataset, it will be more efficient to convert a plink file to GDS
first with the SNPRelate function `snpgdsBED2GDS` (or the SeqArray
function `seqBED2GDS`), then use `coxmeg_gds` instead of `coxmeg_plink`.
