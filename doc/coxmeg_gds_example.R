## ---- message=FALSE-----------------------------------------------------------
library(coxmeg)
library(gdsfmt)
library(SNPRelate)
snpfile <- snpgdsExampleFileName()
snp <- snpgdsOpen(snpfile)

## -----------------------------------------------------------------------------
grm <- snpgdsGRM(snp, method="GCTA", verbose=FALSE)
sigma <- grm$grm
dimnames(sigma) <- list(grm$sample.id, grm$sample.id)
sigma[1:5, 1:5]

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
sex <- read.gdsn(index.gdsn(snp, "sample.annot/sex"))
pop <- read.gdsn(index.gdsn(snp, "sample.annot/pop.group"))
cov <- data.frame(sex, pop)
head(cov)
mm <- model.matrix(~ sex + pop, cov)
cov <- cbind(pheno[,1:2], mm[,-1])
head(cov)

## -----------------------------------------------------------------------------
snp.id <- read.gdsn(index.gdsn(snp, "snp.id"))
re <- coxmeg_gds(snp, pheno, sigma, type='dense', cov=cov, snp.id=snp.id[1:100], spd=FALSE)
head(re$summary)

## -----------------------------------------------------------------------------
snpgdsClose(snp)

## -----------------------------------------------------------------------------
library(SeqArray)
seqfile <- seqExampleFileName()
seq <- seqOpen(seqfile)

## -----------------------------------------------------------------------------
king <- snpgdsIBDKING(seq, verbose=FALSE)
sigma <- GENESIS::kingToMatrix(king, thresh=0.177) * 2
sigma[1:5,1:5]

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
seqSetFilter(seq, variant.sel=1:100)
re <- coxmeg_gds(seq, pheno, sigma, type='bd')
head(re$summary)

## -----------------------------------------------------------------------------
seqClose(seq)

