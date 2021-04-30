context("test gds function")

.testPheno <- function(sample.id, seed=25) {
    n <- length(sample.id)
    set.seed(seed)
    time <- rnorm(n, mean=100, sd=10)
    set.seed(seed)
    status <- rbinom(n, 1, 0.4)
    data.frame(family=sample.id, sample.id,
               time, status, stringsAsFactors=FALSE)
}

.testSparseMatrix <- function(n, n_blocks) {
    n_f <- ceiling(n/n_blocks)
    mat_list <- list()
    size <- rep(n_blocks,n_f)
    offd <- 0.5
    for(i in 1:n_f){
        mat_list[[i]] <- matrix(offd,size[i],size[i])
        diag(mat_list[[i]]) <- 1
    }
    sigma <- as.matrix(bdiag(mat_list))
    sigma <- sigma[1:n, 1:n]
    sigma <- as(sigma,'dgCMatrix')
    return(sigma)
}

test_that("coxmeg_gds matches coxmeg_plink", {
    bed = system.file("extdata", "example_null.bed", package = "coxmeg")
    bed = substr(bed,1,nchar(bed)-4)
    pheno.file = system.file("extdata", "ex_pheno.txt", package = "coxmeg")
    cov.file = system.file("extdata", "ex_cov.txt", package = "coxmeg")
    
    ## building a relatedness matrix
    sigma <- .testSparseMatrix(n=3000, n_blocks=5)
    
    re.plink = coxmeg_plink(pheno.file,sigma,type='bd',bed=bed,tmp_dir=tempdir(),cov_file=cov.file, verbose=FALSE)
    
    gdsfile <- tempfile()
    SNPRelate::snpgdsBED2GDS(bed.fn=paste0(bed,".bed"), fam.fn=paste0(bed,".fam"), bim.fn=paste0(bed,".bim"),
                             out.gdsfn=gdsfile, verbose=FALSE)
    gds <- SNPRelate::snpgdsOpen(gdsfile)
    pheno <- read.table(pheno.file, header=FALSE, as.is=TRUE, na.strings="-9")
    cov <- read.table(cov.file, header=FALSE, as.is=TRUE)
    
    re.gds <- coxmeg_gds(gds,pheno,sigma,type='bd',cov=cov,verbose=FALSE)
    expect_equal(re.plink, re.gds, tolerance=1e-5)
    
    SNPRelate::snpgdsClose(gds)
    unlink(gdsfile)
})


test_that("SNPRelate and SeqArray methods match", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    snpfile <- SNPRelate::snpgdsExampleFileName()
    seqfile <- tempfile()
    SeqArray::seqSNP2GDS(snpfile, seqfile, verbose=FALSE)
    
    snp <- SNPRelate::snpgdsOpen(snpfile)
    seq <- SeqArray::seqOpen(seqfile)
    
    snpsel <- .gdsSelectSNP(snp, maf=0.01, missing.rate=0.01, verbose=FALSE)
    seqsel <- .gdsSelectSNP(seq, maf=0.01, missing.rate=0.01, verbose=FALSE)
    expect_equal(snpsel, seqsel)
    
    seqResetFilter(seq, verbose=FALSE)
    sample.id <- seqGetData(seq, "sample.id")[1:50]
    snp.id <- seqGetData(seq, "variant.id")[1:100]
    snpref <- substr(gdsfmt::read.gdsn(gdsfmt::index.gdsn(snp, "snp.allele")), 1, 1)
    seqref <- seqGetData(seq, "$ref")
    allele.swap <- (snpref != seqref)
    swap.sel <- allele.swap[1:100]
    snpgeno <- .gdsGetGeno(snp, sample.id=sample.id, snp.id=snp.id, verbose=FALSE)
    seqgeno <- .gdsGetGeno(seq, sample.id=sample.id, snp.id=snp.id, verbose=FALSE)
    expect_equivalent(snpgeno[,!swap.sel], seqgeno[,!swap.sel])
    expect_equivalent(snpgeno[,swap.sel], 2-seqgeno[,swap.sel])
    
    snplist <- .gdsSNPList(snp)
    seqlist <- .gdsSNPList(seq)
    snplist$chromosome <- as.character(snplist$chromosome)
    expect_equivalent(snplist[!allele.swap,], seqlist[!allele.swap,])
    expect_equivalent(snplist[allele.swap,1:3], seqlist[allele.swap,1:3])
    expect_equal(snplist$afreq[allele.swap], 1-seqlist$afreq[allele.swap])
    
    SNPRelate::snpgdsClose(snp)
    SeqArray::seqClose(seq)
    unlink(seqfile)
})


test_that("SNPRelate and SeqArray coxmeg_gds match", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    snpfile <- SNPRelate::snpgdsExampleFileName()
    seqfile <- tempfile()
    SeqArray::seqSNP2GDS(snpfile, seqfile, verbose=FALSE)
    
    snp <- SNPRelate::snpgdsOpen(snpfile)
    seq <- SeqArray::seqOpen(seqfile)
    
    sample.id <- seqGetData(seq, "sample.id")
    pheno <- .testPheno(sample.id, seed=5)
    
    # covariance matrix
    #covmat <- SNPRelate::snpgdsGRM(snp, verbose=FALSE)
    #sigma <- covmat$grm
    sigma <- .testSparseMatrix(n=nrow(pheno), n_blocks=5)
    
    # set high MAF so test runs faster
    re.snp <- coxmeg_gds(snp, pheno, sigma, type='bd', maf=0.47, verbose=FALSE)
    re.seq <- coxmeg_gds(seq, pheno, sigma, type='bd', maf=0.47, verbose=FALSE)
    allele.swap <- re.snp$summary$allele != re.seq$summary$allele
    expect_equal(re.snp$summary$beta[!allele.swap], re.seq$summary$beta[!allele.swap], tolerance=1e-4)
    expect_equal(re.snp$summary$beta[allele.swap], -re.seq$summary$beta[allele.swap], tolerance=1e-4)
    expect_equal(re.snp$summary$HR[!allele.swap], re.seq$summary$HR[!allele.swap], tolerance=1e-4)
    expect_equal(re.snp$summary$HR[allele.swap], 1/re.seq$summary$HR[allele.swap], tolerance=1e-4)
    expect_equal(re.snp$summary$sd_beta, re.seq$summary$sd_beta, tolerance=1e-4)
    expect_equal(re.snp$summary$p, re.seq$summary$p, tolerance=1e-4)
    
    SNPRelate::snpgdsClose(snp)
    SeqArray::seqClose(seq)
    unlink(seqfile)
})


test_that("SeqArray utils respect variant filters", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- seqExampleFileName()
    gds <- seqOpen(gdsfile)
    
    x <- .gdsSelectSNP(gds, maf=0.01, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:500, verbose=FALSE)
    x2 <- .gdsSelectSNP(gds, maf=0.01, verbose=FALSE)
    expect_true(length(x) > length(x2))
    
    # with sample.id
    seqResetFilter(gds, verbose=FALSE)
    samp <- seqGetData(gds, "sample.id")[1:50]
    x <- .gdsSelectSNP(gds, sample.id=samp, maf=0.01, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:500, verbose=FALSE)
    x2 <- .gdsSelectSNP(gds, sample.id=samp, maf=0.01, verbose=FALSE)
    expect_true(length(x) > length(x2))
    
    seqClose(gds)  
})


test_that("SeqArray method respects sample filters", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- SeqArray::seqExampleFileName()
    gds <- SeqArray::seqOpen(gdsfile)
    
    sample.id <- SeqArray::seqGetData(gds, "sample.id")
    pheno <- .testPheno(sample.id, seed=7)
    sigma <- .testSparseMatrix(n=nrow(pheno), n_blocks=5)
    
    SeqArray::seqSetFilter(gds, sample.id=sample.id[1:50], variant.sel=1:1000, verbose=FALSE)
    re <- coxmeg_gds(gds, pheno, sigma, type='bd', verbose=FALSE)
    expect_true(re$nsam <= 50)
    
    SeqArray::seqClose(gds)  
})


test_that("SNPRelate method subsets to samples in GDS", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    snpfile <- SNPRelate::snpgdsExampleFileName()
    gds <- SNPRelate::snpgdsOpen(snpfile)
    
    sample.id <- SNPRelate::snpgdsSummary(gds, show=FALSE)$sample.id
    sample.add <- c(sample.id, letters)
    pheno <- .testPheno(sample.add, seed=9)
    sigma <- .testSparseMatrix(n=nrow(pheno), n_blocks=7)
    
    re <- coxmeg_gds(gds, pheno, sigma, type='bd', maf=0.47, verbose=FALSE)
    expect_true(re$nsam <= length(sample.id))
    
    SNPRelate::snpgdsClose(gds)
})


test_that("SeqArray method selects snp.id", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- SeqArray::seqExampleFileName()
    gds <- SeqArray::seqOpen(gdsfile)
    
    sample.id <- SeqArray::seqGetData(gds, "sample.id")
    pheno <- .testPheno(sample.id, seed=7)
    sigma <- .testSparseMatrix(n=nrow(pheno), n_blocks=5)
    
    snp.id <- SeqArray::seqGetData(gds, "variant.id")[1:100]
    re <- coxmeg_gds(gds, pheno, sigma, type='bd', snp.id=snp.id, maf=0, verbose=FALSE)
    expect_true(nrow(re$summary) <= 100)
    
    SeqArray::seqClose(gds)
})


test_that("SNPRelate method selects snp.id", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    snpfile <- SNPRelate::snpgdsExampleFileName()
    gds <- SNPRelate::snpgdsOpen(snpfile)
    
    gdssum <- SNPRelate::snpgdsSummary(gds, show=FALSE)
    pheno <- .testPheno(gdssum$sample.id, seed=9)
    sigma <- .testSparseMatrix(n=nrow(pheno), n_blocks=7)
    
    snp.id <- gdssum$snp.id[1:100]
    re <- coxmeg_gds(gds, pheno, sigma, type='bd', snp.id=snp.id, maf=0, verbose=FALSE)
    expect_true(nrow(re$summary) <= 100)
    
    SNPRelate::snpgdsClose(gds)
})
