setOldClass("SNPGDSFileClass")
setGeneric(".gdsHasSamp", function(gdsobj, ...) standardGeneric(".gdsHasSamp"))
setGeneric(".gdsGetSamp", function(gdsobj, ...) standardGeneric(".gdsGetSamp"))
setGeneric(".gdsSelectSNP", function(gdsobj, ...) standardGeneric(".gdsSelectSNP"))
setGeneric(".gdsGetGeno", function(gdsobj, ...) standardGeneric(".gdsGetGeno"))
setGeneric(".gdsSNPList", function(gdsobj, ...) standardGeneric(".gdsSNPList"))

setMethod(".gdsHasSamp",
          "SNPGDSFileClass",
          function(gdsobj, sample.id) {
              gds.samp <- snpgdsSummary(gdsobj, show=FALSE)$sample.id
              sample.id %in% gds.samp
          })

setMethod(".gdsHasSamp",
          "SeqVarGDSClass",
          function(gdsobj, sample.id) {
              gds.samp <- seqGetData(gdsobj, "sample.id")
              sample.id %in% gds.samp
          })

setMethod(".gdsGetSamp",
          "SNPGDSFileClass",
          function(gdsobj) {
            snpgdsSummary(gdsobj, show=FALSE)$sample.id
          })

setMethod(".gdsGetSamp",
          "SeqVarGDSClass",
          function(gdsobj) {
            seqGetData(gdsobj, "sample.id")
          })

setMethod(".gdsSelectSNP",
          "SNPGDSFileClass",
          function(gdsobj, sample.id=NULL, snp.id=NULL, maf=NaN, 
                   missing.rate=NaN, verbose=TRUE){
              snpgdsSelectSNP(gdsobj, sample.id=sample.id, snp.id=snp.id,
                              maf=maf, missing.rate=missing.rate,
                              verbose=verbose,
                              remove.monosnp=TRUE,
                              autosome.only=FALSE)
          })

setMethod(".gdsSelectSNP",
          "SeqVarGDSClass",
          function(gdsobj, sample.id=NULL, snp.id=NULL, maf=NaN, 
                   missing.rate=NaN, verbose=TRUE){
              if (!is.null(sample.id)) {
                  seqSetFilter(gdsobj, sample.id=sample.id, verbose=verbose)
              }
              if (!is.null(snp.id)) {
                  seqSetFilter(gdsobj, variant.id=snp.id, verbose=verbose)
              }
              seqSetFilterCond(gdsobj, maf=maf, missing.rate=missing.rate, 
                               mac=1L, verbose=verbose)
              seqGetData(gdsobj, "variant.id")
          })


setMethod(".gdsGetGeno",
          "SNPGDSFileClass",
          function(gdsobj, sample.id=NULL, snp.id=NULL, verbose=TRUE){
              snpgdsGetGeno(gdsobj, sample.id=sample.id, snp.id=snp.id,
                            with.id=FALSE, verbose=verbose)
          })

setMethod(".gdsGetGeno",
          "SeqVarGDSClass",
          function(gdsobj, sample.id=NULL, snp.id=NULL, verbose=TRUE){
              seqSetFilter(gdsobj, sample.id=sample.id, variant.id=snp.id, 
                           verbose=verbose)
              seqGetData(gdsobj, "$dosage")
          })


setMethod(".gdsSNPList",
          "SNPGDSFileClass",
          function(gdsobj, sample.id=NULL){
              snpgdsSNPList(gdsobj, sample.id=sample.id)
          })

setMethod(".gdsSNPList",
          "SeqVarGDSClass",
          function(gdsobj, sample.id=NULL){
              seqResetFilter(gdsobj, verbose=FALSE)
              seqSetFilter(gdsobj, sample.id=sample.id, verbose=FALSE)
              snp.id <- seqGetData(gdsobj, "variant.id")
              chromosome <- seqGetData(gdsobj, "chromosome")
              position <- seqGetData(gdsobj, "position")
              allele <- sub(",", "/", seqGetData(gdsobj, "allele"))
              afreq <- seqAlleleFreq(gdsobj)
              data.frame(snp.id, chromosome, position, allele, afreq,
                         stringsAsFactors=FALSE)
          })
