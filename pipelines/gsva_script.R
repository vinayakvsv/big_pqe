prepare_gsva <- function() {
  #gsva
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("GSVA")
  #piano
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("piano")
  #qusage
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("qusage")
  
  #libraries
  library(GSVA)
  library(GSEABase)
  library(KEGG.db)
  library(qusage)
  library(limma)
  library(snow)
  library(rlecuyer)
  
}
#
import_counts <- function(countsfile,imported_counts=TRUE) {
  if (!imported_counts) {
    counts <- read.table(file = countsfile,
                         header = TRUE,
                         sep = "\t",
                         row.names = 1)
  } else {
    counts <- countsfile
  }
  
  gene.names <- sapply(rownames(counts),function(x) {a <- unlist(strsplit(as.character(x),split = "_")); a[2]})
  counts.aggregate <- aggregate(x = counts,by = list(gene.names),FUN = mean)
  rownames(counts.aggregate) <- counts.aggregate$Group.1
  counts.aggregate <- counts.aggregate[,c(2:ncol(counts.aggregate))]
  #rownames(counts.aggregate) <- gene.names
  output <- list(counts.aggregate,gene.names)
  
  names(output) <- c("counts.aggregate","gene.names")
  return(output)
}
#
reformat_counts <- function(counts,capitalize=TRUE) {
  gene.names <- sapply(rownames(counts),function(x) {a <- unlist(strsplit(as.character(x),split = "_")); a[2]})
  counts.aggregate <- aggregate(x = counts,by = list(gene.names),FUN = mean)
  
  if (capitalize) {
    rownames(counts.aggregate) <- toupper(counts.aggregate$Group.1)
  } else {
    rownames(counts.aggregate) <- counts.aggregate$Group.1
  }
  
  counts.aggregate <- counts.aggregate[,c(2:ncol(counts.aggregate))]
  #rownames(counts.aggregate) <- gene.names
  output <- list(counts.aggregate,gene.names)
  
  names(output) <- c("counts.aggregate","gene.names")
  return(output)
}

#
import_gmt <- function(gmtfile) {
  #How do we import a GMT file and convert it to a GeneSet object
  #http://svitsrv25.epfl.ch/R-doc/library/GSEABase/html/getObjects.html
  #BETTER: https://www.rdocumentation.org/packages/qusage/versions/2.4.0/topics/read.gmt
  in.gmt <- read.gmt(file = gmtfile)
  return(in.gmt)
}
#
gsva_enrich <- function(counts,geneset,ikcdf="Poisson",...) {
  enrichresults <- gsva(expr = as.matrix(counts), 
                        gset.idx.list = geneset, 
                        min.sz=2, 
                        max.sz=500,
                        mx.diff=TRUE,
                        verbose=FALSE, 
                        kcdf=ikcdf, 
                        parallel.sz=2,...)
  #issues with bootstrapping addressed here: https://support.bioconductor.org/p/91161/
  #enrichresults$$es.obs
  return(enrichresults)
}
#
calculate_significance <- function(metadata) {
  #this function uses the limma package
  adjPvalueCutoff <- 0.001
  logFCcutoff <- log2(2)
  
  design <- model.matrix(~ factor(metadata$Phenotype))
  colnames(design) <- c("CTRL", "MLLvsALL")
  fit <- lmFit(leukemia_es, design)
  fit <- eBayes(fit)
  allGeneSets <- topTable(fit, coef="MLLvsALL", number=Inf)
  DEgeneSets <- topTable(fit, coef="MLLvsALL", number=Inf,
                         p.value=adjPvalueCutoff, adjust="BH")
  res <- decideTests(fit, p.value=adjPvalueCutoff)
  summary(res)
}
#
#
runSim <- function(p, n, gs.sz, S2N, fracDEgs) {
  #this code was taken from https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf
  sizeDEgs <- round(fracDEgs * gs.sz)
  group.n <- round(n / 2)
  
  sampleEffect <- rnorm(n, mean=0, sd=1)
  sampleEffectDE <- rnorm(n, mean=S2N, sd=0.5)
  probeEffect <- rnorm(p, mean=0, sd=1)
  noise <- matrix(rnorm(p*n, mean=0, sd=1), nrow=p, ncol=n)
  noiseDE <- matrix(rnorm(p*n, mean=0, sd=1), nrow=p, ncol=n)
  M <- outer(probeEffect, sampleEffect, "+") + noise
  M2 <- outer(probeEffect, sampleEffectDE, "+") + noiseDE
  M[1:sizeDEgs, 1:group.n] <- M2[1:sizeDEgs, 1:group.n]
  
  rownames(M) <- paste0("g", 1:nrow(M))
  geneSets <- list(H1GeneSet=paste0("g", 1:(gs.sz)),
                   H0GeneSet=paste0("g", (gs.sz+1):(2*gs.sz)))

  es.gsva <- gsva(M, geneSets, verbose=FALSE, parallel.sz=1)$es.obs
  es.ss <- gsva(M, geneSets, method="ssgsea", verbose=FALSE, parallel.sz=1)
  es.z <- gsva(M, geneSets, method="zscore", verbose=FALSE, parallel.sz=1)
  es.plage <- gsva(M, geneSets, method="plage", verbose=FALSE, parallel.sz=1)
        
  h1.gsva.pval <- t.test(es.gsva["H1GeneSet", 1:group.n],es.gsva["H1GeneSet", (group.n+1):n])$p.value
  h1.ssgsea.pval <- t.test(es.ss["H1GeneSet", 1:group.n],es.ss["H1GeneSet", (group.n+1):n])$p.value
  h1.zscore.pval <- t.test(es.z["H1GeneSet", 1:group.n],es.z["H1GeneSet", (group.n+1):n])$p.value
  h1.plage.pval <- t.test(es.plage["H1GeneSet", 1:group.n],es.plage["H1GeneSet", (group.n+1):n])$p.value

  h0.gsva.pval <- t.test(es.gsva["H0GeneSet", 1:group.n],es.gsva["H0GeneSet", (group.n+1):n])$p.value
  h0.ssgsea.pval <- t.test(es.ss["H0GeneSet", 1:group.n],es.ss["H0GeneSet", (group.n+1):n])$p.value
  h0.zscore.pval <- t.test(es.z["H0GeneSet", 1:group.n],es.z["H0GeneSet", (group.n+1):n])$p.value
  h0.plage.pval <- t.test(es.plage["H0GeneSet", 1:group.n],es.plage["H0GeneSet", (group.n+1):n])$p.value
    
  out <- c(h1.gsva.pval, h1.ssgsea.pval, h1.zscore.pval, h1.plage.pval,
           h0.gsva.pval, h0.ssgsea.pval, h0.zscore.pval, h0.plage.pval)

  return(out)
}

estPwrTypIerr <- function(pvals, alpha=0.05) {
  #this code was taken from https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf
  N <- ncol(pvals)
  out <- c(1 - sum(pvals[1, ] > alpha)/N, 1 - sum(pvals[2, ] > alpha)/N,1 - sum(pvals[3, ] > alpha)/N, 1 - sum(pvals[4, ] > alpha)/N,
    sum(pvals[5, ] <= alpha)/N, sum(pvals[6, ] <= alpha)/N, sum(pvals[7, ] <= alpha)/N, sum(pvals[8, ] <= alpha)/N)
  return(out)
  
}

sum_power_type1error <- function() {  
  #this code was taken from https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf
  set.seed(1)
  exp1 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
                estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
                estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
                estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=0.5, fracDEgs=0.5))))
  
  exp2 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
                estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
                estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
                estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=1.0, fracDEgs=0.5))))
  
  exp3 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
                estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
                estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
                estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=0.5, fracDEgs=0.8))))
  
  exp4 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
                estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
                estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
                estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=1.0, fracDEgs=0.8))))
}

