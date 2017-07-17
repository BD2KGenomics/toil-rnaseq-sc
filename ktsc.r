message("ktsc: Run SC3 on Kallisto equivalence count output (will switch to transcript abundances in future version)")

ktsc <- function(ks = 2:4, idir = ".", odir = ".", debug = FALSE) {
    message("ktsc started.")
    
    message("ktsc loading kallisto output...")
    tsv <- scan(file = getPath(idir, "matrix.tsv"), what = list(row_ec = 0, col_cell = 0, val = 0))
    cells <- scan(file = getPath(idir, "matrix.cells"), what = "string")
    
    message("ktsc creating sparse matrix...")
    sparse <- sparseMatrix(i = tsv$row_ec, j = as.vector(tsv$col_cell), x = tsv$val, index1 = FALSE)
    colnames(sparse) <- cells
    rownames(sparse) <- rownames(sparse, do.NULL = FALSE, prefix = "ec_")
    
    if (debug) {
        message("ktsc DEBUG trimming cells")
        sparse <- sparse[,1:100]
    }
    
    rm(tsv)
    rm(cells)
    
    message("ktsc creating temporary SCESet...")
    tempset <- newSCESet(exprsData = sparse)
    tempset <- calculateQCMetrics(tempset)
    tempset <- sc3_prepare(tempset,ks = ks,n_cores = 1,gene_filter = TRUE)
    
    message("ktsc applying SC3 gene filters...")
    keep <- tempset@featureData@data$sc3_gene_filter
    sparse <- sparse[keep,]
    
    rm(tempset)
    rm(keep)
    
    message("ktsc removing cells with zero variance...")
    vars <- apply(sparse, 2, var)
    varies=(vars != 0)
    sparse <- sparse[,varies]
    
    message("ktsc creating new SCESet...")
    sceset <- newSCESet(exprsData = sparse)
    sceset <- calculateQCMetrics(sceset)
    
    message("ktsc running sc3...")
    sceset <- sc3(sceset, ks = ks, n_cores = 1, biology = TRUE, k_estimator = TRUE,  gene_filter = FALSE)
    
    message("ktsc exporting excel...")
    sc3_export_results_xls(sceset, filename = getPath(odir, "sc3_excel.xls"))
    
    message("ktsc exporting plots...")
    for (k in ks) {
        pngPlot(odir, pngFilename("sc3_consensus_clustering", k), function() sc3_plot_consensus(sceset, k = k))
        pngPlot(odir, pngFilename("sc3_silhouette",k), function() sc3_plot_silhouette(sceset, k = k))
    }
        
    message("ktsc finished.")
    return(sceset)
}

getPath <- function(dir, basename) {
    if (dir != ".") {
        filename <- file.path(dir, basename)
    } else {
        filename <- basename
    }
    return(filename)
}

pngFilename <- function(plotname, k) {
    return(paste(plotname,"_with_k_",k,".png",sep=""))
}

pngPlot <- function(odir, filename, plotf) {
    png(filename = getPath(odir, filename))
    plotf()
    dev.off()
}

args <- commandArgs(TRUE)
if (length(args) != 6) {
    message("Usage: Rscript path/to/script --args kmin kmax idir odir debug")
    message("kmin: int, min # of ks.")
    message("kmax: int, max # of ks.")
    message("idir: input dir (location of kallisto output (use . to export in script dir)")
    message("odir: path to export plots to (use . to export in script dir)")
    message("debug: typically, set this to \"FALSE\". If the dataset has more than 100 cells, set this to \"TRUE\" to only process the first 100 cells. If the dataset has fewer than 100 cells, debug mode will probably crash.")
} else {
    message("Loading dependencies...")
    library(SC3)
    library(scater)
    library(Matrix)
    message(args[1])
    silence_return = ktsc(ks = args[2]:args[3], idir = args[4], odir = args[5], debug = (args[6] == "TRUE"))
}
