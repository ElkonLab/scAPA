require(purrr)

extract_peak_info <- function(x) {
  x <- as.data.frame(x)
  x[,1] <- as.character(x[,1])
  x$UTR_ID <- character(nrow(x))
  x$Peak_Index <- character(nrow(x))
  split_peakID <- function(y) {
    utr.gene <- strsplit(y[1], split = "_")[[1]][1]
    utr.index <- strsplit(y[1], split = "_")[[1]][2]
    y["UTR_ID"] <- paste0(utr.gene,"_", utr.index)
    y["Peak_Index"] <- strsplit(y[1], split = "_")[[1]][3]
    y
  }
  x <- t(apply(X = x, 1, FUN = split_peakID))
  x <- as.data.frame(x)
  x$Peak_Index <- as.numeric(as.character(x$Peak_Index))
  x$UTR_ID <- as.factor(x$UTR_ID)
  x <- dplyr::group_by(x, UTR_ID)
  x <- dplyr::mutate(x, rpeak = dplyr::dense_rank(x = Peak_Index))
  x$Peak_Index <- x$rpeak
  x <- dplyr::select(x, -c(rpeak))
  as.data.frame(x)
}
# Define the class --------------------------------------------------------
#' Single cell RNA-seq 3'UTR or introns peak data.
#' 
#' A list-based S4 class for storing read counts and associated information
#' from single-cell RNA sequencing experiments. 
#' For downstream Alternative Polyadenylation analysis.
#' 
#' @slot .cells.counts a data.frame, first column "Peak_ID", 
#' is the ID of the peak.
#'  Other columns (numeric) are cell barcodes, cells - read counts.
#' @slot .clus.counts If available, the sum of the counts for each 
#' peak ID from the clusters to be analyzed.
#' Data.frame, the first column "Peak_ID", is the ID of the peak.
#' Other columns are the clusters, cells have the sum of the counts 
#' for each cluster for each peak ID.
#' @slot .cluster.anot A Data.frame, the first column is the cell barcode 
#' as it appears in the column names of cells.counts.
#'  Second column is the cluster name corresponding to each cell barcode.
#' @slot .row.Data Data.frame, details regarding the peaks (e.g genomic 
#' location, gene ID)
#' @slot .down.seq Data.frame, the first column is the peak ID, 
#' the second is the genomic sequence 200 nt downstream to it.
#'  For filtering internal priming suspected peaks.
setClass("scAPAList",
         slots = list(cells.counts = "data.frame", clus.counts = "data.frame",
                      cluster.anot = "data.frame",
                      row.Data = "data.frame",
                      norm = "list",
                      down.seq = "data.frame"))


# Setting scAPAList obect -------------------------------------------------
#' scAPAList Constructor
#' 
#' A list-based S4 class for storing read counts and associated information
#' from single-cell RNA sequencing experiments. For downstream Alternative 
#' Polyadenylation analysis.
#'
#' @param .cells.counts a data.frame, first column "Peak_ID", is the 
#' ID of the peak.
#'  Other columns (numeric) are cell barcodes, cells - read counts.
#' @param .clus.counts If available, the sum of the counts for each 
#' peak ID from the clusters to be analyzed.
#' Data.frame, the first column "Peak_ID", is the ID of the peak.
#' Other columns are the clusters, cells have the sum of the counts for 
#' each cluster for each peak ID.
#' @param .cluster.anot A Data.frame, the first column is the cell barcode 
#' as it appears in the column names of cells.counts.
#'  Second column is the cluster name corresponding to each cell barcode.
#' @param .row.Data Data.frame, details regarding the peaks (e.g genomic
#'  location, gene ID)
#' @param .down.seq Data.frame, the first column is the peak ID, the second
#'  is the genomic sequence 200 nt downstream to it.
#'  For filtering internal priming suspected peaks.
set_scAPAList <- function(.cells.counts = data.frame(), 
                          .clus.counts = data.frame(),
                          .cluster.anot= data.frame(), 
                          .row.Data = data.frame(),
                          .down.seq = data.frame()){
  if((nrow(.cells.counts) == 0 | nrow(.cluster.anot) == 0) & nrow(.clus.counts) == 0){
    stop("Provide either nonempty .clus.counts and .cluster.anot,",
         "or a non empty .cells.counts")
  }
  if((nrow(.clus.counts) > 0) & (nrow(.cells.counts) == 0 | nrow(.cluster.anot) == 0)) {
    out <- methods::new("scAPAList", cells.counts = .cells.counts, 
                        clus.counts = .clus.counts,
                        cluster.anot = .cluster.anot, 
                        row.Data = .row.Data,
                        down.seq = .down.seq)
  } else{
.cells.counts$Peak_ID <- as.character(.cells.counts$Peak_ID)
.row.Data_t <- .row.Data[,c(4,1,2,3,5,6)]
list.tointersect <- list(.cells.counts, .clus.counts, .row.Data_t,
                         .down.seq)

existing.data <- lapply(list.tointersect, function(x){nrow(x) > 0})
list.tointersect <- list.tointersect[unlist(existing.data)]
list.tointersect <- lapply(list.tointersect, function(x){x[,1]})
peaks <- purrr::reduce(.f = intersect, .x = list.tointersect)
if(nrow(.cells.counts) > 0 & nrow(.cluster.anot) > 0){
cells <- intersect(colnames(.cells.counts)[-1], .cluster.anot[,1])
.cells.counts<- .cells.counts[,c("Peak_ID", cells)]
.cluster.anot <- .cluster.anot[match(cells, .cluster.anot[,1]),]
}

match_peaks <- function(y) {
  if(ncol(y) > 0) y[match(x = peaks, table = y[,1]),]
  else y
}
.cells.counts <- match_peaks(.cells.counts)
.clus.counts <- match_peaks(.clus.counts)
.row.Data_t <- match_peaks(.row.Data_t)
.down.seq <- match_peaks(.down.seq)
.row.Data_t <- .row.Data_t[,c(2:4,1,5,6)]

out <- methods::new("scAPAList", cells.counts = .cells.counts,
                    clus.counts = .clus.counts,
      cluster.anot = .cluster.anot, row.Data = .row.Data, 
      down.seq = .down.seq)
}
out
}

# subsetting --------------------------------------------------------------

setMethod("[",
          c(x = "scAPAList"),
          function(x,i,j,drop="missing") {
            if(is.character(i)){
             i <- ifelse(test = nrow(x@cells.counts) > 0,
                         yes = which(x@cells.counts[,1] %in% i),
                         no =  which(x@clus.counts[,1] %in% i))}
            .sbcells.counts <- x@cells.counts[i, ]
            .sbclus.counts <- x@clus.counts[i, ]
            .subrow.Data <- x@row.Data[i,]
            .subdown.seq <- x@down.seq[i,]
            .subnorm <- lapply(x@norm, FUN = function(x){x[i,]})
            methods::new("scAPAList", cells.counts =.sbcells.counts, 
                         clus.counts = .sbclus.counts,
                row.Data = .subrow.Data, down.seq = .subdown.seq, 
                norm = .subnorm, cluster.anot= x@cluster.anot)
          })


# Viewing functions -------------------------------------------------------

setMethod("show",
          c(object = "scAPAList"),
          function(object) {
            if(nrow(object@clus.counts) > 0){
            cat("an scAPA object\n\n# peaks:\t", nrow(object@clus.counts),
                "\n# cells:\t", ncol(object@cells.counts) -1,"\nclusters:\t",
                colnames(object@clus.counts)[-1], "\n\n")
            print(head(object@clus.counts))
            } else {
              cat("an scAPA object\n\n# peaks:\t", nrow(object@cells.counts),
                  "\n# cells:\t", ncol(object@cells.counts) -1,"\nclusters:\t",
                  unique(as.character(object@cluster.anot[,2])), "\n\n",
                  "Clusters counts were not calculated yet.\n")
            }
          })

setMethod("dim",
          c(x = "scAPAList"),
          function(x) {
            dim(x@clus.counts)
          })

setMethod("head",
          c(x = "scAPAList"),
          function(x) {
            show(object = x)
          })


# Calculating clusters' counts --------------------------------------------
#' Calculates clusters' peaks counts 
#' 
#' Sums the reads from each peak from all the cells belonging to the same cluster.
#'Returns an scAPAList where the cluster counts are calculated and assigned to the clus.counts slot.
#'
#'@param exp.mat an scAPAList object \linkS4class{scAPAList}. Must have nonempty .cluster.anot and .cells.counts slots.
setGeneric("calc_clusters_counts", function(exp.mat,cluster.anot){

  tidy.exp.mat <- t(exp.mat[,-1])
  peak.id <- exp.mat[,1]
  colnames(tidy.exp.mat) <- exp.mat$Peak_ID
  tidy.exp.mat <- cbind.data.frame(Cell_BC = rownames(tidy.exp.mat),
                                   tidy.exp.mat)
  colnames(cluster.anot) <- c("cells", "cluster")
  tidy.exp.mat <- merge(x = cluster.anot, y= tidy.exp.mat, by.x = 1,
                        by.y  = 1)
  clusters.counts.tidy <- dplyr::group_by(tidy.exp.mat[,-1], cluster )
  clusters.counts.tidy <- dplyr::summarise_all(clusters.counts.tidy, sum)
  clusters.counts <- t(clusters.counts.tidy[,-1])
  colnames(clusters.counts) <- as.data.frame(clusters.counts.tidy[,1])[,1]
  clusters.counts <- cbind.data.frame(Peak_ID = peak.id,
                                     clusters.counts)
  rownames(clusters.counts) <- NULL
  clusters.counts

})


setMethod("calc_clusters_counts", c(exp.mat = "scAPAList"), function(exp.mat){
            x <- exp.mat
            .exp.mat <- x@cells.counts
            .cluster.anot <- x@cluster.anot
            x@clus.counts <- calc_clusters_counts(exp.mat = .exp.mat,
                                                  cluster.anot=.cluster.anot)
            x
          })


# CPM normalization -------------------------------------------------------
#'Calculates the CPMs of the cluster's peaks.
#'
#'Calculates the Counts per Milion (CPMs) of the closter's peaks. 
#'
#'@param x An scAPAList with a nonempty clust.count slot
#'
#'Normalized matrix is assigned to x@norm$CPM
setGeneric("calc_cpm", function(x){
  lib_size = colSums(x)
  norm <- t((t(x) * 10^6)/lib_size)
})

setMethod("calc_cpm",
          c(x = "scAPAList"),
          function(x){
            .x <- data.matrix(x@clus.counts[,-1])
            CPMmat <- calc_cpm(x = .x)
            x@norm$CPM <- CPMmat
            x
          })

# Internal priming -------------------------------------------------------
#'Filter peaks that may stem from internal priming.
#'
#'Filter out peaks that may stem from internal priming.
#'If a sequance specified by int.priming.seq, is found in an interval
#' btween 'left' and 'right' nt 
#'downstream the peak's 3' edge, the peak will be filtered out.
#'Returns a filtered scAPAList.
#'@param x An scAPAList with a nonempty down.seq slot
#'@param int.priming.seq Charechter, the sequance to filter for. 
#'Defult for 3'UTRs "AAAAAAAA".
#'@param left numeric, defult is 10
#'@param right numeric defult is 140

setGeneric("filter_IP", function(x, int.priming.seq, left, right){
  standardGeneric("filter_IP")
})

# Annotate ----------------------------------------------------------------
setGeneric("annotate", function(x, org){
  standardGeneric("annotate")
})
setMethod("annotate",
          c(x = "scAPAList"),
          function(x, org){
  if(org == "Mm") annot <- mmanots
  if(org == "Hs") annot <- hsanots
  order <- x@row.Data[,4]
  x@row.Data$Ensemble_ID <- gsub(x = x@row.Data[,4], 
                                 pattern = "_\\d*_\\d*$",
                                 replacement = "")
  x@row.Data <- merge(x@row.Data, annot, by = "Ensemble_ID")
  x@row.Data <- x@row.Data[,c(8,1,5,2:4,7)]
  x@row.Data <- x@row.Data[match(order, x@row.Data[,3]),]
  colnames(x@row.Data)[3:6] <- c("Peak_ID", "Chromosome",
                                 "Peak_Start", "Peak_End")
  x
})

# Density plot ------------------------------------------------------------
setGeneric("plot_seq_pos_density", 
           function(x, int.priming.seq = "AAAAAAAA"){
  standardGeneric("plot_seq_pos_density")
})
setMethod("plot_seq_pos_density",
          c(x = "scAPAList"),
          function(x, int.priming.seq){
            .down.seq.data <- split(x = as.character(x@down.seq[,2]),
                                    f = x@down.seq[,1])
            title <- paste0("Histogram of ",int.priming.seq, " positions")
            int.priming.seq.neg <- gsub(x = int.priming.seq, pattern = "A",
                                        replacement = "T")
            int.priming.seq.pos <- paste0(int.priming.seq, "|",
                                          int.priming.seq.neg)
            list_pos <- function(x, .int.priming.seq = int.priming.seq){
              gregexpr(pattern  = .int.priming.seq, text = x)
            }
            pos <- lapply(FUN = list_pos, X = .down.seq.data)
            pos <- lapply(X = pos, FUN = unlist)
            keep <- lapply(X = pos, FUN = function(x){all(x > 0)})
            keep <- unlist(keep)
            pos <- pos[keep]
            pos <- unlist(pos)
            pos <- data.frame(id = names(pos), pos = pos)
            hist(pos[,2], breaks = 10, 
                 xlab = "Distance to peaks' 3' edge (nt)", main = title )
            # g <- ggplot2::ggplot(pos, ggplot2::aes(x = pos, color = pos))
            # g <- g + ggplot2::geom_density()
            # g <- g + ggplot2::xlab("Distance to peaks' 3' edge (nt)")
            # g <- g + ggplot2::theme_bw()
            # g <- g + ggplot2::xlim(c(1,150))
            #  print(g)
            })

# results ----------------------------------------------------------------
#' Single cell RNA-seq 3'UTR or introns analysis results data.
#' 
#' A list-based S4 class for storing the results of alternative polyadenilation (APA) analysis of
#' single-cell RNA sequencing data.
#' 
#' @slot .clus.counts a list. Each element is a data.frame corresponding to a 3'UTR/Inton ID. Columns are cell clusters, raws are
#' Peaks from the 3'UTR (or intronic peak and 3'UTR counts for introns)
#'  Other columns (numeric) are cell barcodes, cells - read counts.
#' @slot .cells.counts  Same as .clus.counts, but columns correspond to cells.
#' @slot .pvalues A List containing matricis. Matricis contain p-values of dynamical (APA) events.
#' Each matrix correspond to a comperison btwen clusters. 
#' Rawnames are 3'UTR/Intron ID. First column is p-value, second column is FDR orreted p-value.
#' @slot peak.pvaues A List containing matricis. Matricis contain p-values of peaks of dynamical (APA) events with multiple peaks (>3).
#' Each matrix correspond to a comperison btwen clusters. 
#' Rawnames are the peak ID. First column is p-value, second column is FDR orreted p-value.
#' @slot .pAi.clus a matrix of pA index. rawnames are UTR IDs, columns are clusters.
#' @slot pAi.cells Same as .pAi.clus, but columns are cells.
#' @slot ppui.clus A matrix of proximal usage index (or )
# Class definition --------------------------------------------------------
setClass("scAPAresults",
         slots = list(clus.counts = "list", cells.counts = "list",
                      pvalues = "list",
                      peak.pvaues = "list",
                      pAi.clus = "matrix",
                      pAi.cells = "matrix",
                      ppui.clus = "matrix",
                      ppui.cells = c("matrix","numeric"),
                      down.seq = "data.frame",
                      metadata = "data.frame"))

setGeneric("set_scAPAresults", function(x, int = F, cpm = NULL, counts){
  standardGeneric("set_scAPAresults")
})

# Subsetting --------------------------------------------------------------
setMethod("[",
          c(x = "scAPAresults"),
          function(x,i,j,drop="missing") {
            if(is.character(i)) i <- which(names(results@clus.counts) %in% i)
            if(length(x@cells.counts) > 0) {
              .subcells.counts <- x@cells.counts[i]
            } else .subcells.counts <- list()
            .subclus.counts <- x@clus.counts[i]
            .subpAi.clus <- x@pAi.clus
            if(nrow(x@pAi.clus) > 0) .subpAi.clus <- x@pAi.clus[i, ]
            .subpAi.cells <- x@pAi.cells
            if(nrow(x@pAi.cells) > 0) .subpAi.cells <- x@pAi.cells[i, ]
            .subppui.clus <- x@ppui.clus
            if(nrow(x@ppui.clus) > 0) .subppui.clus <- x@ppui.clus[i,]
            .subppui.cells <- x@ppui.cells
            if(nrow(x@ppui.cells) > 0) {
              ii <- i
              names(ii) <- NULL
              .subppui.cells <- as.matrix(x@ppui.cells[ii,])
              } else subppui.cells <- matrix()
            if(nrow(x@down.seq) > 0) {
              .subdown.seq <- x@down.seq[i,]
            } else .subdown.seq <- data.frame()
            peakstokeep <- lapply(.subclus.counts, row.names)
            peakstokeep <- as.character(unlist(peakstokeep))
            .submetadata <- x@metadata[which(x@metadata[,4] %in% peakstokeep),]
            .subpvalues <- lapply(x@pvalues, FUN = function(x){x[i,]})
            methods::new("scAPAresults", cells.counts = .subcells.counts,
                         clus.counts = .subclus.counts,
                         pAi.clus = .subpAi.clus,
                         pAi.cells = .subpAi.cells, ppui.clus = .subppui.clus,
                         ppui.cells = .subppui.cells,
                         metadata = .submetadata, pvalues = .subpvalues)
          })

# Creat the object --------------------------------------------------------
setMethod(f = "set_scAPAresults", signature =  c(x = "scAPAList"),
          definition = function(x, int, cpm, counts){
            if(!int){
              peakID <- x@clus.counts[,1]
              peak.info <- extract_peak_info(x = peakID)
              to.split <- cbind.data.frame(peak.info[,-1],x@clus.counts[,-1])
              rownames(to.split) <- peak.info[,1]
              .clus.counts <- split(x = to.split[, -1], f = to.split[,1],
                                    drop = T)
              row.n <- lapply(.clus.counts, nrow)
              f.list.utr <- .clus.counts[row.n > 1]
              clus.counts <- f.list.utr
              peakstokeep <- lapply(clus.counts, row.names)
              peakstokeep <- as.character(unlist(peakstokeep))
              .metadata <- x@row.Data[which(x@row.Data[,4] %in% peakstokeep),]
              if(ncol(x@cells.counts) == 0){
                out <- methods::new("scAPAresults", clus.counts = clus.counts,
                                    cells.counts= list(), metadata = .metadata,
                                    down.seq = data.frame())
              } else{
                to.split <- cbind.data.frame(peak.info[,-1],
                                             x@cells.counts[,-1])
                rownames(to.split) <- peak.info[,1]
                .clus.counts <- split(x = to.split[, -1], f = to.split[,1],
                                      drop = T)
                row.n <- lapply(.clus.counts, nrow)
                f.list.utr <- .clus.counts[row.n > 1]
                .cells.counts <- f.list.utr
                .peak.pvaues = list()
                out <- methods::new("scAPAresults", clus.counts = clus.counts,
                                    cells.counts= .cells.counts, 
                                    metadata = .metadata,
                                    down.seq = x@down.seq, 
                                    peak.pvaues = .peak.pvaues)
              }
            }
            if(int){
              cells <- x@cells.counts
              x <- calc_cpm(x)
              total <- (rowSums(x@norm$CPM) >=1 ) & (rowSums(x@clus.counts[,-1]) > counts)
              y <- x[which((grepl(x = x@clus.counts[,1], pattern = "_Int_")) & (total)),]
              y <- calc_cpm(y)
              total <- rowSums(y@norm$CPM)
              y <- y[total > cpm,]
              introns <- cbind.data.frame(y@row.Data[,c(1:3,6)], y@clus.counts)
              introns$gene <- gsub(pattern = "_Int_.*$", replacement = "",
                                   x = introns[,5])
              if(nrow(cells) > 0){
                introns.cells <- cells[which(cells$Peak_ID %in% introns$Peak_ID),]
                introns.cells <- cbind.data.frame(introns[,c(1:4, ncol(introns))],
                                                  introns.cells)
              }
              z <- x[which(!grepl(x = x@clus.counts[,1], pattern = "_Int_")),]
              utrs <- cbind.data.frame(z@row.Data[,c(1:3,6)], z@clus.counts)
              utrs$gene <- gsub(pattern = "_.*$", replacement = "", x = utrs[,5])
              if(nrow(cells) > 0){
                utrs.cells <- cells[which(cells$Peak_ID %in% utrs$Peak_ID),]
                utrs.cells <- cbind.data.frame(utrs[,c(1:4, ncol(utrs))],
                                               utrs.cells)
              }
              int_utr <- merge(introns, utrs, by = "gene")
              if(nrow(cells) > 0){
                int_utr.cells <- merge(introns.cells, utrs.cells, 
                                       by = "gene")
              }
              int_utr <- dplyr::filter(int_utr, ifelse(Strand.x == "+",
                                                       End.x < End.y,
                                                       End.x > End.y))
              if(nrow(cells) > 0){
                int_utr.cells <- dplyr::filter(int_utr.cells, 
                                               ifelse(Strand.x == "+",
                                                                     End.x < End.y,
                                                                     End.x > End.y))
              }
              
              f <- 6 + ncol(x@clus.counts) - 1
              group_col <- c("gene", "Chr.x",  "Start.x", "End.x", "Strand.x",
                             "Peak_ID.x", 
                             paste0(colnames(y@clus.counts[-1]),".x"))
              int_utr <-  dplyr::group_by(int_utr, .dots=group_col)
              if(nrow(cells) > 0){
                group_col.cells <- c("gene", "Chr.x",  "Start.x", "End.x", 
                                     "Strand.x",
                                     "Peak_ID.x", 
                                     paste0(colnames(x@cells.counts[-1]),".x"))
                int_utr.cells <-  dplyr::group_by_at(int_utr.cells, 
                                                     .vars = group_col.cells)
              }
              
              tosum <- paste0(colnames(y@clus.counts[-1]),
                              ".y")
              int_utr <- dplyr::summarise_at(.tbl = int_utr, 
                                             .vars =tosum ,
                                             .funs = sum)
              int_utr <- as.data.frame(int_utr)
              if(nrow(cells) > 0){
                tosum.cells <- paste0(colnames(x@cells.counts[-1]),
                                      ".y")
                int_utr.cells <- dplyr::summarise_at(.tbl = int_utr.cells, 
                                                     .vars =tosum.cells, 
                                                     .funs = sum)
                int_utr.cells <- as.data.frame(int_utr.cells)
              }
              p <- split(int_utr, f = int_utr$Peak_ID.x, drop = T)
              clustnu <- ncol(x@clus.counts) -1
              dat.fan <- function(x){
                feature <- c(1,2)
                introns <- as.vector(x[,c(7:(6+clustnu))])
                utr <- as.vector(x[,c((6+clustnu+1):(6+clustnu+clustnu))]) 
                names(introns) <- colnames(y@clus.counts[-1])
                names(utr) <- colnames(y@clus.counts[-1])
                out <- rbind.data.frame(introns, utr)
                out <- cbind.data.frame(feature, out)
              }
              .clus.counts <- lapply(p, dat.fan)
              if(nrow(cells) > 0){
                p.cells <- split(int_utr.cells, f = int_utr.cells$Peak_ID.x,
                                 drop = T)
                cellsnu <- ncol(x@cells.counts) -1
                dat.fan.cells <- function(y){
                  feature <- c(1,2)
                  introns <- as.vector(y[,c(7:(6+cellsnu))])
                  utr <- as.vector(y[,c((6+cellsnu+1):(6+2*cellsnu))]) 
                  names(introns) <- colnames(x@cells.counts[-1])
                  names(utr) <- colnames(x@cells.counts[-1])
                  out <- rbind.data.frame(introns, utr)
                  out <- cbind.data.frame(feature, out)
                }
                .cells.counts <- lapply(p.cells, dat.fan.cells)
              }
              else{.cells.counts <- list()
              }
              .metadata <- int_utr[,c(2,3,4,6,1,5)]
              
              colnames(.metadata) <- c("Chr", "Start", "End", "Intron_ID",
                                       "Gene", "Strand")
              .metadata <- .metadata[.metadata$Intron_ID %in% names(.clus.counts),]
              .metadata <- .metadata[match(names(.clus.counts), 
                                           .metadata$Intron_ID),]
              .down.seq <- y@down.seq
              .down.seq <- .down.seq[.down.seq$Peak_ID %in% names(.clus.counts),]
              .down.seq <- .down.seq[match(names(.clus.counts), 
                                           .down.seq$Peak_ID),]
              out <- methods::new("scAPAresults",
                                  clus.counts = .clus.counts,
                                  cells.counts= .cells.counts,
                                  metadata = .metadata,
                                  down.seq=.down.seq)
            }
            out
          })

# Show --------------------------------------------------------------------
setMethod("show",
          c(object = "scAPAresults"),
          function(object) {
            if(length(object@cells.counts) > 0){
              cat("an scAPA results object\n\n# 3' UTR with APA: ",
                  length(object@clus.counts), "\n# cells: ",
                  ncol(object@cells.counts[[1]]),"\n clusters: ",
                  colnames(object@clus.counts[[1]][-1]), "\n",
                  "metadata: ", colnames(object@metadata), "\n")
            } else {
              cat("an scAPA results object\n\n# 3' UTR with APA:\t",
                  length(object@clus.counts), "\no cell expression data",
                  "\nclusters:\t",
                  colnames(object@clus.counts[[1]][-1]), "\n")
            }
            if(length(object@pvalues) > 0) {
              cat("# significant (FDR < 5%) APA events: ",
                  sum(object@pvalues[[1]][,2] < 0.05, na.rm = T), "\n")
              
            }
            if(nrow(object@pAi.clus) > 0) {
              cat("\nwighted mean of pA index summary:\n")
              print(summary(object@pAi.clus))
            }
            if(nrow(object@ppui.clus) > 0) {
              cat("\nproximal peak usage index summary:\n")
              print(summary(object@ppui.clus))
            }
          })

setMethod("head",
          c(x = "scAPAresults"),
          function(x) {
            show(object = x)
          })
setMethod("summary",
          c(object = "scAPAresults"),
          function(object) {
            show(object = object)
          })

# Test APA ----------------------------------------------------------------
setGeneric("test_APA", function(x, clus = "all"){
  standardGeneric("test_APA")
})
setMethod("test_APA",
          c(x = "scAPAresults"),
          function(x, clus) {
            names.pval <- paste0(clus, collapse = "_")
            x@pvalues <- c(x@pvalues, names.pval = matrix())
            test_utr <- function(z, .clus) {
              .clus = clus
              y <- z[, -c(1)]
              if(length(.clus) > 1) y <- y[, .clus]
              u <- suppressWarnings(chisq.test(y)$p.value)
            }
            ptable <- lapply(X = x@clus.counts, FUN = test_utr, 
                             .clus = clus)
            ptable <- do.call(what = rbind, args = ptable)
            ptable <- as.matrix(ptable)
            colnames(ptable) <- "pval"
            qval <- p.adjust(ptable[,1])
            ptable <- cbind(ptable, qval)
            x@pvalues$names.pval <- ptable
            names(x@pvalues) <- c(names(x@pvalues)[-length(x@pvalues)],
                                  names.pval)
            x
          })

setGeneric("test_peaks", function(x, clus = "all", sig.level = 0.05){
  standardGeneric("test_peaks")
})
setMethod(f = "test_peaks", c(x = "scAPAresults"),
          function(x, clus, sig.level) {
            if(is.null(x@pvalues[[clus]])){
              mes <- paste0("P-values for the clusters were not calculated",
                            " yet. Use test_APA first")
              stop(mes)
            }
            apapval <- x@pvalues$all
            
            names.pval <- paste0(clus, collapse = "_")
            names.pval <- paste0("Peaks_", names.pval)
            x@peak.pvaues <- c(x@peak.pvaues, names.pval = data.frame())
            test_inpeaks <- function(z, .clus) {
              .clus = clus
              y <- z[, -c(1)]
              if(length(.clus) > 1) y <- y[, .clus]
              expected <- colSums(y)/sum(y)
              calc_pui <- function(y, .psudo = 1){
                peak.counts <- y
                peak.counts[,(colSums(peak.counts) == 0)] <-  NA
                peak.counts <- peak.counts + .psudo
                utr.geo.mean = apply(peak.counts, 2, EnvStats::geoMean)
                i.peak.counts <- y + .psudo
                out <- log2(t(t(i.peak.counts)/utr.geo.mean))
                colnames(out) <- paste0(colnames(out),"_PUI" )
                out
              }
              chi.fit <- function(x) {
                suppressWarnings(chisq.test(p = expected, x = x)$p.value)
              }
              
              pval <- apply(y, 1, chi.fit)
              pui <- calc_pui(y)
              cbind(pui, pval)
            }
            sig <- x@clus.counts[which(apapval[,2] < sig.level)]
            moretwo <- lapply(sig, function(x){nrow(x) > 2})
            moretwo <- unlist(moretwo, use.names = F)
            totest <- sig[which(moretwo)]
            ptable <- lapply(X = totest, FUN = test_inpeaks, .clus = clus)
            ptable <- do.call(what = rbind, args = ptable)
            qval <- p.adjust(ptable[,ncol(ptable)])
            ptable <- cbind(ptable, qval)
            ptable <- as.data.frame(ptable)
            Peak_ID <- rownames(ptable)
            UTR_ID <- gsub(x = Peak_ID, pattern = "_\\d$", replacement = "")
            Gene_ID <- gsub(x = UTR_ID, pattern = "_\\d$", replacement = "")
            ptable <- cbind.data.frame(Gene_ID, UTR_ID, Peak_ID, ptable)
            rownames(ptable) <- NULL
            x@peak.pvaues$names.pval <- ptable
            names(x@peak.pvaues) <- c(names(x@peak.pvaues)[-length(x@peak.pvaues)],
                                      names.pval)
            x
          })
# Calculate pA index ------------------------------------------------------
setGeneric("calc_pAi_mat", function(x){
  calc_mean_pAi <- function(z) {
    peak.counts <- z[,-1]
    peak.index <- z[,1]
    utr.total.couns = colSums(peak.counts)
    frac.counts <- t((t(peak.counts))/utr.total.couns)
    colSums(frac.counts * peak.index)
  }
  mean.pA.mat <- lapply(X = x, FUN = calc_mean_pAi)
  mean.pA.mat <- do.call(what = rbind, args = mean.pA.mat)
  mean.pA.mat <- as.data.frame(mean.pA.mat)
  colnames(mean.pA.mat) <- paste0(colnames(mean.pA.mat), "_<pA>")
  data.matrix(mean.pA.mat)
})

setMethod("calc_pAi_mat",
          c(x = "scAPAresults"),
          function(x) {
            x@pAi.clus <- calc_pAi_mat(x = x@clus.counts)
            if(length(x@cells.counts) > 0){
              x@pAi.cells <- calc_pAi_mat(x = x@cells.counts)
            }
            x
          })
# Calculate PUI -----------------------------------------------------------
setGeneric("calc_p_pui_mat", function(x, psudo =1, int = FALSE){
  calc_prox_pui <- function(z, .psudo = psudo){
    peak.counts <- z[,-1]
    peak.counts[,(colSums(peak.counts) == 0)] <-  NA
    peak.counts <- peak.counts + .psudo
    utr.geo.mean = apply(peak.counts, 2, EnvStats::geoMean)
    prox.peak.counts <- z[z[, 1] == 1, -1] + .psudo
    prox.peak.counts <- as.numeric(prox.peak.counts)
    out <- log2(prox.peak.counts/utr.geo.mean)
    out
  }
  ppui.mat <- lapply(X = x, FUN = calc_prox_pui)
  ppui.mat <- do.call(what = rbind, args = ppui.mat)
  if(!int){
    colnames(ppui.mat) <- paste0(colnames(ppui.mat), "_proximal_PUI")
  }
  if(int){
    colnames(ppui.mat) <- paste0(colnames(ppui.mat), "_intronic_PUI")
  }
  ppui.mat
})
setMethod("calc_p_pui_mat",
          c(x = "scAPAresults"),
          function(x, psudo, int) {
            x@ppui.clus <- calc_p_pui_mat(x = x@clus.counts,
                                          psudo = psudo)
            if(length(x@cells.counts) > 0){
              # x@ppui.cells <- calc_p_pui_mat(x = x@cells.counts,
              #                                psudo = psudo)
              ppui.cells <- calc_p_pui_mat(x = x@cells.counts,
                                           psudo = psudo)
              t <-  apply(X = ppui.cells, MARGIN = 2, FUN = mean, na.rm = T)
              t <- as.matrix(t)
              if(!int){
                colnames(t) <- "Mean_proximal_PUI"
              }
              if(int){
                colnames(t) <- "Mean_intronic_PUI"
              }
              x@ppui.cells <- t
            }
            x
          })


# Summarise results ------------------------------------------------------

setGeneric("disply_results", function(x, org, int = FALSE){
  standardGeneric("disply_results")
})
setMethod("disply_results",
          c(x = "scAPAresults"),
          function(x, org, int){
            if(org == "Mm") annot <- mmanots
            if(org == "Hs") annot <- hsanots
            if(!int){
              out <- data.frame(UTR_ID = row.names(x@pvalues[[1]]),
                                x@pvalues[[1]])
              if(nrow(x@pAi.clus) > 0){
                out <- merge(out, x@pAi.clus, by = "row.names")
              }
              out <- merge(out, x@ppui.clus, by.x = "UTR_ID", 
                           by.y = "row.names")
              out <- merge(annot[,c(1,2,6)], out, by = "UTR_ID")
              out <- out[,-4]
              out <- out[,c(3,2,1,4:11)]
              colnames(out)[4:5] <- c("p-value", "q-value") 
            }
            if(int){
              out <- data.frame(Intron_ID = row.names(x@pvalues[[1]]), 
                                x@pvalues[[1]])
              if(nrow(x@pAi.clus) > 0){
                out <- merge(out, x@pAi.clus, by = "row.names")
              }
              out <- merge(out, x@ppui.clus, by.x = "Intron_ID", 
                           by.y = "row.names")
              out <- merge(y = out, x = x@metadata[,c(5,4,1:3,6)],
                           by = "Intron_ID")
              out <- merge(annot[,c(1,2)], out, by.x = "Ensemble_ID",
                           by.y = "Gene")
              out <- out[,-3]
              colnames(out)[7:8] <- c("p-value", "q-value") 
            }
          })

# Internal priming --------------------------------------------------------


setGeneric("find_internal_prim_seq", function(x, int.priming.seq, left, right){
  standardGeneric("find_internal_prim_seq")
})
setMethod("find_internal_prim_seq",
          c(x = "scAPAresults"),
          function(x, int.priming.seq, left, right){
            .down.seq.data <- split(x = as.character(x@down.seq[,2]),
                                    f = x@down.seq[,1])
            list_pos <- function(x, .int.priming.seq = int.priming.seq){
              gregexpr(pattern  = .int.priming.seq, text = x)
            }
            pos <- lapply(FUN = list_pos, X = .down.seq.data)
            pos <- lapply(X = pos, FUN = unlist)
            tofilter <- lapply(X = pos, FUN = function(x){any((x < right) & (x > left))})
            tofilter <- unlist(tofilter)
            tofilter <- data.frame(Peak_ID = names(tofilter), has_seq = tofilter)
            tofilter <- tofilter[match(names(x@clus.counts), tofilter[,1]),]
            rownames(tofilter) <- NULL
            tofilter
          })
# Internal primming -------------------------------------------------------
setMethod("filter_IP",
          c(x = "scAPAresults"),
          function(x, int.priming.seq, left, right){
            .down.seq.data <- split(x = as.character(x@down.seq[,2]),
                                    f = x@down.seq[,1])
            list_pos <- function(x, .int.priming.seq = int.priming.seq){
              gregexpr(pattern  = .int.priming.seq, text = x)
            }
            pos <- lapply(FUN = list_pos, X = .down.seq.data)
            pos <- lapply(X = pos, FUN = unlist)
            tofilter <- lapply(X = pos, FUN = function(x){any((x < right) & (x > left))})
            tofilter <- unlist(tofilter)
            tofilter <- data.frame(Peak_ID = names(tofilter), has_seq = tofilter)
            tofilter <- tofilter[match(names(x@clus.counts), tofilter[,1]),]
            rownames(tofilter) <- NULL
            keep.int.pr <- as.vector(!tofilter[,2])
            x <- x[which(keep.int.pr),]
            x
          })
setMethod("filter_IP",
          c(x = "scAPAList"),
          function(x, int.priming.seq, left, right){
            .down.seq.data <- split(x = as.character(x@down.seq[,2]),
                                    f = x@down.seq[,1])
            list_pos <- function(x, .int.priming.seq = int.priming.seq){
              gregexpr(pattern  = .int.priming.seq, text = x)
            }
            pos <- lapply(FUN = list_pos, X = .down.seq.data)
            pos <- lapply(X = pos, FUN = unlist)
            tofilter <- lapply(X = pos, 
                               FUN = function(x){
                                 any((x < right) & (x > left))})
            tofilter <- unlist(tofilter)
            tofilter <- data.frame(Peak_ID = names(tofilter), 
                                   has_seq = tofilter)
            tofilter <- tofilter[match(x@clus.counts[,1], 
                                       tofilter[,1]),]
            rownames(tofilter) <- NULL
            keep.int.pr <- as.vector(!tofilter[,2])
            x <- x[which(keep.int.pr),]
            x
          })
# Write peaks pvalues -----------------------------------------------------
setGeneric("write.peaks", function(.x, .f){
  standardGeneric("write.peaks")
})

setMethod("write.peaks",
          c(.x = "scAPAresults"),
          function(.x, .f){
            if (org == "Mm") ano <- mmanots
            if(org == "Hs") ano <- hsanots
            peakspval <- merge(.x@peak.pvaues$Peaks_all, ano[,c(1,2)],
                               by.x = "Gene_ID",
                               by.y = "Ensemble_ID")
            peakspval <- unique(peakspval)
            peakspval$di <- paste0("Gene Symbol: ", peakspval$Gene_Symbol,
                                   ", Gene ID ",
                                   peakspval$Gene_ID, ", UTR ID ", 
                                   peakspval$UTR_ID)
            colnames(peakspval) <- gsub(pattern = "pval", replacement = "p-value",
                                        x = colnames(peakspval))
            colnames(peakspval) <- gsub(pattern = "qval",
                                        replacement = "q-value",
                                        x = colnames(peakspval))
            peakspval <- split(peakspval[,c(3:8)], f = peakspval$di, drop = T)
            for(i in 1:length(peakspval)){
              write(x = names(peakspval)[i], file = .f, append = T)
              write.table(x = peakspval[[i]], file = .f, append = T, quote = F,
                          sep = "\t", 
                          row.names = F, col.names = T)
            }
          })

# anotate -----------------------------------------------------------------

setGeneric("annotate_results", function(x, org){
  standardGeneric("annotate_results")
})
setMethod("annotate_results",
          c(x = "scAPAresults"),
          function(x, org){
            if(org == "Mm") annot <- mmanots
            if(org == "Hs") annot <- hsanots
            order <- x@metadata[,4]
            x@metadata$Ensemble_ID <- gsub(x = x@metadata[,4], 
                                           pattern = "_Int_\\d*_\\d*$", 
                                           replacement = "")
            x@metadata <- merge(y = x@metadata, x = annot[,1:2], 
                                by = "Ensemble_ID")
            x@metadata  <- x@metadata[,c(2,1,6,3:5,8)]
            x
          })

