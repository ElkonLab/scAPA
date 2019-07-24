Output
================

List of outputs
===============

The files are created in the directory scAPA/outs in the directory of the analyzed BAM files. If intronic peaks are analyzed as well, then for most files, an analogous file is created as well for the intron analysis. E.g, summary.UTR.txt is created for summarising UTR analysis, and summary.Int.txt is created for summarising intronic analysis.

Text files
----------

-   **summary.UTR.txt** summary of the 3'UTR peak analysis.
-   **ThreeUTR.peaks.txt (Intron.peaks.tt)** The 3'UTR peaks that passed filtering. The file contains the peak ID, gene symbol, ensemble ID, and genomic location.
-   **APA.events.txt** For 3'UTRs with more than one peak, gives the p-value, FDR corrected q-value of APA event tested across the clusters. Also given is the proximal PUI index for each cell cluster for each 3'UTR.
-   **UTRs.with.multiple.peaks.txt** a file containing the results of testing the differential usage of individual peaks across clusters. Only peaks from 3'UTR that came up significant in the analysis and had more than 2 peaks are analyzed. The file contains tables for each 3'UTR. The p-value and FDR q-value for each peak (from chi-square test for goodness of fit) are given. The peak usage index (PUI) for each cluster is given. Higher PUI in a cluster means higher usage of the peak in the cluster.
-   **Mean.Cell.PPUI.txt** gives the mean proximal PUI index for the cells analyzed (cells in the cluster annotation files.)

R files
-------

These are R object to be viewed with scAPA package can be used for viewing the results or for downstream analysis. \* **Peaks.RDS** an R scAPAList object. The object contains data regarding the peaks identified before filtering such as counts in the different clusters. See for more details. To read in R, load scAPA package and use R function readRDS(object = "Peaks.RDS")

-   **Results.RDS (Results.int.RDS)** an R scAPAList object. The object contains data regarding APA analysis such as p-values, proximal (intronic) PUI, cell's mean proximal PUI. To read in R, load scAPA package and use R function readRDS(object = "Reasults.RDS")

Plots
-----

-   **Proximal.PUI.ECDF.pdf** a plot of the empirical cumulative distribution of the proximal peak usage index for all clusters analyzed. The distribution is calculated for the 3'UTR with significant APA events.
