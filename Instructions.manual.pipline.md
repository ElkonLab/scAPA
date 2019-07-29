Instructions for Manual Implementation of Our Pipeline
================

This is a file containing instruction for manually implementing the first two steps of our pipeline (see for descriptions of the steps). **Note that scAPA.shell.script.R automatically performs all of the five steps of the pipline**.

Prerequisites
-------------

1.  Fasta and chromosome length files for human (hg19) or (and) mouse (mm10). Can be downloaded from [UCSC website](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/)
2.  [Samtools](http://www.htslib.org/download/) version 0.1.19-96b5f2294a or above.
3.  [Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) version v2.27.1-17-gf59480f or above.
4.  [Homer](http://homer.ucsd.edu/homer/introduction/install.html) version v4.8.3 or above. Uses homer's "MakeTagdirectory" and "findPeaks".
5.  [UMI\_tools](https://github.com/CGATOxford/UMI-tools/blob/master/doc/QUICK_START.md)
6.  [Drop-seq\_tools-2.2.0](https://github.com/broadinstitute/Drop-seq/releases/tag/v2.2.0)
7.  R version 3.5.3, or above.

Downloading the example files
-----------------------------

The files used for this example are from [Lukassen et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6132189/). The original fastq files can be obtained from [GSE104556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556). **Please notice that the BAM files here have been downsampled, to contain 25% of the reads.** The result will differ from the results obtained from analyzing the full BAM files. Downlad the folder down.sampled.spermatogenesis from [this link](https://drive.google.com/open?id=1xK7lR2ECfJ-Cjb1f4bYnaA5JjYtdqrGA).

Install scAPA package
---------------------

Here we make use of scAPA package. Install and load:

``` r
install.packages("devtools")
require(devtools)
devtools::install_github("ElkonLab/scAPA/scAPA") 
require(scAPA)
```

### 1. Defining 3'UTR peaks

#### a. PCR duplicates removal

-   PCR duplicates are removed using UMI-tools dedup. As UMI tools dedup requires that each line in the BAM file has a molecular barcode tag, Drop-seq tools is first used to filter the BAMs, leaving only reads for which cell ranger counts produced corrected molecular barcode tag.

<!-- -->

    sh -c 'for sample in down.sampled.SRR6129050 down.sampled.SRR6129051 ; do FilterBam TAG_RETAIN=UB I=${sample}.bam O=UB.${sample}.bam; done'

-   UMI tools requires indexed bams:

<!-- -->

    sh -c 'for sample in down.sampled.SRR6129050 down.sampled.SRR6129051 ; do samtools index UB.${sample}.bam; done'

-   Then UMI tools is ran with “method=unique” so that cellranger corrected molecular barcodes are used:

<!-- -->

    sh -c 'for sample in down.sampled.SRR6129050 down.sampled.SRR6129051 ;do umi_tools dedup -I UB.${sample}.bam -S dedup.${sample}.bam --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB; done'

#### b. Peak detection

-   Homer is used to create a tag directory (Tagdirectory) from the PCR duplicate removed BAMs:

<!-- -->

    makeTagDirectory Tagdirectory dedup.down.sampled.SRR6129050.bamdedup.down.sampled.SRR6129051.bam

-   Homer findPeaks is used to identify peaks. By default, findPeaks adjusts reads to the center of their fragment. To avoid this, fragLength is set to the average read length. To find peaks of variable width, Homer is set to find peaks of width 50nt and a minimum distance of 1 nt between peaks.

<!-- -->

    findPeaks Tagdirectory -size 50 -fragLength 100 -minDist 1 -strand separate -o Peakfile

-   The function merge\_peaks from scAPA package uses bedtools to merge peaks less than 100 nt apart:

``` r
merge_peaks(bedtools.path = bedtools.path, peaks.file = "Peakfile", path = "./")
```

bedtools.path is the path to bedtools

-   scAPA function intersect\_peaks uses bedtools to intersect peaks file with a 3' UTR bed file to create a bed file of the 3’ UTR peaks. The peaks are annotated according to their 3’ UTR and their position within it.

``` r
peaks.bed <- intersect_peaks(org = "Mm", bed.name = "./merge.peakfile.bed",
                             path = "", bedtools.path = bedtools.path, 
                             introns = F)
write.bed(.x = peaks.bed, f = "./peaks.bed")
```

#### c. Separating Peaks with bimodal UMI counts distribution

Adjacent pA sites may result in a single peak. To detect and separate such peaks the R package mclust is used to fit a Gaussian finite mixture model with two components to the UMI counts distribution in the interval of each peak. The input to mclust, the UMI counts distribution, is prepared as follows: \* To detect reads from the union of all reads from all samples, the duplicate-removed BAM files are merged.

    samtools merge merged.bam dedup.down.sampled.SRR6129050.bam dedup.down.sampled.SRR6129051.bam

-   Two BEDGRAPHS files are produced from this BAM, one for the plus strand and one for the minus strand, using bedtools genomcove

<!-- -->

    bedtools genomecov -ibam merged.bam -bg -strand + | awk 'BEGIN {OFS = " "}{print $1, $2, $3, $4, ".", "+"}' > merged.wig

    bedtools genomecov -ibam merged.bam -bg -strand - | awk 'BEGIN {OFS = " "}{print $1, $2, $3, $4, ".", "-"}' > merged.wig

-   The BEDGRAPHS files are converted into a BED format and bedtools intersect is used to intersect them with the peaks' BED file.

<!-- -->

    bedtools intersect -s -wb -b peaks.bed -a merged.wig > intersected.wig

-   The intersected file is read in R and converted to a list such that each element of the list, corresponding to a specific peak, is a numeric vector whose values represent read coverage observed across the peak.

``` r
peaks.wig <- read.delim(file = "intersected.wig", header = F)
peaks.wig <- split(x = peaks.wig, f = peaks.wig$V10, drop = T)
```

-   mclust is used to fit an equal variance Gaussian finite mixture model with two components to (G=2, modelNames="E") to each list element. If the predicted means of the two fitted Gaussian components are separated by more than three standard deviations and at least 75 nt, the peak is split into two, according to mclust classification. The peak's bed is edited accordingly (Correct peak index). All this is done using the function creat\_mclus, from scAPA package.

``` r
bed <- plyr::rbind.fill(lapply(1:length(peaks.wig),
                                           FUN = scAPA::creat_mclus))
```

### 2. Quantifying the usage of each peak in each cell cluster

This step uses featureCounts to count the reads that overlap each peak in each cluster ("cell type"). A separate BAM file for each cell cluster is generated. First, Drop-seq tools is used to split the reads in each sample BAM into separate BAMs that correspond to the different clusters. This is done using FilterBAMByTag together with text files, where each file contains the cell barcodes that belong to a certain cell cluster, from a certain sample.

-   The following is used to write the text files.

``` r
    annotation.files <- c("clusters_down.sampled.SRR6129050.txt",
                          "clusters_down.sampled.SRR6129051.txt")
    for (j in annotation.files)) {
    annotations <- read.delim(file = j, header = F)
    list.cells[, 1] <- as.character(annotations[, 1])
    list.cells <- split(x = list.cells, f = list.cells[, 2], drop = T)
     for (i in 1:length(list.cells)) {
       cluster <- names(list.cells)[i]
            cluster.cells <- as.data.frame(list.cells[[i]][, 1])
            tagValuefile <- paste0(cluster, "_", sample)
            write.bed(.x = cluster.cells, f = tagValuefile)
            scAPA::write.bed(.x = cluster.cells, f = tagValuefile)
    }
```

The object annotation.files contains the paths to the cell cluster annotation files.

-   Now use FilterBamByTag:

<!-- -->


    FilterBamByTag TAG=CB I=dedup.down.sampled.SRR6129050.bam O=down.sampled.SRR6129050_ES.bam TAG_VALUES_FILE=tag_value_file.down.sampled.SRR6129050_ES.txt
     
    FilterBamByTag TAG=CB I=dedup.down.sampled.SRR6129050.bam O=down.sampled.SRR6129050_RS.bam TAG_VALUES_FILE=tag_value_file.down.sampled.SRR6129050_RS.txt
     
    FilterBamByTag TAG=CB I=dedup.down.sampled.SRR6129050.bam O=down.sampled.SRR6129050_SC.bam TAG_VALUES_FILE=tag_value_file.down.sampled.SRR6129050_SC.txt
     
    FilterBamByTag TAG=CB I=dedup.down.sampled.SRR6129051.bam O=down.sampled.SRR6129051_ES.bam TAG_VALUES_FILE=tag_value_file.down.sampled.SRR6129051_ES.txt
     
    FilterBamByTag TAG=CB I=dedup.down.sampled.SRR6129051.bam O=down.sampled.SRR6129051_RS.bam TAG_VALUES_FILE=tag_value_file.down.sampled.SRR6129051_RS.txt
     
    FilterBamByTag TAG=CB I=dedup.down.sampled.SRR6129051.bam O=down.sampled.SRR6129051_SC.bam TAG_VALUES_FILE=tag_value_file.down.sampled.SRR6129051_SC.txt

-   Next, bams of the same cluster are merged:

<!-- -->



    samtools merge ES.bam down.sampled.SRR6129050_ES.bam down.sampled.SRR6129051_ES.bam

    samtools merge RS.bam down.sampled.SRR6129050_RS.bam down.sampled.SRR6129051_RS.bam

    samtools merge SC.bam down.sampled.SRR6129050_SC.bam down.sampled.SRR6129051_SC.bam

We end up with 3 BAMS: ES.bam, RS.bam, and SC.bam.

Rsubread package featureCounts function is used.

-   The annotation file is a SAF data.frame edited from the peaks bed file.

``` r
utr.saf <- bed[, c(4, 1, 2, 3, 6)]
colnames(utr.saf) <- c("GeneID", "Chr", "Start", "End", "Strand")
```

largestOverlap = True is specified so that reads spanning two peaks are counted according to their largest overlap.

``` r
  counts <- Rsubread::featureCounts(files = c("ES.bam", "RS.bam", "CS.bam"), isGTFAnnotationFile = F,
                          strandSpecific = 1, annot.ext = utr.saf,
                          largestOverlap = T)
```

-   The object *counts* is edited for downstream analysis:

``` r
  co <- cbind.data.frame(rownames(counts$counts), counts$counts)
  colnames(co) <- c("Peak_ID", "ES", "RS", "SC")
```

-   We also produce an object storing the 200 nt downstream sequences for each peak, using scAPA function read\_down.seq.

``` r
  aseq <- read_down.seq(saf = utr.saf, 
                            char.length.path = char.length.path,
                            fasta.path = fasta.path, chr.modify = T)
  aseq <- aseq[,c(4, 6)]
```

Here char.length.path is the path to the chromosomes length file, and fasta.path is the path to the fasta file. We also produce an object with information regarding the genomic location of the peaks:

``` r
row.Data <- counts$annotation[,c(2,3,4,1,6,5)]
```

-   Finnaly, we use set\_scAPAList from scAPA to creat a scAPAList object for downstream analysis:

``` r
a <- set_scAPAList(.clus.counts = co, .row.Data = row.Data, 
                     .down.seq = aseq)
```

There is no need to match the order of the peaks in the three object as the function set\_scAPAList does so automatically. For the next steps of the analysis, see scAPA package [vignette](scAPA_vignette.md).
