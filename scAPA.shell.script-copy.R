c = 30
CPM.cuttoff = 10
A.number = 8
filter.border.left = 10
filter.border.right = 140
ICPM.cuttoff = 10
IA.number = 7
Ifilter.border.left = 1
Ifilter.border.right = 200
IC.cuttoff = 50
org = "mm10"
drop.seq.tools.path = "/home/ubuntu/Drop-seq_tools-2.2.0"
char.length.path = paste0("/data/", org, ".chrom.sizes")
fasta.path = paste0("/data/", org, ".fa")
path.to.files = "/data/down.sampled.spermatogenesis"

# Load R packages and function scripts------------------------------------------
# This function load packages and install them in case they are not installed
require(package = "dplyr", warn.conflicts = F)
require(package = "tidyr", warn.conflicts = F)
require(package = "ggplot2", warn.conflicts = F)
require(package = "EnvStats", warn.conflicts = F)
if(c > 1){
require(package = "parallel", warn.conflicts = F)
}
require(package = "mclust", warn.conflicts = F)
require("Rsubread")
require("scAPA")

# Input ---------------------------
setwd(path.to.files)
bams.files <- list.files(pattern = ".bam$")
samples.vector <- gsub(x = bams.files, pattern = ".bam", replacement = "")
bams <- paste(bams.files, collapse = " ")
samples <- gsub(x = bams, pattern = ".bam", replacement = "")

if (!dir.exists("scAPA")) dir.create("scAPA")
if (!dir.exists("scAPA/Log.files")) dir.create("scAPA/Log.files")
if (!dir.exists("scAPA/temp")) dir.create("scAPA/temp")
if (!dir.exists("scAPA/outs")) dir.create("scAPA/outs")

setwd(paste0(path.to.files, "/scAPA"))

script.start.message <- paste0(Sys.time(),
                                "\t Started scAPA.script.R.\n",
                                "-p = ", path.to.files, "\n",
                                "-org = ", org, "\n",
                                "-c = ", c, "\n",
                                "-cpm = ", CPM.cuttoff, "\n",
                                "-a = ", A.number, "\n",
                                "-u = ", filter.border.left, "\n",
                                "-d = ", filter.border.right, "\n")

write(x = script.start.message, file = "scAPA.script.log", append = F)

# Our pipeline consists of the following 5 steps:
# 1. Defining 3'UTR peaks
#2. Quantifying the usage of each peak in each cell cluster
#3. Filtering peaks
#4. Detecting dynamical APA events
#5. Inferring global trends of APA modulation

step = 1

# 1.Defining 3'UTR peaks PCR duplicates removal ---------------------------
# FilterBAM ---------------------------------------------------------------
write_log_start("Stage 1: Defining 3'UTR peaks\n\n", command = NA)
write_log_start("Stage 1a: PCR duplicates removal\n\n", command = NA)

command <- paste0("for sample in ", samples, " ;", " do ",
                            drop.seq.tools.path, "FilterBam TAG_RETAIN=UB",
                            " I=../${sample}.bam O=./temp/UB.${sample}.bam &> ",
                            "./Log.files/FilterBAM.${sample}.bam.out; done")

write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1
  
write_log_start("FilterBAM", command = command)
write_log(stage = "FilterBAM")

# ii. Use UMI tools with “method=unique” so that cellranger corrected molecular barcodes are used

# samtools index ----------------------------------------------------------
command <- paste0("for sample in ", samples, " ; do ",
                        "samtools index temp/UB.${sample}.bam",
                        " &> ./Log.files/Index.${sample}.out; done")

write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1
  
write_log_start("samtools index", command = command)

# umi_tools ---------------------------------------------------------------
command <- paste0("for sample in ", samples, " ;do ",
                  "umi_tools dedup -I ",
                  "temp/UB.${sample}.bam -S temp/dedup.${sample}.bam",
                   " --method=unique --extract-umi-method=tag ",
                   "--umi-tag=UB --cell-tag=CB > ",
                    "./Log.files/umi_tools.${sample}.out;",
                    " done")

write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

write_log_start("umi_tools,", command = command)

command = "rm temp/UB.*.bam*"
write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

# Peak detection -------------------------------------------------------
write_log_start("Stage 1b: Peak detection\n\n", command = NA)

# makeTagDirectory --------------------------------------------------------
command <- paste0("makeTagDirectory ./temp/Tagdirectory ",
                  "./temp/dedup.* &> ./Log.files/",
                  "makeTagDirectory.out")

write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

write_log_start("makeTagDirectory", command = command)

# # findPeaks ---------------------------------------------------------------
command <- paste0("findPeaks ./temp/Tagdirectory -size 50 -frag",
                  "Length 100 -minDist 1 -strand separate -o ",
                  "./temp/Peakfile &> ./Log.files/findPeaks.out")

write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

write_log_start("findPeaks,", command = command)
write_log(stage = "findPeaks")

#RUN SHELL SCRIPTS BEFORE RUNNING NEXT STEPS

# Peak analysis -----------------------------------------------------------
setwd("./temp")
merge_peaks(bedtools.path = bedtools.path, peaks.file = "Peakfile", path = "./")
peaks.bed <- intersect_peaks(org = org, bed.name = "./merge.peakfile.bed",
                             path = "", bedtools.path = bedtools.path,
                             introns = F)
write.bed(.x = peaks.bed, f = "./peaks.bed")

# Spliting bimodal peaks --------------------------------------------------
write_log_start(paste0("Stage 1c: Separating Peaks with bimodal UMI",
                       " counts distribution\n\n"), command = NA,
                f = "../scAPA.script.log")

# Samtools merge
dedup.bams <- paste0("dedup.", bams.files, collapse = " ")
command <- paste0(samtools.path, "samtools merge merged.bam ",
                                dedup.bams,
                                " &> ../Log.files/samtoolsmerge.out")

write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

write_log_start(f = "../scAPA.script.log", "samtools merge", command=command)
write_log(f = "../scAPA.script.log", stage = "samtools merge")

# Bedgraphs
command <- paste0(bedtools.path,
                           "bedtools genomecov -ibam merged.bam ",
                           "-bg -strand + | awk 'BEGIN ",
                           "{OFS = \"\t\"}{print $1",
                           ", $2, $3, $4, \".\", \"+\"}' > merged.wig")

write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

command <- paste0(bedtools.path, "bedtools genomecov -ibam merged.bam ",
                            "-bg -strand - | awk 'BEGIN {OFS = \"\t\"}{print $1",
                            ", $2, $3, $4, \".\", \"-\"}' >> merged.wig")
write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

command = "rm merged.bam"
write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

command <- paste0(bedtools.path, "bedtools intersect -s -wb ",
                                "-b peaks.bed -a merged.wig > intersected.wig")
write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

command = "rm merged.wig"
write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

#RUN BASH SCRIPTS BEFORE RUNNING NEXT STEPS

peaks.wig <- read.delim(file = "intersected.wig", header = F)
peaks.wig <- split(x = peaks.wig, f = peaks.wig$V10, drop = T)

# mclust
mclust.command <- paste0("Mclust()")
write_log_start(f = "../scAPA.script.log", "Mclust,", command = mclust.command)

if(c > 1){
bed <- plyr::rbind.fill(parallel::mclapply(1:length(peaks.wig),
                                           FUN = creat_mclus,
                                           mc.cores = c,
                                           mc.preschedule = T))
}
if(c == 1){
bed <- plyr::rbind.fill(lapply(X = 1:length(peaks.wig),
                                           FUN = creat_mclus))
}

rm("peaks.wig")

if (!exists("bed")) {
    write_log(f = "../scAPA.script.log", stage = "Mclust", success = F)
    stop("Mclust did not end, shell script did NOT end!")
}

write_log(f = "../scAPA.script.log", stage = "Mclust")

command = "rm intersected.wig merged.wig"
write(command, paste0(path.to.files, "/scAPA/", "step", step, ".sh"))
step = step+1

# Split bams ----------------------------------------

allclusters <- character()
for (j in 1:length(samples.vector)) {
    sample <- samples.vector[j]
    list.file <- paste0("../../clusters_", sample, ".txt")
    list.cells <- read.delim(file = list.file, header = F)
    list.cells[, 1] <- as.character(list.cells[, 1])
    list.cells <- split(x = list.cells, f = list.cells[, 2], drop = T)
    allclusters <- c(allclusters, names(list.cells))
    for (i in 1:length(list.cells)) {
        cluster <- names(list.cells)[i]
        cluster.cells <- as.data.frame(list.cells[[i]][, 1])
        tagValuefile <- paste0(cluster, "_", sample)
        write.bed(.x = cluster.cells, f = tagValuefile)
        command <- paste0(drop.seq.tools.path,
                                         "FilterBamByTag ",
                                         "TAG=CB I=dedup.",
                                         sample, ".bam O=",
                                         sample, "_", cluster, ".bam ",
                                         "TAG_VALUES_FILE=", tagValuefile,
                                         " &> ../Log.files/",
                                         "FilterBamByTag_",
            tagValuefile)
        
        write(command, paste0(path.to.files, "/scAPA/", "step", step, "_", tagValuefile, ".sh"))
        step = step+1
        
        write_log_start(f = "../scAPA.script.log",
                        paste0("FilterBAMbyTag for ", tagValuefile),
                        command = command)

        write_log(f = "../scAPA.script.log",
                  stage = paste0("FilterBAMbyTag.command: ",
                                 sample, "_", cluster, ".bam"))
    }
}

#RUN SHELL SCRIPTS BEFORE RUNNING NEXT STEPS

clusters <- unique(allclusters)
for (i in 1:length(clusters)) {
    bams.to.merge <- paste0(samples.vector, "_", clusters[i], ".bam",
                            collapse = " ")
    command <- paste0("samtools merge ", clusters[i],
                                     ".bam ", bams.to.merge)
    write(command, paste0(path.to.files, "/scAPA/", "step", step, "_", clusters[i], ".sh"))
    step = step+1
    
    write_log_start(f = "../scAPA.script.log",
                    paste0("samtools merge for ", clusters[i]),
                    command = command)

    write_log(f = "../scAPA.script.log",
              stage = paste0("samtools merge for ", clusters[i]))
}

#RUN SHELL SCRIPTS BEFORE RUNNING NEXT STEPS

write_log_start(f = "../scAPA.script.log", "featureCounts for 3UTRs",
                command = NA)
bam.cluster.files <- paste0(clusters, ".bam")
counts <- Rsubread::featureCounts(files = bam.cluster.files, isGTFAnnotationFile = F,
                        strandSpecific = 1, annot.ext = utr.saf,
                        largestOverlap = T, nthreads = c)
co <- cbind.data.frame(rownames(counts$counts), counts$counts)
colnames(co) <- c("Peak_ID", clusters)
meta <- counts$annotation
meta <- meta[, c(2, 3, 4, 1, 6, 5)]
metadata <- read_down.seq(saf = utr.saf,
                          char.length.path = char.length.path,
                          fasta.path = fasta.path, chr.modify = T)
aseq <- metadata[, c(4, 6)]
a <- set_scAPAList(.clus.counts = co, .row.Data = meta, .down.seq = aseq)
saveRDS(object = a, file = "../outs/Peaks.RDS")
write_log(f = "../scAPA.script.log", stage = "featureCounts for 3UTRs")

# Peak filtering ----------------------------------------------------------
write_log_start("Stage 3: Peak filtering (3'UTRs)\n\n",
                command = NA, f = "../scAPA.script.log")

a <- calc_cpm(a)
keep.cpm <- rowSums(a@norm$CPM) > CPM.cuttoff
a.fil <- a[keep.cpm, ]

# Internal priming filtering
IP.seq <- paste0(rep("A", times = A.number), collapse = "")
a.fil <- filter_IP(x = a.fil, int.priming.seq = IP.seq,
                   left = filter.border.left, right = filter.border.right)

# Statistical analysis ----------------------------------------------------
write_log_start("Stage 4: Statistical analysis (3'UTRs)\n\n",
                command = NA, f = "../scAPA.script.log")
results <- set_scAPAresults(a.fil)
# Chi-squre test for APA
results <- test_APA(results, clus = "all")
# For APA with more than one peak and whose p-value is < sig.level,
# perform chi-squered test for goodness of fit
results <- test_peaks(results, clus = "all", sig.level = 0.05)

# Inferring global trends -------------------------------------------------
write_log_start("Stage 5: Inferring global trends  (3'UTRs)\n\n", command = NA,
                f = "../scAPA.script.log")
results <- calc_pAi_mat(results)
results <- calc_p_pui_mat(results)
saveRDS(object = results, file = "../outs/results.RDS")

# Write final outputs -----------------------------------------------------
setwd("../")
results.output <- disply_results(x = results, org = org)
a.fil <- annotate(a.fil, org = org)

write.table(x = results.output, file = "./outs/APA.events.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

write.table(x = a.fil@row.Data, file = "./outs/ThreeUTR.peaks.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

write.peaks(.x = results, .f = "./outs/UTRs.with.multiple.peaks.txt")

file.create("./outs/summary.UTR.txt")
cat("Number of peaks passed filter:\t", nrow(a.fil@clus.counts),
    "\nNumber of significant (FDR < 5%) APA events:",
    sum(results@pvalues[[1]][, 2] < 0.05, na.rm = T), "\n",
    file = "./outs/summary.UTR.txt")

# Plot --------------------------------------------------------------------
sig.utrs <- as.vector(results@pvalues$all[,2] < 0.05)
sig <- results[which(sig.utrs),]
tidy.ppui <- as.data.frame(sig@ppui.clus)
colnames(tidy.ppui) <- gsub(x = colnames(tidy.ppui),
                            pattern = "_proximal_PUI", replacement = "")
tidy.ppui <- tidyr::gather(data = tidy.ppui)
colnames(tidy.ppui) <- c("Cluster", "value")
pdf("./outs/Proximal.PUI.ECDF.pdf")
p <- ggplot2::ggplot(data = tidy.ppui,
                     ggplot2::aes(x = value, color = Cluster))
p <- p + ggplot2::stat_ecdf(size = 1)
p <- p + ggplot2::theme_bw()
p <- p + ggplot2::xlab("Proximal peak usage index")
p <- p + ggplot2::ylab("Cumulative fraction")
print(p)
dev.off()