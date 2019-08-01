# #!/usr/bin/env Rscript
script.args = commandArgs(trailingOnly = TRUE)
if (((script.args[1] == "-h") | (script.args[1] == "--help") | (is.na(script.args[1])))) {
    cat("\n Usage: scAPA.script.R -p <path.to.files> -org <organism> -sp ",
        "<path.to.script.Dir> [options]\n\n",
        "Basic Options:\n\n\t-p \tThe path to the folder with input bam files,"
        , "and cell cluster annotatios files",
        " and cell.cluster.list.txt\n\t\t/path/to/files\n\n\t-org \tOrganism,",
        " either Mm for mouse, or Hs for human\n\n\t-c\tThe number of cores to",
        "use.\n\t-sp The path to the directory of ",
        "scAPA.shell.script.\n\n\t-loc\tPath to users' library trees within",
        "which packages are installed and looked for.\n\nOutput Options",
        "\n\n\t-wig\tWeather to generate cluster wig files.",
        " Defult value: false\n\nAnalysis ", "\n\n\t-ChangePoint\t weather to",
        " use Change Point. Defult is false.",
        "Options\n\n\t-sc\tIf true, (default) counts read from individual",
        " cells (slower). Otherwise, counts reads from ",
        "clusters.\n\n\t-int\tIf true, ",
        "(default) perform intronic APA analysis as well as 3'UTR",
        " analysis.\n\nFiltering",
        "Options 3'UTRS\n\n\t-cpm\tConsider only peaks with more than",
        " a total sum of CPMs over all cell clusters larger than ",
        "-cpm.\n\t\tDefult",
        " value: 10\n\n\t-a\tNumeric value, filter out peaks with",
        " -a consecutive As", "in the region -u to -d downstream thire 3' ",
        "edge.\n\t\tDefult value: 8\n\n\t-u",
        "\tNumeric value, filter out peaks with -a consecutive As in the ",
        "region -u to -d downstream thire 3' edge.\n\t\tDefult value: ",
        " 10\n\n\t-d\tNumeric ",
        "value, filter out peaks with -a consecutive As in the region -u to ",
        "-d downstream thire 3' edge.", "\n\t\tDefult value: 140",
        "\n\nFiltering Options introns\n\n\t-Icpm\tConsider only peaks ",
        "with more than",
        " a total sum of CPMs over all cell clusters larger than -Icpm.",
        "\n\n\t-Ico\tConsider only peaks with more than",
        " a total sum of counts over all cell clusters larger than",
        " -Ico\n\t\tDefult",
        " value: 50\n\n\t-Ia\tNumeric value, filter out peaks with -Ia",
        " consecutive As", "in the region -u to -d downstream thire 3' ",
        "edge.\n\t\tDefult value: 7\n\n\t-Iu",
        "\tNumeric value, filter out peaks with -Ia consecutive As in the ",
        "region -Iu to -Id downstream thire 3' edge.\n\t\tDefult value: ",
        "1\n\n\t-Id\tNumeric ",
        "value, filter out peaks with -Ia consecutive As in the region",
        " -Iu to -Id ", "downstream thire 3' edge.",
        "\n\t\tDefult value: 200\n\n")
    q(save = "no", status = 1, runLast = FALSE)
}

# Read in user's arguments ------------------------------------------------
read_args <- function(arg = script.args, arg.string, defult = NULL) {
    nchr.arg <- nchar(arg.string)
    arg.string <- paste0(arg.string, " .*")
    pos <- regexpr(pattern = arg.string, text = arg)
    if (pos > 0) {
        pos <- pos + nchr.arg + 1
        out <- substring(arg, first = pos)
        out <- substring(out, first = 1, last = (regexpr(pattern = " ",
                                                         text = out) - 1))
    } else {
        out <- defult
    }
    out
}

# test arguments:
if (length(script.args) == 0) {
    stop(paste0("scAPA.script.R -p <path.to.files> -org <organism>.\n",
                "The following arguments are optional: -cpm <CPM.cuttoff>",
                "  -a <A.number> -u <filter.border.left> -d ",
                "<filter.border.right>"),
        call. = FALSE)
}
if (length(script.args) < 3) {
    stop(paste0("At least three argument must be supplied (path to files ,",
                "organism \"Hs\"/\"Mm\" and path to ",
                "scAPA.shell.script folder)"),
         call. = FALSE)
}
script.args <- paste(script.args, collapse = " ")
script.args <- paste0(script.args, " ")
# Read arguments
path.to.files <- read_args(arg.string = "-p", arg = script.args)
org <- read_args(arg.string = "-org", arg = script.args)
if (!(org %in% c("Mm", "Hs"))) {
    stop(paste0("The second argument (organism) must be either Mm ",
                "(for mouse) or Hs (for human)"), call. = FALSE)
}
path.to.config <- read_args(arg.string = "-sp", arg = script.args)
c <- read_args(arg.string = "-c", defult = 30, arg = script.args)
c <- as.numeric(c)
# The following arguments are optional. If not provided, defult values are set.
sc <- read_args(arg.string = "-sc", defult = "true", arg = script.args)
sc <- ifelse(sc == "true", TRUE, FALSE)
int <- read_args(arg.string = "-int", defult = "true", arg = script.args)
int <- ifelse(int == "true", TRUE, FALSE)
CPM.cuttoff <- read_args(arg.string = "-cpm", defult = 10,
                         arg = script.args)
A.number <- read_args(arg.string = "-a", defult = 8, arg = script.args)
filter.border.left <- read_args(arg.string = "-u", defult = 10,
                                arg = script.args)
filter.border.right <- read_args(arg.string = "-d", defult = 140,
                                 arg = script.args)
ICPM.cuttoff <- read_args(arg.string = "-Icpm", defult = 10,
                          arg = script.args)
IA.number <- read_args(arg.string = "-Ia", defult = 7,
                       arg = script.args)
Ifilter.border.left <- read_args(arg.string = "-Iu", defult = 1,
                                 arg = script.args)
Ifilter.border.right <- read_args(arg.string = "-Id", defult = 200,
                                  arg = script.args)
IC.cuttoff <- read_args(arg.string = "-Ico", defult = 50, arg = script.args)
ChangePoint <- read_args(arg.string = "-ChangePoint ", defult = FALSE,
                         arg = script.args)
ChangePoint <- ifelse(ChangePoint == "true", TRUE, FALSE)
wig <- read_args(arg.string = "-wig", defult = "false",
                 arg = script.args)
wig <- ifelse(wig == "true", TRUE, FALSE)
loc <- read_args(arg.string = "-loc", defult = "false", arg = script.args)
loc <- ifelse(wig == "false", NULL, loc)

# Read in configuration file ----------------------------------------------
setwd(path.to.config)
configfile <- as.character(read.delim("configfile.txt", header = F)[, 1])
configfile[configfile == "PATH"] <- ""
if (org == "Mm") {
    char.length.path <- configfile[2]
    fasta.path <- configfile[4]
}
if (org == "Hs") {
    char.length.path <- configfile[6]
    fasta.path <- configfile[8]
}
drop.seq.tools.path <- configfile[10]
samtools.path <- configfile[12]
umi_tools.path <- configfile[14]
homer.path <- configfile[16]
bedtools.path <- configfile[18]
path.to.chang.point <- configfile[20]
# Load R packeges and function scripts------------------------------------------
# This function load packeges and install them in case they are not installed
if(!is.null(loc)) .libPaths(loc)
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
log.files.messege <- paste0("\nThe script's log file is written into: ",
                            path.to.files, "/scAPA/scAPA.script.log.\n",
                            "The log files of the tools used are ",
                            "written to the directory: ", path.to.files,
    "/scAPA/Log.files.\n\n")
cat(x = log.files.messege)
if (!int) {
    script.start.messsege <- paste0(Sys.time(),
                                    "\t Started scAPA.script.R.\n",
                                    "-p = ", path.to.files, "\n",
                                    "-org = ", org, "\n",
                                    "-c = ", c, "\n",
                                    "-sp = ", path.to.config, "\n",
                                    "-sc = ", sc, "\n",
                                    "-int = ", int, "\n",
                                    "-wig = ", wig, "\n",
                                    "-ChangePoint = ", ChangePoint, "\n",
                                    "-cpm = ", CPM.cuttoff, "\n",
                                    "-a = ", A.number, "\n",
                                    "-u = ", filter.border.left, "\n",
                                    "-d = ", filter.border.right, "\n")
} else {
    script.start.messsege <- paste0(Sys.time(),
                                    "\t Started scAPA.script.R.\n",
                                    "-p = ", path.to.files, "\n",
                                    "-org = ", org, "\n",
                                    "-c = ", c, "\n",
                                    "-sp = ", path.to.config, "\n",
                                    "-sc = ", sc, "\n",
                                    "-int = ", int, "\n",
                                    "-wig = ", wig, "\n",
                                    "-ChangePoint = ", ChangePoint, "\n",
                                    "-cpm = ", CPM.cuttoff, "\n",
                                    "-a = ", A.number, "\n",
                                    "-u = ", filter.border.left, "\n",
                                    "-d = ", filter.border.right, "\n",
                                    "-Icpm = ", ICPM.cuttoff, "\n",
                                    "-Ia = ", IA.number, "\n",
                                    "-Iu = ", Ifilter.border.left, "\n",
                                    "-Id = ", Ifilter.border.right, "\n")
}
write(x = script.start.messsege, file = "scAPA.script.log", append = F)

# Our pipeline consists of the following 5 steps:
# 1. Defining 3'UTR peaks
#2. Quantifying the usage of each peak in each cell cluster
#3. Filtering peaks
#4. Detecting dynamical APA events
#5. Inferring global trends of APA modulation

# 1.Defining 3'UTR peaks PCR duplicates removal ---------------------------
# FilterBAM ---------------------------------------------------------------
write_log_start("Stage 1: Defining 3'UTR peaks\n\n", command = NA)
write_log_start("Stage 1a: PCR duplicates removal\n\n", command = NA)
FilterBAM.command <- paste0("sh -c 'for sample in ", samples, " ;", " do ",
                            drop.seq.tools.path, "FilterBam TAG_RETAIN=UB",
                            " I=../${sample}.bam O=./temp/UB.${sample}.bam &> ",
                            "./Log.files/FilterBAM.${sample}.bam.out; done'")
write_log_start("FilterBAM", command = FilterBAM.command)
system(command = FilterBAM.command, wait = T)
if (!all(file.exists(paste0("./temp/UB.", samples.vector, ".bam")))) {
    write_log(stage = "FilterBAM", success = F)
    stop("FilterBAM did not end, shell script did NOT end!")
}
write_log(stage = "FilterBAM")

# ii. Use UMI tools with “method=unique” so that cellranger corrected molecular barcodes are used

# samtools index ----------------------------------------------------------
index.command <- paste0("sh -c 'for sample in ", samples, " ; do ",
                        samtools.path,
                        "samtools index temp/UB.${sample}.bam",
                        " &> ./Log.files/Index.${sample}.out; done'")
write_log_start("samtools index", command = index.command)
print(index.command)
system(command = index.command, wait = T)
if (!all(file.exists(paste0("./temp/UB.", samples.vector, ".bam.bai")))) {
    write_log(stage = "samtools index", success = F)
    stop("samtools index did not end, shell script did NOT end!")
}
write_log(stage = "samtools index")

# umi_tools ---------------------------------------------------------------
umi_tools.command <- paste0("sh -c 'for sample in ", samples, " ;do ",
                            umi_tools.path, "umi_tools dedup -I ",
                            "temp/UB.${sample}.bam -S temp/dedup.${sample}.bam",
                            " --method=unique --extract-umi-method=tag ",
                            "--umi-tag=UB --cell-tag=CB &> ",
                            "./Log.files/umi_tools.${sample}.out;",
                            " done'")
write_log_start("umi_tools,", command = umi_tools.command)
print(umi_tools.command)
system(command = umi_tools.command, wait = T)
if (!all(file.exists(paste0("temp/dedup.", samples.vector, ".bam")))) {
    write_log(stage = "umi_tools", success = F)
    stop("umi_tools did not end, shell script did NOT end!")
}
write_log(stage = "umi_tools")
system(command = "rm temp/UB.*.bam*", wait = T)

# Peak detection -------------------------------------------------------
write_log_start("Stage 1b: Peak detection\n\n", command = NA)

# makeTagDirectory --------------------------------------------------------
makeTagDirectory.command <- paste0(homer.path,
                                   "makeTagDirectory ./temp/Tagdirectory ",
                                   "./temp/dedup.* &> ./Log.files/",
                                   "makeTagDirectory.out")
write_log_start("makeTagDirectory", command = makeTagDirectory.command)
print(makeTagDirectory.command)
system(command = makeTagDirectory.command, wait = T)
if (!dir.exists(paste0("./temp/Tagdirectory"))) {
    write_log(stage = "makeTagDirectory", success = F)
    stop("makeTagDirectory did not end, shell script did NOT end!")
}
write_log(stage = "makeTagDirectory")
# # findPeaks ---------------------------------------------------------------
findPeaks.command <- paste0(homer.path,
                            "findPeaks ./temp/Tagdirectory -size 50 -frag",
                            "Length 100 -minDist 1 -strand separate -o ",
                            "./temp/Peakfile &> ./Log.files/findPeaks.out")
write_log_start("findPeaks,", command = findPeaks.command)
system(command = findPeaks.command, wait = T)
if (!file.exists(paste0("./temp/Peakfile"))) {
    write_log(stage = "findPeaks", success = F)
    stop("findPeaks did not end, shell script did NOT end!")
}
write_log(stage = "findPeaks")


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
samtools.merge.commad <- paste0(samtools.path, "samtools merge merged.bam ",
                                dedup.bams,
                                " &> ../Log.files/samtoolsmerge.out")
write_log_start(f = "../scAPA.script.log", "samtools merge",
                command = samtools.merge.commad)
system(command = samtools.merge.commad, wait = T)
if (!file.exists(paste0("merged.bam"))) {
    write_log(f = "../scAPA.script.log", stage = "samtools merge", success = F)
    stop("samtools merge did not end, shell script did NOT end!")
}
write_log(f = "../scAPA.script.log", stage = "samtools merge")

# Bedgraphs
wig.plus.command <- paste0(bedtools.path,
                           "bedtools genomecov -ibam merged.bam ",
                           "-bg -strand + | awk 'BEGIN ",
                           "{OFS = \"\t\"}{print $1",
                           ", $2, $3, $4, \".\", \"+\"}' > merged.wig")
system(command = wig.plus.command, wait = T)
wig.minus.command <- paste0(bedtools.path, "bedtools genomecov -ibam merged.bam ",
                            "-bg -strand - | awk 'BEGIN {OFS = \"\t\"}{print $1",
                            ", $2, $3, $4, \".\", \"-\"}' >> merged.wig")
system(command = wig.minus.command, wait = T)
system("rm merged.bam", wait = T)
intersect.wig.command <- paste0(bedtools.path, "bedtools intersect -s -wb ",
                                "-b peaks.bed -a merged.wig > intersected.wig")
system(intersect.wig.command, wait = T)
peaks.wig <- read.delim(file = "intersected.wig", header = F)
peaks.wig <- split(x = peaks.wig, f = peaks.wig$V10, drop = T)
system("rm intersected.wig merged.wig", wait = T)
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
if (!exists("bed")) {
    write_log(f = "../scAPA.script.log", stage = "Mclust", success = F)
    stop("Mclust did not end, shell script did NOT end!")
}
write_log(f = "../scAPA.script.log", stage = "Mclust")
rm("peaks.wig")
# SAFs --------------------------------------------------------------------
utr.saf <- bed[, c(4, 1, 2, 3, 6)]
rm("bed")
colnames(utr.saf) <- c("GeneID", "Chr", "Start", "End", "Strand")
if (int) {
    merge_peaks(bedtools.path = bedtools.path, peaks.file = "Peakfile", path = "./")
    int.peaks.bed <- intersect_peaks(org = org,
                                     bed.name = "./merge.peakfile.bed",
                                     path = "", bedtools.path = bedtools.path,
                                     introns = T)
    int.peaks.bed <- as.data.frame(int.peaks.bed)
    int.saf <- int.peaks.bed[, c(4, 1, 2, 3, 6)]
    rm("int.peaks.bed")
    colnames(int.saf) <- c("GeneID", "Chr", "Start", "End", "Strand")
}

# Split bams for wig or if sc = false----------------------------------------
if (wig | !sc | ChangePoint) {
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
            FilterBAMbyTag.command <- paste0(drop.seq.tools.path,
                                             "FilterBamByTag ",
                                             "TAG=CB I=dedup.",
                                             sample, ".bam O=",
                                             sample, "_", cluster, ".bam ",
                                             "TAG_VALUES_FILE=", tagValuefile,
                                             " &> ../Log.files/",
                                             "FilterBamByTag_",
                tagValuefile)
            write_log_start(f = "../scAPA.script.log",
                            paste0("FilterBAMbyTag for ", tagValuefile),
                            command = FilterBAMbyTag.command)
            system(command = FilterBAMbyTag.command, wait = T)
            if (!file.exists(paste0(sample, "_", cluster, ".bam"))) {
                write_log(f = "../scAPA.script.log",
                          stage = paste0("FilterBAMbyTag for ", tagValuefile),
                          success = F)
                stop("FilterBAMbyTag.command did not end, shell script did NOT end!")
            }
            write_log(f = "../scAPA.script.log",
                      stage = paste0("FilterBAMbyTag.command: ",
                                     sample, "_", cluster, ".bam"))
        }
    }
    clusters <- unique(allclusters)
    for (i in 1:length(clusters)) {
        bams.to.merge <- paste0(samples.vector, "_", clusters[i], ".bam",
                                collapse = " ")
        samtools.merge.command <- paste0("samtools merge ", clusters[i],
                                         ".bam ", bams.to.merge)
        print(samtools.merge.command)
        write_log_start(f = "../scAPA.script.log",
                        paste0("samtools merge for ", clusters[i]),
                        command = FilterBAMbyTag.command)
        system(command = samtools.merge.command, wait = T)
        write_log(f = "../scAPA.script.log",
                  stage = paste0("samtools merge for ", clusters[i]))
        if (wig) {
            wig.command <- paste0(bedtools.path,
                                  "bedtools genomecov -trackline -bg ",
                                  "-trackopts \"name=", clusters[i],
                                  "-ibam ", clusters[i],
                                  ".bam > ../outs",
                                  clusters[i], ".wig")
            print(wig.command)
            write_log_start(f = "../scAPA.script.log",
                            paste0("bedtools genomecov for ", clusters[i]),
                            command = wig.command)
            system(command = wig.command, wait = T)
            write_log(f = "../scAPA.script.log",
                      stage = paste0("bedtools genomecov for ", clusters[i]))
        }
    }

}

if (!sc) {
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
    if (int) {
        write_log_start(f = "../scAPA.script.log",
                        "featureCounts for introns", command = NA)
        counts_int <- Rsubread::featureCounts(files = bam.cluster.files,
                                    isGTFAnnotationFile = F,
                                    annot.ext = int.saf, strandSpecific = 1,
                                    largestOverlap = T, nthreads = c)
        co_int <- cbind.data.frame(rownames(counts_int$counts),
                                   counts_int$counts)
        colnames(co_int) <- c("Peak_ID", clusters)
        meta_int <- counts_int$annotation
        meta_int <- meta_int[, c(2, 3, 4, 1, 6, 5)]
        metadata_int <- read_down.seq(saf = int.saf,
                                      char.length.path = char.length.path,
                                      fasta.path = fasta.path, chr.modify = T)
        aseq_int <- metadata_int[, c(4, 6)]
        a.int <- set_scAPAList(.clus.counts = co_int, .row.Data = meta_int,
                               .down.seq = aseq_int)
        saveRDS(object = a.int, file = "../outs/Peaks.int.RDS")
        write_log(f = "../scAPA.script.log",
                  stage = "featureCounts for introns")
    }
}

# ChangePoint -------------------------------------------------------------
if (ChangePoint) {
    allcombintations <- paste0(" -t ", bam.cluster.files, " -c ",
                               rep(bam.cluster.files,
                                   each = length(bam.cluster.files)))
    toexclude <- seq(from = 1, to = length(allcombintations),
                     by = (length(bam.cluster.files) + 1))
    allcombintations <- allcombintations[-toexclude]
    if (org == "Mm")
        write.bed(.x = mm10.utr, f = "3UTR.BED")
    if (org == "Hs")
        write.bed(.x = hg19.utr, f = "3UTR.BED")
    allcombintations_names <- gsub(pattern = " -t ", replacement = "",
                                   x = allcombintations)
    allcombintations_names <- gsub(pattern = " -c ", replacement = "_vs_",
                                   x = allcombintations_names)
    allcombintations_names <- paste0("ChangePoint.", allcombintations_names)
    ChangePoint.command <- paste0("perl ", path.to.chang.point,
                                  "change_point.pl -g 3UTR.BED",
                                  allcombintations, " -d s -o ",
                                  allcombintations_names,
                                  " &> ../Log.files/", allcombintations_names,
        ".out")
    write_log_start(f = "../scAPA.script.log", "Change Point", command = NA)
    for (i in 1:length(ChangePoint.command)) {
        system(command = ChangePoint.command[i], wait = T)
    }
    system("mv ChangePoint.* ../outs", wait = T)
}

# Sc ----------------------------------------------------------------------
if (sc) {
    write_log_start(f = "../scAPA.script.log", "Splitting bams", command = NA)
    if (!dir.exists("CellBams"))
        dir.create("CellBams")
    full.list.cells <- data.frame()
    for (j in 1:length(samples.vector)) {
        sample <- samples.vector[j]
        list.file <- paste0("../../clusters_", sample, ".txt")
        list.cells <- read.delim(file = list.file, header = F)
        list.cells[, 1] <- as.character(list.cells[, 1])
        list.cells <- as.data.frame(list.cells)
        list.cells$Cell <- paste0(sample, "_", list.cells$V1)
        full.list.cells <- rbind.data.frame(full.list.cells,
                                            list.cells[, c(3, 2)])
        list.cells <- split(x = list.cells, f = list.cells$V1, drop = T)
        split_bams <- function(x) {
            cell <- x[, 1]
            bam <- paste0(x[, 3], ".bam")
            FilterBAMbyTag.command <- paste0(drop.seq.tools.path,
                                             "FilterBamByTag ",
                                             "TAG=CB I=dedup.",
                                             sample, ".bam O=",
                                             "./CellBams/", bam,
                                             " TAG_VALUE=", cell,
                                             " &>> ../Log.files/",
                                             "FillterBam_forcells.out")
            system(command = FilterBAMbyTag.command, wait = T)
        }
        if(c > 1){
        parallel::mclapply(X = list.cells, FUN = split_bams, mc.cores = c,
                           mc.preschedule = T)
        }
        if (c ==1){
        lapply(X = list.cells, FUN = split_bams)
        }
    }
    write_log(f = "../scAPA.script.log", stage = "Splitting bams")

    # Count reads --------------------------------------------------------------
    write_log_start("Stage 2: Quantifying peak usage\n\n", command = NA,
                    f = "../scAPA.script.log")
    write_log_start(f = "../scAPA.script.log", "featureCounts for 3UTRs",
                    command = NA)
    bam.cluster.files <- list.files(path = "./CellBams", pattern = ".bam$",
                                    full.names = F)
    cellnames <- gsub(x = bam.cluster.files, pattern = ".bam", replacement = "")
    bam.cluster.files <- list.files(path = "./CellBams", pattern = ".bam$",
                                    full.names = T)
    counts <- Rsubread::featureCounts(files = bam.cluster.files, isGTFAnnotationFile = F,
                            strandSpecific = 1, annot.ext = utr.saf,
                            largestOverlap = T, nthreads = c)
    co <- cbind.data.frame(rownames(counts$counts), counts$counts)
    colnames(co) <- c("Peak_ID", cellnames)
    meta <- counts$annotation
    meta <- meta[, c(2, 3, 4, 1, 6, 5)]
    metadata <- read_down.seq(saf = utr.saf, char.length.path = char.length.path,
                              fasta.path = fasta.path, chr.modify = T)
    aseq <- metadata[, c(4, 6)]
    a <- set_scAPAList(.cells.counts = co, .row.Data = meta, .down.seq = aseq,
                       .cluster.anot = full.list.cells)
    rm(list = c("co", "meta", "aseq"))
    saveRDS(object = a, file = "../outs/Peaks.RDS")
    write_log(f = "../scAPA.script.log", stage = "featureCounts for 3UTRs")
    if (int) {
        write_log_start(f = "../scAPA.script.log",
                        "featureCounts for introns", command = NA)
        counts_int <- Rsubread::featureCounts(files = bam.cluster.files,
                                    isGTFAnnotationFile = F,
                                    strandSpecific = 1,
                                    annot.ext = int.saf,
                                    largestOverlap = T,
                                    nthreads = c)
        co_int <- cbind.data.frame(rownames(counts_int$counts),
                                   counts_int$counts)
        colnames(co_int) <- c("Peak_ID", cellnames)
        meta_int <- counts_int$annotation
        meta_int <- meta_int[, c(2, 3, 4, 1, 6, 5)]
        metadata_int <- read_down.seq(saf = int.saf,
                                      char.length.path = char.length.path,
                                      fasta.path = fasta.path,
                                      chr.modify = T)
        aseq_int <- metadata_int[, c(4, 6)]
        a.int <- set_scAPAList(.cells.counts = co_int, .row.Data = meta_int,
                               .down.seq = aseq_int,
                               .cluster.anot = full.list.cells)
        rm(list = c("co_int", "meta_int", "aseq_int"))
        saveRDS(object = a.int, file = "../outs/Peaks.int.RDS")
        write_log(f = "../scAPA.script.log",
                  stage = "featureCounts for introns")
    }
    rm(list = c("counts"))
    system(command = "rm -r ./CellBams")
}

# Peak filtering ----------------------------------------------------------
write_log_start("Stage 3: Peak filtering (3'UTRs)\n\n",
                command = NA, f = "../scAPA.script.log")
if (sc) a <- calc_clusters_counts(a)
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
if(sc){
cells_pui <- as.data.frame(results@ppui.cells)
cells_pui$cell <- gsub(pattern = "_proximal_PUI", replacement = "",
                       x = rownames(cells_pui))
cells_pui$sample <- gsub(pattern = "_.*", replacement = "",
                         x = cells_pui$cell)
cells_pui$cell <- gsub(pattern = ".*_", replacement = "",
                       x = cells_pui$cell)
cells_pui <- cells_pui[, c(3, 2, 1)]
colnames(cells_pui) <- c("Sample", "Cell_BC", "Mean_Proximal_PUI")
write.table(x = cells_pui, file = "./outs/Mean.Cell.PPUI.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)
}
write.table(x = results.output, file = "./outs/APA.events.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

write.table(x = a.fil@row.Data, file = "./outs/ThreeUTR.peaks.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

write.peaks(.x = results, .f = "./outs/UTRs.with.multiple.peaks.txt")

file.create("./outs/summary.UTR.txt")
cat("Number of peaks passed fillter:\t", nrow(a.fil@clus.counts),
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

# Introns -----------------------------------------------------------------
# Peak filtering ----------------------------------------------------------
write_log_start("Stage 3: Peak filtering (introns)\n\n",
                command = NA, f = "./scAPA.script.log")
if (int) {
    if (sc) a.int <- calc_clusters_counts(a.int)
    results.int <- set_scAPAresults(x = a.int, int = T, cpm = ICPM.cuttoff,
                                  counts = IC.cuttoff)

    # Internal priming filtering
    IP.seq.I <- paste0(rep("A", times = IA.number), collapse = "")
    results.int <- filter_IP(x = results.int, int.priming.seq = IP.seq.I,
                          left = Ifilter.border.left,
                          right = Ifilter.border.right)

    # Statistical analysis ----------------------------------------------------
    write_log_start("Stage 4: Statistical analysis (introns)\n\n",
                    command = NA, f = "./scAPA.script.log")
    # Chi-squre test for APA
    results.int <- test_APA(results.int, clus = "all")
    # Inferring global trends -------------------------------------------------
    write_log_start("Stage 5: Inferring global trends (introns)\n\n",
                    command = NA, f = "./scAPA.script.log")
    results.int <- calc_p_pui_mat(results.int)
    saveRDS(object = results.int, file = "./outs/results_introns.RDS")

    # Write final outputs -----------------------------------------------------
    results.int.output <- disply_results(x = results.int, org = org, int = T)
    results.int <- annotate_results(results.int, org = org)
    if(sc){
    cells_pui <- as.data.frame(results.int@ppui.cells)
    cells_pui$cell <- gsub(pattern = "_proximal_PUI", replacement = "",
                           x = rownames(cells_pui))
    cells_pui$sample <- gsub(pattern = "_.*", replacement = "",
                             x = cells_pui$cell)
    cells_pui$cell <- gsub(pattern = ".*_", replacement = "",
                           x = cells_pui$cell)
    cells_pui <- cells_pui[, c(3, 2, 1)]
    colnames(cells_pui) <- c("Sample", "Cell_BC", "Mean_intronic_PUI")
    write.table(x = cells_pui, file = "./outs/Intronic.Mean.Cell.IPUI.txt",
                quote = F, sep = "\t", col.names = T, row.names = F)
    }
    write.table(x = results.int.output, file = "./outs/Intronic.APA.events.txt",
                quote = F, sep = "\t", col.names = T, row.names = F)
    write.table(x = results.int@metadata, file = "./outs/Intronic.peaks.txt",
                quote = F, sep = "\t", col.names = T, row.names = F)
    file.create("./outs/summary.Introns.txt")
    cat("Number of peaks passed fillter:\t", length(results.int@clus.counts),
        "\nNumber of significant (FDR < 5%) APA events:",
        sum(results.int@pvalues[[1]][, 2] < 0.05, na.rm = T),
        "\n", file = "./outs/summary.Introns.txt")

    # Plot --------------------------------------------------------------------
    sig.utrs <- as.vector(results.int@pvalues$all[,2] < 0.05)
    sig <- results.int[which(sig.utrs),]
    tidy.ppui <- as.data.frame(sig@ppui.clus)
    colnames(tidy.ppui) <- gsub(x = colnames(tidy.ppui),
                                pattern = "_proximal_PUI", replacement = "")
    tidy.ppui <- tidyr::gather(data = tidy.ppui)
    colnames(tidy.ppui) <- c("Cluster", "value")
    pdf("./outs/Intronic.PUI.ECDF.pdf")
    p <- ggplot2::ggplot(data = tidy.ppui,
                         ggplot2::aes(x = value, color = Cluster))
    p <- p + ggplot2::stat_ecdf(size = 1)
    p <- p + ggplot2::theme_bw()
    p <- p + ggplot2::xlab("Intronic peak usage index")
    p <- p + ggplot2::ylab("Cumulative fraction")
    print(p)
    dev.off()
}
# End of script -----------------------------------------------------------
# Removing remining temporary files
write(x = cat(Sys.time(), "\tscAPA.script.R Finished.\nOutputs are in the ",
              "directory:\n", path.to.files, "/scAPA/outs"),
      file = "./scAPA.script.log", append = T)
cat("scAPA.script.R Finished.\nOutputs are in the directory:\n",
    path.to.files, "./scAPA/outs")
