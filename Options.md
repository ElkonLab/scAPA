Options
================

### Basic Options

-   **-p** The path to the directory with the bam files and cell cluster annotations files to be analyzed (required)

-   **-org** Organism, either Mm for mouse, or Hs for human (required).

-   **-scp** The path to the directory of scAPA.shell.script (required). Required for reading the configuration file.

-   **-c** The number of cores to use. The default is 30.

### Output Optios

-   **-wig** Weather to generate cluster wig files. The default is false.

-   **-ChangePoint** weather to use Change Point. The default is false.

### Analysis Options

-   **-sc** If true, (default) counts read from individual cells (slower). Otherwise, counts read from clusters.

-   **-int** If true (default) performs intronic APA analysis as well as 3'UTR analysis

### Filtering Options for 3'UTRs analysis

-   **-cpm** Consider only peaks with more than a total sum of CPMs over all cell clusters larger than -cpm. Default value: 10

-   **-a** Numeric value, filter out peaks with -a consecutive As in the region -u to -d downstream their 3' edge. Default value: 8

-   **-u** Numeric value, filter out peaks with -a consecutive As in the region -u to -d downstream their 3' edge. Default value: 10

-   **-d** Numeric value, filter out peaks with -a consecutive As in the region -u to -d downstream their 3' edge. Default value: 140

### Filtering Options for intron analysis

-   **-Icpm** Consider only peaks with more than a total sum of CPMs over all cell clusters larger than -Icpm.

-   **-Ico** Consider only peaks with more than a total sum of counts over all cell clusters larger than -Ico. Default value: 50

-   **-Ia** Numeric value, filter out peaks with -Ia consecutive As in the region -Iu to -Id downstream thire 3' edge Default value: 7

-   **-Iu** Numeric value, filter out peaks with -a consecutive As in the region -Iu to -Id downstream their 3' edge. Default value: 1

-   **-Id** Numeric value, filter out peaks with -Ia consecutive As in the region -Iu to -Id downstream their 3' edge. Default value: 200
