
## Files

step0_process_atacseq.sh : A shell file that launches the ATAC-seq pipeline to call peaks on raw, bulk human and mouse retina ATAC-seq data

step1_filter_peaks.sh : A shell file that performs the following filters on the retina optimal IDR peaks: removes peaks within 20 kB of a transcription start site, removes peaks overlapping ENCODE blacklist v2 regions, removes peaks overlapping protein coding exons, and removes peaks longer than 1 kB in length.
