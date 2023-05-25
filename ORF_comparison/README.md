# ORF_comparison

The script `orf_callers_comparison.Rmd` can be used to recreate the comparisons of the ORF callers by number of ORFs detected, ORF lengths and ORF categories. The sessioninfo associated with this R script is stored in `comparison_sessioninfo.txt`.

Before usage, the file path strings and file path vectors under the heading `# Define directories, files and colors` have to be initialized as documented. Afterwards, consecutive execution of all code blocks should recreate the analyses. Please note that for PRICE, it is expected that the .tsv output file is present in the same directory as the .cit.bed output file.