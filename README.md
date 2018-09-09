# Core-genome-size-estimation
Estimate the core genome size for a group of isolates 

A Perl script that estimates the cluster core genome size across a group of isolates, for a given depth of coverage, by assessing BAM files. An optional parameter in which the user also provides a BED file of masked regions on the reference genome will result in disregarding regions that do not harbor phylogenetically informative hqSNPs (i.e., regions of homologous recombination). This Perl script can be used with output from any SNP pipeline that generates BAM files. The script will operate post-SNP calling to support isolate relatedness conclusions. The script needs [samtools](https://github.com/samtools/), [bedtools](http://bedtools.readthedocs.io/en/latest/), and awk.

### Usage
    perl estimate_core_genome_from_bam.pl -bam /path/to/bam/files -genome mapping_reference.fasta -depth 10 -out output_folder
    
    OR
    
    perl estimate_core_genome_from_bam.pl -bam /path/to/bam/files -genome mapping_reference.fasta -depth 10 -out output_folder -masked masked_regions.bed 

