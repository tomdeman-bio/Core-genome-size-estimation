# Core-genome-size-estimation
Estimate the core genome size for a group of isolates 

A Perl script that estimates the core genome size of a group of isolates, for a given depth of coverage, by assessing read pileups in BAM files. Optionally, a BED file containing masked regions located on the reference genome can be used as input to predict a more accurate size of the core genome analyzed during identification of core SNPs. The script needs [samtools](https://github.com/samtools/), [bedtools](http://bedtools.readthedocs.io/en/latest/), and awk.

### Usage
    perl estimate_core_genome_from_bam.pl -bam /path/to/bam/files -genome mapping_reference.fasta -depth 10 > output.txt
    
    OR
    
    perl estimate_core_genome_from_bam.pl -bam /path/to/bam/files -genome mapping_reference.fasta -depth 10 -masked masked_regions.bed > output.txt
---------------------------------------------------------------------------------------------------------------------------------------
