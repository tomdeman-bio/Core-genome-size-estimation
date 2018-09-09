#! /usr/bin/perl

#Written by Tom de Man
#script needs sorted bam files 
#script also needs samtools, bedtools, and awk in order to operate
#bed files containing masked regions on the reference genome are optional

use strict;
use warnings;
use Getopt::Long;
use String::ShellQuote qw(shell_quote);
use Array::Utils qw(:all);

#global variables
my $sortbam_path;
my $genome;
my $depth;
my $out;
my $bed;

my @bams;
my $genome_size;
my $genomes;
my $coregenome_size;
my $relative;
my $core_pos;
my $masked_pos;
my $masked_regions;
my @core_masked;
my $handle;

my $starttime = localtime;
my $version = "1.2";

GetOptions(	"bam=s"		=> \$sortbam_path,
		"genome=s"	=> \$genome,
		"depth=s"	=> \$depth,
		"out=s"		=> \$out,
		"masked=s"	=> \$bed,
		"help|?"	=> sub {Usage()}
		);

if (($sortbam_path) && ($genome) && ($depth) && ($out)) {
	@bams = &get_files("bam");
	print STDERR "Hi $ENV{USER}, you are now running $0 version: $version on $starttime \n\n";
	
	#create output folder
	mkdir "$out";
	open ($handle, ">$out/" . "Core_genome_size_stats.txt");
	
	print $handle "$0 version: $version \n\n";
	print $handle "BAM files included in core genome estimate: \n";
	foreach (@bams) {
		print $handle "$_\n";
	}
	
	$genome_size = &genome_size_calc;
	
	($genomes, $coregenome_size, $relative, $core_pos) = &estimate_core(@bams);
	print $handle "Core genome size for $genomes genomes is: $coregenome_size base pairs, which equals $relative% of the mapping reference genome\n\n";
} else {
	&Usage;
}

#parsing optional parameters	
if ($bed) {
	($masked_pos, $masked_regions) = &subtract_masked_regions($bed);
	#masked positions from the BED file, often phage DNA or read mapping cliffs
	@core_masked = intersect(@$masked_pos, @$core_pos);
	my $core_masked_len = scalar(@core_masked);
	print $handle "Total masked regions: $masked_regions bp\n"; 
	print $handle "Total masked regions in core genome: $core_masked_len bp \n\n";
	my $coregenome_nocoremasked = $coregenome_size - $core_masked_len;
	
	my $percentage = ($coregenome_nocoremasked/$genome_size)* 100;
	my $rounded = sprintf "%.2f", $percentage;
	
	print $handle "Core genome size excluding masked regions for $genomes genomes is: $coregenome_nocoremasked base pairs, which equals $rounded% of the mapping reference genome\n";
}

sub genome_size_calc {
	print STDERR "\n";
	print $handle "your mapping reference is $genome \n\n";
	open FASTA, "$genome" or die "cannot open $genome for reading \n";
	my $total_bases = 0;
	while (<FASTA>) {
		if (!(/^>/)) {
			chomp;
			s/\r//g;
			$total_bases += length;
		}
	}
	return $total_bases;
	close FASTA;
}

sub estimate_core {
	#create fasta index
	print STDERR "running samtools..... generating a contig list file\n";
	system("samtools faidx $genome");
	open INDEX, "$genome.fai" or die "cannot open $genome.fai for reading \n";
	open (my $fh, '>', "$genome.contig");
	while (<INDEX>) {
		chomp;
		my @split = split ("\t", $_);
		print $fh "$split[0]\t$split[1]\n";
	}
	close $fh;
	close INDEX;

	my $cnt = 0;
	#calculate genome wide coverage for each nucleotide position, for each sorted BAM file
	foreach my $file (@_) {
		$cnt += 1;
		my $bam = "$sortbam_path/$file";
		my $con_len = "$genome.contig";
		my $out = "$sortbam_path/$file.$depth.cov";
		print STDERR "running bedtools for sample $cnt..... generating genome coverage data \n";
		system("bedtools genomecov -ibam ".shell_quote($bam)." -g ".shell_quote($con_len)." -d | awk '{if(\$3>=$depth){ print \$0}}' > ".shell_quote($out)."");
	}
	
	#get the .cov files	
	my @covs = &get_files("cov");
	open (my $fh2, ">$out/" . "Pairwise_core_genome_sizes.txt");
	my @cov2d;
	foreach (@covs) {
		my $countl = 0;
		my @n;
		my $c = "$sortbam_path/$_";
		my @cov_split = split ("\\.", $_);
		open COV, $c or die "cannot open $c for reading \n";
		while (<COV>) {
			chomp;
			$countl +=1;
			my @split = split ("\t", $_);
			my $contig_position = $split[0]."_".$split[1];
			push @n, $contig_position;
		}
		close COV;
		push @cov2d, \@n;
		print $fh2 "reference"."\t"."$cov_split[0]"."\t"."$countl"."\n";
 
	}
	close $fh2;
	
	my $rows = scalar @cov2d;
	print STDERR "You are going to estimate a core genome for $rows isolates ..... how exciting!!! \n";	

	my $start_ref = shift @cov2d;
	my @overlap = @$start_ref;
	my $remainder = $rows - 1;
	#pairwise comparison via intersect
	for (my $i=0; $i < $remainder; $i++) {
		my $comparison = shift @cov2d;
		@overlap = intersect(@$comparison, @overlap);
	}
	my $core = scalar @overlap;
	my $percentage = ($core/$genome_size)* 100;
	my $rounded = sprintf "%.2f", $percentage;
	
	return $rows, $core, $rounded, \@overlap;

}

sub subtract_masked_regions {
	my $bed_file = shift;
	my $bp_regions = 0;
	my @bed_pos;
	open BED, "$bed_file" or die "cannot open $bed_file for reading \n";
	while (<BED>) {
		chomp;
		my @split = split ("\t", $_);
		my $size = $split[2]-$split[1];
		for (my $i=$split[1]; $i <= $split[2]; $i++) {
			my $pos = "$split[0]"."_"."$i";
			push @bed_pos, $pos;
		}
		$bp_regions += $size;
	}
	
	return \@bed_pos, $bp_regions;
	close BED;
}

sub get_files {
	my $ext = qr/$_[0]/;
	my @bamfiles;
	opendir(DIR, $sortbam_path) or die "cannot open $sortbam_path \n";
	my @files = readdir(DIR);
	close DIR;

	foreach my $file (@files){
		next if (!($file =~ /\.$ext$/));
		push @bamfiles, $file;
	}
	return @bamfiles;
}

sub Usage {
	print STDERR "\n Please provide mandatory input files and values (-bam, -genome, -depth)!!!\n\n";
	print STDERR "\n Usage:  perl $0 -bam <BAM file path> -genome <genome FASTA file> -depth <minimum depth of coverage to include in output> -out <output folder>\n\n";
	print STDERR "\n optional input: -masked <bed file> \n";
	exit;
}