#!/usr/bin/perl -w

#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
#
# grab_1kg_exome_pop_MAFs.pl - make a table of MAFs by 1000Genomes popuation
# by Jason Corneveaux, Bioinformatician, TGen Neurogenomics
# 12/30/2011
#
#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==

$pops = "1kg.samples";

$REGION = $ARGV[0];

# We need to strip chr prefix if its there to match 1kg
# We then check if a whole chromosome or region is specifed
# and split up as needed. Need $CHR for the download link
# and $REGION  for naming

$REGION =~ s/chr//gi; 
if ($REGION =~ /:/g) {
	$CHR = $` if $REGION =~ /:/;
} else {
	$CHR = $REGION;
}

# Build a hash of populations and genders for each sample
open(POPS,"<$pops");
while(<POPS>) {
	chomp($_);
	@data = split("\t",$_);
	$pop_hash{$data[1]} = $data[0];
	$sex = $data[4];
	$sex =~ s/female/2/g;
	$sex =~ s/male/1/g;
	$sex_hash{$data[1]} = $sex;
}
close(POPS);

# Get a list of the uniq population codes.
@uniq = keys %pop_hash;
foreach $pop (@uniq) {
	push @pops, $pop_hash{$pop};
}
%seen = ();
@uniqpops = grep { ! $seen{$_} ++ } @pops;

# The workhors routine
get_data();

sub get_data {
	system("tabix -fh ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.chr$CHR.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.gz $REGION > ALL.chr$CHR.phase1_integrated_calls.20101123.$REGION.vcf");
	system("vcftools --vcf ALL.chr$CHR.phase1_integrated_calls.20101123.$REGION.vcf --plink-tped --out 1000G_integrated_calls.20101123.$REGION");
	
	system("rm ALL.chr$CHR.phase1_integrated_calls.20101123.$REGION.vcf*");
	system("cp 1000G_integrated_calls.20101123.$REGION.tfam 1000G_integrated_calls.20101123.$REGION.tfam.orig");
	
	# Add population codes and genders to the tfam file
	open(FAMo,"<1000G_integrated_calls.20101123.$REGION.tfam.orig");
	open(FAMn,">1000G_integrated_calls.20101123.$REGION.tfam");
	while(<FAMo>) {
		chomp($_);
		@data = split("\t",$_);
		$fam = $pop_hash{$data[1]};
		$sex = $sex_hash{$data[1]};
		print FAMn "$fam\t$data[1]\t$data[2]\t$data[3]\t$sex\t$data[5]\n";
	
	}
	close(FAMo);
	close(FAMn);

	# Use PLINK to get the MAFs within each population
	foreach $pop (@uniqpops) {
		print "$pop\n";
		system("grep $pop 1000G_integrated_calls.20101123.$REGION.tfam | cut -f1,2 > $pop.samples");
		system("plink --tfile 1000G_integrated_calls.20101123.$REGION --noweb --keep $pop.samples --freq --out $pop.$REGION");
		system("rm $pop.$REGION.log $pop.samples");
	}

	# Use R out of laziness to make a master table for the MAFs by populations 
	system("Rscript --vanilla make_MAF_table.R 1000G_integrated_calls.20101123.$REGION.tped");
	system("rm 1000G_integrated_calls.20101123.$REGION.t* 1000G_integrated_calls.20101123.$REGION.log");
	system("rm *.$REGION.frq");
}

