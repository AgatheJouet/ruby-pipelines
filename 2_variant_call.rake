#!/usr/bin/env rake
#Agathe Jouet - 14/11/2017 - modified 09/07/2018
#usage: source ruby-2.3.1 && rake -f scripts/2_variant_call.rake [:task] fasta=path/to/fasta vcf_basename=vcf_basename loci=loci_type [--trace]

ENV["fasta"] ? @fasta = ENV["fasta"] : nil
ENV["vcf_basename"] ? @vcf_basename = ENV["vcf_basename"] : nil
ENV["loci"] ? @loci = ENV["loci"] : nil
READ1 = FileList['reads/*_1.fq.gz']
BAM = READ1.pathmap("bam/raw/#{@loci}/%n")
FILTERED_BAM = BAM.pathmap("%{raw,filtered}p_best_lte0.25-clipped.bam")


file "fasta/#{@loci}_list.txt" => ["#{@fasta}"] do
        sh "grep '>' #{@fasta} | sed 's/>//g' - > fasta/#{@loci}_list.txt"
end

desc "Creating list of loci"
task :loci_list => ["fasta/#{@loci}_list.txt"] do
        puts "List of loci created"
end


directory "vcf/#{@loci}"

file "vcf/#{@loci}/#{@vcf_basename}_all.vcf" => ["vcf/#{@loci}", "fasta/#{@loci}_list.txt"] do
        bam_to_merge = FILTERED_BAM.join(' ')
        sh "source samtools-1.3.1 && source bcftools-1.3.1 && samtools mpileup -Buf #{@fasta} #{bam_to_merge} --output-tags DP | bcftools call -mv - | bcftools view - -o vcf/#{@loci}/#{@vcf_basename}_all.vcf"
end

desc "Merging BAMs into VCF"
task :merge_bam => ["vcf/#{@loci}/#{@vcf_basename}_all.vcf"] do
        puts "BAMs merged into VCF"
end


file "vcf/#{@loci}/#{@vcf_basename}_all_d10.recode.vcf" => ["vcf/#{@loci}/#{@vcf_basename}_all.vcf"] do
	sh "source vcftools-0.1.14 && vcftools --vcf vcf/#{@loci}/#{@vcf_basename}_all.vcf --minDP 10 --out vcf/#{@loci}/#{@vcf_basename}_all_d10 --recode"
end

desc "Changing sites with read depth lower than 10x to missing"
task :filter_low_depth => ["vcf/#{@loci}/#{@vcf_basename}_all_d10.recode.vcf"] do
	puts "Sites with low depth changed to missing"
end


file "vcf/#{@loci}/#{@vcf_basename}_all_d10_no-indels.recode.vcf" => ["vcf/#{@loci}/#{@vcf_basename}_all_d10.recode.vcf"] do
	sh "source vcftools-0.1.14 && vcftools --vcf vcf/#{@loci}/#{@vcf_basename}_all_d10.recode.vcf --recode --out vcf/#{@loci}/#{@vcf_basename}_all_d10_no-indels --remove-indels"
end

desc "Filtering indels"
task :filter_indels => ["vcf/#{@loci}/#{@vcf_basename}_all_d10_no-indels.recode.vcf"] do
	puts "Indels filtered"
end


file "vcf/#{@loci}/#{@vcf_basename}_all_d10_no-indels_no-missing.recode.vcf"  => ["vcf/#{@loci}/#{@vcf_basename}_all_d10_no-indels.recode.vcf"] do
        sh "source vcftools-0.1.14 && vcftools --vcf vcf/#{@loci}/#{@vcf_basename}_all_d10_no-indels.recode.vcf --max-missing 1 --recode --out vcf/#{@loci}/#{@vcf_basename}_all_d10_no-indels_no-missing"
end

desc "Removing sites with missing genotypes in at least one individual"
task :filter_missing => ["vcf/#{@loci}/#{@vcf_basename}_all_d10_no-indels_no-missing.recode.vcf"] do
        puts "Missing values filtered out from VCFs"
end


desc "default"
task :default => :filter_missing
