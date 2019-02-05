#!/usr/bin/env rake
#Agathe Jouet - 29/12/2017 - modified 21/06/2018 for rice project
#usage: source ruby-2.3.1 && rake [:task] -f scripts/3b_phased_vcf-1.rake fasta=path/to/fasta vcf=path/to/vcf loci=loci_type nb_loci=nb_of_loci [--trace]
##vcf=path/to/vcf without ".recode.vcf" extension - just to have simpler names

ENV["fasta"] ? @fasta = ENV["fasta"] : nil
ENV["vcf"] ? @vcf = ENV["vcf"] : nil
ENV["loci"] ? @loci = ENV["loci"] : nil
ENV["nb_loci"] ? @nb_loci = ENV["nb_loci"] : nil


file "#{@vcf}.recode.vcf.gz" => ["#{@vcf}.recode.vcf"] do
        sh "source tabix-0.2.6 && bgzip -c #{@vcf}.recode.vcf > #{@vcf}.recode.vcf.gz && tabix -p vcf #{@vcf}.recode.vcf.gz"
end

desc "Indexing VCF"
task :vcf_index => ["#{@vcf}.recode.vcf.gz"] do
        puts "VCF indexed"
end


file "fasta/#{@loci}_list.txt" => ["#{@fasta}", "#{@vcf}.recode.vcf.gz"] do
        sh "grep '>' #{@fasta} | sed 's/>//g' - > fasta/#{@loci}_list.txt"
end

desc "Creating list of loci"
task :loci_list => ["fasta/#{@loci}_list.txt"] do
        puts "List of loci created"
end


file "#{@vcf}_biallelic_only.vcf" => [:loci_list, "#{@vcf}.recode.vcf"] do
        sh "source bcftools-1.3.1 && bcftools view --max-alleles 2 #{@vcf}.recode.vcf > #{@vcf}_biallelic_only.vcf"
end

desc "Removing non-biallelic SNPs before running Shapeit"
task :remove_non_biallelic => ["#{@vcf}_biallelic_only.vcf"] do
        puts "Non biallelic SNPs removed from VCF"
end


file "#{@vcf}_biallelic_only_no-fully-missing-snps.recode.vcf" => [:remove_non_biallelic] do
        sh "source vcftools-0.1.14 && vcftools --vcf #{@vcf}_biallelic_only.vcf --max-missing 0.02 --recode --out #{@vcf}_biallelic_only_no-fully-missing-snps"
end

desc "Removing fully missing SNPs before running Shapeit"
task :remove_missing_snps => ["#{@vcf}_biallelic_only_no-fully-missing-snps.recode.vcf"] do
        puts "Fully missing SNPs removed from VCF"
end


desc "Calling rakefile for phasing preparation (checking loci with variants and indexing bam files)"
task :phasing_prep_rakefile => [:remove_missing_snps] do
	sh "source ruby-2.3.1 && rake -f scripts/3b_phased_vcf-2.rake fasta=#{@fasta} vcf=#{@vcf} loci=#{@loci} -j #{@nb_loci} --trace"
end


desc "default task"
task :default => :phasing_prep_rakefile
