#!/usr/bin/env rake
#Agathe Jouet - 29/12/2017 - modified 21/06/2018 (removed "albugo instaces")
#3b_phased_vcf-2.rake

ENV["fasta"] ? @fasta = ENV["fasta"] : nil
ENV["vcf"] ? @vcf = ENV["vcf"] : nil
ENV["loci"] ? @loci = ENV["loci"] : nil
ENV["nb_loci"] ? @nb_loci = ENV["nb_loci"] : nil
LOCI=File.readlines("fasta/#{@loci}_list.txt").each {|locus| locus.chomp!}
BAM=FileList["bam/filtered/#{@loci}/*.bam"]

directory("vcf/#{@loci}/single_loci")

LOCI.each do |locus|
	file "vcf/#{@loci}/single_loci/#{locus}.recode.vcf" => ["vcf/#{@loci}/single_loci"] do
		sh "source vcftools-0.1.14 && vcftools --vcf #{@vcf}_biallelic_only_no-fully-missing-snps.recode.vcf --chr #{locus} --recode --out vcf/#{@loci}/single_loci/#{locus}"
	end
end


desc "Splitting VCF"
task :split_vcf => LOCI.map {|x| "vcf/#{@loci}/single_loci/" + x + ".recode.vcf"} do
	puts "VCF split"
end


file "vcf/#{@loci}/single_loci/#{@loci}_1_snp.txt" => ["#{@vcf}_biallelic_only_no-fully-missing-snps.recode.vcf", :split_vcf] do
        sh "grep '^[^#]' #{@vcf}_biallelic_only_no-fully-missing-snps.recode.vcf | cut -f 1 | sort -u > vcf/#{@loci}/single_loci/#{@loci}_1_snp.txt"
end


desc "Print kept loci (loci that have variants)"
task :print_kept_loci => ["vcf/#{@loci}/single_loci/#{@loci}_1_snp.txt"] do
        puts "Loci with at least one variant printed"
end


BAM.each do |bam|
        file "#{bam}.bai" => ["vcf/#{@loci}/single_loci/#{@loci}_1_snp.txt"] do
                sh "source samtools-1.3.1 && samtools index #{bam}"
        end
end

desc "Indexing BAMs"
multitask :index_bam => BAM.map{|x| x + ".bai"} do
        puts "BAMs indexed"
end


desc "Calling rakefile for treating missing individuals"
task :missing_indiv_rakefile => [:index_bam] do
        sh "source ruby-2.3.1 && rake -f scripts/3b_phased_vcf-3.rake fasta=#{@fasta} vcf=#{@vcf} loci=#{@loci} -j #{@nb_loci} --trace"
end


desc "default task"
task :default => :missing_indiv_rakefile
