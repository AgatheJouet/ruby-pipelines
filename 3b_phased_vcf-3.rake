#!/usr/bin/env rake
#Agathe Jouet - 29/12/2017
#3b_phased_vcf-3.rake

ENV["fasta"] ? @fasta = ENV["fasta"] : nil
ENV["vcf"] ? @vcf = ENV["vcf"] : nil
ENV["loci"] ? @loci = ENV["loci"] : nil
ENV["nb_loci"] ? @nb_loci = ENV["nb_loci"] : nil
BAM=FileList["bam/filtered/#{@loci}/*.bam"]
BAM_BASENAME = BAM.pathmap("%f")
KEPT_LOCI = File.readlines("vcf/#{@loci}/single_loci/#{@loci}_1_snp.txt").each {|kept_locus| kept_locus.chomp!}


KEPT_LOCI.each do |kept_locus|
        file "vcf/#{@loci}/single_loci/#{kept_locus}.imiss" => ["vcf/#{@loci}/single_loci/#{@loci}_1_snp.txt"] do
                sh "source vcftools-0.1.14 && vcftools --vcf vcf/#{@loci}/single_loci/#{kept_locus}.recode.vcf --missing-indv --out vcf/#{@loci}/single_loci/#{kept_locus}"
        end
end

desc "Writing missing stats per locus"
multitask :missing_stats => KEPT_LOCI.map{|x| "vcf/#{@loci}/single_loci/" + x + ".imiss"} do
        puts "Missing stats reported"
end


KEPT_LOCI.each do |kept_locus|
	file "vcf/#{@loci}/single_loci/#{@loci}_10_indiv.txt" => ["vcf/#{@loci}/single_loci/#{kept_locus}.imiss"] do
		sh "bash scripts/sed_imiss.sh vcf/#{@loci}/single_loci/#{kept_locus}.imiss >> vcf/#{@loci}/single_loci/#{@loci}_10_indiv.txt"
	end
end

desc "Writing VCF filenames with at least 10 individuals"
task :ten_indiv_vcf => ["vcf/#{@loci}/single_loci/#{@loci}_10_indiv.txt"] do
        puts "VCF filenames with at least 10 individuals printed"
end


desc "Calling rakefile for phasing (phasing data for all loci)"
task :phasing_rakefile => [:ten_indiv_vcf] do
        sh "source ruby-2.3.1 && rake -f scripts/3b_phased_vcf-4.rake fasta=#{@fasta} vcf=#{@vcf} loci=#{@loci} -j #{@nb_loci} --trace"
end


desc "default task"
task :default => :phasing_rakefile
