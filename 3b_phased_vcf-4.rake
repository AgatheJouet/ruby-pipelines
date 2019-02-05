#!/usr/bin/env rake
#Agathe Jouet - 29/12/2017
#3b_phased_vcf-4.rake

ENV["fasta"] ? @fasta = ENV["fasta"] : nil
ENV["vcf"] ? @vcf = ENV["vcf"] : nil
ENV["loci"] ? @loci = ENV["loci"] : nil
ENV["nb_loci"] ? @nb_loci = ENV["nb_loci"] : nil
BAM=FileList["bam/filtered/#{@loci}/*.bam"]
BAM_BASENAME = BAM.pathmap("%f")
KEPT_LOCI = File.readlines("vcf/#{@loci}/single_loci/#{@loci}_10_indiv.txt").each {|kept_locus| kept_locus.chomp!}


KEPT_LOCI.each do |kept_locus|
        file "vcf/#{@loci}/single_loci/#{kept_locus}.imiss-indv" => ["vcf/#{@loci}/single_loci/#{kept_locus}.imiss"] do
                sh "awk '$5 == 1 {print $1}' vcf/#{@loci}/single_loci/#{kept_locus}.imiss | sed 's/^bam/--remove-indv bam/g' | awk '$1=$1' ORS=' ' > vcf/#{@loci}/single_loci/#{kept_locus}.imiss-indv"
        end
end

desc "Printing remove-indv for vcftools"
multitask :remove_indv_print => KEPT_LOCI.map{|x| "vcf/#{@loci}/single_loci/" + x + ".imiss-indv"} do
        puts "Remove-indv printed"
end


KEPT_LOCI.each do |kept_locus|
        file "vcf/#{@loci}/single_loci/#{kept_locus}_no_miss.recode.vcf" => ["vcf/#{@loci}/single_loci/#{kept_locus}.imiss-indv"] do
                remove_indv = File.open("vcf/#{@loci}/single_loci/#{kept_locus}.imiss-indv").map(&:chomp).join(" ")
                sh "source vcftools-0.1.14 && vcftools --vcf vcf/#{@loci}/single_loci/#{kept_locus}.recode.vcf #{remove_indv} --out vcf/#{@loci}/single_loci/#{kept_locus}_no_miss --recode"
        end
end

desc "Removing completely missing individuals from nlr vcf"
multitask :remove_missing_indv => KEPT_LOCI.map{|x| "vcf/#{@loci}/single_loci/" + x + "_no_miss.recode.vcf"} do
        puts "Missing individuals removed from nlr vcfs"
end


directory("shapeit/bam_list/#{@loci}")

KEPT_LOCI.each do |kept_locus|
        BAM_BASENAME.each do |bam|
                file "shapeit/bam_list/#{@loci}/#{kept_locus}.txt" => ["shapeit/bam_list/#{@loci}", KEPT_LOCI.map{|x| "vcf/#{@loci}/single_loci/" + x + "_no_miss.recode.vcf"}].flatten do
                        sh "echo -e 'bam/filtered/#{@loci}/#{bam} bam/filtered/#{@loci}/#{bam} #{kept_locus}' >> shapeit/bam_list/#{@loci}/#{kept_locus}.txt"
                end
        end
end


desc "Creating one bam list per locus"
multitask :bam_list => KEPT_LOCI.map{|x| "shapeit/bam_list/#{@loci}/" + x + ".txt"} do
        puts "bam list created for each locus"
end


directory("shapeit/extract_pirs/#{@loci}")

KEPT_LOCI.each do |kept_locus|
        file "shapeit/extract_pirs/#{@loci}/#{kept_locus}_extract-pirs" => ["shapeit/extract_pirs/#{@loci}", KEPT_LOCI.map{|x| "shapeit/bam_list/#{@loci}/" + x + ".txt"}].flatten do
                sh "source shapeit-2.20 && source samtools-1.3.1 && extractPIRs --bam shapeit/bam_list/#{@loci}/#{kept_locus}.txt --vcf vcf/#{@loci}/single_loci/#{kept_locus}_no_miss.recode.vcf --out shapeit/extract_pirs/#{@loci}/#{kept_locus}_extract-pirs"
        end
end

desc "Extracting PIRs with Shapeit"
multitask :extract_pirs => KEPT_LOCI.map{|x| "shapeit/extract_pirs/#{@loci}/" + x + "_extract-pirs"} do
        puts "PIRs extracted"
end


directory("shapeit/phased_files/#{@loci}")

KEPT_LOCI.each do |kept_locus|
        file "shapeit/phased_files/#{@loci}/#{kept_locus}.hap" => ["shapeit/phased_files/#{@loci}", KEPT_LOCI.map{|x| "shapeit/extract_pirs/#{@loci}/" + x + "_extract-pirs"}].flatten do
                sh "source samtools-1.3.1 && source shapeit-2.20 && shapeit -assemble --input-vcf vcf/#{@loci}/single_loci/#{kept_locus}_no_miss.recode.vcf --input-pir shapeit/extract_pirs/#{@loci}/#{kept_locus}_extract-pirs --states 1000 --burn 60 --prune 60 --main 300 --effective-size 1000000 --window 0.5 -O shapeit/phased_files/#{@loci}/#{kept_locus}.hap --output-grap shapeit/phased_files/#{@loci}/#{kept_locus}.hap.graph"
        end
end

desc "Running ShapeIt"
multitask :run_shapeit => KEPT_LOCI.map{|x| "shapeit/phased_files/#{@loci}/" + x + ".hap"} do
        puts "ShapeIt has run - prepare for conversion to vcf"
end


KEPT_LOCI.each do |kept_locus|
        file "shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf" => ["shapeit/phased_files/#{@loci}/#{kept_locus}.hap", KEPT_LOCI.map{|x| "shapeit/phased_files/#{@loci}/" + x + ".hap"}].flatten do
                sh "source samtools-1.3.1 && source shapeit-2.20 && shapeit -convert --input-haps shapeit/phased_files/#{@loci}/#{kept_locus}.hap --output-vcf shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf"
        end
end

desc "ShapeIt outfiles to vcf conversion"
multitask :shapeit2vcf => KEPT_LOCI.map{|x| "shapeit/phased_files/#{@loci}/" + x + ".hap.vcf"} do
        puts "ShapeIt output files converted back to vcf"
end


desc "default task"
task :default => :shapeit2vcf
