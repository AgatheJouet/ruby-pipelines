#!/usr/bin/env rake
#Agathe Jouet - 13/11/2017 - modified 21/06/2018 adaptation for rice
#usage: source ruby-2.3.1 && rake [:task] -f scripts/1_mapping_to_targets.rake fasta=path/to/fasta loci=loci_type -j $nb_bam [--trace]

ENV["fasta"] ? @fasta = ENV["fasta"] : nil
ENV["loci"] ? @loci = ENV["loci"] : nil
READ1 = FileList['reads/*_1.fq.gz']
READ2 = FileList['reads/*_2.fq.gz']
BAM = READ1.pathmap("bam/raw/#{@loci}/%n")
FILTERED_BAM = BAM.pathmap("%{raw,filtered}p_best_lte0.25-clipped.bam")

file "#{@fasta}.amb" => ["#{@fasta}"] do
        sh "source bwa-0.7.4; bwa index #{@fasta}"
end                                                                                                                                                  

file "#{@fasta}.ann" => ["#{@fasta}"] do
end

file "#{@fasta}.bwt" => ["#{@fasta}"] do
end

file "#{@fasta}.pac" => ["#{@fasta}"] do
end

desc "Indexing the gene fasta file"
task :index_fasta => ["#{@fasta}.amb", "#{@fasta}.ann", "#{@fasta}.bwt", "#{@fasta}.pac"] do
        puts "Fasta indexed"
end


directory "bam/raw/#{@loci}"

[BAM, READ1, READ2].transpose.each do |bam, read1, read2|
        file bam => ["bam/raw/#{@loci}", read1, read2, :index_fasta] do
                sh "source bwa-0.7.4; source samtools-1.3.1; bwa mem -T 0 #{@fasta} #{read1} #{read2} | samtools view -hbS - | samtools sort -n - -o #{bam}"
        end
end

desc "mapping reads to reference genes"
multitask :mapping => BAM do
        puts "Mapping completed"
end


directory "bam/filtered/#{@loci}"

FILTERED_BAM.zip(BAM).each do |f_bam, bam|
        file f_bam => ["bam/filtered/#{@loci}", :mapping] do
                sh "source ngsutils-0.5.7 && source samtools-1.3.1 && bamutils best #{bam} #{bam}_best.bam && bamutils removeclipping #{bam}_best.bam #{bam}_best_clip-info.bam && rm #{bam}_best.bam && bamutils filter #{bam}_best_clip-info.bam #{bam}_best_lte0.25.bam -lte ZC 0.25 && rm #{bam}_best_clip-info.bam && samtools sort #{bam}_best_lte0.25.bam > #{f_bam} && rm #{bam}_best_lte0.25.bam"
        end
end                                                                                                                                                  

desc "Filtering BAM"
multitask :filter_bam => FILTERED_BAM do
        puts "BAM filtered"
end


desc "default"
task :default => :filter_bam
