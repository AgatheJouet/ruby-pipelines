#!/usr/bin/env rake
#Agathe Jouet - 06/02/2018 - modified 11/07/2018
#usage: source ruby-2.3.1 && rake [:task] -f scripts/4b_phased_consensus.rake fasta=path/to/reference/fasta loci=loci_type -j nb_of_loci [--trace]
#Before using this script, you need to create a file that is necessary to order the loci fasta files into pops. See POPS below.

ENV["fasta"] ? @fasta = ENV["fasta"] : nil
ENV["loci"] ? @loci = ENV["loci"] : nil
ENV["nb_loci"] ? @nb_loci = ENV["nb_loci"] : nil
KEPT_LOCI = File.readlines("vcf/#{@loci}/single_loci/#{@loci}_10_indiv_no_LOC02g26500.txt").each {|kept_locus| kept_locus.chomp!}
INDEXED_VCF = KEPT_LOCI.map {|x| "shapeit/phased_files/#{@loci}/" + x + ".hap.vcf.gz"}
LOCI_REF_FASTA = KEPT_LOCI.map {|x| "fasta/per_gene/" + x + ".fas"}
HAP1 = KEPT_LOCI.map {|x| "fasta/#{@loci}/hap1/" + x + "_hap1.fas"}
HAP2 = KEPT_LOCI.map {|x| "fasta/#{@loci}/hap2/" + x + "_hap2.fas"}
HAPS = KEPT_LOCI.map {|x| "fasta/#{@loci}/haps/" + x + "_haps.fas"}
HAPS_UNWRAPPED = KEPT_LOCI.map {|x| "fasta/#{@loci}/haps/" + x + "_haps_unwrapped.fas"}
POPS = File.readlines("scripts/rice_pops.txt").map {|x| x.chomp!}.map {|x| x.split(',')}.group_by(&:first).map{|key, value| [key, value.map(&:last)]}.to_h         #Needs to be created as a comma-separated table with two columns: the first with the name of the pop, the second with the name of the sample followed by an underscore "_"
POPS_KEY = POPS.map {|k, v| k}


KEPT_LOCI.each do |kept_locus|
	file "shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf.gz" => ["shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf"] do
	        sh "source tabix-0.2.6 && bgzip -c shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf > shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf.gz && tabix -p vcf shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf.gz"
	end
end

desc "Indexing VCFs"
multitask :vcf_index => INDEXED_VCF do
        puts "VCFs indexed"
end


directory "fasta/per_gene/"

KEPT_LOCI.each do |kept_locus|
	file "fasta/per_gene/#{kept_locus}.fas" => ["#{@fasta}", "fasta/per_gene/", "shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf.gz"] do
		sh "grep -A 1 '>#{kept_locus}' #{@fasta} > fasta/per_gene/#{kept_locus}.fas"
	end
end

desc "Creating one reference fasta per locus"
multitask :vcf_index => LOCI_REF_FASTA do
        puts "Per gene reference fastas created"
end


directory "fasta/#{@loci}/hap1/"

KEPT_LOCI.each do |kept_locus|
	SAMPLES = File.readlines("shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf").select {|x| x =~ /#CHROM/}.map {|x| x.gsub(/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t/,'').gsub(/.bam\t/,'.bam ')}.map {|x| x.chomp!}.join(' ').split(' ')
	SAMPLES.each do |sample|
		file "fasta/#{@loci}/hap1/#{kept_locus}_hap1.fas" => ["fasta/#{@loci}/hap1/", "fasta/per_gene/#{kept_locus}.fas"] do
				sh "source tabix-0.2.6 && source bcftools-1.3.1 && echo #{sample} >> fasta/#{@loci}/hap1/#{kept_locus}_hap1.fas && bcftools consensus -f fasta/per_gene/#{kept_locus}.fas -s #{sample} -H 1 shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf.gz >> fasta/#{@loci}/hap1/#{kept_locus}_hap1.fas"
		end
	end
end

desc "Calling haplotypes 1 from phased VCFs"
multitask :hap1_call => HAP1 do
	puts "Haplotypes 1 called"
end


directory "fasta/#{@loci}/hap2/"

KEPT_LOCI.each do |kept_locus|
        SAMPLES = File.readlines("shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf").select {|x| x =~ /#CHROM/}.map {|x| x.gsub(/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t/,'').gsub(/.bam\t/,'.bam ')}.map {|x| x.chomp!}.join(' ').split(' ')
        SAMPLES.each do	|sample|
	        file "fasta/#{@loci}/hap2/#{kept_locus}_hap2.fas" => ["fasta/#{@loci}/hap2/", "fasta/per_gene/#{kept_locus}.fas", "fasta/#{@loci}/hap1/#{kept_locus}_hap1.fas"] do
                        sh "source tabix-0.2.6 && source bcftools-1.3.1 && echo #{sample} >> fasta/#{@loci}/hap2/#{kept_locus}_hap2.fas && bcftools consensus -f fasta/per_gene/#{kept_locus}.fas -s #{sample} -H 2 shapeit/phased_files/#{@loci}/#{kept_locus}.hap.vcf.gz >> fasta/#{@loci}/hap2/#{kept_locus}_hap2.fas"
        	end
	end
end

desc "Calling haplotypes 2 from phased VCFs"
multitask :hap2_call =>	HAP2 do
        puts "Haplotypes 2 called"
end


directory "fasta/#{@loci}/haps/"

[KEPT_LOCI, HAP1, HAP2].transpose.each do |kept_locus, hap1, hap2|
	file "fasta/#{@loci}/haps/#{kept_locus}_haps.fas" => ["fasta/#{@loci}/haps/", "#{hap1}", "#{hap2}"] do
		sh "bash scripts/sed_haps.sh #{hap1} 1 > #{kept_locus}_hap1.txt && bash scripts/sed_haps.sh #{hap2} 2 > #{kept_locus}_hap2.txt && cat #{kept_locus}_hap1.txt #{kept_locus}_hap2.txt > fasta/#{@loci}/haps/#{kept_locus}_haps.fas && rm #{kept_locus}_hap1.txt && rm #{kept_locus}_hap2.txt"
	end
end

desc "Creating one fasta per gene, including both haplotypes"
multitask :paste_haps => HAPS do
	puts "Haplotype 1 or 2 info added and files formatted to fasta, haps pasted in one fasta"
end


KEPT_LOCI.each do |kept_locus|
        file "fasta/#{@loci}/haps/#{kept_locus}_haps_unwrapped.fas" => ["fasta/#{@loci}/haps/#{kept_locus}_haps.fas"] do
                sh "source perl-5.16.3 && perl scripts/fast1line.pl fasta/#{@loci}/haps/#{kept_locus}_haps.fas fasta/#{@loci}/haps/#{kept_locus}_haps_unwrapped.fas"
        end
end

desc "Unwrapping fasta"
multitask :unwrap_fasta => HAPS_UNWRAPPED do
        puts "Fasta unwrapped"
end


directory "fasta/#{@loci}/haps/per_pop/"

KEPT_LOCI.each do |kept_locus|
        POPS.each_pair do |key, values|
		values.each do |value|
			file "fasta/#{@loci}/haps/per_pop/#{kept_locus}_#{key}.fas" => ["fasta/#{@loci}/haps/per_pop/", "scripts/rice_pops.txt", "fasta/#{@loci}/haps/#{kept_locus}_haps_unwrapped.fas"] do
				sh "if grep -q '#{value}' fasta/#{@loci}/haps/#{kept_locus}_haps_unwrapped.fas; then grep -A 1 '#{value}' fasta/#{@loci}/haps/#{kept_locus}_haps_unwrapped.fas >> fasta/#{@loci}/haps/per_pop/#{kept_locus}_#{key}.fas; else echo '--'; fi"
	                end
	        end
	end
end

desc "Creating sub-fasta files per race"
multitask :per_race_fasta => KEPT_LOCI.map{|x| "fasta/#{@loci}/haps/per_pop/" + x + "_"}.product(POPS_KEY).map(&:join).map{|y| y + ".fas"} do
        puts "Sub-fasta per pop created"
end


KEPT_LOCI.each do |kept_locus|
        POPS.each_key do |key|
		file "fasta/#{@loci}/haps/per_pop/#{kept_locus}-#{key}.fas" => ["fasta/#{@loci}/haps/per_pop/#{kept_locus}_#{key}.fas"] do
			sh "if [ -e fasta/#{@loci}/haps/per_pop/#{kept_locus}_#{key}.fas ]; then sed '/--/d' fasta/#{@loci}/haps/per_pop/#{kept_locus}_#{key}.fas > fasta/#{@loci}/haps/per_pop/#{kept_locus}-#{key}.fas; else echo 'no locus for this population'; fi && if [ -e fasta/#{@loci}/haps/per_pop/#{kept_locus}_#{key}.fas ]; then rm fasta/#{@loci}/haps/per_pop/#{kept_locus}_#{key}.fas; else echo 'nothing to delete'; fi"
		end
	end
end

desc "Formatting to final per locus, per race fasta files"
multitask :final_fasta_format => KEPT_LOCI.map{|x| "fasta/#{@loci}/haps/per_pop/" + x + "-"}.product(POPS_KEY).map(&:join).map{|y| y + ".fas"} do
	puts "Final fasta files formatted"
end
		 

desc "default task"
task :default => :final_fasta_format
