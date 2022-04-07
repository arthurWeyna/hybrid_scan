
import os 
import numpy as np
import re


##PARAMETERS (feel free to edit)
#Path to directory for data files
DATADIR="/home/arthur/Bureau/postdocs/isem/hybscan/data"
#Path to tools directory
TOOLDIR="/home/arthur/Bureau/postdocs/isem/hybscan/hybrid_scan/tools"
#Path to file with ids and probes to use for each id (UCE runs)
UCETARGETS="/home/arthur/Bureau/postdocs/isem/hybscan/hybrid_scan/uce_list"
#Path to file with ids and database/augustus_species to use for each id (BUSCO runs)
BUSCOTARGETS="/home/arthur/Bureau/postdocs/isem/hybscan/hybrid_scan/busco_list"


FASTP_THREADS=4
FASTP_OPT="-q 20 -u 70 -n 40 -l 40 -w 1"
MEGAHIT_THREADS=4
MEGAHIT_OPT="--k-min 31 --k-max 101 --k-step 10"
BUSCO_THREADS=4
BUSCO_OPT="-m genome"
BWA_THREADS=4
BWA_OPT="-k 19"
FREEBAYES_OPT="--min-alternate-count 1 -z 0.05"
VCFTOOLS_FILTER_OPT="--remove-indels --minQ 30 --minDP 5"
CONTAM_FILTER_ALT_PROB="0.05 0.2" #give a whitespace-seperated list of parameter (e) values to run
STAN_THREADS=2
STAN_OPT="2 1000 10000" #1: chains number #2:warmup iterations #3: total iteration number per chain


##construct dic with probes to use for each id (UCE runs) 
#Manually 
UCEDIC={"SRR5437981":"formicidae", "SRR5406035":"hymenoptera"}
#Or from UCETARGETS, a file with two columns: id probes
UCEDIC={l.strip().split(" ")[0]:l.strip().split(" ")[1] for l in open(UCETARGETS, "r") if not re.search("^#", l)}

##construct dic with database/augustus species to use for each id (BUSCO runs)
#Manually 
BUSCODIC={"SRR1325015":["hymenoptera_odb10", "camponotus_floridanus"], "SRR4292931":["hymenoptera_odb10", "camponotus_floridanus"], "SRR4292935":["hymenoptera_odb10", "camponotus_floridanus"]}
#Or from BUSCOTARGETS, a file with three columns: id database species
BUSCODIC={l.strip().split(" ")[0]:[l.strip().split(" ")[1], l.strip().split(" ")[2]] for l in open(BUSCOTARGETS, "r") if not re.search("^#", l)}


##CONSTRAINTS
wildcard_constraints:
	ID="[^/\.]*",
	EXT="[^/]*",
	EXT2="[^/]*",
	NB="[0-9]*",

ruleorder:
	fastp_pese_run > fastp_pe_run > fastp_se_run
ruleorder:
	megahit_pese_run > megahit_pe_run > megahit_se_run
ruleorder:
	bwa_pese_run > bwa_pe_run > bwa_se_run


#######################################################
#####TARGETS###########################################
#######################################################
#Require final output. Comment out unwanted output.
rule all: 
	input: 
		######UCE######
		####hybrid scan results for all ids in UCEDIC
		##using ANGSD
		expand("{PATH}/alluce.fastp.megahit.phyluce.bwa.fastp.angsd.divestim.gathered.txt", PATH=DATADIR, ID=list(UCEDIC.keys())),
		##using snp calling and contamination filtering
		expand("{PATH}/alluce.fastp.megahit.phyluce.bwa.fastp.freebayes.filter.contamfilter{PROB}.het.divestim.gathered.txt", PATH=DATADIR, ID=list(UCEDIC.keys()), PROB=CONTAM_FILTER_ALT_PROB.split(" ")),
		####plots to vizualize the effect of filtering for contamination
		##before
		expand("{PATH}/{ID}.fastp.megahit.phyluce.bwa.fastp.freebayes.filter.refalt.png", PATH=DATADIR, ID=list(UCEDIC.keys())),
		##after
		expand("{PATH}/{ID}.fastp.megahit.phyluce.bwa.fastp.freebayes.filter.contamfilter{PROB}.refalt.png", PATH=DATADIR, ID=list(UCEDIC.keys()), PROB=CONTAM_FILTER_ALT_PROB.split(" ")),
		######BUSCO######
		####hybrid scan results for all ids in BUSCODIC
		##using ANGSD
		expand("{PATH}/allbusco.fastp.megahit.busco.bwa.fastp.angsd.divestim.gathered.txt", PATH=DATADIR, ID=list(BUSCODIC.keys())),
		##using snp calling and contamination filtering
		expand("{PATH}/allbusco.fastp.megahit.busco.bwa.fastp.freebayes.filter.contamfilter{PROB}.het.divestim.gathered.txt", PATH=DATADIR, ID=list(BUSCODIC.keys()), PROB=CONTAM_FILTER_ALT_PROB.split(" ")),
		####plots to vizualize the effect of filtering for contamination
		##before
		expand("{PATH}/{ID}.fastp.megahit.busco.bwa.fastp.freebayes.filter.refalt.png", PATH=DATADIR, ID=list(BUSCODIC.keys())),
		##after
		expand("{PATH}/{ID}.fastp.megahit.busco.bwa.fastp.freebayes.filter.contamfilter{PROB}.refalt.png", PATH=DATADIR, ID=list(BUSCODIC.keys()), PROB=CONTAM_FILTER_ALT_PROB.split(" ")),



#######################################################
#####RULES#############################################
#######################################################

#######################################################
#run fastp on single-end data
rule fastp_se_run:
	input:	
		single="{PATH}/{ID}{EXT}fastq.gz"
	output:
		single=temp("{PATH}/{ID}{EXT}fastp.fastq.gz"),
		json=temp("{PATH}/{ID}{EXT}fastp.json"),
		html=temp("{PATH}/{ID}{EXT}fastp.html")
	threads:
		int(FASTP_THREADS)
	shell:
        	"fastp -j {output.json} -h {output.html} {FASTP_OPT} -i {input.single} -o {output.single} -w {threads};"

#######################################################
#run fastp on paired-end data
rule fastp_pe_run:
	input:	
		forward="{PATH}/{ID}_1{EXT}fastq.gz",
		reverse="{PATH}/{ID}_2{EXT}fastq.gz"
	output:
		forward=temp("{PATH}/{ID}_1{EXT}fastp.fastq.gz"),
		reverse=temp("{PATH}/{ID}_2{EXT}fastp.fastq.gz"),
		single=temp("{PATH}/{ID}{EXT}fastp.fastq.gz"),
		json=temp("{PATH}/{ID}{EXT}fastp.json"),
		html=temp("{PATH}/{ID}{EXT}fastp.html")
	threads:
		int(FASTP_THREADS)
	shell:
        	"fastp -j {output.json} -h {output.html} {FASTP_OPT} --detect_adapter_for_pe -c -i {input.forward} -I {input.reverse} -o {output.forward} -O {output.reverse} --unpaired1 {output.single} --unpaired2 {output.single} -w {threads};"

#######################################################
#run fastp when there is both se and pe data
rule fastp_pese_run:
	input:	
		forward="{PATH}/{ID}_1{EXT}fastq.gz",
		reverse="{PATH}/{ID}_2{EXT}fastq.gz",
		single="{PATH}/{ID}{EXT}fastq.gz"
	output:
		forward=temp("{PATH}/{ID}_1{EXT}fastp.fastq.gz"),
		reverse=temp("{PATH}/{ID}_2{EXT}fastp.fastq.gz"),
		singlep=temp("{PATH}/{ID}{EXT}fastp.pe.fastq.gz"),
		singles=temp("{PATH}/{ID}{EXT}fastp.se.fastq.gz"),
		single=temp("{PATH}/{ID}{EXT}fastp.fastq.gz"),
		jsonp=temp("{PATH}/{ID}{EXT}fastp.pe.json"),
		jsons=temp("{PATH}/{ID}{EXT}fastp.se.json"),
		htmlp=temp("{PATH}/{ID}{EXT}fastp.pe.html"),
		htmls=temp("{PATH}/{ID}{EXT}fastp.se.html")
	threads:
		int(FASTP_THREADS)
	shell:
        	"fastp -j {output.jsonp} -h {output.htmlp} {FASTP_OPT} --detect_adapter_for_pe -c -i {input.forward} -I {input.reverse} -o {output.forward} -O {output.reverse} --unpaired1 {output.singlep} --unpaired2 {output.singlep} -w {threads};"
        	"fastp -j {output.jsons} -h {output.htmls} {FASTP_OPT} -i {input.single} -o {output.singles} -w {threads};"
		#Put se data and orphan reads from pe data in the same file
		"nonemptysingle=`ls -l {output.singles} {output.singlep} | awk '{{if ($5 > 0) print $9}}' | tr '\n' ' '`;"
		"echo $nonemptysingle | xargs zcat | pigz -p{threads} > {output.single}"
		
#######################################################
#run megahit (assembly) on single-end data
rule megahit_se_run:
	input:
		single="{PATH}/{ID}{EXT}fastq.gz"
	output: 
		cont="{PATH}/{ID}{EXT}megahit.fasta",
		log="{PATH}/{ID}{EXT}megahit.log",
		dir=temp(directory("{PATH}/{ID}{EXT}megahit")),
	threads: int(MEGAHIT_THREADS)
	shell:
        	"megahit -t {threads} -r {input.single} -o {output.dir} {MEGAHIT_OPT};"
		"mv {output.dir}/log {output.log};"
        	""" awk '!/>/{{print}}; />/{{cnt++; print ">contig"cnt}}' {output.dir}/final.contigs.fa > {output.cont}; """

#######################################################
#run megahit (assembly) on paired-end data
rule megahit_pe_run:
	input:
		forward="{PATH}/{ID}_1{EXT}fastq.gz",
		reverse="{PATH}/{ID}_2{EXT}fastq.gz",
	output: 
		cont="{PATH}/{ID}{EXT}megahit.fasta",
		log="{PATH}/{ID}{EXT}megahit.log",
		dir=temp(directory("{PATH}/{ID}{EXT}megahit")),
	threads: int(MEGAHIT_THREADS)
	shell:
        	"megahit -t {threads} -1 {input.forward} -2 {input.reverse} -o {output.dir} {MEGAHIT_OPT};"
		"mv {output.dir}/log {output.log};"
        	""" awk '!/>/{{print}}; />/{{cnt++; print ">contig"cnt}}' {output.dir}/final.contigs.fa > {output.cont}; """


#######################################################
#run megahit (assembly) when there is both se and pe data 
rule megahit_pese_run:
	input:
		forward="{PATH}/{ID}_1{EXT}fastq.gz",
		reverse="{PATH}/{ID}_2{EXT}fastq.gz",
		single="{PATH}/{ID}{EXT}fastq.gz"
	output: 
		cont="{PATH}/{ID}{EXT}megahit.fasta",
		log="{PATH}/{ID}{EXT}megahit.log",
		dir=temp(directory("{PATH}/{ID}{EXT}megahit")),
	threads: int(MEGAHIT_THREADS)
	shell:
        	"megahit -t {threads} -1 {input.forward} -2 {input.reverse} -r {input.single} -o {output.dir} {MEGAHIT_OPT};"
		"mv {output.dir}/log {output.log};"
        	""" awk '!/>/{{print}}; />/{{cnt++; print ">contig"cnt}}' {output.dir}/final.contigs.fa > {output.cont}; """

#######################################################
#run phyluce (uce identification and extraction into fasta) using specified set of probes
rule phyluce_run:
	input: 
		cont="{PATH}/{ID}{EXT}fasta",
		probes=lambda wildcards: expand("{TOOLS}/uce_probes/{PR}.fasta", TOOLS=TOOLDIR, PR=UCEDIC[wildcards.ID])
	output:
		dir=temp(directory("{PATH}/{ID}{EXT}phyluce")),
		uce="{PATH}/{ID}{EXT}phyluce.fasta",
		log="{PATH}/{ID}{EXT}phyluce.log",
	params:
		env="phyluce_env"
	shell:
		"set +eu;"
        	". $(conda info --base)/etc/profile.d/conda.sh;"
        	"conda activate {params.env};"
		"mkdir {output.dir} {output.dir}/input;"
		""" sed 's/[Y|S|W|K|M|R|y|s|w|k|m|r]/N/g' {input.cont} | awk '/>/{{sub(">contig", ">NODE_", $0); sub("$", "_length_100_cov_10", $0); print $0}}; !/>/{{print}}' > {output.dir}/input/{wildcards.ID}.contigs.fasta; """
		"phyluce_assembly_match_contigs_to_probes --contigs {output.dir}/input --probes {input.probes} --output {output.dir}/output --log-path {output.dir};"
		"echo '[all]\n{wildcards.ID}' > {output.dir}/output/{wildcards.ID}.conf;"
		"phyluce_assembly_get_match_counts --locus-db {output.dir}/output/probe.matches.sqlite --taxon-list-config {output.dir}/output/{wildcards.ID}.conf --taxon-group 'all' --incomplete-matrix --output {output.dir}/output/{wildcards.ID}.match.conf --log-path {output.dir};"
		"phyluce_assembly_get_fastas_from_match_counts --contigs {output.dir}/input --locus-db {output.dir}/output/probe.matches.sqlite --match-count-output {output.dir}/output/{wildcards.ID}.match.conf --output {output.dir}/output/{wildcards.ID}.uce.fasta --log-path {output.dir};"
		""" sed 's/_'{wildcards.ID}'.*$//g' {output.dir}/output/{wildcards.ID}.uce.fasta | awk '/^>/{{head=$0}}; !/^>/{{seq[head]=seq[head]""$0}}; END{{PROCINFO["sorted_in"]="@ind_str_asc"; for(h in seq){{print h; print seq[h]}}}}' > {output.uce}; """
		"ls -1rt {output.dir}/*log | xargs cat > {output.log};"

#######################################################
#download database for busco (should run only once per unique database)
rule busco_database_download:
	output:
		directory("{PATH}/lineages/{DB}_odb{NB}")
	shell:
		"busco --download_path {wildcards.PATH} --download {wildcards.DB}_odb{wildcards.NB};"

#######################################################
#run busco (gene identification and extraction into fasta) using specified database
rule busco_run:
	input: 
		cont="{PATH}/{ID}{EXT}fasta",
		db=lambda wildcards: expand("{DATA}/busco_database/lineages/{DB}", DATA=DATADIR, DB=BUSCODIC[wildcards.ID][0])
	output:
		dir=directory("{PATH}/{ID}{EXT}busco"),
		fasta="{PATH}/{ID}{EXT}busco.fasta",
		log="{PATH}/{ID}{EXT}busco.log",
	threads:
		int(BUSCO_THREADS)
	params:
		db=lambda wildcards: expand("{DB}", DB=BUSCODIC[wildcards.ID][0]),
		species=lambda wildcards: expand("{SP}", SP=BUSCODIC[wildcards.ID][1])
	shell:
		"busco -f --download_path {wildcards.PATH} -i {input.cont} --out_path {wildcards.PATH} -o {wildcards.ID}{wildcards.EXT}busco -l {input.db} -c {threads} --augustus --augustus_species {params.species} {BUSCO_OPT};"
		"cp {output.dir}/logs/busco.log {output.log};"
		""" awk 'FNR==1{{file=FILENAME; sub(".fna", "", file); sub("^.*\\/", "", file)}}; FNR!=1{{seq[file]=seq[file]""$0}}; END{{PROCINFO["sorted_in"]="@ind_str_asc"; for(f in seq){{print ">busco_"f; print seq[f]}}}}' {output.dir}/run_{params.db}/busco_sequences/single_copy_busco_sequences/*fna > {output.fasta}; """
		
	
#######################################################
#index fasta file with bwa
rule bwa_indexing:
	input:
        	cont="{PATH}.fasta"
	output:
		amb=temp("{PATH}.fasta.amb"),
		ann=temp("{PATH}.fasta.ann"),
		pac=temp("{PATH}.fasta.pac"),
		bwt=temp("{PATH}.fasta.bwt.2bit.64"),
		num=temp("{PATH}.fasta.0123"),
	shell:
		"bwa-mem2 index {input.cont}"	

#######################################################
#align single-end reads using bwa 
rule bwa_se_run:
	input:
		index=expand("{{PATH}}/{{ID}}{{EXT}}fasta.{EXT2}", EXT2="amb ann pac bwt.2bit.64 0123".split(" ")),
		single="{PATH}/{ID}{EXT2}fastq.gz",
		cont="{PATH}/{ID}{EXT}fasta",
	output:
		bam=temp("{PATH}/{ID}{EXT}bwa{EXT2}bam"),
	threads: 
		int(BWA_THREADS)
	shell:
		"bwa-mem2 mem {BWA_OPT} -t {threads} {input.cont} {input.single} | samtools view -bu -F 260 -@ {threads} | samtools sort -m 1G -@ {threads} > {output.bam}; """

#######################################################
#align paired-end reads using bwa 
rule bwa_pe_run:
	input:
		index=expand("{{PATH}}/{{ID}}{{EXT}}fasta.{EXT2}", EXT2="amb ann pac bwt.2bit.64 0123".split(" ")),
		forward="{PATH}/{ID}_1{EXT2}fastq.gz",
		reverse="{PATH}/{ID}_2{EXT2}fastq.gz",
		cont="{PATH}/{ID}{EXT}fasta",
	output:
		bam=temp("{PATH}/{ID}{EXT}bwa{EXT2}bam"),
	threads: 
		int(BWA_THREADS)
	shell:
		" bwa-mem2 mem {BWA_OPT} -t {threads} {input.cont} {input.forward} {input.reverse} | samtools view -bu -F 260 -@ {threads} | samtools sort -m 1G -@ {threads} > {output.bam}; """

#######################################################
#align reads using bwa when there is both se and pe data
rule bwa_pese_run:
	input:
		index=expand("{{PATH}}/{{ID}}{{EXT}}fasta.{EXT2}", EXT2="amb ann pac bwt.2bit.64 0123".split(" ")),
		forward="{PATH}/{ID}_1{EXT2}fastq.gz",
		reverse="{PATH}/{ID}_2{EXT2}fastq.gz",
		single="{PATH}/{ID}{EXT2}fastq.gz",
		cont="{PATH}/{ID}{EXT}fasta",
	output:
		bamse=temp("{PATH}/{ID}{EXT}bwa{EXT2}se.bam"),
		bampe=temp("{PATH}/{ID}{EXT}bwa{EXT2}pe.bam"),
		bam=temp("{PATH}/{ID}{EXT}bwa{EXT2}bam"),
	threads: 
		int(BWA_THREADS)
	shell:
		"bwa-mem2 mem {BWA_OPT} -t {threads} {input.cont} {input.single} | samtools view -bu -F 260 -@ {threads} | samtools sort -m 1G -@ {threads} > {output.bamse}; """
		" bwa-mem2 mem {BWA_OPT} -t {threads} {input.cont} {input.forward} {input.reverse} | samtools view -bu -F 260 -@ {threads} | samtools sort -m 1G -@ {threads} > {output.bampe}; """
		#merge pe and se alignments
		"samtools merge -@ {threads} -o {output.bam} {output.bamse} {output.bampe};"

#######################################################
#index bam file with samtools
rule samtools_bam_index:
	input:
		"{PATH}.bam"
	output:
		temp("{PATH}.bam.bai")
	shell:
		"samtools index {input};"

#######################################################
#index fasta file with samtools
rule samtools_fasta_index:
	input:
		"{PATH}.fasta"
	output:
		temp("{PATH}.fasta.fai")
	shell:
		"samtools faidx {input};"

#######################################################
#estimate number of heterozygous sites from a bam file per contig using ANGSD and bam file
rule angsd_run:
	input:
		bam="{PATH}/{ID}{EXT}bwa{EXT2}bam",
		bamindex="{PATH}/{ID}{EXT}bwa{EXT2}bam.bai",
		cont="{PATH}/{ID}{EXT}fasta",
	output:
		bamlist=temp("{PATH}/{ID}{EXT}bwa{EXT2}angsd.bamlist"),
		safidx=temp("{PATH}/{ID}{EXT}bwa{EXT2}angsd.saf.idx"),
		safpos=temp("{PATH}/{ID}{EXT}bwa{EXT2}angsd.saf.pos.gz"),
		saf=temp("{PATH}/{ID}{EXT}bwa{EXT2}angsd.saf.gz"),
		sfs=temp("{PATH}/{ID}{EXT}bwa{EXT2}angsd.sfs"),
		thetaidx=temp("{PATH}/{ID}{EXT}bwa{EXT2}angsd.thetas.idx"),
		theta=temp("{PATH}/{ID}{EXT}bwa{EXT2}angsd.thetas.gz"),
		pestpg=temp("{PATH}/{ID}{EXT}bwa{EXT2}angsd.thetas.idx.pestPG"),
		arg=temp("{PATH}/{ID}{EXT}bwa{EXT2}angsd.arg"),
		stats="{PATH}/{ID}{EXT}bwa{EXT2}angsd.stats.txt",
		counts="{PATH}/{ID}{EXT}bwa{EXT2}angsd.counts.txt",
	params:
		prefix="{PATH}/{ID}{EXT}bwa{EXT2}angsd",
		env="angsd_env"
	shell:
		"samtools faidx {input.cont};" 
		"set +eu;"
        	". $(conda info --base)/etc/profile.d/conda.sh;"
        	"conda activate {params.env};"
		"echo {input.bam} > {output.bamlist};"
		"angsd -anc {input.cont} -doSaf 1 -GL 1 -bam {output.bamlist} -out {params.prefix};"
		"realSFS {output.safidx} > {output.sfs};"
		"angsd -anc {input.cont} -doThetas 1 -doSaf 1 -pest {output.sfs} -bam {output.bamlist} -out {params.prefix} -GL 1;"
		"thetaStat do_stat {output.thetaidx};"
		"cp {output.pestpg} {output.stats};"
		""" awk 'NR!=1{{print $2" "$14+1" "$4}}' {output.stats} | sort -k 1 > {output.counts}; """

#######################################################
#call snp using freebayes 
rule freebayes_individual_run:
	input:
		bam="{PATH}/{ID}{EXT}bwa{EXT2}bam",
		cont="{PATH}/{ID}{EXT}fasta",
		contindex="{PATH}/{ID}{EXT}fasta.fai",
	output:
		vcf="{PATH}/{ID}{EXT}bwa{EXT2}freebayes.vcf",
	shell:
		"freebayes -f {input.cont} {FREEBAYES_OPT} {input.bam} > {output}"

#######################################################
#filter snp using vcftools
rule vcftools_filter_run:
	input:
		"{PATH}.vcf"
	output:
		"{PATH}.filter.vcf"
	shell:
		"vcftools --vcf {input} {VCFTOOLS_FILTER_OPT} --recode --out {wildcards.PATH}.filter;"
		"mv {wildcards.PATH}.filter.recode.vcf {output};"
 
#######################################################
#house-made filter for contamination
rule contamfilter_run:
	input:
		"{PATH}.vcf"
	output:
		"{PATH}.contamfilter{ALTPROB}.vcf"
	shell:
		"python {TOOLDIR}/scripts/vcf_contam_filter.py {input} {wildcards.ALTPROB} > {output};"

#######################################################
#count number of snp (assumed to be heterozygous sites, all filters should be done before this step) per contigs from vcf 
rule het_sites_count_run:
	input:
		"{PATH}.vcf"
	output:
		"{PATH}.het.counts.txt"
	shell:
		"python {TOOLDIR}/scripts/vcf_het_sites_count.py {input} | sort -k 1 > {output};"

#######################################################
#compile divergence model for use by stan (should run only once)
rule stan_model_train:
	input:
		"{PATH}.stan"
	output:
		"{PATH}.rds"
	shell:
		"Rscript {TOOLDIR}/stan/stan_train_model.R {input};"

#######################################################
#estimate divergence using stan
rule divergence_estimations_run:
	input:
		rds=expand("{TOOLS}/stan/divergence_model.rds", TOOLS=TOOLDIR),
		counts="{PATH}.counts.txt",
	output:
		estimates="{PATH}.divestim.txt",
		posterior="{PATH}.divposterior.txt",
	threads:
		int(STAN_THREADS)
	shell: 
		"Rscript {TOOLDIR}/scripts/stan_divergence_estimation.R {input.counts} {output.estimates} {output.posterior} {input.rds} {STAN_OPT} {threads};"

#######################################################
#gather results for all id
rule divergence_estimations_uce_gather:
	input:
		lambda wildcards: expand("{PATH}/{ID}{EXT}txt", PATH=wildcards.PATH, ID=UCEDIC.keys(), EXT=wildcards.EXT)
	output:
		"{PATH}/alluce{EXT}gathered.txt",	
	shell:
		""" awk 'NR==1{{print}}; FNR!=1{{print}}' {input} > {output}; """

#######################################################
#gather results for all id
rule divergence_estimations_busco_gather:
	input:
		lambda wildcards: expand("{PATH}/{ID}{EXT}txt", PATH=wildcards.PATH, ID=BUSCODIC.keys(), EXT=wildcards.EXT)
	output:
		"{PATH}/allbusco{EXT}gathered.txt",	
	shell:
		""" awk 'NR==1{{print}}; FNR!=1{{print}}' {input} > {output}; """


#######################################################
#produce plots of the ref and alt depths in a vcf (useful to check the effect of the contamination filter)
rule vcf_refalt_plot:
	input:
		"{PATH}.vcf"
	output:
		refalt="{PATH}.refalt.txt",
		plot="{PATH}.refalt.png",
	shell:
		"python {TOOLDIR}/scripts/vcf_sites_refalt_counts.py {input} > {output.refalt};"
		"Rscript {TOOLDIR}/scripts/snp_refalt_plot.R {output.refalt};"


