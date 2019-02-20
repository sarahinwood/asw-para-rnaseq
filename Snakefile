#!/usr/bin/env python3
import pathlib2
import os
import pandas

#############
# FUNCTIONS #
#############

def resolve_path(x):
	return(str(pathlib2.Path(x).resolve(strict=False)))

def find_read_files(read_dir):
#Make list of files
	path_generator = os.walk(read_dir, followlinks = True)
	my_files = list((dirpath, filenames)
		for (dirpath, dirname, filenames)
		in path_generator)
#Make new dictionary & populate with files (flowcell = key)
	my_fastq_files = {}
	for dirpath, filenames in my_files:
		for filename in filenames:
			if filename.endswith('.fastq.gz'):
				my_flowcell = pathlib2.Path(dirpath).name
				my_fastq = str(pathlib2.Path(dirpath,filename))
				if my_flowcell in my_fastq_files:
					my_fastq_files[my_flowcell].append(my_fastq)
				else:
					my_fastq_files[my_flowcell]= []
					my_fastq_files[my_flowcell].append(my_fastq)
	return(my_fastq_files)

def sample_name_to_fastq(wildcards):
	sample_row = sample_key[sample_key['Sample_name'] == wildcards.sample]
	sample_id = sample_row.iloc[-1]['OGF_sample_ID']
	sample_flowcell = sample_row.iloc[-1]['Flow_cell']
	sample_all_fastq = [x for x in all_fastq[sample_flowcell]
						if '-{}-'.format(sample_id) in x]
	sample_r1 = sorted(list(x for x in sample_all_fastq
							if '_R1_' in os.path.basename(x)))
	sample_r2 = sorted(list(x for x in sample_all_fastq
							if '_R2_' in os.path.basename(x)))
	return({'r1': sample_r1, 'r2': sample_r2})

###########
# GLOBALS #
###########

read_dir = 'data/reads'

sample_key_file = 'data/sample_key.csv'

bbduk_adapters = '/adapters.fa'

star_reference_folder = 'output/star/star_reference'

#containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
salmon_container = 'shub://TomHarrop/singularity-containers:salmon_0.11.1'
star_container = 'shub://TomHarrop/singularity-containers:star_2.7.0c'

#########
# SETUP #
#########
# generate name to filename dictionary
all_fastq = find_read_files(read_dir)

sample_key = pandas.read_csv(sample_key_file)

all_samples = sorted(set(sample_key['Sample_name']))

#########
# RULES #
#########

rule target:
	input:
	 expand('output/asw_salmon/{sample}_quant/quant.sf', sample = all_samples)

rule asw_salmon_quant:
	input:
		index_output = 'output/asw_salmon/transcripts_index/hash.bin',
		left = 'output/renamed/{sample}_Unmapped_mate1.fq.gz',
		right = 'output/renamed/{sample}_Unmapped_mate2.fq.gz'
	output:
		quant = 'output/asw_salmon/{sample}_quant/quant.sf',
		eq = 'output/asw_salmon/{sample}_quant/aux_info/eq_classes.txt',
	params:
		index_outdir = 'output/asw_salmon/transcripts_index',
		outdir = 'output/asw_salmon/{sample}_quant'
	threads:
		20
	singularity:
		salmon_container
	log:
		'output/logs/salmon/asw_salmon_quant_{sample}.log'
	shell:
		'salmon quant '
		'-i {params.index_outdir} '
		'-l ISR '
		'--dumpEq '
		'-1 {input.left} '
		'-2 {input.right} '
		'-o {params.outdir} '
		'--writeUnmappedNames '
		'-p {threads} '
		'&> {log}'

rule asw_salmon_index:
	input:
		transcriptome_length_filtered = 'data/asw_transcriptome/isoforms_by_length.fasta'
	output:
		'output/asw_salmon/transcripts_index/hash.bin'
	params:
		outdir = 'output/asw_salmon/transcripts_index'
	threads:
		20
	singularity:
		salmon_container
	log:
		'output/logs/asw_salmon_index.log'
	shell:
		'salmon index '
		'-t {input.transcriptome_length_filtered} '
		'-i {params.outdir} '
		'-p {threads} '
		'&> {log}'

rule rename_unmapped:
	input:
		unmapped_left = 'output/star/star_pass2/{sample}.Unmapped.out.mate1',
		unmapped_right = 'output/star/star_pass2/{sample}.Unmapped.out.mate2'
	output:
		fq_left = temp('output/renamed/{sample}_Unmapped_mate1.fq.gz'),
		fq_right = temp('output/renamed/{sample}_Unmapped_mate2.fq.gz')
	shell:
		'cp {input.unmapped_left} {output.fq_left} & '
		'cp {input.unmapped_right} {output.fq_right} & '
		'wait '

rule star_second_pass:
	input:
		left = 'output/bbduk_trim/{sample}_r1.fq.gz',
		right = 'output/bbduk_trim/{sample}_r2.fq.gz',
		junctions = expand('output/star/star_pass1/{sample}.SJ.out.tab', sample=all_samples)
	output:
		bam = ('output/star/star_pass2/{sample}.Aligned.sortedByCoord.out.bam'),
		unmapped_left = 'output/star/star_pass2/{sample}.Unmapped.out.mate1',
		unmapped_right = 'output/star/star_pass2/{sample}.Unmapped.out.mate2'
	threads:
		30
	params:
		genome_dir = star_reference_folder,
		prefix = 'output/star/star_pass2/{sample}.'
	log:
		'output/logs/star/star_pass2_{sample}.log'
	singularity:
		star_container
	shell:
		'STAR '
		'--runThreadN {threads} '
		'--genomeDir {params.genome_dir} '
		'--sjdbFileChrStartEnd {input.junctions} '
		'--outSAMtype BAM SortedByCoordinate '
		'--outBAMcompression 10 '
		'--outReadsUnmapped Fastx '
		'--readFilesIn {input.left} {input.right} '
		'--readFilesCommand zcat '
		'--outFileNamePrefix {params.prefix} '
		'&> {log}'

		##sjdbFileChrStartEnd option? /path/to/sjdbFile.txt --> need all samples at once?

rule star_first_pass:
	input:
		left = 'output/bbduk_trim/{sample}_r1.fq.gz',
		right = 'output/bbduk_trim/{sample}_r2.fq.gz',
		star_reference = 'output/star/star_reference/Genome'
	output:
		sjdb = 'output/star/star_pass1/{sample}.SJ.out.tab'
	params:
		genome_dir = star_reference_folder,
		prefix = 'output/star/star_pass1/{sample}.'
	threads:
		30
	log:
		'output/logs/star/star_pass1_{sample}.log'
	singularity:
		star_container
	shell:
		'STAR '
		'--runThreadN {threads} '
		'--genomeDir {params.genome_dir} '
		'--outSJfilterReads Unique '
		'--outSAMtype None '
		'--readFilesIn {input.left} {input.right} '
		'--readFilesCommand zcat '
		'--outFileNamePrefix {params.prefix} '
		'&> {log}'

rule star_reference:
	input:
		fasta = 'data/Mh_assembly.fa'
	output:
		'output/star/star_reference/Genome'
	params:
		genome_dir = star_reference_folder
	threads:
		20
	log:
		'output/logs/star/star_reference.log'
	singularity:
		star_container
	shell:
		'STAR '
		'--runThreadN {threads} '
		'--runMode genomeGenerate '
		'--genomeDir {params.genome_dir} '
		'--genomeFastaFiles {input.fasta} '
		'&> {log} '

rule bbduk_trim:
	input:
		r1 = 'output/joined/{sample}_r1.fq.gz',
		r2 = 'output/joined/{sample}_r2.fq.gz'
	output:
		r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
		r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
	params:
		adapters = bbduk_adapters
	log:
		'output/logs/bbduk_trim/{sample}.log'
	threads:
		20
	singularity:
		bbduk_container
	shell:
		'bbduk.sh '
		'in={input.r1} '
		'in2={input.r2} '
		'out={output.r1} '
		'out2={output.r2} '
		'ref={params.adapters} '
		'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
		'&> {log}'

rule cat_reads:
	input:
		unpack(sample_name_to_fastq)
	output: 
		r1 = temp('output/joined/{sample}_r1.fq.gz'),
		r2 = temp('output/joined/{sample}_r2.fq.gz')
	threads:
		1
	shell:
		'cat {input.r1} > {output.r1} & '
		'cat {input.r2} > {output.r2} & '
		'wait'



