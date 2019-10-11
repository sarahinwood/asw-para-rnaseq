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