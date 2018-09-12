#!/usr/bin/env python3

import pandas

###########
# GLOBALS #
###########

read_dir = 'data/reads'

full_sample_key_file = 'data/full_sample_key.csv'

bbduk_adapters = '/adapters.fa'

#containers

salmon_container = 'shub://TomHarrop/singularity-containers:salmon_0.11.1'

#########
# SETUP #
#########

full_sample_key = pandas.read_csv(full_sample_key_file)

all_samples = sorted(set(full_sample_key['Sample_name']))

#########
# RULES #
#########

rule target:
    input:
        expand('output/salmon/asw/{sample}_quant/quant.sf',
            sample=all_samples),
        expand('output/salmon/mh/{sample}_quant/quant.sf',
            sample=all_samples)

rule asw_salmon_quant:
    input:
        index_output = 'output/salmon/asw/transcripts_index/hash.bin',
        left = 'data/bbduk_trim/{sample}_r1.fq.gz',
        right = 'data/bbduk_trim/{sample}_r2.fq.gz'
    output:
        'output/salmon/asw/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/salmon/asw/transcripts_index',
        outdir = 'output/salmon/asw/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/asw_salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.left} '
        '-2 {input.right} '
        '-o {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule asw_salmon_index:
    input:
        transcriptome_length_filtered = 'data/sorted_fasta/asw_isoforms_by_length.fasta'
    output:
        'output/salmon/asw/transcripts_index/hash.bin'
    params:
        outdir = 'output/salmon/asw/transcripts_index'
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

rule mh_salmon_quant:
    input:
        index_output = 'output/salmon/mh/transcripts_index/hash.bin',
        left = 'data/bbduk_trim/{sample}_r1.fq.gz',
        right = 'data/bbduk_trim/{sample}_r2.fq.gz'
    output:
        'output/salmon/mh/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/salmon/mh/transcripts_index',
        outdir = 'output/salmon/mh/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/mh_salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.left} '
        '-2 {input.right} '
        '-o {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule mh_salmon_index:
    input:
        transcriptome_length_filtered = 'data/sorted_fasta/mh_isoforms_by_length.fasta'
    output:
        'output/salmon/mh/transcripts_index/hash.bin'
    params:
        outdir = 'output/salmon/mh/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/mh_salmon_index.log'
    shell:
        'salmon index '
        '-t {input.transcriptome_length_filtered} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'





