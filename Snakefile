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

rule target:
    input:
        expand('output/salmon/{sample}_quant/quant.sf',
            sample=all_samples)

rule salmon_quant:
    input:
        index_output = 'output/salmon/transcripts_index/hash.bin',
        left = 'data/bbduk_trim/{sample}_r1.fq.gz',
        right = 'data/bbduk_trim/{sample}_r2.fq.gz'
    output:
        'output/salmon/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/salmon/transcripts_index',
        outdir = 'output/salmon/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.left} '
        '-2 {input.right} '
        '-o {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule salmon_index:
    input:
        transcriptome_length_filtered = 'data/isoforms_by_length.fasta'
    output:
        'output/salmon/transcripts_index/hash.bin'
    params:
        outdir = 'output/salmon/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/salmon_index.log'
    shell:
        'salmon index '
        '-t {input.transcriptome_length_filtered} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'









