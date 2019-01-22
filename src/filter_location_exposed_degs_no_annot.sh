#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=data/asw_transcriptome/Trinity.fasta \
include=t \
names=output/exposed/deseq2_ru_v_linc/degs_with_no_annot.txt \
out=output/exposed/location_no_annot/degs_no_annot.fasta