#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=data/Trinity.fasta \
include=t \
names=output/asw_timecourse/deseq2/degs_with_no_annot.txt \
out=output/asw_timecourse/no_annot/degs_no_annot.fasta