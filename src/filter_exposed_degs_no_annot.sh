#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=data/asw_Trinity.fasta \
include=t \
names=output/exposed/deseq2/degs_with_no_annot.txt \
out=output/exposed/no_annot/degs_no_annot.fasta