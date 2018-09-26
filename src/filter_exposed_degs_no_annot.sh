#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=data/sorted_fasta/asw_isoforms_by_length.fasta \
include=t \
names=output/exposed/dese2/degs_with_no_annot.txt \
out=output/exposed/no_annot/degs_no_annot.fasta