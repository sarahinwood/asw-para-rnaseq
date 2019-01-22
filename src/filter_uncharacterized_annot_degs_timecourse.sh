#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=data/asw_transcriptome/Trinity.fasta \
include=t \
names=output/asw_timecourse/interproscan/interproscan_ids.txt \
out=output/asw_timecourse/interproscan/interproscan_degs.fasta