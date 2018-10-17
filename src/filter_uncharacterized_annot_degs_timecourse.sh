#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=data/asw_Trinity.fasta \
include=t \
names=output/asw_timecourse/interproscan/unchar_hypo_annot_ids.txt \
out=output/asw_timecourse/interproscan/unchar_hypo_annot_degs.fasta