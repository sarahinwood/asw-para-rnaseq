#!/usr/bin/env bash

set -eu

bin/interproscan-5.31-70.0/interproscan.sh \
--input output/asw_timecourse/interproscan/unchar_hypo_annot_degs.fasta \
--seqtype n \
--output-dir output/asw_timecourse/interproscan \
--cpu 50 \
--goterms