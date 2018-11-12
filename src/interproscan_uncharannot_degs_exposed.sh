#!/usr/bin/env bash

set -eu

bin/interproscan-5.31-70.0/interproscan.sh \
--input output/exposed/interproscan/interproscan_degs.fasta \
--seqtype n \
--output-dir output/exposed/interproscan \
--cpu 50 \
--goterms