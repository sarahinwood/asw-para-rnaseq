#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=data/asw_Trinity.fasta \
include=t \
names=output/asw_timecourse/interproscan/interproscan_ids.txt \
out=output/asw_timecourse/interproscan/interproscan_degs.fasta