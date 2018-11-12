#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=data/Trinity.fasta \
include=t \
names=output/exposed/interproscan/interproscan_ids.txt \
out=output/exposed/interproscan/interproscan_degs.fasta