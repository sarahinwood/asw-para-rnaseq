#!/usr/bin/env bash

set -eu

blastx \
	-query output/asw_timecourse/no_annot/degs_no_annot.fasta \
	-db bin/db/blastdb/nr/nr \
	-num_threads 50 \
	-outfmt 6 > output/exposed/no_annot/blastx.outfmt6