#!/usr/bin/env bash

set -eu

blastx \
	-query output/exposed/no_annot/degs_no_annot.fasta \
	-db bin/db/blastdb/refseq_protein/refseq_protein.pal \
	-num_threads 50 \
	-outfmt 6 > output/exposed/no_annot/blastx.outfmt6