#!/usr/bin/env bash

set -eu

blastx \
	-query output/exposed/noannot/degs_no_annot.fasta \
	-db bin/db/blastdb/refseq_protein/refseq_protein.pal \
	-num_threads 50 \
	-outfmt 6 > output/exposed/noannot/blastx.outfmt6