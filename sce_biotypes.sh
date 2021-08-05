#!/bin/bash

wget -N ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.36.gtf.gz

zcat Saccharomyces_cerevisiae.R64-1-1.36.gtf.gz \
| grep -w gene | cut -f9 | awk '{print $2,$(NF)}' \
| tr -d '";,' | sort -u \
| tr ' ' '\t' > Saccharomyces_cerevisiae.R64-1-1.36.gtf.biotypes.tsv
