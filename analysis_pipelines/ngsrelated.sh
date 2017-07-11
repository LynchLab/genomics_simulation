#!/bin/bash
zcat ../sequences/angsd_calls.vcf.mafs.gz | cut -f -d '	' > ../sequences/freq
ngsrelate -g ../sequences/angsd_calls.vcf.glf.gz -n 6 -f ../sequences/freq -s 1
