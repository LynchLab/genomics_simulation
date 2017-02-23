vg construct -f -S -a -r ref.fasta -v vars.vcf.gz > graph.vg
vg index -x graph.xg -g graph.gcsa -k 11 graph.vg
vg map -x graph.xg -g graph.gcsa -f reads.fq > mapped.gam
vg index -d mapped.gam.index -N mapped.gam
vg srpe -S vars.vcf -r ref.fa -I insertion_seqs.fa mapped.gam mapped.gam.index/ graph.vg
