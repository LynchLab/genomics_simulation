This the pipeline we use to test MAPGD.

* Reference simulation: By in large genomes do not have a random structure, but instead transpositions, duplications, selective constraint etc. create substantial sequence similarity between different regions of  the genome. This may create systematic errors when attempting to align reads to a reference genome/population, and we would like to capture these potential errors in our simulation. Ideally the de novo assembly of the genome should be simulated, so that the effects of the artifacts created by the assembly process can be correctly modeled, but we do not currently model this process.
	NEEDS: NOTHING
	MAKES: REFERENCE FILE
	STATUS: DONE

* Population simulation: Individuals within a population are related to each other by a pedigree structure. For instance, in finite populations we cannot model the genotypes of two individuals as being independent samples. Nor can we assume that the coalescent trees between freely recombining sites are independent of each other, because the underlying pedigree which relates individuals influences the coalescent tree. We explicitly model this pedigree, but currently can only model freely recombining loci.

	NEEDS: NOTHING
	MAKES: VARIANT PRESENCE/ABSENSE, PHENOTYPIC SCORES
	TODO: ADD LINKAGE

* Variant simulation: Variants themselves are not random, but occur at different rates depending on the type of variant and the local sequences context. 

	NEEDS: REFERENCE, VARIANT PRESENCE/ABSENCE
	MAKES: INDIVIDUAL FASTAS
	TODO: ADD CONTEXT MODIFIRES

* Sequencing simulation: Sequencers do not just make random  errors when sequencing DNA. Luckily many programs have been developed to simulate sequencing, and we can just make use of the programs that already exist.

	NEEDS: INDIVIDUAL FASTAS, REFERENCE
	MAKES: BAMFILES
	STATUS: DONE

* Analysis pipelines:
	Labeled frequency calling:	
		NEEDS: BAMFILES, REFFERENCE
		STATUS: MAPGD, GATK, BCFTOOLS, ANGSD
	Relatedness:
		NEEDS: BAMFILES, REFERENCE, PHENOTYPIC SCORES
		STATUS: MAPGD, GCTA, PLINK,
	Pooled frequency calling:
		NEEDS: BAMFILE, REFERENCE
		STATUS: MAPGD, BRESEQ	
