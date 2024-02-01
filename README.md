# Venu_SticklebackRecombination

Fine-scale contemporary recombination variation and its fitness consequences in adaptively diverging stickleback fish

Vrinda Venu1,*, Enni Harjunmaa1, Andreea Dreau1, Shannon Brady2, Devin Absher3, David M. Kingsley2,4, Felicity C. Jones1,5,†

Affiliations
1.	Friedrich Miescher Laboratory of the Max Planck Society, Max Planck Ring 9, 72076 Tübingen, Germany
2.	Dept. of Developmental Biology, Stanford University School of Medicine, Stanford, CA 94305-5329, USA
3.	HudsonAlpha Institute for Biotechnology, 601 Genome Way, Huntsville, AL 35806, USA
4.	Howard Hughes Medical Institute, Chevy Chase, MD 20815, USA
5.	Groningen Institute for Evolutionary Life Sciences (GELIFES), University of Groningen, Nijenborgh 7, Groningen, The Netherlands

\* vrinda@lanl.gov

† corresponding author: fcjones@tuebingen.mpg.de

Abstract

Despite deep evolutionary conservation, recombination varies greatly across the genome, among individuals, sexes and populations and can be a major evolutionary force in the wild. Yet this variation in recombination and its impact on adaptively diverging populations is not well understood. To elucidate the nature and potential consequences of recombination rate variation, we characterized fine-scale recombination landscapes by combining pedigrees, functional genomics and field fitness measurements in an adaptively divergent pair of marine and freshwater threespine stickleback populations from River Tyne, Scotland. Through whole-genome sequencing of large nuclear families, we identified the genomic location of almost 50,000 crossovers and built recombination maps for 36 marine, freshwater, and hybrid individuals at 3.8 kilobase resolution. Using these maps, we quantified the factors driving variation in recombination rate: we find strong heterochiasmy between sexes (68% of the variation) but also differences among ecotypes (21.8%). Hybrids show evidence of significant recombination suppression, both in overall map length and in individual loci. We further tested and found reduced recombination rates both within single marine–freshwater adaptive loci and between loci on the same chromosome, suggestive of selection on linked ‘cassettes’. We tested theory supporting the evolution of linked selection using temporal sampling along a natural hybrid zone, and found that recombinants with shuffled alleles across loci show traits associated with reduced fitness. Our results support predictions that divergence in cis-acting recombination modifiers whose mechanisms are disrupted in hybrids, may have an important role to play in the maintenance of differences among adaptively diverging populations.

____

*Code*

This Github provides the code for 

1) crossover detection from whole genome sequending of nuclear families

2) population genetic simulations exploring the effects of recombination heterochiasmy and recombination suppression in cis(coupled) heterozygotes in the context of populations undergoing a adaptive divergence with gene flow with a selection-migration balance.  




**1) Crossover Detection**

This script detect crossovers for each chromosome for a given parent (maternal or paternal) and outputs a bed file with chromosome, individual_id, and CO boundaries. A .txt file for each chromosome with SNP filtering information is also produced.

The script requires scaffold boundary coordinates, duohmm corrected haplotype file, genotyping error file, and duohmm corrected sample file as input. 

*To invoke:*

	Rscript crossover_detection.R \
	        ${chromosome}_scaffold_boundary_file.txt \
	        ${chromosome}.maternal.shapeit.duohmm-corrected.haps \
	        ${chromosome}.maternal.shapeit.duohmm.generr \
	        ${chromosome}.maternal.shapeit.duohmm-corrected.sample \
	        mum \
	        ${chromosome} \
	        <output_prefix>;

____


**2) Population Genetic simulations of a selection migration balance (hybrid zone) with heterochiasmy and recombination suppression in cis(coupled) heterozygotes.**

We investigated how heterochiasmy and recombination suppression influence the spread of freshwater adaptive alleles to neighbouring freshwater populations via hybrid zones. Specifically we used a five population model to simulate a selection-migration balance maintaining divergence among a freshwater (population 1), intertidal (population 2), and marine population (population 3). A second hybrid zone with intertidal (population 4) and freshwater population (population 5) was included with connection to the first via marine population (population 3), resulting in parallel hybrid zones with access to the same marine population.


		Population1 <-> Population2 <-> Population3 <-> Population4 <-> Population5
		(freshwater)    (intertidal)      (marine)      (intertidal)    (freshwater)
	Ne:        200             200              600            200             200


	chromosome homolog1:          Locus1 ----- Locus2 ----- Locus3
	 				       X            X  sex-specific recombination rate 
	chromosome homolog2:          Locus1 ----- Locus2 ----- Locus3	  (rr.female, rr.male)


Using diploid forward simulations to track alleles at three linked autosomal loci of equal recombination distance apart, we simulated for each of 120 generations:
	1. the migration of individuals into and out of each population (m=0.2*Ne individuals per generation per year with a female-biased migration ratio of 3females:1male mimicking ratios found in migrating marine fish);
	2. for each individual, the production of haploid gametes via meiosis with sex-specific recombination rate among neighboring loci specified by probabilities (rr.male and rr.female)
	3. the pairing of male and female gametes to produce 100 diploid offspring per breeding pair (clutch size = 100).
	4. strong selection on offspring with survival dependent on the additive effects of the multilocus haplotypes carried by each individual and the habitat in which they are found (e.g. in freshwater and marine habitats selection favors the additive number of freshwater and marine alleles respectively; in intertidal populations selection favors individuals with cis-linked haplotypes of marine and/or freshwater alleles on both chromosome homologs; see Table S9 for more information)
	5. finally to maintain equal population sizes, we imposed a population carrying capacity of 200 by randomly selecting 200 of the surviving offspring from all breeding pairs in the current generation (g) to continue as breeders in the subsequent generation (g+1).

We started each simulation with a two-locus secondary contact scenario where Population 1 carried freshwater adaptive alleles at Locus 1 and 3, while all other populations carried marine alleles. (Locus 2 is initially invariant among all populations, with no effect on offspring survival). We allowed migration, breeding (with recombination) and selection to occur for 30 generations during which a selection-migration balance is quickly established and maintains divergent allele frequencies in Population 1 and 3.  Even with the maintenance of this divergence, the freshwater adaptive allele spreads from population 1 through the hybrid zone into population 3 and subsequently invades the parallel hybrid zone (populations 4 and 5) with a second selection-migration balance maintaining divergence in allele frequencies between populations 3 and 5. Then, after 30 generations, a new freshwater-beneficial mutation was introduced at Locus 2, by the arrival into freshwater Population 1 of N=0.01*ne migrants homozygous for freshwater alleles at all three loci. With locus 2 no longer neutral, from this point on the survival probabilities of offspring now depended on their three locus haplotypes. 

We ran 8 different simulations that differed in the extent of heterochiasmy and recombination suppression in cis-(coupled)-heterozygotes (see Table S8), and for each simulation tracked allele frequencies at each locus in each population.

We were specifically interested in quantifying differences in how quickly the new freshwater adaptive mutation at Locus 2 was able to establish and increase to high frequency in Population 1, the speed at which this new allele was able to spread and establish in Population 5, and the maintenance and magnitude of allele frequency differences between populations at either end of the hybrid zones (pop1 vs pop3; and pop5 vs pop3) .

*To invoke:* 

	Recombination_SelMigBalance.Rscript

____

Further details and description for both analyses are outlined in the supplementary methods of the paper.
