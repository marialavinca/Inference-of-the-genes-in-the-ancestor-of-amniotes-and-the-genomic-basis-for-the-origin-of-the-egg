# Inference of the genes in the ancestor of amniotes and the genomic basis for the origin of the egg
Amniotes are the first fully terrestrial vertebrate animals with several evolutionary innovations in their common ancestor that allowed them to become fully independent of the aquatic environment, including a more complex egg with shell and additional structures. During evolution, organisms’ form and functions evolve, and so do their genomes, which are the ultimately responsible for the observed changes. In fact, genomes are evolutionarily labile and experience changes in gene content and structure. This project aims to investigate the genomic basis for the origin of amniotes and their evolutionary innovations. From a bioinformatic point of view, this work involves: (i). choosing the highest quality vertebrate genomes to use in our analysis; (ii). estimating the genes that originated in the common ancestor of reptiles, birds and mammals by searching for sequence similarity and clustering of homologous genes; (iii). functionally characterizing the novel genes that originated in the common ancestor of amniotes and identifying any relationship with the origin of the amniote egg. 

## Keywords: Amniote, Comparative genomics, Egg, Gene gain and loss, Genomes, Evolution, Homology, Gene ontology.

DOI data: https://doi.org/10.5281/zenodo.3385935

## Structure followed of the scripts and data published:

-Selection_Of_Species

     *data_Selection_Of_Species.txt (txt file that points to Zenodo repository where the main data for this method part is contained.)
     
-Genome_Quality_Assessment_With_BUSCO    
     
     *BUSCO_scores.xlsx (Excel sheet containing the detailed scores of the genome quality assessment with BUSCO.)

-Phylogenetic_Aware_Parsing_Script

     *Species_Labels

           +Labels95.png (Table of the species of the Ensembl release 95 and the labels used after filtering the species with BUSCO.)
           
           +Labels96.png (Table of the species of the Ensembl release 96 and the labels used after filtering the species with BUSCO.)

     *Searching_Protein_Sequence_Similarities_With_DIAMOND

           +data_Searching_Protein_Sequence_Similarities_With_DIAMOND.txt (txt file that points to Zenodo repository where the main data for this method part is contained.)

     *Gene_Clustering_In_HGs_With_MCL

           +MCL_row_counter.pl (perl script that parses the output of MCL to produce a taxonomic occupancy table. The result file has a HG in each row and one species per column.)

     *Preparation_Of_PAPS_and_Inferring_Ancestral_And_Novel_Genes

           +Ancestral_Novel_HGs_searches.txt (HG searches done by clades.)
           
           +Ancestral_Novel_HGs_searches_results.png (Result graph of the HG searches done by clades.)
           
           +Ancestral_Novel_HGs_searches_results.xlsx (Excel sheet containing the detailed number of HGs by clades.)
           
           +Create_DBs.pl (perl script that speeds up the subsequent steps in the PAPS.)
           
           +Phylogenetic Aware Parsing Script.pl (perl script that allows for interactive searches for specific patterns of HG distribution across the phylogeny with the hash of hashes modified for this work.) 

-Removal_Of_False_Positives

     *search_HG_with_Fase_Pos.pl (perl script that identifies the sequence of significant hits to identify false positives in the HGs.)

-Annotation_With_GO_terms.tar.gz

     *Result_graphs (interactive map, revigo map and scatterplot result graph of the GO terms and enriched functions.)
     
           +Interactive_graph.png
           
           +revigo_treemap.pdf
           
           +scatterplot_from_R.pdf
     
     *Amniote_novel_human_9HGs_GOannot.csv (Csv sheet containing the details of the functions (GO: biological process) of the 10 human genes found in the HGs of the amniote.)
     
     *query_GO_terms_from_ensembl.sh (bash script that contains the queries to obtain GO annotations.)
     
     *topGO_fisherClassic_results.txt (Results of topGO for the Fisher's test.)
