#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(BiocManager)
options(repos = BiocManager::repositories())

library(shiny)
library(shinythemes)
library(WGCNA)
library(ggplot2)
library(DT)
#library(easyalluvial)
#library(parcats)
library(readr)
library(readxl)
library(RColorBrewer)
library(heatmaply)
library(plotly)
library(shinycssloaders)
library(vembedr)

options(stringsAsFactors = FALSE);
##############EColi_EPEC#######
#Reading the data
Escherichia_coli_EPEC_0127_H6_E2348_69 = read.csv("Escherichia_coli_EPEC_0127_H6_E2348_69_R_input.csv", sep=";")
dim(Escherichia_coli_EPEC_0127_H6_E2348_69) #analising input

#getting TPM values
Escherichia_coli_EPEC_0127_H6_E2348_69_expr =as.data.frame(Escherichia_coli_EPEC_0127_H6_E2348_69[,2:ncol(Escherichia_coli_EPEC_0127_H6_E2348_69)])
dim(Escherichia_coli_EPEC_0127_H6_E2348_69_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Escherichia_coli_EPEC_0127_H6_E2348_69_expr) = Escherichia_coli_EPEC_0127_H6_E2348_69[,1]



#########################ETEC#############################
#Reading the data
Escherichia_coli_ETEC_H10407 = read.csv("Escherichia_coli_ETEC_H10407_R_input.csv", sep=";")
dim(Escherichia_coli_ETEC_H10407) #analising input

#getting TPM values
Escherichia_coli_ETEC_H10407_expr =as.data.frame(Escherichia_coli_ETEC_H10407[,2:ncol(Escherichia_coli_ETEC_H10407)])
dim(Escherichia_coli_ETEC_H10407_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Escherichia_coli_ETEC_H10407_expr) = Escherichia_coli_ETEC_H10407[,1]



#####################UPEC################################
#Reading the data
Escherichia_coli_UPEC_536 = read.csv("Escherichia_coli_UPEC_536_R_input.csv", sep=";")
dim(Escherichia_coli_UPEC_536) #analising input

#getting TPM values
Escherichia_coli_UPEC_536_expr =as.data.frame(Escherichia_coli_UPEC_536[,2:ncol(Escherichia_coli_UPEC_536)])
dim(Escherichia_coli_UPEC_536_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Escherichia_coli_UPEC_536_expr) = Escherichia_coli_UPEC_536[,1]

#filtering out missing values




##############################################################


#######################HP_G27####################

#Reading the data
Helicobacter_pylori_G27 = read.csv("Helicobacter_pylori_G27_R_input2.csv", sep=";")
dim(Helicobacter_pylori_G27) #analising input

#getting TPM values
Helicobacter_pylori_G27_expr =as.data.frame(Helicobacter_pylori_G27[,2:ncol(Helicobacter_pylori_G27)])
dim(Helicobacter_pylori_G27_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Helicobacter_pylori_G27_expr) = Helicobacter_pylori_G27[,1]

#########################Salmonella#############################
Salmonella = read.csv("Salmonella_enterica_R_input.csv", sep=";")
dim(Salmonella) #analising input

#getting TPM values
Salmonella_expr =as.data.frame(Salmonella[,2:ncol(Salmonella)])
dim(Salmonella_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Salmonella_expr) = Salmonella[,1]

############################################


#######################HP_J99####################

#Reading the data
Helicobacter_pylori_J99 = read.csv("Helicobacter_pylori_J99_R_input2.csv", sep=";")
dim(Helicobacter_pylori_J99) #analising input

#getting TPM values
Helicobacter_pylori_J99_expr =as.data.frame(Helicobacter_pylori_J99[,2:ncol(Helicobacter_pylori_J99)])
dim(Helicobacter_pylori_J99_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Helicobacter_pylori_J99_expr) = Helicobacter_pylori_J99[,1]

#############################################################
#######################Acinetobacter_baumannii####################

#Reading the data
Acinetobacter_baumannii = read.csv("Acinetobacter_baumannii_R_input.csv", sep=";")
dim(Acinetobacter_baumannii) #analising input

#getting TPM values
Acinetobacter_baumannii_expr =as.data.frame(Acinetobacter_baumannii[,2:ncol(Acinetobacter_baumannii)])
dim(Acinetobacter_baumannii_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Acinetobacter_baumannii_expr) = Acinetobacter_baumannii[,1]

#############################################################
#######################Achromobacter_xylosoxidans####################

#Reading the data
Achromobacter_xylosoxidans = read.csv("Achromobacter_xylosoxidans_SOLR_R_input.csv", sep=";")
dim(Achromobacter_xylosoxidans) #analising input

#getting TPM values
Achromobacter_xylosoxidans_expr =as.data.frame(Achromobacter_xylosoxidans[,2:ncol(Achromobacter_xylosoxidans)])
dim(Achromobacter_xylosoxidans_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Achromobacter_xylosoxidans_expr) = Achromobacter_xylosoxidans[,1]

#############################################################
#######################Aggregatibacter_actinomycetemcomitans####################

#Reading the data
Aggregatibacter_actinomycetemcomitans= read.csv("Aggregatibacter_actinomycetemcomitans_R_input.csv", sep=";")
dim(Aggregatibacter_actinomycetemcomitans) #analising input

#getting TPM values
Aggregatibacter_actinomycetemcomitans_expr =as.data.frame(Aggregatibacter_actinomycetemcomitans[,2:ncol(Aggregatibacter_actinomycetemcomitans)])
dim(Aggregatibacter_actinomycetemcomitans_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Aggregatibacter_actinomycetemcomitans_expr) = Aggregatibacter_actinomycetemcomitans[,1]

#######################Borrelia_burgdorferi####################

#Reading the data
Borrelia_burgdorferi= read.csv("Borrelia_burgdorferi_B31_R_input2.csv", sep=";")
dim(Borrelia_burgdorferi) #analising input

#getting TPM values
Borrelia_burgdorferi_expr =as.data.frame(Borrelia_burgdorferi[,2:ncol(Borrelia_burgdorferi)])
dim(Borrelia_burgdorferi_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Borrelia_burgdorferi_expr) = Borrelia_burgdorferi[,1]


#######################Burkholderia_pseudomallei####################

#Reading the data
Burkholderia_pseudomallei = read.csv("Burkholderia_pseudomallei_R_input2.csv", sep=";")
dim(Burkholderia_pseudomallei) #analising input

#getting TPM values
Burkholderia_pseudomallei_expr =as.data.frame(Burkholderia_pseudomallei[,2:ncol(Burkholderia_pseudomallei)])
dim(Burkholderia_pseudomallei_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Burkholderia_pseudomallei_expr) = Burkholderia_pseudomallei[,1]


#######################Campylobacter_jejuni####################

#Reading the data
Campylobacter_jejuni = read.csv("Campylobacter_jejuni_R_input.csv", sep=";")
dim(Campylobacter_jejuni) #analising input

#getting TPM values
Campylobacter_jejuni_expr =as.data.frame(Campylobacter_jejuni[,2:ncol(Campylobacter_jejuni)])
dim(Campylobacter_jejuni_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Campylobacter_jejuni_expr) = Campylobacter_jejuni[,1]


###########################Francisella_tularensis################

#Reading the data
Francisella_tularensis = read.csv("Francisella_tularensis_R_input.csv", sep=";")
dim(Francisella_tularensis) #analising input

#getting TPM values
Francisella_tularensis_expr =as.data.frame(Francisella_tularensis[,2:ncol(Francisella_tularensis)])
dim(Francisella_tularensis_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Francisella_tularensis_expr) = Francisella_tularensis[,1]



###########################Haemophilus_influenzae################

#Reading the data
Haemophilus_influenzae = read.csv("Haemophilus_influenzae_R_input2.csv", sep=";")
dim(Haemophilus_influenzae) #analising input

#getting TPM values
Haemophilus_influenzae_expr =as.data.frame(Haemophilus_influenzae[,2:ncol(Haemophilus_influenzae)])
dim(Haemophilus_influenzae_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Haemophilus_influenzae_expr) = Haemophilus_influenzae[,1]


###########################Klebsiella_pneumoniae################

#Reading the data
Klebsiella_pneumoniae= read.csv("Klebsiella_pneumoniae_R_input.csv", sep=";")
dim(Klebsiella_pneumoniae) #analising input

#getting TPM values
Klebsiella_pneumoniae_expr =as.data.frame(Klebsiella_pneumoniae[,2:ncol(Klebsiella_pneumoniae)])
dim(Klebsiella_pneumoniae_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Klebsiella_pneumoniae_expr) = Klebsiella_pneumoniae[,1]

###########################Legionella_pneumophila################

#Reading the data
Legionella_pneumophila= read.csv("Legionella_pneumophila_R_input2.csv", sep=";")
dim(Legionella_pneumophila) #analising input

#getting TPM values
Legionella_pneumophila_expr =as.data.frame(Legionella_pneumophila[,2:ncol(Legionella_pneumophila)])
dim(Legionella_pneumophila_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Legionella_pneumophila_expr) = Legionella_pneumophila[,1]

###########################Listeria_monocytogenes################

#Reading the data
Listeria_monocytogenes= read.csv("Listeria_monocytogenes_R_input.csv", sep=";")
dim(Listeria_monocytogenes) #analising input

#getting TPM values
Listeria_monocytogenes_expr =as.data.frame(Listeria_monocytogenes[,2:ncol(Listeria_monocytogenes)])
dim(Listeria_monocytogenes_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Listeria_monocytogenes_expr) = Listeria_monocytogenes[,1]


###########################Mycobacterium_tuberculosis################

#Reading the data
Mycobacterium_tuberculosis= read.csv("Mycobacterium_tuberculosis_R_input.csv", sep=";")
dim(Mycobacterium_tuberculosis) #analising input

#getting TPM values
Mycobacterium_tuberculosis_expr =as.data.frame(Mycobacterium_tuberculosis[,2:ncol(Mycobacterium_tuberculosis)])
dim(Mycobacterium_tuberculosis_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Mycobacterium_tuberculosis_expr) = Mycobacterium_tuberculosis[,1]


###########################Neisseria_gonorrhoeae################

#Reading the data
Neisseria_gonorrhoeae= read.csv("Neisseria_gonorrhoeae_R_input.csv", sep=";")
dim(Neisseria_gonorrhoeae) #analising input

#getting TPM values
Neisseria_gonorrhoeae_expr =as.data.frame(Neisseria_gonorrhoeae[,2:ncol(Neisseria_gonorrhoeae)])
dim(Neisseria_gonorrhoeae_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Neisseria_gonorrhoeae_expr) = Neisseria_gonorrhoeae[,1]

###########################Neisseria_meningitidis################

#Reading the data
Neisseria_meningitidis= read.csv("Neisseria_meningitidis_R_input.csv", sep=";")
dim(Neisseria_meningitidis) #analising input

#getting TPM values
Neisseria_meningitidis_expr =as.data.frame(Neisseria_meningitidis[,2:ncol(Neisseria_meningitidis)])
dim(Neisseria_meningitidis_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Neisseria_meningitidis_expr) = Neisseria_meningitidis[,1]

###########################Pseudomonas_aeruginosa################

#Reading the data
Pseudomonas_aeruginosa= read.csv("Pseudomonas_aeruginosa_R_input.csv", sep=";")
dim(Pseudomonas_aeruginosa) #analising input

#getting TPM values
Pseudomonas_aeruginosa_expr =as.data.frame(Pseudomonas_aeruginosa[,2:ncol(Pseudomonas_aeruginosa)])
dim(Pseudomonas_aeruginosa_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Pseudomonas_aeruginosa_expr) = Pseudomonas_aeruginosa[,1]


###########################Staphylococcus_aureus_MRSA252################

#Reading the data
Staphylococcus_aureus_MRSA252= read.csv("Staphylococcus_aureus_MRSA252_R_input.csv", sep=";")
dim(Staphylococcus_aureus_MRSA252) #analising input

#getting TPM values
Staphylococcus_aureus_MRSA252_expr =as.data.frame(Staphylococcus_aureus_MRSA252[,2:ncol(Staphylococcus_aureus_MRSA252)])
dim(Staphylococcus_aureus_MRSA252_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Staphylococcus_aureus_MRSA252_expr) = Staphylococcus_aureus_MRSA252[,1]


###########################Staphylococcus_aureus_MSSA476################

#Reading the data
Staphylococcus_aureus_MSSA476= read.csv("Staphylococcus_aureus_MSSA476_R_input.csv", sep=";")
dim(Staphylococcus_aureus_MSSA476) #analising input

#getting TPM values
Staphylococcus_aureus_MSSA476_expr =as.data.frame(Staphylococcus_aureus_MSSA476[,2:ncol(Staphylococcus_aureus_MSSA476)])
dim(Staphylococcus_aureus_MSSA476_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Staphylococcus_aureus_MSSA476_expr) = Staphylococcus_aureus_MSSA476[,1]


###########################Staphylococcus_epidermidis################

#Reading the data
Staphylococcus_epidermidis= read.csv("Staphylococcus_epidermidis_R_input.csv", sep=";")
dim(Staphylococcus_epidermidis) #analising input

#getting TPM values
Staphylococcus_epidermidis_expr =as.data.frame(Staphylococcus_epidermidis[,2:ncol(Staphylococcus_epidermidis)])
dim(Staphylococcus_epidermidis_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Staphylococcus_epidermidis_expr) = Staphylococcus_epidermidis[,1]

###########################Streptococcus_agalactiae################

#Reading the data
Streptococcus_agalactiae= read.csv("Streptococcus_agalactiae_R_input.csv", sep=";")
dim(Streptococcus_agalactiae) #analising input

#getting TPM values
Streptococcus_agalactiae_expr =as.data.frame(Streptococcus_agalactiae[,2:ncol(Streptococcus_agalactiae)])
dim(Streptococcus_agalactiae_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Streptococcus_agalactiae_expr) = Streptococcus_agalactiae[,1]


###########################Streptococcus_pneumoniae################

#Reading the data
Streptococcus_pneumoniae= read.csv("Streptococcus_pneumoniae_R_input.csv", sep=";")
dim(Streptococcus_pneumoniae) #analising input

#getting TPM values
Streptococcus_pneumoniae_expr =as.data.frame(Streptococcus_pneumoniae[,2:ncol(Streptococcus_pneumoniae)])
dim(Streptococcus_pneumoniae_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Streptococcus_pneumoniae_expr) = Streptococcus_pneumoniae[,1]

###########################Streptococcus_pyogenes################

#Reading the data
Streptococcus_pyogenes= read.csv("Streptococcus_pyogenes_R_input.csv", sep=";")
dim(Streptococcus_pyogenes) #analising input

#getting TPM values
Streptococcus_pyogenes_expr =as.data.frame(Streptococcus_pyogenes[,2:ncol(Streptococcus_pyogenes)])
dim(Streptococcus_pyogenes_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Streptococcus_pyogenes_expr) = Streptococcus_pyogenes[,1]


###########################Streptococcus_suis################

#Reading the data
Streptococcus_suis= read.csv("Streptococcus_suis_R_input.csv", sep=";")
dim(Streptococcus_suis) #analising input

#getting TPM values
Streptococcus_suis_expr =as.data.frame(Streptococcus_suis[,2:ncol(Streptococcus_suis)])
dim(Streptococcus_suis_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Streptococcus_suis_expr) = Streptococcus_suis[,1]

#############################Yersinia###############################

#Reading the data
Yersinia_pseudot_all = read.csv("Early_readcounts_with_PATHOgenex_T3SS_invivo_duplicates.csv")
Yersinia_pseudot_TPM = Yersinia_pseudot_all[,c(1,4:36,49:51)]
dim(Yersinia_pseudot_TPM)
Yersinia_pseudot_TPM_expr=as.data.frame(t(Yersinia_pseudot_TPM[,2:37]))
names(Yersinia_pseudot_TPM_expr)=Yersinia_pseudot_TPM[,1]


###########################Enterococcus_faecalis################

#Reading the data
Enterococcus_faecalis= read.csv("Enterococcus_faecalis_R_input.csv", sep=";")
dim(Enterococcus_faecalis) #analising input

#getting TPM values
Enterococcus_faecalis_expr =as.data.frame(Enterococcus_faecalis[,2:ncol(Enterococcus_faecalis)])
dim(Enterococcus_faecalis_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Enterococcus_faecalis_expr) = Enterococcus_faecalis[,1]


###########################Shigella_flexneri################

#Reading the data
Shigella_flexneri= read.csv("Shigella_flexneri_R_input.csv", sep=";")
dim(Shigella_flexneri) #analising input

#getting TPM values
Shigella_flexneri_expr =as.data.frame(Shigella_flexneri[,2:ncol(Shigella_flexneri)])
dim(Shigella_flexneri_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Shigella_flexneri_expr) = Shigella_flexneri[,1]

#################Vibrio_cholerae#####################

#Reading the data
Vibrio_cholerae= read.csv("Vibrio_cholerae_R_input.csv", sep=";")
dim(Vibrio_cholerae) #analising input

#getting TPM values
Vibrio_cholerae_expr =as.data.frame(Vibrio_cholerae[,2:ncol(Vibrio_cholerae)])
dim(Vibrio_cholerae_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Vibrio_cholerae_expr) = Vibrio_cholerae[,1]


##################PGFams##########################
PGFams_locustags13df64a3fa583c4875d73be0d98473ba3c3d01cbe802357f34ddfce7e49976a1 <- read_excel("PGFams_locustags13df64a3fa583c4875d73be0d98473ba3c3d01cbe802357f34ddfce7e49976a1.xlsx")

Annotations_all_PGFAM <- PGFams_locustags13df64a3fa583c4875d73be0d98473ba3c3d01cbe802357f34ddfce7e49976a1

Annotations_all_PGFAM_new <- Annotations_all_PGFAM[ which(Annotations_all_PGFAM$comp_genome_1_percent_identity > 0.95  & Annotations_all_PGFAM$comp_genome_1_seq_coverage > 0.90), ]
Annotations_core <- read_csv("Annotations_core.csv")
###################E_Coli_EPEC_UPEC_ETEC################################
Locus_Tag_UPEC= colnames(Escherichia_coli_UPEC_536_expr)
probes_UPEC = Locus_Tag_UPEC
Genes=Annotations_all_PGFAM_new$ref_genome_pgfam_id[match(probes_UPEC,Annotations_all_PGFAM_new$comp_genome_1_patric_id)]
UPEC_STRAIN_expr = as.data.frame(cbind(Genes,t(Escherichia_coli_UPEC_536_expr)))
#attach(UPEC_STRAIN_expr)
UPEC_STRAIN_expr_df = UPEC_STRAIN_expr[!is.na(UPEC_STRAIN_expr$Genes), ]
UPEC_STRAIN_expr_df2= UPEC_STRAIN_expr_df[, 2:ncol(UPEC_STRAIN_expr_df)]
UPEC_STRAIN_expr_df2[is.na(UPEC_STRAIN_expr_df2)] <- 0
rownames(UPEC_STRAIN_expr_df2)= NULL
rownames(UPEC_STRAIN_expr_df2)= make.names(UPEC_STRAIN_expr_df$Genes, unique=TRUE)




Locus_Tag_EPEC= colnames(Escherichia_coli_EPEC_0127_H6_E2348_69_expr)
probes_EPEC = Locus_Tag_EPEC
Genes=Annotations_all_PGFAM_new$ref_genome_pgfam_id[match(probes_EPEC,Annotations_all_PGFAM_new$comp_genome_1_patric_id)]
EPEC_STRAIN_expr = as.data.frame(cbind(Genes,t(Escherichia_coli_EPEC_0127_H6_E2348_69_expr)))
EPEC_STRAIN_expr_df = EPEC_STRAIN_expr[!is.na(EPEC_STRAIN_expr$Genes), ]
EPEC_STRAIN_expr_df2= EPEC_STRAIN_expr_df[, 2:ncol(EPEC_STRAIN_expr_df)]
EPEC_STRAIN_expr_df2[is.na(EPEC_STRAIN_expr_df2)] <- 0
rownames(EPEC_STRAIN_expr_df2)= NULL
rownames(EPEC_STRAIN_expr_df2)= make.names(EPEC_STRAIN_expr_df$Genes, unique=TRUE)


Locus_Tag_ETEC= colnames(Escherichia_coli_ETEC_H10407_expr)
probes_ETEC = Locus_Tag_ETEC
Genes=Annotations_all_PGFAM_new$ref_genome_pgfam_id[match(probes_ETEC,Annotations_all_PGFAM_new$comp_genome_1_patric_id)]
ETEC_STRAIN_expr = as.data.frame(cbind(Genes,t(Escherichia_coli_ETEC_H10407_expr)))
#attach(ETEC_STRAIN_expr)
ETEC_STRAIN_expr_df = ETEC_STRAIN_expr[!is.na(ETEC_STRAIN_expr$Genes), ]
ETEC_STRAIN_expr_df2= ETEC_STRAIN_expr_df[, 2:ncol(ETEC_STRAIN_expr_df)]
ETEC_STRAIN_expr_df2[is.na(ETEC_STRAIN_expr_df2)] <- 0
rownames(ETEC_STRAIN_expr_df2)= NULL
rownames(ETEC_STRAIN_expr_df2)= make.names(ETEC_STRAIN_expr_df$Genes, unique=TRUE)

#####UPECandEPEC#######

UPEC_STRAIN_expr_t= t(UPEC_STRAIN_expr)
EPEC_STRAIN_expr_t= t(EPEC_STRAIN_expr)
list_of_data = list(UPEC_STRAIN_expr_df2,EPEC_STRAIN_expr_df2)
common_names = Reduce(intersect, lapply(list_of_data, row.names))
list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,] })

UPEC_commongenes_1 = as.data.frame(list_of_data[1])
EPEC_commongenes_1 = as.data.frame(list_of_data[2])

UPEC_commongenes=as.data.frame(sapply(UPEC_commongenes_1, as.numeric))
rownames(UPEC_commongenes)= rownames(UPEC_commongenes_1)
EPEC_commongenes=as.data.frame(sapply(EPEC_commongenes_1, as.numeric))
rownames(EPEC_commongenes)= rownames(EPEC_commongenes_1)

nSets=2
setLabels_analysis=NULL
setLabels_analysis["UPEC vs EPEC"] = list(c("Strain UPEC", "Strain EPEC"))
shortLabels = c("UPEC", "EPEC")
UPEC_EPEC_multiExpr = vector(mode = "list", length = nSets)

UPEC_commongenes_t = t(UPEC_commongenes)
EPEC_commongenes_t = t(EPEC_commongenes)
UPEC_commongenes_t_ord = UPEC_commongenes_t[ , order(colnames(UPEC_commongenes_t))]
EPEC_commongenes_t_ord = EPEC_commongenes_t[ , order(colnames(EPEC_commongenes_t))]

UPEC_EPEC_multiExpr[[1]] = list(data = UPEC_commongenes_t_ord)
UPEC_EPEC_multiExpr[[2]] = list(data = EPEC_commongenes_t_ord)
names(UPEC_EPEC_multiExpr[[1]]$data) = colnames(UPEC_commongenes_t_ord)
names(UPEC_EPEC_multiExpr[[2]]$data) = colnames(EPEC_commongenes_t_ord)



#####UPECandETEC#######

UPEC_STRAIN_expr_t= t(UPEC_STRAIN_expr)
ETEC_STRAIN_expr_t= t(ETEC_STRAIN_expr)
list_of_data = list(UPEC_STRAIN_expr_df2,ETEC_STRAIN_expr_df2)
common_names = Reduce(intersect, lapply(list_of_data, row.names))
list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,] })

UPEC_commongenes_1 = as.data.frame(list_of_data[1])
ETEC_commongenes_1 = as.data.frame(list_of_data[2])

UPEC_commongenes=as.data.frame(sapply(UPEC_commongenes_1, as.numeric))
rownames(UPEC_commongenes)= rownames(UPEC_commongenes_1)
ETEC_commongenes=as.data.frame(sapply(ETEC_commongenes_1, as.numeric))
rownames(ETEC_commongenes)= rownames(ETEC_commongenes_1)

nSets=2
setLabels_analysis["UPEC vs ETEC"]= list(c("Strain UPEC", "Strain ETEC"))
shortLabels = c("UPEC", "ETEC")
UPEC_ETEC_multiExpr = vector(mode = "list", length = nSets)

UPEC_commongenes_t = t(UPEC_commongenes)
ETEC_commongenes_t = t(ETEC_commongenes)
UPEC_commongenes_t_ord = UPEC_commongenes_t[ , order(colnames(UPEC_commongenes_t))]
ETEC_commongenes_t_ord = ETEC_commongenes_t[ , order(colnames(ETEC_commongenes_t))]

UPEC_ETEC_multiExpr[[1]] = list(data = UPEC_commongenes_t_ord)
UPEC_ETEC_multiExpr[[2]] = list(data = ETEC_commongenes_t_ord)

names(UPEC_ETEC_multiExpr[[1]]$data) = colnames(UPEC_commongenes_t_ord)
names(UPEC_ETEC_multiExpr[[2]]$data) = colnames(ETEC_commongenes_t_ord)

#####EPECandETEC#######

EPEC_STRAIN_expr_t= t(EPEC_STRAIN_expr)
ETEC_STRAIN_expr_t= t(ETEC_STRAIN_expr)
list_of_data = list(EPEC_STRAIN_expr_df2,ETEC_STRAIN_expr_df2)
common_names = Reduce(intersect, lapply(list_of_data, row.names))
list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,] })

EPEC_commongenes_1 = as.data.frame(list_of_data[1])
ETEC_commongenes_1 = as.data.frame(list_of_data[2])

EPEC_commongenes=as.data.frame(sapply(EPEC_commongenes_1, as.numeric))
rownames(EPEC_commongenes)= rownames(EPEC_commongenes_1)
ETEC_commongenes=as.data.frame(sapply(ETEC_commongenes_1, as.numeric))
rownames(ETEC_commongenes)= rownames(ETEC_commongenes_1)

nSets=2
setLabels_analysis["EPEC vs ETEC"]= list(c("Strain EPEC", "Strain ETEC"))
shortLabels = c("EPEC", "ETEC")
EPEC_ETEC_multiExpr = vector(mode = "list", length = nSets)

EPEC_commongenes_t = t(EPEC_commongenes)
ETEC_commongenes_t = t(ETEC_commongenes)
EPEC_commongenes_t_ord = EPEC_commongenes_t[ , order(colnames(EPEC_commongenes_t))]
ETEC_commongenes_t_ord = ETEC_commongenes_t[ , order(colnames(ETEC_commongenes_t))]

EPEC_ETEC_multiExpr[[1]] = list(data = EPEC_commongenes_t_ord)
EPEC_ETEC_multiExpr[[2]] = list(data = ETEC_commongenes_t_ord)

names(EPEC_ETEC_multiExpr[[1]]$data) = colnames(EPEC_commongenes_t_ord)
names(EPEC_ETEC_multiExpr[[2]]$data) = colnames(ETEC_commongenes_t_ord)

##The three strains together#####

EPEC_STRAIN_expr_t= t(EPEC_STRAIN_expr)
ETEC_STRAIN_expr_t= t(ETEC_STRAIN_expr)
UPEC_STRAIN_expr_t= t(UPEC_STRAIN_expr)
list_of_data = list(EPEC_STRAIN_expr_df2,ETEC_STRAIN_expr_df2,UPEC_STRAIN_expr_df2)
common_names = Reduce(intersect, lapply(list_of_data, row.names))
list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,] })

EPEC_commongenes_1 = as.data.frame(list_of_data[1])
ETEC_commongenes_1= as.data.frame(list_of_data[2])
UPEC_commongenes_1= as.data.frame(list_of_data[3])

EPEC_commongenes=as.data.frame(sapply(EPEC_commongenes_1, as.numeric))
rownames(EPEC_commongenes)= rownames(EPEC_commongenes_1)
ETEC_commongenes=as.data.frame(sapply(ETEC_commongenes_1, as.numeric))
rownames(ETEC_commongenes)= rownames(ETEC_commongenes_1)
UPEC_commongenes=as.data.frame(sapply(UPEC_commongenes_1, as.numeric))
rownames(UPEC_commongenes)= rownames(UPEC_commongenes_1)


nSets=3
setLabels_analysis["EPEC vs ETEC vs UPEC"]= list(c("Strain EPEC", "Strain ETEC","Strain UPEC" ))
setLabels = c("EPEC", "ETEC", "UPEC")
EPEC_ETEC_UPEC_multiExpr = vector(mode = "list", length = nSets)

EPEC_commongenes_t = t(EPEC_commongenes)
ETEC_commongenes_t = t(ETEC_commongenes)

UPEC_commongenes_t = t(UPEC_commongenes)
UPEC_commongenes_t_ord = UPEC_commongenes_t[ , order(colnames(UPEC_commongenes_t))]


EPEC_commongenes_t_ord = EPEC_commongenes_t[ , order(colnames(EPEC_commongenes_t))]
ETEC_commongenes_t_ord = ETEC_commongenes_t[ , order(colnames(ETEC_commongenes_t))]

EPEC_ETEC_UPEC_multiExpr[[1]] = list(data = EPEC_commongenes_t_ord)
EPEC_ETEC_UPEC_multiExpr[[2]] = list(data = ETEC_commongenes_t_ord)
EPEC_ETEC_UPEC_multiExpr[[3]] = list(data = UPEC_commongenes_t_ord)

names(EPEC_ETEC_UPEC_multiExpr[[1]]$data) = colnames(EPEC_commongenes_t_ord)
names(EPEC_ETEC_UPEC_multiExpr[[2]]$data) = colnames(ETEC_commongenes_t_ord)
names(EPEC_ETEC_UPEC_multiExpr[[3]]$data) = colnames(UPEC_commongenes_t_ord)




###################HP_G27_HP_J99################################
Locus_Tag_HP_G27= colnames(Helicobacter_pylori_G27_expr)
probes_HP_G27 = Locus_Tag_HP_G27
Genes=Annotations_all_PGFAM_new$ref_genome_pgfam_id[match(probes_HP_G27,Annotations_all_PGFAM_new$comp_genome_1_patric_id)]
HP_G27_STRAIN_expr = as.data.frame(cbind(Genes,t(Helicobacter_pylori_G27_expr)))
#attach(HP_G27_STRAIN_expr)
HP_G27_STRAIN_expr_df = HP_G27_STRAIN_expr[!is.na(HP_G27_STRAIN_expr$Genes), ]
HP_G27_STRAIN_expr_df2= HP_G27_STRAIN_expr_df[, 2:ncol(HP_G27_STRAIN_expr_df)]
HP_G27_STRAIN_expr_df2[is.na(HP_G27_STRAIN_expr_df2)] <- 0
rownames(HP_G27_STRAIN_expr_df2)= NULL
rownames(HP_G27_STRAIN_expr_df2)= make.names(HP_G27_STRAIN_expr_df$Genes, unique=TRUE)


Locus_Tag_HP_J99= colnames(Helicobacter_pylori_J99_expr)
probes_HP_J99 = Locus_Tag_HP_J99
Genes=Annotations_all_PGFAM_new$ref_genome_pgfam_id[match(probes_HP_J99,Annotations_all_PGFAM_new$comp_genome_1_patric_id)]
HP_J99_STRAIN_expr = as.data.frame(cbind(Genes,t(Helicobacter_pylori_J99_expr)))
#attach(HP_J99_STRAIN_expr)
HP_J99_STRAIN_expr_df = HP_J99_STRAIN_expr[!is.na(HP_J99_STRAIN_expr$Genes), ]
HP_J99_STRAIN_expr_df2= HP_J99_STRAIN_expr_df[, 2:ncol(HP_J99_STRAIN_expr_df)]
HP_J99_STRAIN_expr_df2[is.na(HP_J99_STRAIN_expr_df2)] <- 0
rownames(HP_J99_STRAIN_expr_df2)= NULL
rownames(HP_J99_STRAIN_expr_df2)= make.names(HP_J99_STRAIN_expr_df$Genes, unique=TRUE)

HP_J99_STRAIN_expr_t= t(HP_J99_STRAIN_expr)
HP_G27_STRAIN_expr_t= t(HP_G27_STRAIN_expr)
list_of_data = list(HP_J99_STRAIN_expr_df2,HP_G27_STRAIN_expr_df2)
common_names = Reduce(intersect, lapply(list_of_data, row.names))
list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,] })

HP_J99_commongenes_1 = as.data.frame(list_of_data[1])
HP_G27_commongenes_1 = as.data.frame(list_of_data[2])

HP_J99_commongenes=as.data.frame(sapply(HP_J99_commongenes_1, as.numeric))
rownames(HP_J99_commongenes)= rownames(HP_J99_commongenes_1)
HP_G27_commongenes=as.data.frame(sapply(HP_G27_commongenes_1, as.numeric))
rownames(HP_G27_commongenes)= rownames(HP_G27_commongenes_1)

nSets=2
setLabels_analysis["HP_J99 vs HP_G27"]= list(c("Strain HP_J99", "Strain HP_G27"))
shortLabels = c("HP_J99", "HP_G27")
HP_J99_HP_G27_multiExpr = vector(mode = "list", length = nSets)

HP_J99_commongenes_t = t(HP_J99_commongenes)
HP_G27_commongenes_t = t(HP_G27_commongenes)
HP_J99_commongenes_t_ord = HP_J99_commongenes_t[ , order(colnames(HP_J99_commongenes_t))]
HP_G27_commongenes_t_ord = HP_G27_commongenes_t[ , order(colnames(HP_G27_commongenes_t))]

HP_J99_HP_G27_multiExpr[[1]] = list(data = HP_J99_commongenes_t_ord)
HP_J99_HP_G27_multiExpr[[2]] = list(data = HP_G27_commongenes_t_ord)

names(HP_J99_HP_G27_multiExpr[[1]]$data) = colnames(HP_J99_commongenes_t_ord)
names(HP_J99_HP_G27_multiExpr[[2]]$data) = colnames(HP_G27_commongenes_t_ord)


###################MRSA252_MSSA476################################
Locus_Tag_MSSA476= colnames(Staphylococcus_aureus_MSSA476_expr)
probes_MSSA476 = Locus_Tag_MSSA476
Genes=Annotations_all_PGFAM_new$ref_genome_pgfam_id[match(probes_MSSA476,Annotations_all_PGFAM_new$comp_genome_1_patric_id)]
MSSA476_STRAIN_expr = as.data.frame(cbind(Genes,t(Staphylococcus_aureus_MSSA476_expr)))
#attach(MSSA476_STRAIN_expr)
MSSA476_STRAIN_expr_df = MSSA476_STRAIN_expr[!is.na(MSSA476_STRAIN_expr$Genes), ]
MSSA476_STRAIN_expr_df2= MSSA476_STRAIN_expr_df[, 2:ncol(MSSA476_STRAIN_expr_df)]
MSSA476_STRAIN_expr_df2[is.na(MSSA476_STRAIN_expr_df2)] <- 0
rownames(MSSA476_STRAIN_expr_df2)= NULL
rownames(MSSA476_STRAIN_expr_df2)= make.names(MSSA476_STRAIN_expr_df$Genes, unique=TRUE)

Locus_Tag_MRSA252= colnames(Staphylococcus_aureus_MRSA252_expr)
probes_MRSA252 = Locus_Tag_MRSA252
Genes=Annotations_all_PGFAM_new$ref_genome_pgfam_id[match(probes_MRSA252,Annotations_all_PGFAM_new$comp_genome_1_patric_id)]
MRSA252_STRAIN_expr = as.data.frame(cbind(Genes,t(Staphylococcus_aureus_MRSA252_expr)))
#attach(MRSA252_STRAIN_expr)
MRSA252_STRAIN_expr_df = MRSA252_STRAIN_expr[!is.na(MRSA252_STRAIN_expr$Genes), ]
MRSA252_STRAIN_expr_df2= MRSA252_STRAIN_expr_df[, 2:ncol(MRSA252_STRAIN_expr_df)]
MRSA252_STRAIN_expr_df2[is.na(MRSA252_STRAIN_expr_df2)] <- 0
rownames(MRSA252_STRAIN_expr_df2)= NULL
rownames(MRSA252_STRAIN_expr_df2)= make.names(MRSA252_STRAIN_expr_df$Genes, unique=TRUE)


MSSA476_STRAIN_expr_t= t(MSSA476_STRAIN_expr)
MRSA252_STRAIN_expr_t= t(MRSA252_STRAIN_expr)
list_of_data = list(MSSA476_STRAIN_expr_df2,MRSA252_STRAIN_expr_df2)
common_names = Reduce(intersect, lapply(list_of_data, row.names))
list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,] })

MSSA476_commongenes_1 = as.data.frame(list_of_data[1])
MRSA252_commongenes_1 = as.data.frame(list_of_data[2])

MSSA476_commongenes=as.data.frame(sapply(MSSA476_commongenes_1, as.numeric))
rownames(MSSA476_commongenes)= rownames(MSSA476_commongenes_1)
MRSA252_commongenes=as.data.frame(sapply(MRSA252_commongenes_1, as.numeric))
rownames(MRSA252_commongenes)= rownames(MRSA252_commongenes_1)

nSets=2
setLabels_analysis["MSSA476 vs MRSA252"]= list(c("Strain MSSA476", "Strain MRSA252"))
shortLabels = c("MSSA476", "MRSA252")
MSSA476_MRSA252_multiExpr = vector(mode = "list", length = nSets)

MSSA476_commongenes_t = t(MSSA476_commongenes)
MRSA252_commongenes_t = t(MRSA252_commongenes)
MSSA476_commongenes_t_ord = MSSA476_commongenes_t[ , order(colnames(MSSA476_commongenes_t))]
MRSA252_commongenes_t_ord = MRSA252_commongenes_t[ , order(colnames(MRSA252_commongenes_t))]

MSSA476_MRSA252_multiExpr[[1]] = list(data = MSSA476_commongenes_t_ord)
MSSA476_MRSA252_multiExpr[[2]] = list(data = MRSA252_commongenes_t_ord)

names(MSSA476_MRSA252_multiExpr[[1]]$data) = colnames(MSSA476_commongenes_t_ord)
names(MSSA476_MRSA252_multiExpr[[2]]$data) = colnames(MRSA252_commongenes_t_ord)

#################################ACID_STRESS###################
########UPEC_As_Sp###############
Escherichia_coli_UPEC_536_As = read.csv("UPEC_As_1As_2As_3Oss_1_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Escherichia_coli_UPEC_536_expr_As =as.data.frame(Escherichia_coli_UPEC_536_As[,2:ncol(Escherichia_coli_UPEC_536_As)])
rownames(Escherichia_coli_UPEC_536_expr_As) = Escherichia_coli_UPEC_536_As[,1]

########EPEC_As1###############
Escherichia_coli_EPEC_0127_H6_E2348_69_As = read.csv("EPEC_As_1As_2As_3_UP.csv", sep=",")
Escherichia_coli_EPEC_0127_H6_E2348_69_expr_As =as.data.frame(Escherichia_coli_EPEC_0127_H6_E2348_69_As[,2:ncol(Escherichia_coli_EPEC_0127_H6_E2348_69_As)])
rownames(Escherichia_coli_EPEC_0127_H6_E2348_69_expr_As) = Escherichia_coli_EPEC_0127_H6_E2348_69_As[,1]

########EPEC_As2###############
Escherichia_coli_EPEC_0127_H6_E2348_69_As2= read.csv("EPEC_As_1As_2As_3Mig_1Mig_2_UP.csv", sep=",")
Escherichia_coli_EPEC_0127_H6_E2348_69_expr_As2 =as.data.frame(Escherichia_coli_EPEC_0127_H6_E2348_69_As2[,2:ncol(Escherichia_coli_EPEC_0127_H6_E2348_69_As2)])
rownames(Escherichia_coli_EPEC_0127_H6_E2348_69_expr_As2) = Escherichia_coli_EPEC_0127_H6_E2348_69_As2[,1]

##########ENTFA_As##############
Enterococcus_faecalis_As = read.csv("ENTFA_As_1As_2As_3_UP.csv", sep=",")
Enterococcus_faecalis_expr_As =as.data.frame(Enterococcus_faecalis_As[,2:ncol(Enterococcus_faecalis_As)])
rownames(Enterococcus_faecalis_expr_As) = Enterococcus_faecalis_As[,1]

#########Listeria_As#########
Listeria_monocytogenes_As = read.csv("Listeria_As_1As_2As_3_UP.csv", sep=",")
Listeria_monocytogenes_expr_As =as.data.frame(Listeria_monocytogenes_As[,2:ncol(Listeria_monocytogenes_As)])
rownames(Listeria_monocytogenes_expr_As) = Listeria_monocytogenes_As[,1]

#########PSEUDO_As#########
Pseudomonas_aeruginosa_As = read.csv("PSEUDO_As_1As_2As_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_As =as.data.frame(Pseudomonas_aeruginosa_As[,2:ncol(Pseudomonas_aeruginosa_As)])
rownames(Pseudomonas_aeruginosa_expr_As) = Pseudomonas_aeruginosa_As[,1]

#########PSEUDO_As_Tm#########
Pseudomonas_aeruginosa_As_Tm = read.csv("PSEUDO_As_1As_2As_3Tm_1Tm_2Tm_3_UP_Nd_1Nd_2_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_As_Tm =as.data.frame(Pseudomonas_aeruginosa_As_Tm[,2:ncol(Pseudomonas_aeruginosa_As_Tm)])
rownames(Pseudomonas_aeruginosa_expr_As_Tm) = Pseudomonas_aeruginosa_As_Tm[,1]

###########Escherichia_coli_ETEC_H10407_As############
Escherichia_coli_ETEC_H10407_As =read.csv("ETEC_As_1As_2As_3_UP.csv",  sep=",")
Escherichia_coli_ETEC_H10407_expr_As =as.data.frame(Escherichia_coli_ETEC_H10407_As[,2:ncol(Escherichia_coli_ETEC_H10407_As)])
rownames(Escherichia_coli_ETEC_H10407_expr_As) = Escherichia_coli_ETEC_H10407_As[,1]

###########Salmonella_As############
Salmonella_As =read.csv("SALMT_As_1As_2As_3_UP.csv",  sep=",")
Salmonella_expr_As =as.data.frame(Salmonella_As[,2:ncol(Salmonella_As)])
rownames(Salmonella_expr_As) = Salmonella_As[,1]

###########MRSA252_As############
Staphylococcus_aureus_MRSA252_As =read.csv("MRSA252_As_1As_2As_3_UP_Ns_2_DOWN.csv",  sep=",")
Staphylococcus_aureus_MRSA252_expr_As =as.data.frame(Staphylococcus_aureus_MRSA252_As[,2:ncol(Staphylococcus_aureus_MRSA252_As)])
rownames(Staphylococcus_aureus_MRSA252_expr_As) = Staphylococcus_aureus_MRSA252_As[,1]

###########MSSA476_As############
Staphylococcus_aureus_MSSA476_As =read.csv("MSSA476_As_1As_2As_3Bs_3_UP_Nd_1Nd_3_DOWN.csv",  sep=",")
Staphylococcus_aureus_MSSA476_expr_As =as.data.frame(Staphylococcus_aureus_MSSA476_As[,2:ncol(Staphylococcus_aureus_MSSA476_As)])
rownames(Staphylococcus_aureus_MSSA476_expr_As) = Staphylococcus_aureus_MSSA476_As[,1]

###########MSSA476_As_Ns############
Staphylococcus_aureus_MSSA476_As_Ns =read.csv("MSSA476_As_1As_2As_3_UP_Ns_1Ns_2Ns_3_DOWN.csv",  sep=",")
Staphylococcus_aureus_MSSA476_expr_As_Ns =as.data.frame(Staphylococcus_aureus_MSSA476_As_Ns[,2:ncol(Staphylococcus_aureus_MSSA476_As_Ns)])
rownames(Staphylococcus_aureus_MSSA476_expr_As_Ns) = Staphylococcus_aureus_MSSA476_As_Ns[,1]

###########KLEBS_As_Sp############
Klebsiella_pneumoniae_As_Sp =read.csv("KLEBS_As_1As_2As_3_UP_Sp_1Sp_2Sp_3_DOWN.csv",  sep=",")
Klebsiella_pneumoniae_expr_As_Sp =as.data.frame(Klebsiella_pneumoniae_As_Sp[,2:ncol(Klebsiella_pneumoniae_As_Sp)])
rownames(Klebsiella_pneumoniae_expr_As_Sp) = Klebsiella_pneumoniae_As_Sp[,1]

###########SHIF_As############
Shigella_flexneri_As =read.csv("SHIF_As_1As_2As_3Ctrl_2Ctrl_3_UP_Sp_2Sp_3_DOWN.csv",  sep=",")
Shigella_flexneri_expr_As =as.data.frame(Shigella_flexneri_As[,2:ncol(Shigella_flexneri_As)])
rownames(Shigella_flexneri_expr_As) = Shigella_flexneri_As[,1]

###########SEPI_As############
Staphylococcus_epidermidis_As =read.csv("SEPI_As_1As_2As_3_UP_Tm_1_DOWN.csv",  sep=",")
Staphylococcus_epidermidis_expr_As =as.data.frame(Staphylococcus_epidermidis_As[,2:ncol(Staphylococcus_epidermidis_As)])
rownames(Staphylococcus_epidermidis_expr_As) = Staphylococcus_epidermidis_As[,1]

###########BURK_As############
Burkholderia_pseudomallei_As =read.csv("BURK_As_1As_2As_3_UP_Sp_1_DOWN.csv",  sep=",")
Burkholderia_pseudomallei_expr_As =as.data.frame(Burkholderia_pseudomallei_As[,2:ncol(Burkholderia_pseudomallei_As)])
rownames(Burkholderia_pseudomallei_expr_As) = Burkholderia_pseudomallei_As[,1]

##################SEPI_As_up_Mig_up_##########################
Staphylococcus_epidermidis_Mig_up_As_up=read.csv("SEPI_As_1As_2As_3Mig_1Mig_2Mig_3_UP.csv", sep=",")
Staphylococcus_epidermidis_expr_Mig_up_As_up=as.data.frame(Staphylococcus_epidermidis_Mig_up_As_up[,2:ncol(Staphylococcus_epidermidis_Mig_up_As_up)])
rownames(Staphylococcus_epidermidis_expr_Mig_up_As_up) =Staphylococcus_epidermidis_Mig_up_As_up[,1]

################Listeria_As_up_Oxs_up###############
Listeria_monocytogenes_Oxs_up_As_up=read.csv("Listeria_As_1As_2As_3Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Listeria_monocytogenes_expr_Oxs_up_As_up=as.data.frame(Listeria_monocytogenes_Oxs_up_As_up[,2:ncol(Listeria_monocytogenes_Oxs_up_As_up)])
rownames(Listeria_monocytogenes_expr_Oxs_up_As_up) =Listeria_monocytogenes_Oxs_up_As_up[,1]

#########Listeria_As_Up_Nd_down#########
Listeria_monocytogenes_As_up_Nd_down = read.csv("Listeria_As_1As_2As_3Mig_1Mig_2_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Listeria_monocytogenes_expr_As_up_Nd_down =as.data.frame(Listeria_monocytogenes_As_up_Nd_down[,2:ncol(Listeria_monocytogenes_As_up_Nd_down)])
rownames(Listeria_monocytogenes_expr_As_up_Nd_down) = Listeria_monocytogenes_As_up_Nd_down[,1]

##################SEPI_As_up_Tm_down##########################
Staphylococcus_epidermidis_As_up_Tm_down=read.csv("SEPI_As_1As_2As_3Li_1Li_3_UP_Sp_2Tm_1Tm_2Tm_3_DOWN.csv", sep=",")
Staphylococcus_epidermidis_expr_As_up_Tm_down=as.data.frame(Staphylococcus_epidermidis_As_up_Tm_down[,2:ncol(Staphylococcus_epidermidis_As_up_Tm_down)])
rownames(Staphylococcus_epidermidis_expr_As_up_Tm_down) =Staphylococcus_epidermidis_As_up_Tm_down[,1]


##########SP_STRESS################
###########NGON_Sp############
Neisseria_gonorrhoeae_Sp =read.csv("NGON_As_1As_2Ctrl_2_UP_Li_2Sp_1Sp_2Sp_3_DOWN.csv",  sep=",")
Neisseria_gonorrhoeae_expr_Sp =as.data.frame(Neisseria_gonorrhoeae_Sp[,2:ncol(Neisseria_gonorrhoeae_Sp)])
rownames(Neisseria_gonorrhoeae_expr_Sp) = Neisseria_gonorrhoeae_Sp[,1]

########EPEC_Sp###############
Escherichia_coli_EPEC_0127_H6_E2348_69_Sp = read.csv("EPEC_As_1As_2Vic_3_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Escherichia_coli_EPEC_0127_H6_E2348_69_expr_Sp =as.data.frame(Escherichia_coli_EPEC_0127_H6_E2348_69_Sp[,2:ncol(Escherichia_coli_EPEC_0127_H6_E2348_69_Sp)])
rownames(Escherichia_coli_EPEC_0127_H6_E2348_69_expr_Sp) = Escherichia_coli_EPEC_0127_H6_E2348_69_Sp[,1]

########ACIB_Sp1###############
Acinetobacter_baumannii_Sp1 = read.csv("ACIB_Oxs_2_UP_Mig_1Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Acinetobacter_baumannii_expr_Sp1 =as.data.frame(Acinetobacter_baumannii_Sp1[,2:ncol(Acinetobacter_baumannii_Sp1)])
rownames(Acinetobacter_baumannii_expr_Sp1) = Acinetobacter_baumannii_Sp1[,1]

########ACIB_Sp2###############
Acinetobacter_baumannii_Sp2 = read.csv("ACIB_Mig_2Mig_3_UP_Mig_1Nd_1Nd_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Acinetobacter_baumannii_expr_Sp2 =as.data.frame(Acinetobacter_baumannii_Sp2[,2:ncol(Acinetobacter_baumannii_Sp2)])
rownames(Acinetobacter_baumannii_expr_Sp2) = Acinetobacter_baumannii_Sp2[,1]

########ACIB_Sp3###############
Acinetobacter_baumannii_Sp3 = read.csv("ACIB_Li_1_UP_Mig_1Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Acinetobacter_baumannii_expr_Sp3 =as.data.frame(Acinetobacter_baumannii_Sp3[,2:ncol(Acinetobacter_baumannii_Sp3)])
rownames(Acinetobacter_baumannii_expr_Sp3) = Acinetobacter_baumannii_Sp3[,1]


########ACIB_Sp4###############
Acinetobacter_baumannii_Sp4 = read.csv("ACIB_Vic_1Vic_2_UP_Mig_1Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Acinetobacter_baumannii_expr_Sp4 =as.data.frame(Acinetobacter_baumannii_Sp4[,2:ncol(Acinetobacter_baumannii_Sp4)])
rownames(Acinetobacter_baumannii_expr_Sp4) = Acinetobacter_baumannii_Sp4[,1]

########SPYO_Sp_Mig_down1###############
Streptococcus_pyogenes_Sp_Mg = read.csv("SPYO_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Streptococcus_pyogenes_expr_Sp_Mg =as.data.frame(Streptococcus_pyogenes_Sp_Mg[,2:ncol(Streptococcus_pyogenes_Sp_Mg)])
rownames(Streptococcus_pyogenes_expr_Sp_Mg) = Streptococcus_pyogenes_Sp_Mg[,1]

################SPYO_Sp_Mig_down2###############
Streptococcus_pyogenes_Mig2_down_Sp_down=read.csv("SPYO_Vic_1_UP_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Streptococcus_pyogenes_expr_Mig2_down_Sp_down =as.data.frame(Streptococcus_pyogenes_Mig2_down_Sp_down[,2:ncol(Streptococcus_pyogenes_Mig2_down_Sp_down)])
rownames(Streptococcus_pyogenes_expr_Mig2_down_Sp_down) =Streptococcus_pyogenes_Mig2_down_Sp_down[,1]

########BURK_Sp_down###############
Burkholderia_pseudomallei_Sp_down = read.csv("BURK_As_2Oxs_3_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Burkholderia_pseudomallei_expr_Sp_down =as.data.frame(Burkholderia_pseudomallei_Sp_down[,2:ncol(Burkholderia_pseudomallei_Sp_down)])
rownames(Burkholderia_pseudomallei_expr_Sp_down) = Burkholderia_pseudomallei_Sp_down[,1]

########BURK_Sp_up_createdbyCansu###############
Burkholderia_pseudomallei_Sp_up=read.csv("BURK_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Sp_up =as.data.frame(Burkholderia_pseudomallei_Sp_up[,2:ncol(Burkholderia_pseudomallei_Sp_up)])
rownames(Burkholderia_pseudomallei_expr_Sp_up) = Burkholderia_pseudomallei_Sp_up[,1]

########ENTFA_Sp_up_createdbyCansu###############
Enterococcus_faecalis_Sp_up=read.csv("ENTFA_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Enterococcus_faecalis_expr_Sp_up =as.data.frame(Enterococcus_faecalis_Sp_up[,2:ncol(Enterococcus_faecalis_Sp_up)])
rownames(Enterococcus_faecalis_expr_Sp_up) =Enterococcus_faecalis_Sp_up[,1]

################ACHX_Sp_createdbyCansu###############
Achromobacter_xylosoxidans_Sp=read.csv("ACHX_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Achromobacter_xylosoxidans_expr_Sp =as.data.frame(Achromobacter_xylosoxidans_Sp[,2:ncol(Achromobacter_xylosoxidans_Sp)])
rownames(Achromobacter_xylosoxidans_expr_Sp) =Achromobacter_xylosoxidans_Sp[,1]

################BBURG_Sp_createdbyCansu###############
Borrelia_burgdorferi_Sp=read.csv("BBURG_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Borrelia_burgdorferi_expr_Sp =as.data.frame(Borrelia_burgdorferi_Sp[,2:ncol(Borrelia_burgdorferi_Sp)])
rownames(Borrelia_burgdorferi_expr_Sp) =Borrelia_burgdorferi_Sp[,1]

################BBURG_Sp2_createdbyCansu###############
Borrelia_burgdorferi_Sp2=read.csv("BBURG_Oss_1Oss_2Sp_1Sp_2Sp_3_UP.csv", sep=",")
Borrelia_burgdorferi_expr_Sp2 =as.data.frame(Borrelia_burgdorferi_Sp2[,2:ncol(Borrelia_burgdorferi_Sp2)])
rownames(Borrelia_burgdorferi_expr_Sp2) =Borrelia_burgdorferi_Sp2[,1]

################BURK_Sp_Tm_createdbyCansu###############
Burkholderia_pseudomallei_Sp_Tm=read.csv("BURK_Sp_1Sp_2Sp_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Sp_Tm =as.data.frame(Burkholderia_pseudomallei_Sp_Tm[,2:ncol(Burkholderia_pseudomallei_Sp_Tm)])
rownames(Burkholderia_pseudomallei_expr_Sp_Tm) =Burkholderia_pseudomallei_Sp_Tm[,1]

################CAMPY_Sp_createdbyCansu###############
Campylobacter_jejuni_Sp=read.csv("Campy_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Campylobacter_jejuni_expr_Sp =as.data.frame(Campylobacter_jejuni_Sp[,2:ncol(Campylobacter_jejuni_Sp)])
rownames(Campylobacter_jejuni_expr_Sp) =Campylobacter_jejuni_Sp[,1]


################EPEC_Sp_createdbyCansu###############
Escherichia_coli_EPEC_Sp=read.csv("EPEC_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Escherichia_coli_EPEC_expr_Sp =as.data.frame(Escherichia_coli_EPEC_Sp[,2:ncol(Escherichia_coli_EPEC_Sp)])
rownames(Escherichia_coli_EPEC_expr_Sp) =Escherichia_coli_EPEC_Sp[,1]

################ETEC_Sp_createdbyCansu###############
Escherichia_coli_ETEC_Sp=read.csv("ETEC_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Sp =as.data.frame(Escherichia_coli_ETEC_Sp[,2:ncol(Escherichia_coli_ETEC_Sp)])
rownames(Escherichia_coli_ETEC_expr_Sp) =Escherichia_coli_ETEC_Sp[,1]

################FRAT_Sp_createdbyCansu###############
Francisella_tularensis_Sp=read.csv("FRAT_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Francisella_tularensis_expr_Sp =as.data.frame(Francisella_tularensis_Sp[,2:ncol(Francisella_tularensis_Sp)])
rownames(Francisella_tularensis_expr_Sp) =Francisella_tularensis_Sp[,1]

################HINF_Sp_createdbyCansu###############
Haemophilus_influenzae_Sp=read.csv("HINF_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Haemophilus_influenzae_expr_Sp =as.data.frame(Haemophilus_influenzae_Sp[,2:ncol(Haemophilus_influenzae_Sp)])
rownames(Haemophilus_influenzae_expr_Sp) =Haemophilus_influenzae_Sp[,1]

################HINF_Sp2_createdbyCansu###############
Haemophilus_influenzae_Sp2=read.csv("HINF_Ns_2Ns_3Sp_1Sp_2Sp_3_UP_Mig_1_DOWN.csv", sep=",")
Haemophilus_influenzae_expr_Sp2 =as.data.frame(Haemophilus_influenzae_Sp2[,2:ncol(Haemophilus_influenzae_Sp2)])
rownames(Haemophilus_influenzae_expr_Sp2) =Haemophilus_influenzae_Sp2[,1]

################HP_G27_Sp_createdbyCansu###############
Helicobacter_pylori_G27_Sp=read.csv("HP_G27_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Helicobacter_pylori_G27_expr_Sp =as.data.frame(Helicobacter_pylori_G27_Sp[,2:ncol(Helicobacter_pylori_G27_Sp)])
rownames(Helicobacter_pylori_G27_expr_Sp) =Helicobacter_pylori_G27_Sp[,1]

################Klebs_Sp_createdbyCansu###############
Klebsiella_pneumoniae_Sp=read.csv("KLEBS_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Sp=as.data.frame(Klebsiella_pneumoniae_Sp[,2:ncol(Klebsiella_pneumoniae_Sp)])
rownames(Klebsiella_pneumoniae_expr_Sp) =Klebsiella_pneumoniae_Sp[,1]

################Listeria_Sp_Vic_up_createdbyCansu###############
Listeria_monocytogenes_Sp_Vic_up=read.csv("Listeria_Sp_1Sp_2Sp_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Listeria_monocytogenes_expr_Sp_Vic_up=as.data.frame(Listeria_monocytogenes_Sp_Vic_up[,2:ncol(Listeria_monocytogenes_Sp_Vic_up)])
rownames(Listeria_monocytogenes_expr_Sp_Vic_up) =Listeria_monocytogenes_Sp_Vic_up[,1]

################MRSA252_Sp_up_createdbyCansu###############
Staphylococcus_aureus_MRSA252_Sp_up=read.csv("MRSA252_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Sp_up=as.data.frame(Staphylococcus_aureus_MRSA252_Sp_up[,2:ncol(Staphylococcus_aureus_MRSA252_Sp_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Sp_up) =Staphylococcus_aureus_MRSA252_Sp_up[,1]

################MSSA476_Sp_Nd_up_createdbyCansu###############
Staphylococcus_aureus_MSSA476_Sp_up_Nd_up=read.csv("MSSA476_Nd_1Nd_2Nd_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Staphylococcus_aureus_MSSA476_expr_Sp_up_Nd_up=as.data.frame(Staphylococcus_aureus_MSSA476_Sp_up_Nd_up[,2:ncol(Staphylococcus_aureus_MSSA476_Sp_up_Nd_up)])
rownames(Staphylococcus_aureus_MSSA476_expr_Sp_up_Nd_up) =Staphylococcus_aureus_MSSA476_Sp_up_Nd_up[,1]

################MSSA476_Sp_up_createdbyCansu###############
Staphylococcus_aureus_MSSA476_Sp_up=read.csv("MSSA476_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Staphylococcus_aureus_MSSA476_expr_Sp_up=as.data.frame(Staphylococcus_aureus_MSSA476_Sp_up[,2:ncol(Staphylococcus_aureus_MSSA476_Sp_up)])
rownames(Staphylococcus_aureus_MSSA476_expr_Sp_up) =Staphylococcus_aureus_MSSA476_Sp_up[,1]


################NGON_Sp1_up_createdbyCansu###############
Neisseria_gonorrhoeae_Sp1_up=read.csv("NGON_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Neisseria_gonorrhoeae_expr_Sp1_up=as.data.frame(Neisseria_gonorrhoeae_Sp1_up[,2:ncol(Neisseria_gonorrhoeae_Sp1_up)])
rownames(Neisseria_gonorrhoeae_expr_Sp1_up) =Neisseria_gonorrhoeae_Sp1_up[,1]

################NGON_Sp2_up_createdbyCansu###############
Neisseria_gonorrhoeae_Sp2_up=read.csv("NGON_Sp_1Sp_2Sp_3_UPyellow.csv", sep=",")
Neisseria_gonorrhoeae_expr_Sp2_up=as.data.frame(Neisseria_gonorrhoeae_Sp2_up[,2:ncol(Neisseria_gonorrhoeae_Sp2_up)])
rownames(Neisseria_gonorrhoeae_expr_Sp2_up) =Neisseria_gonorrhoeae_Sp2_up[,1]

################NMEN_Sp_up_Li_up_createdbyCansu###############
Neisseria_meningitidis_Sp_up_Li_up=read.csv("NMEN_Li_1Li_2Li_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Neisseria_meningitidis_expr_Sp_up_Li_up=as.data.frame(Neisseria_meningitidis_Sp_up_Li_up[,2:ncol(Neisseria_meningitidis_Sp_up_Li_up)])
rownames(Neisseria_meningitidis_expr_Sp_up_Li_up) =Neisseria_meningitidis_Sp_up_Li_up[,1]

################NMEN_Sp_up_createdbyCansu###############
Neisseria_meningitidis_Sp_up=read.csv("NMEN_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Neisseria_meningitidis_expr_Sp_up=as.data.frame(Neisseria_meningitidis_Sp_up[,2:ncol(Neisseria_meningitidis_Sp_up)])
rownames(Neisseria_meningitidis_expr_Sp_up) =Neisseria_meningitidis_Sp_up[,1]

################PSEUDO_Sp_up_Nd_up_createdbyCansu###############
Pseudomonas_aeruginosa_Sp_up_Nd_up=read.csv("PSEUDO_Nd_1Nd_2Nd_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Sp_up_Nd_up=as.data.frame(Pseudomonas_aeruginosa_Sp_up_Nd_up[,2:ncol(Pseudomonas_aeruginosa_Sp_up_Nd_up)])
rownames(Pseudomonas_aeruginosa_expr_Sp_up_Nd_up) =Pseudomonas_aeruginosa_Sp_up_Nd_up[,1]

################PSEUDO_Sp_up_createdbyCansu###############
Pseudomonas_aeruginosa_Sp_up=read.csv("PSEUDO_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Sp_up=as.data.frame(Pseudomonas_aeruginosa_Sp_up[,2:ncol(Pseudomonas_aeruginosa_Sp_up)])
rownames(Pseudomonas_aeruginosa_expr_Sp_up) =Pseudomonas_aeruginosa_Sp_up[,1]

################SALMT_Sp1_up_createdbyCansu###############
Salmonella_enterica_Sp1_up=read.csv("SALMT_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Salmonella_enterica_expr_Sp1_up=as.data.frame(Salmonella_enterica_Sp1_up[,2:ncol(Salmonella_enterica_Sp1_up)])
rownames(Salmonella_enterica_expr_Sp1_up) =Salmonella_enterica_Sp1_up[,1]

################SALMTSp2_up_createdbyCansu###############
Salmonella_enterica_Sp2_up=read.csv("SALMT_Sp_1Sp_2Sp_3Vic_3_UP.csv", sep=",")
Salmonella_enterica_expr_Sp2_up=as.data.frame(Salmonella_enterica_Sp2_up[,2:ncol(Salmonella_enterica_Sp2_up)])
rownames(Salmonella_enterica_expr_Sp2_up) =Salmonella_enterica_Sp2_up[,1]

################SEPI_Sp_up_Tm_up_createdbyCansu###############
Staphylococcus_epidermidis_Sp_up_Tm_up=read.csv("SEPI_Sp_1Sp_2Sp_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Staphylococcus_epidermidis_expr_Sp_up_Tm_up=as.data.frame(Staphylococcus_epidermidis_Sp_up_Tm_up[,2:ncol(Staphylococcus_epidermidis_Sp_up_Tm_up)])
rownames(Staphylococcus_epidermidis_expr_Sp_up_Tm_up) =Staphylococcus_epidermidis_Sp_up_Tm_up[,1]

################SEPI_Sp_up_createdbyCansu###############
Staphylococcus_epidermidis_Sp_up=read.csv("SEPI_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Staphylococcus_epidermidis_expr_Sp_up=as.data.frame(Staphylococcus_epidermidis_Sp_up[,2:ncol(Staphylococcus_epidermidis_Sp_up)])
rownames(Staphylococcus_epidermidis_expr_Sp_up) =Staphylococcus_epidermidis_Sp_up[,1]

################SHIF_Sp_up_createdbyCansu###############
Shigella_flexneri_Sp_up=read.csv("SHIF_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Shigella_flexneri_expr_Sp_up=as.data.frame(Shigella_flexneri_Sp_up[,2:ncol(Shigella_flexneri_Sp_up)])
rownames(Shigella_flexneri_expr_Sp_up) =Shigella_flexneri_Sp_up[,1]

################STRPN_Sp1_up_createdbyCansu###############
Streptococcus_pneumoniae_Sp1_up=read.csv("STRPN_Sp_1Sp_2Sp_3Tm_1Tm_3_UP.csv", sep=",")
Streptococcus_pneumoniae_expr_Sp1_up=as.data.frame(Streptococcus_pneumoniae_Sp1_up[,2:ncol(Streptococcus_pneumoniae_Sp1_up)])
rownames(Streptococcus_pneumoniae_expr_Sp1_up) =Streptococcus_pneumoniae_Sp1_up[,1]

################STRPN_Sp2_up_createdbyCansu###############
Streptococcus_pneumoniae_Sp2_up=read.csv("STRPN_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Streptococcus_pneumoniae_expr_Sp2_up=as.data.frame(Streptococcus_pneumoniae_Sp2_up[,2:ncol(Streptococcus_pneumoniae_Sp2_up)])
rownames(Streptococcus_pneumoniae_expr_Sp2_up) =Streptococcus_pneumoniae_Sp2_up[,1]

################UPEC_Sp_up_createdbyCansu###############
Escherichia_coli_UPEC_Sp_up=read.csv("UPEC_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Escherichia_coli_UPEC_expr_Sp_up=as.data.frame(Escherichia_coli_UPEC_Sp_up[,2:ncol(Escherichia_coli_UPEC_Sp_up)])
rownames(Escherichia_coli_UPEC_expr_Sp_up) =Escherichia_coli_UPEC_Sp_up[,1]

################Vibrio_Sp_up_Li_up_createdbyCansu###############
Vibrio_cholerae_Sp_up_Li_up=read.csv("Vibrio_Li_1Li_2Li_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Sp_up_Li_up=as.data.frame(Vibrio_cholerae_Sp_up_Li_up[,2:ncol(Vibrio_cholerae_Sp_up_Li_up)])
rownames(Vibrio_cholerae_expr_Sp_up_Li_up) =Vibrio_cholerae_Sp_up_Li_up[,1]

################Vibrio_Sp_up_createdbyCansu###############
Vibrio_cholerae_Sp_up=read.csv("Vibrio_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Sp_up=as.data.frame(Vibrio_cholerae_Sp_up[,2:ncol(Vibrio_cholerae_Sp_up)])
rownames(Vibrio_cholerae_expr_Sp_up) =Vibrio_cholerae_Sp_up[,1]


################Ypstb_Sp1_up_createdbyCansu###############
Yersinia_pseudotuberculosis_Sp1_up=read.csv("YPSTB_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Sp1_up=as.data.frame(Yersinia_pseudotuberculosis_Sp1_up[,2:ncol(Yersinia_pseudotuberculosis_Sp1_up)])
rownames(Yersinia_pseudotuberculosis_expr_Sp1_up) =Yersinia_pseudotuberculosis_Sp1_up[,1]

################Ypstb_Sp2_up_createdbyCansu###############
Yersinia_pseudotuberculosis_Sp2_up=read.csv("YPSTB_Sp_1Sp_2Sp_3_UPpurple.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Sp2_up=as.data.frame(Yersinia_pseudotuberculosis_Sp2_up[,2:ncol(Yersinia_pseudotuberculosis_Sp2_up)])
rownames(Yersinia_pseudotuberculosis_expr_Sp2_up) =Yersinia_pseudotuberculosis_Sp2_up[,1]

################Ypstb_Sp_up_Mig_down_createdbyCansu###############
Yersinia_pseudotuberculosis_Sp_up_Mig_down=read.csv("YPSTB_Sp_1Sp_2Sp_3Vic_1Vic_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Sp_up_Mig_down=as.data.frame(Yersinia_pseudotuberculosis_Sp_up_Mig_down[,2:ncol(Yersinia_pseudotuberculosis_Sp_up_Mig_down)])
rownames(Yersinia_pseudotuberculosis_expr_Sp_up_Mig_down) =Yersinia_pseudotuberculosis_Sp_up_Mig_down[,1]

################AGGA_Sp_up_createdbyCansu###############
Aggregatibacter_actinomycetemcomitans_Sp_up=read.csv("AGGA_Mig_2Mig_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Sp_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Sp_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Sp_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Sp_up) =Aggregatibacter_actinomycetemcomitans_Sp_up[,1]

########BURK_Sp_up_Mig_up_createdbyCansu###############
Burkholderia_pseudomallei_Sp_up_Mig_up=read.csv("BURK_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Sp_up_Mig_up =as.data.frame(Burkholderia_pseudomallei_Sp_up_Mig_up[,2:ncol(Burkholderia_pseudomallei_Sp_up_Mig_up)])
rownames(Burkholderia_pseudomallei_expr_Sp_up_Mig_up) = Burkholderia_pseudomallei_Sp_up_Mig_up[,1]

########Entfa_Sp_up_Mig_up_createdbyCansu###############
Enterococcus_faecalis_Sp_up_Mig_up=read.csv("ENTFA_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Enterococcus_faecalis_expr_Sp_up_Mig_up =as.data.frame(Enterococcus_faecalis_Sp_up_Mig_up[,2:ncol(Enterococcus_faecalis_Sp_up_Mig_up)])
rownames(Enterococcus_faecalis_expr_Sp_up_Mig_up) =Enterococcus_faecalis_Sp_up_Mig_up[,1]

########Spyo_Sp_up_Mig_up_createdbyCansu###############
Streptococcus_pyogenes_Sp_up_Mig_up=read.csv("SPYO_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Streptococcus_pyogenes_expr_Sp_up_Mig_up =as.data.frame(Streptococcus_pyogenes_Sp_up_Mig_up[,2:ncol(Streptococcus_pyogenes_Sp_up_Mig_up)])
rownames(Streptococcus_pyogenes_expr_Sp_up_Mig_up) =Streptococcus_pyogenes_Sp_up_Mig_up[,1]



########STAGA_Mig_up_Sp_up_createdbyCansu###############
Streptococcus_agalactiae_Sp_up_Mig_up=read.csv("STAGA_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Streptococcus_agalactiae_expr_Sp_up_Mig_up =as.data.frame(Streptococcus_agalactiae_Sp_up_Mig_up[,2:ncol(Streptococcus_agalactiae_Sp_up_Mig_up)])
rownames(Streptococcus_agalactiae_expr_Sp_up_Mig_up) =Streptococcus_agalactiae_Sp_up_Mig_up[,1]



########STAGA_Sp_up_createdbyCansu###############
Streptococcus_agalactiae_Sp_up=read.csv("STAGA_Mig_1Mig_3Sp_1Sp_2Sp_3Vic_1_UP.csv", sep=",")
Streptococcus_agalactiae_expr_Sp_up =as.data.frame(Streptococcus_agalactiae_Sp_up[,2:ncol(Streptococcus_agalactiae_Sp_up)])
rownames(Streptococcus_agalactiae_expr_Sp_up) =Streptococcus_agalactiae_Sp_up[,1]

########UPEC_Sp_As###############
Escherichia_coli_UPEC_536_As = read.csv("UPEC_As_1As_2As_3Oss_1_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Escherichia_coli_UPEC_536_expr_As =as.data.frame(Escherichia_coli_UPEC_536_As[,2:ncol(Escherichia_coli_UPEC_536_As)])
rownames(Escherichia_coli_UPEC_536_expr_As) = Escherichia_coli_UPEC_536_As[,1]

###########KLEBS_Sp_As############
Klebsiella_pneumoniae_As_Sp =read.csv("KLEBS_As_1As_2As_3_UP_Sp_1Sp_2Sp_3_DOWN.csv",  sep=",")
Klebsiella_pneumoniae_expr_As_Sp =as.data.frame(Klebsiella_pneumoniae_As_Sp[,2:ncol(Klebsiella_pneumoniae_As_Sp)])
rownames(Klebsiella_pneumoniae_expr_As_Sp) = Klebsiella_pneumoniae_As_Sp[,1]

##################KLEBS_SP_Bs##########################
Klebsiella_pneumoniae_Bs_Sp= read.csv("KLEBS_Bs_1Bs_2Bs_3Ctrl_1Ctrl_2Ctrl_3_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Klebsiella_pneumoniae_expr_Bs_Sp=as.data.frame(Klebsiella_pneumoniae_Bs_Sp[,2:ncol(Klebsiella_pneumoniae_Bs_Sp)])
rownames(Klebsiella_pneumoniae_expr_Bs_Sp) = Klebsiella_pneumoniae_Bs_Sp[,1]

##################Pseudo_Sp1_down_Ns_up##########################
Pseudomonas_aeruginosa_Ns1_up_Sp_down=read.csv("PSEUDO_Li_2Ns_1Ns_2Ns_3_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_Ns1_up_Sp_down=as.data.frame(Pseudomonas_aeruginosa_Ns1_up_Sp_down[,2:ncol(Pseudomonas_aeruginosa_Ns1_up_Sp_down)])
rownames(Pseudomonas_aeruginosa_expr_Ns1_up_Sp_down) =Pseudomonas_aeruginosa_Ns1_up_Sp_down[,1]

##################Pseudo_Sp2_down_Ns_up##########################
Pseudomonas_aeruginosa_Ns2_up_Sp_down=read.csv("PSEUDO_Ns_1Ns_2Ns_3_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_Ns2_up_Sp_down=as.data.frame(Pseudomonas_aeruginosa_Ns2_up_Sp_down[,2:ncol(Pseudomonas_aeruginosa_Ns2_up_Sp_down)])
rownames(Pseudomonas_aeruginosa_expr_Ns2_up_Sp_down) =Pseudomonas_aeruginosa_Ns2_up_Sp_down[,1]

##################KLEBS_Sp_down##########################
Klebsiella_pneumoniae_Sp_down=read.csv("KLEBS_Li_1Li_2Nd_2Nd_3_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Klebsiella_pneumoniae_expr_Sp_down=as.data.frame(Klebsiella_pneumoniae_Sp_down[,2:ncol(Klebsiella_pneumoniae_Sp_down)])
rownames(Klebsiella_pneumoniae_expr_Sp_down) =Klebsiella_pneumoniae_Sp_down[,1]

################PSEUDO_Sp1_down_createdbyCansu###############
Pseudomonas_aeruginosa_Sp_down=read.csv("PSEUDO_As_2As_3Oss_1Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_Sp_down=as.data.frame(Pseudomonas_aeruginosa_Sp_down[,2:ncol(Pseudomonas_aeruginosa_Sp_down)])
rownames(Pseudomonas_aeruginosa_expr_Sp_down) =Pseudomonas_aeruginosa_Sp_down[,1]

################PSEUDO_Sp2_down_createdbyCansu###############
Pseudomonas_aeruginosa_Sp2_down=read.csv("PSEUDO_Li_2_UP_Oss_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_Sp2_down=as.data.frame(Pseudomonas_aeruginosa_Sp2_down[,2:ncol(Pseudomonas_aeruginosa_Sp2_down)])
rownames(Pseudomonas_aeruginosa_expr_Sp2_down) =Pseudomonas_aeruginosa_Sp2_down[,1]




########SPYO_Sp_down_Li_up_Mig_down###############
Streptococcus_pyogenes_Mig_Sp_Li = read.csv("SPYO_Li_1Li_2Li_3_UP_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Streptococcus_pyogenes_expr_Mig_Sp_Li =as.data.frame(Streptococcus_pyogenes_Mig_Sp_Li[,2:ncol(Streptococcus_pyogenes_Mig_Sp_Li)])
rownames(Streptococcus_pyogenes_expr_Mig_Sp_Li) = Streptococcus_pyogenes_Mig_Sp_Li[,1]

#################UPEC_Sp_Mig_down##########################
Escherichia_coli_UPEC_Mig_Sp_down=read.csv("UPEC_Li_1Li_2Ns_2_UP_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Escherichia_coli_UPEC_expr_Mig_Sp_down=as.data.frame(Escherichia_coli_UPEC_Mig_Sp_down[,2:ncol(Escherichia_coli_UPEC_Mig_Sp_down)])
rownames(Escherichia_coli_UPEC_expr_Mig_Sp_down) =Escherichia_coli_UPEC_Mig_Sp_down[,1]

##################SHIF_Sp_Mig_down##########################
Shigella_flexneri_Mig_Sp_down=read.csv("SHIF_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Shigella_flexneri_expr_Mig_Sp_down=as.data.frame(Shigella_flexneri_Mig_Sp_down[,2:ncol(Shigella_flexneri_Mig_Sp_down)])
rownames(Shigella_flexneri_expr_Mig_Sp_down) =Shigella_flexneri_Mig_Sp_down[,1]

##########ND_STRESS################
###########ENTFA_Nd_down############
Enterococcus_faecalis_Nd =read.csv("ENTFA_Nd_1Nd_2Nd_3_DOWN.csv",  sep=",")
Enterococcus_faecalis_expr_Nd =as.data.frame(Enterococcus_faecalis_Nd[,2:ncol(Enterococcus_faecalis_Nd)])
rownames(Enterococcus_faecalis_expr_Nd) = Enterococcus_faecalis_Nd[,1]



###########MSSA476_Nd_Mig############
Staphylococcus_aureus_MSSA476_Nd_Mig =read.csv("MSSA476_Ns_3_UP_Mig_1Mig_2Mig_3Nd_1Nd_2Nd_3Sp_3_DOWN.csv",  sep=",")
Staphylococcus_aureus_MSSA476_expr_Nd_Mig =as.data.frame(Staphylococcus_aureus_MSSA476_Nd_Mig[,2:ncol(Staphylococcus_aureus_MSSA476_Nd_Mig)])
rownames(Staphylococcus_aureus_MSSA476_expr_Nd_Mig) = Staphylococcus_aureus_MSSA476_Nd_Mig[,1]

###########MRSA252_Nd_Mig############
Staphylococcus_aureus_MRSA252_Nd_Mig =read.csv("MRSA252_Mig_1Mig_2Mig_3Nd_1Nd_2Nd_3_DOWN.csv",  sep=",")
Staphylococcus_aureus_MRSA252_expr_Nd_Mig =as.data.frame(Staphylococcus_aureus_MRSA252_Nd_Mig[,2:ncol(Staphylococcus_aureus_MRSA252_Nd_Mig)])
rownames(Staphylococcus_aureus_MRSA252_expr_Nd_Mig) = Staphylococcus_aureus_MRSA252_Nd_Mig[,1]

###########STAGA_Nd1_down############
Streptococcus_agalactiae_Nd =read.csv("STAGA_Mig_1Nd_1Nd_2Nd_3_DOWN.csv",  sep=",")
Streptococcus_agalactiae_expr_Nd =as.data.frame(Streptococcus_agalactiae_Nd[,2:ncol(Streptococcus_agalactiae_Nd)])
rownames(Streptococcus_agalactiae_expr_Nd) = Streptococcus_agalactiae_Nd[,1]

#################STAGA_Nd2_down##########################
Streptococcus_agalactiae_Nd2_down=read.csv("STAGA_Tm_1Vic_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Nd2_down=as.data.frame(Streptococcus_agalactiae_Nd2_down[,2:ncol(Streptococcus_agalactiae_Nd2_down)])
rownames(Streptococcus_agalactiae_expr_Nd2_down) =Streptococcus_agalactiae_Nd2_down[,1]

###########HP_G27_Nd1_don############
Helicobacter_pylori_G27_Nd =read.csv("HP_G27_As_1As_2_UP_Nd_1Nd_2Nd_3_DOWN.csv",  sep=",")
Helicobacter_pylori_G27_expr_Nd =as.data.frame(Helicobacter_pylori_G27_Nd[,2:ncol(Helicobacter_pylori_G27_Nd)])
rownames(Helicobacter_pylori_G27_expr_Nd) = Helicobacter_pylori_G27_Nd[,1]

###########HPG27_Nd2_down_createdbyCansu############
Helicobacter_pylori_G27_Nd2_down=read.csv("HP_G27_Sp_1_UP_Nd_1Nd_2Nd_3_DOWN.csv",  sep=",")
Helicobacter_pylori_G27_expr_Nd2_down=as.data.frame(Helicobacter_pylori_G27_Nd2_down[,2:ncol(Helicobacter_pylori_G27_Nd2_down)])
rownames(Helicobacter_pylori_G27_expr_Nd2_down) =Helicobacter_pylori_G27_Nd2_down[,1]


###########HP_J99_Nd############
Helicobacter_pylori_J99_Nd =read.csv("HP_J99_Ctrl_2_UP_Nd_1Nd_2Nd_3Vic_1Vic_2_DOWN.csv",  sep=",")
Helicobacter_pylori_J99_expr_Nd =as.data.frame(Helicobacter_pylori_J99_Nd[,2:ncol(Helicobacter_pylori_J99_Nd)])
rownames(Helicobacter_pylori_J99_expr_Nd) = Helicobacter_pylori_J99_Nd[,1]

###########ACHX_Nd_up_createdbyCansu############
Achromobacter_xylosoxidans_Nd_up=read.csv("ACHX_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Achromobacter_xylosoxidans_expr_Nd_up =as.data.frame(Achromobacter_xylosoxidans_Nd_up[,2:ncol(Achromobacter_xylosoxidans_Nd_up)])
rownames(Achromobacter_xylosoxidans_expr_Nd_up) = Achromobacter_xylosoxidans_Nd_up[,1]

###########ACIB_Nd_up_createdbyCansu############
Acinetobacter_baumannii_Nd_up=read.csv("ACIB_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Acinetobacter_baumannii_expr_Nd_up =as.data.frame(Acinetobacter_baumannii_Nd_up[,2:ncol(Acinetobacter_baumannii_Nd_up)])
rownames(Acinetobacter_baumannii_expr_Nd_up) = Acinetobacter_baumannii_Nd_up[,1]

###########AGGA_Nd1_up_createdbyCansu############
Aggregatibacter_actinomycetemcomitans_Nd_up=read.csv("AGGA_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Nd_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Nd_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Nd_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Nd_up) =Aggregatibacter_actinomycetemcomitans_Nd_up[,1]

###########AGGA_Nd2_up_createdbyCansu############
Aggregatibacter_actinomycetemcomitans_Nd2_up=read.csv("AGGA_Nd_1Nd_2Nd_3_UP_Ns_2_DOWN.csv",  sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Nd2_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Nd2_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Nd2_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Nd2_up) =Aggregatibacter_actinomycetemcomitans_Nd2_up[,1]

###########BBURG_Nd_up_createdbyCansu############
Borrelia_burgdorferi_Nd_up=read.csv("BBURG_Nd_1Nd_2Nd_3_UP_Sp_1Sp_3_DOWN.csv",  sep=",")
Borrelia_burgdorferi_expr_Nd_up=as.data.frame(Borrelia_burgdorferi_Nd_up[,2:ncol(Borrelia_burgdorferi_Nd_up)])
rownames(Borrelia_burgdorferi_expr_Nd_up) =Borrelia_burgdorferi_Nd_up[,1]

###########BURK_Nd1_up_createdbyCansu############
Burkholderia_pseudomallei_Nd1_up=read.csv("BURK_Nd_1Nd_2Nd_3Tm_2Tm_3_UP.csv",  sep=",")
Burkholderia_pseudomallei_expr_Nd1_up=as.data.frame(Burkholderia_pseudomallei_Nd1_up[,2:ncol(Burkholderia_pseudomallei_Nd1_up)])
rownames(Burkholderia_pseudomallei_expr_Nd1_up) =Burkholderia_pseudomallei_Nd1_up[,1]

###########BURK_Nd2_up_createdbyCansu############
Burkholderia_pseudomallei_Nd2_up=read.csv("BURK_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Burkholderia_pseudomallei_expr_Nd2_up=as.data.frame(Burkholderia_pseudomallei_Nd2_up[,2:ncol(Burkholderia_pseudomallei_Nd2_up)])
rownames(Burkholderia_pseudomallei_expr_Nd2_up) =Burkholderia_pseudomallei_Nd2_up[,1]


###########BURK_Nd3_up_createdbyCansu############
Burkholderia_pseudomallei_Nd3_up=read.csv("BURK_Nd_1Nd_2Nd_3_UPsalmon.csv",  sep=",")
Burkholderia_pseudomallei_expr_Nd3_up=as.data.frame(Burkholderia_pseudomallei_Nd3_up[,2:ncol(Burkholderia_pseudomallei_Nd3_up)])
rownames(Burkholderia_pseudomallei_expr_Nd3_up) =Burkholderia_pseudomallei_Nd3_up[,1]

###########Campy_Nd_up_createdbyCansu############
Campylobacter_jejuni_Nd_up=read.csv("Campy_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Campylobacter_jejuni_expr_Nd_up=as.data.frame(Campylobacter_jejuni_Nd_up[,2:ncol(Campylobacter_jejuni_Nd_up)])
rownames(Campylobacter_jejuni_expr_Nd_up) =Campylobacter_jejuni_Nd_up[,1]

###########ENFTA_Nd1_up_createdbyCansu############
Enterococcus_faecalis_Nd1_up=read.csv("ENTFA_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Enterococcus_faecalis_expr_Nd1_up=as.data.frame(Enterococcus_faecalis_Nd1_up[,2:ncol(Enterococcus_faecalis_Nd1_up)])
rownames(Enterococcus_faecalis_expr_Nd1_up) =Enterococcus_faecalis_Nd1_up[,1]

###########ENFTA_Nd2_up_createdbyCansu############
Enterococcus_faecalis_Nd2_up=read.csv("ENTFA_Nd_1Nd_2Nd_3_UP_Mig_2Mig_3_DOWN.csv",  sep=",")
Enterococcus_faecalis_expr_Nd2_up=as.data.frame(Enterococcus_faecalis_Nd2_up[,2:ncol(Enterococcus_faecalis_Nd2_up)])
rownames(Enterococcus_faecalis_expr_Nd2_up) =Enterococcus_faecalis_Nd2_up[,1]


###########EPEC_Nd_up_Oss_up_createdbyCansu############
Escherichia_coli_EPEC_Nd_up_Oss_up=read.csv("EPEC_Nd_1Nd_2Nd_3Oss_1Oss_2Oss_3_UP.csv",  sep=",")
Escherichia_coli_EPEC_expr_Nd_up_Oss_up=as.data.frame(Escherichia_coli_EPEC_Nd_up_Oss_up[,2:ncol(Escherichia_coli_EPEC_Nd_up_Oss_up)])
rownames(Escherichia_coli_EPEC_expr_Nd_up_Oss_up) =Escherichia_coli_EPEC_Nd_up_Oss_up[,1]


###########EPEC_Nd_up_createdbyCansu############
Escherichia_coli_EPEC_Nd_up=read.csv("EPEC_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Escherichia_coli_EPEC_expr_Nd_up=as.data.frame(Escherichia_coli_EPEC_Nd_up[,2:ncol(Escherichia_coli_EPEC_Nd_up)])
rownames(Escherichia_coli_EPEC_expr_Nd_up) =Escherichia_coli_EPEC_Nd_up[,1]

###########ETEC_Nd1_up_createdbyCansu############
Escherichia_coli_ETEC_Nd1_up=read.csv("ETEC_As_2Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Escherichia_coli_ETEC_expr_Nd1_up=as.data.frame(Escherichia_coli_ETEC_Nd1_up[,2:ncol(Escherichia_coli_ETEC_Nd1_up)])
rownames(Escherichia_coli_ETEC_expr_Nd1_up) =Escherichia_coli_ETEC_Nd1_up[,1]


###########ETEC_Nd2_up_createdbyCansu############
Escherichia_coli_ETEC_Nd2_up=read.csv("ETEC_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Escherichia_coli_ETEC_expr_Nd2_up=as.data.frame(Escherichia_coli_ETEC_Nd2_up[,2:ncol(Escherichia_coli_ETEC_Nd2_up)])
rownames(Escherichia_coli_ETEC_expr_Nd2_up) =Escherichia_coli_ETEC_Nd2_up[,1]

###########ETEC_Nd_up_Mig_upcreatedbyCansu############
Escherichia_coli_ETEC_Nd_up_Mig_up=read.csv("ETEC_Mig_1Mig_2Mig_3Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Escherichia_coli_ETEC_expr_Nd_up_Mig_up=as.data.frame(Escherichia_coli_ETEC_Nd_up_Mig_up[,2:ncol(Escherichia_coli_ETEC_Nd_up_Mig_up)])
rownames(Escherichia_coli_ETEC_expr_Nd_up_Mig_up) =Escherichia_coli_ETEC_Nd_up_Mig_up[,1]


###########FRAT_Nd1_up_createdbyCansu############
Francisella_tularensis_Nd1_up=read.csv("FRAT_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Francisella_tularensis_expr_Nd1_up=as.data.frame(Francisella_tularensis_Nd1_up[,2:ncol(Francisella_tularensis_Nd1_up)])
rownames(Francisella_tularensis_expr_Nd1_up) =Francisella_tularensis_Nd1_up[,1]

###########FRAT_Nd2_up_createdbyCansu############
Francisella_tularensis_Nd2_up=read.csv("FRAT_Nd_1Nd_2Nd_3_UP_Li_1Li_3_DOWN.csv",  sep=",")
Francisella_tularensis_expr_Nd2_up=as.data.frame(Francisella_tularensis_Nd2_up[,2:ncol(Francisella_tularensis_Nd2_up)])
rownames(Francisella_tularensis_expr_Nd2_up) =Francisella_tularensis_Nd2_up[,1]

###########FRAT_Nd1_down############
Francisella_tularensis_Nd =read.csv("FRAT_Nd_1Nd_2Nd_3_DOWN.csv",  sep=",")
Francisella_tularensis_expr_Nd =as.data.frame(Francisella_tularensis_Nd[,2:ncol(Francisella_tularensis_Nd)])
rownames(Francisella_tularensis_expr_Nd) = Francisella_tularensis_Nd[,1]

###########FRAT_Nd2_down_createdbyCansu############
Francisella_tularensis_Nd2_down=read.csv("FRAT_Li_1Li_2Nd_1Nd_2Nd_3_DOWN.csv",  sep=",")
Francisella_tularensis_expr_Nd2_down=as.data.frame(Francisella_tularensis_Nd2_down[,2:ncol(Francisella_tularensis_Nd2_down)])
rownames(Francisella_tularensis_expr_Nd2_down) =Francisella_tularensis_Nd2_down[,1]

###########HINF_Nd_up_createdbyCansu############
Haemophilus_influenzae_Nd_up=read.csv("HINF_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Haemophilus_influenzae_expr_Nd_up=as.data.frame(Haemophilus_influenzae_Nd_up[,2:ncol(Haemophilus_influenzae_Nd_up)])
rownames(Haemophilus_influenzae_expr_Nd_up) =Haemophilus_influenzae_Nd_up[,1]

###########HPG27_Nd_up_createdbyCansu############
Helicobacter_pylori_G27_Nd_up=read.csv("HP_G27_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Helicobacter_pylori_G27_expr_Nd_up=as.data.frame(Helicobacter_pylori_G27_Nd_up[,2:ncol(Helicobacter_pylori_G27_Nd_up)])
rownames(Helicobacter_pylori_G27_expr_Nd_up) =Helicobacter_pylori_G27_Nd_up[,1]


###########KLEBS_Nd1_up_createdbyCansu############
Klebsiella_pneumoniae_Nd_up=read.csv("KLEBS_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Klebsiella_pneumoniae_expr_Nd_up=as.data.frame(Klebsiella_pneumoniae_Nd_up[,2:ncol(Klebsiella_pneumoniae_Nd_up)])
rownames(Klebsiella_pneumoniae_expr_Nd_up) =Klebsiella_pneumoniae_Nd_up[,1]

###########KLEBS_Nd2_up_createdbyCansu############
Klebsiella_pneumoniae_Nd2_up=read.csv("KLEBS_Li_2Nd_1Nd_2Nd_3Tm_1Tm_3_UP.csv",  sep=",")
Klebsiella_pneumoniae_expr_Nd2_up=as.data.frame(Klebsiella_pneumoniae_Nd2_up[,2:ncol(Klebsiella_pneumoniae_Nd2_up)])
rownames(Klebsiella_pneumoniae_expr_Nd2_up) =Klebsiella_pneumoniae_Nd2_up[,1]

###########Listeria_Nd1_up_createdbyCansu############
Listeria_monocytogenes_Nd1_up=read.csv("Listeria_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Listeria_monocytogenes_expr_Nd1_up=as.data.frame(Listeria_monocytogenes_Nd1_up[,2:ncol(Listeria_monocytogenes_Nd1_up)])
rownames(Listeria_monocytogenes_expr_Nd1_up) =Listeria_monocytogenes_Nd1_up[,1]

###########Listeria_Nd2_up_createdbyCansu############
Listeria_monocytogenes_Nd2_up=read.csv("Listeria_Nd_1Nd_2Nd_3Oss_2_UP.csv",  sep=",")
Listeria_monocytogenes_expr_Nd2_up=as.data.frame(Listeria_monocytogenes_Nd2_up[,2:ncol(Listeria_monocytogenes_Nd2_up)])
rownames(Listeria_monocytogenes_expr_Nd2_up) =Listeria_monocytogenes_Nd2_up[,1]

###########Listeria_Nd3_up_createdbyCansu############
Listeria_monocytogenes_Nd3_up=read.csv("Listeria_Nd_1Nd_2Nd_3Sp_3Vic_1Vic_2_UP.csv",  sep=",")
Listeria_monocytogenes_expr_Nd3_up=as.data.frame(Listeria_monocytogenes_Nd3_up[,2:ncol(Listeria_monocytogenes_Nd3_up)])
rownames(Listeria_monocytogenes_expr_Nd3_up) =Listeria_monocytogenes_Nd3_up[,1]



###########MRSA252_Nd_up_Ns_down_createdbyCansu############
Staphylococcus_aureus_MRSA252_Nd_up_Ns_down=read.csv("MRSA252_As_2As_3Nd_1Nd_2Nd_3_UP_Ns_1Ns_2Ns_3_DOWN.csv",  sep=",")
Staphylococcus_aureus_MRSA252_expr_Nd_up_Ns_down=as.data.frame(Staphylococcus_aureus_MRSA252_Nd_up_Ns_down[,2:ncol(Staphylococcus_aureus_MRSA252_Nd_up_Ns_down)])
rownames(Staphylococcus_aureus_MRSA252_expr_Nd_up_Ns_down) =Staphylococcus_aureus_MRSA252_Nd_up_Ns_down[,1]

###########MRSA252_Nd_up_createdbyCansu############
Staphylococcus_aureus_MRSA252_Nd_up=read.csv("MRSA252_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Staphylococcus_aureus_MRSA252_expr_Nd_up=as.data.frame(Staphylococcus_aureus_MRSA252_Nd_up[,2:ncol(Staphylococcus_aureus_MRSA252_Nd_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Nd_up) =Staphylococcus_aureus_MRSA252_Nd_up[,1]

###########MSSA476_Nd_up_createdbyCansu############
Staphylococcus_aureus_MSSA476_Nd_up=read.csv("MSSA476_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Staphylococcus_aureus_MSSA476_expr_Nd_up=as.data.frame(Staphylococcus_aureus_MSSA476_Nd_up[,2:ncol(Staphylococcus_aureus_MSSA476_Nd_up)])
rownames(Staphylococcus_aureus_MSSA476_expr_Nd_up) =Staphylococcus_aureus_MSSA476_Nd_up[,1]

###########NMEN_Nd_up_createdbyCansu############
Neisseria_meningitidis_Nd_up=read.csv("NMEN_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Neisseria_meningitidis_expr_Nd_up=as.data.frame(Neisseria_meningitidis_Nd_up[,2:ncol(Neisseria_meningitidis_Nd_up)])
rownames(Neisseria_meningitidis_expr_Nd_up) =Neisseria_meningitidis_Nd_up[,1]

###########Pseudo_Nd_up_Ns_up_createdbyCansu############
Pseudomonas_aeruginosa_Nd_up_Ns_up=read.csv("PSEUDO_Nd_1Nd_2Nd_3Ns_1Ns_2Ns_3_UP.csv",  sep=",")
Pseudomonas_aeruginosa_expr_Nd_up_Ns_up=as.data.frame(Pseudomonas_aeruginosa_Nd_up_Ns_up[,2:ncol(Pseudomonas_aeruginosa_Nd_up_Ns_up)])
rownames(Pseudomonas_aeruginosa_expr_Nd_up_Ns_up) =Pseudomonas_aeruginosa_Nd_up_Ns_up[,1]

###########Pseudo_Nd_up_Oss_up_createdbyCansu############
Pseudomonas_aeruginosa_Nd_up_Oss_up=read.csv("PSEUDO_Nd_1Nd_2Nd_3Oss_1Oss_2Oss_3_UP.csv",  sep=",")
Pseudomonas_aeruginosa_expr_Nd_up_Oss_up=as.data.frame(Pseudomonas_aeruginosa_Nd_up_Oss_up[,2:ncol(Pseudomonas_aeruginosa_Nd_up_Oss_up)])
rownames(Pseudomonas_aeruginosa_expr_Nd_up_Oss_up) =Pseudomonas_aeruginosa_Nd_up_Oss_up[,1]

###########Pseudo_Nd_up_createdbyCansu############
Pseudomonas_aeruginosa_Nd_up=read.csv("PSEUDO_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Pseudomonas_aeruginosa_expr_Nd_up=as.data.frame(Pseudomonas_aeruginosa_Nd_up[,2:ncol(Pseudomonas_aeruginosa_Nd_up)])
rownames(Pseudomonas_aeruginosa_expr_Nd_up) =Pseudomonas_aeruginosa_Nd_up[,1]

###########SALMT_Nd_up_createdbyCansu############
Salmonella_enterica_Nd_up=read.csv("SALMT_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Salmonella_enterica_expr_Nd_up=as.data.frame(Salmonella_enterica_Nd_up[,2:ncol(Salmonella_enterica_Nd_up)])
rownames(Salmonella_enterica_expr_Nd_up) =Salmonella_enterica_Nd_up[,1]

###########SEPI_Nd1_up_createdbyCansu############
Staphylococcus_epidermidis_Nd1_up=read.csv("SEPI_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Staphylococcus_epidermidis_expr_Nd1_up=as.data.frame(Staphylococcus_epidermidis_Nd1_up[,2:ncol(Staphylococcus_epidermidis_Nd1_up)])
rownames(Staphylococcus_epidermidis_expr_Nd1_up) =Staphylococcus_epidermidis_Nd1_up[,1]

###########SEPI_Nd2_up_createdbyCansu############
Staphylococcus_epidermidis_Nd2_up=read.csv("SEPI_Nd_1Nd_2Nd_3_UP_Sp_1Sp_3_DOWN.csv",  sep=",")
Staphylococcus_epidermidis_expr_Nd2_up=as.data.frame(Staphylococcus_epidermidis_Nd2_up[,2:ncol(Staphylococcus_epidermidis_Nd2_up)])
rownames(Staphylococcus_epidermidis_expr_Nd2_up) =Staphylococcus_epidermidis_Nd2_up[,1]

###########SEPI_Nd_up_Tm_down_createdbyCansu############
Staphylococcus_epidermidis_Nd_up_Tm_down=read.csv("SEPI_Nd_1Nd_2Nd_3_UP_Tm_1Tm_2Tm_3_DOWN.csv",  sep=",")
Staphylococcus_epidermidis_expr_Nd_up_Tm_down=as.data.frame(Staphylococcus_epidermidis_Nd_up_Tm_down[,2:ncol(Staphylococcus_epidermidis_Nd_up_Tm_down)])
rownames(Staphylococcus_epidermidis_expr_Nd_up_Tm_down) =Staphylococcus_epidermidis_Nd_up_Tm_down[,1]

###########SHIF_Nd_up_createdbyCansu############
Shigella_flexneri_Nd_up=read.csv("SHIF_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Shigella_flexneri_expr_Nd_up=as.data.frame(Shigella_flexneri_Nd_up[,2:ncol(Shigella_flexneri_Nd_up)])
rownames(Shigella_flexneri_expr_Nd_up) =Shigella_flexneri_Nd_up[,1]

###########SPYO_Nd_up_createdbyCansu############
Streptococcus_pyogenes_Nd_up=read.csv("SPYO_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Streptococcus_pyogenes_expr_Nd_up=as.data.frame(Streptococcus_pyogenes_Nd_up[,2:ncol(Streptococcus_pyogenes_Nd_up)])
rownames(Streptococcus_pyogenes_expr_Nd_up) =Streptococcus_pyogenes_Nd_up[,1]

###########STAGA_Nd_up_createdbyCansu############
Streptococcus_agalactiae_Nd_up=read.csv("STAGA_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Streptococcus_agalactiae_expr_Nd_up=as.data.frame(Streptococcus_agalactiae_Nd_up[,2:ncol(Streptococcus_agalactiae_Nd_up)])
rownames(Streptococcus_agalactiae_expr_Nd_up) =Streptococcus_agalactiae_Nd_up[,1]

###########UPEC_Nd1_up_createdbyCansu############
Escherichia_coli_UPEC_Nd1_up=read.csv("UPEC_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Escherichia_coli_UPEC_expr_Nd1_up=as.data.frame(Escherichia_coli_UPEC_Nd1_up[,2:ncol(Escherichia_coli_UPEC_Nd1_up)])
rownames(Escherichia_coli_UPEC_expr_Nd1_up) =Escherichia_coli_UPEC_Nd1_up[,1]

###########UPEC_Nd2_up_createdbyCansu############
Escherichia_coli_UPEC_Nd2_up=read.csv("UPEC_Nd_1Nd_2Nd_3Oxs_2_UP.csv",  sep=",")
Escherichia_coli_UPEC_expr_Nd2_up=as.data.frame(Escherichia_coli_UPEC_Nd2_up[,2:ncol(Escherichia_coli_UPEC_Nd2_up)])
rownames(Escherichia_coli_UPEC_expr_Nd2_up) =Escherichia_coli_UPEC_Nd2_up[,1]

##################UPEC_Nd3_up##########################
Escherichia_coli_UPEC_Nd3_up=read.csv("UPEC_Mig_3Nd_1Nd_2Nd_3Sp_3_UP.csv", sep=",")
Escherichia_coli_UPEC_expr_Nd3_up=as.data.frame(Escherichia_coli_UPEC_Nd3_up[,2:ncol(Escherichia_coli_UPEC_Nd3_up)])
rownames(Escherichia_coli_UPEC_expr_Nd3_up) =Escherichia_coli_UPEC_Nd3_up[,1]

###########Vibrio_Nd_up_createdbyCansu############
Vibrio_cholerae_Nd_up=read.csv("Vibrio_Nd_1Nd_2Nd_3Sp_1_UP.csv",  sep=",")
Vibrio_cholerae_expr_Nd_up=as.data.frame(Vibrio_cholerae_Nd_up[,2:ncol(Vibrio_cholerae_Nd_up)])
rownames(Vibrio_cholerae_expr_Nd_up) =Vibrio_cholerae_Nd_up[,1]

###########YPSTB_Nd_up_createdbyCansu############
Yersinia_pseudotuberculosis_Nd_up=read.csv("YPSTB_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Yersinia_pseudotuberculosis_expr_Nd_up=as.data.frame(Yersinia_pseudotuberculosis_Nd_up[,2:ncol(Yersinia_pseudotuberculosis_Nd_up)])
rownames(Yersinia_pseudotuberculosis_expr_Nd_up) =Yersinia_pseudotuberculosis_Nd_up[,1]

###########YPSTB_Nd_up_Vic_down_createdbyCansu############
Yersinia_pseudotuberculosis_Nd_up_Vic_down=read.csv("YPSTB_Nd_1Nd_2Nd_3_UP_Mig_1Mig_3Vic_1Vic_2Vic_3_DOWN.csv",  sep=",")
Yersinia_pseudotuberculosis_expr_Nd_up_Vic_down=as.data.frame(Yersinia_pseudotuberculosis_Nd_up_Vic_down[,2:ncol(Yersinia_pseudotuberculosis_Nd_up_Vic_down)])
rownames(Yersinia_pseudotuberculosis_expr_Nd_up_Vic_down) =Yersinia_pseudotuberculosis_Nd_up_Vic_down[,1]


###########ENTFA_Nd_up_Oxs_up_Vic_down_createdbyCansu############
Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down=read.csv("ENTFA_Nd_1Nd_2Nd_3Oxs_1Oxs_2Oxs_3_UP_Vic_1Vic_2Vic_3_DOWN.csv",  sep=",")
Enterococcus_faecalis_expr_Nd_up_Oxs_up_Vic_down=as.data.frame(Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down[,2:ncol(Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down)])
rownames(Enterococcus_faecalis_expr_Nd_up_Oxs_up_Vic_down) =Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down[,1]

###########KLEBS_Nd_up_Oxs_up_createdbyCansu############
Klebsiella_pneumoniae_Nd_up_Oxs_up=read.csv("KLEBS_Nd_1Nd_2Nd_3Oxs_1Oxs_2Oxs_3_UP.csv",  sep=",")
Klebsiella_pneumoniae_expr_Nd_up_Oxs_up=as.data.frame(Klebsiella_pneumoniae_Nd_up_Oxs_up[,2:ncol(Klebsiella_pneumoniae_Nd_up_Oxs_up)])
rownames(Klebsiella_pneumoniae_expr_Nd_up_Oxs_up) =Klebsiella_pneumoniae_Nd_up_Oxs_up[,1]

###########SALMT_Nd_up_Oxs_up_createdbyCansu############
Salmonella_enterica_Nd_up_Oxs_up=read.csv("SALMT_Nd_1Nd_2Nd_3Oxs_1Oxs_2Oxs_3_UP.csv",  sep=",")
Salmonella_enterica_expr_Nd_up_Oxs_up=as.data.frame(Salmonella_enterica_Nd_up_Oxs_up[,2:ncol(Salmonella_enterica_Nd_up_Oxs_up)])
rownames(Salmonella_enterica_expr_Nd_up_Oxs_up) =Salmonella_enterica_Nd_up_Oxs_up[,1]

###########SHIF_Nd_up_Oxs_up_createdbyCansu############
Shigella_flexneri_Nd_up_Oxs_up=read.csv("SHIF_Nd_1Nd_2Nd_3Oxs_1Oxs_2Oxs_3_UP.csv",  sep=",")
Shigella_flexneri_expr_Nd_up_Oxs_up=as.data.frame(Shigella_flexneri_Nd_up_Oxs_up[,2:ncol(Shigella_flexneri_Nd_up_Oxs_up)])
rownames(Shigella_flexneri_expr_Nd_up_Oxs_up) =Shigella_flexneri_Nd_up_Oxs_up[,1]

################MSSA476_Nd_Sp_up_createdbyCansu###############
Staphylococcus_aureus_MSSA476_Sp_up_Nd_up=read.csv("MSSA476_Nd_1Nd_2Nd_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Staphylococcus_aureus_MSSA476_expr_Sp_up_Nd_up=as.data.frame(Staphylococcus_aureus_MSSA476_Sp_up_Nd_up[,2:ncol(Staphylococcus_aureus_MSSA476_Sp_up_Nd_up)])
rownames(Staphylococcus_aureus_MSSA476_expr_Sp_up_Nd_up) =Staphylococcus_aureus_MSSA476_Sp_up_Nd_up[,1]

################Pseudo_Nd_up_Sp_up_createdbyCansu###############
Pseudomonas_aeruginosa_Sp_up_Nd_up=read.csv("PSEUDO_Nd_1Nd_2Nd_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Sp_up_Nd_up=as.data.frame(Pseudomonas_aeruginosa_Sp_up_Nd_up[,2:ncol(Pseudomonas_aeruginosa_Sp_up_Nd_up)])
rownames(Pseudomonas_aeruginosa_expr_Sp_up_Nd_up) =Pseudomonas_aeruginosa_Sp_up_Nd_up[,1]

##################Listeria_Nd_down_Li_up_##########################
Listeria_monocytogenes_Li_up_Nd_down=read.csv("Listeria_Li_1Li_2Li_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Listeria_monocytogenes_expr_Li_up_Nd_down=as.data.frame(Listeria_monocytogenes_Li_up_Nd_down[,2:ncol(Listeria_monocytogenes_Li_up_Nd_down)])
rownames(Listeria_monocytogenes_expr_Li_up_Nd_down) =Listeria_monocytogenes_Li_up_Nd_down[,1]

##################Vibrio_Nd_up_Mig_up##########################
Vibrio_cholerae_Mig_up_Nd_up=read.csv("Vibrio_Mig_1Mig_2Mig_3Nd_1Nd_2Nd_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Mig_up_Nd_up=as.data.frame(Vibrio_cholerae_Mig_up_Nd_up[,2:ncol(Vibrio_cholerae_Mig_up_Nd_up)])
rownames(Vibrio_cholerae_expr_Mig_up_Nd_up) =Vibrio_cholerae_Mig_up_Nd_up[,1]

##################AGGA_Nd_down_Mig_up##########################
Aggregatibacter_actinomycetemcomitans_Mig_up_Nd_down=read.csv("AGGA_Mig_1Mig_2Mig_3_UP_Nd_1Nd_2Nd_3Tm_2_DOWN.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Mig_up_Nd_down=as.data.frame(Aggregatibacter_actinomycetemcomitans_Mig_up_Nd_down[,2:ncol(Aggregatibacter_actinomycetemcomitans_Mig_up_Nd_down)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Mig_up_Nd_down) =Aggregatibacter_actinomycetemcomitans_Mig_up_Nd_down[,1]

##################MRSA252_Nd_up_Mig_up##########################
Staphylococcus_aureus_MRSA252_Mig_up_Nd_up=read.csv("MRSA252_Mig_1Mig_2Mig_3Nd_1Nd_2Nd_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Mig_up_Nd_up=as.data.frame(Staphylococcus_aureus_MRSA252_Mig_up_Nd_up[,2:ncol(Staphylococcus_aureus_MRSA252_Mig_up_Nd_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Mig_up_Nd_up) =Staphylococcus_aureus_MRSA252_Mig_up_Nd_up[,1]

##################STAGA_Nd_down_Mig_up##########################
Streptococcus_agalactiae_Mig_up_Nd_down=read.csv("STAGA_Mig_1Mig_2Mig_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Mig_up_Nd_down=as.data.frame(Streptococcus_agalactiae_Mig_up_Nd_down[,2:ncol(Streptococcus_agalactiae_Mig_up_Nd_down)])
rownames(Streptococcus_agalactiae_expr_Mig_up_Nd_down) =Streptococcus_agalactiae_Mig_up_Nd_down[,1]


################Listeria_Nd_down_Oss_up###############
Listeria_monocytogenes_Oss_up_Nd_down=read.csv("Listeria_Oss_1Oss_2Oss_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Listeria_monocytogenes_expr_Oss_up_Nd_down=as.data.frame(Listeria_monocytogenes_Oss_up_Nd_down[,2:ncol(Listeria_monocytogenes_Oss_up_Nd_down)])
rownames(Listeria_monocytogenes_expr_Oss_up_Nd_down) =Listeria_monocytogenes_Oss_up_Nd_down[,1]

################STAGA_Nd_down_Oss_up_###############
Streptococcus_agalactiae_Oss_up_Nd_down=read.csv("STAGA_Oss_1Oss_2Oss_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Oss_up_Nd_down =as.data.frame(Streptococcus_agalactiae_Oss_up_Nd_down[,2:ncol(Streptococcus_agalactiae_Oss_up_Nd_down)])
rownames(Streptococcus_agalactiae_expr_Oss_up_Nd_down) =Streptococcus_agalactiae_Oss_up_Nd_down[,1]

################SPYO_Nd_up_Tm_up###############
Streptococcus_pyogenes_Tm_up_Nd_up=read.csv("SPYO_Nd_1Nd_2Nd_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Streptococcus_pyogenes_expr_Tm_up_Nd_up =as.data.frame(Streptococcus_pyogenes_Tm_up_Nd_up[,2:ncol(Streptococcus_pyogenes_Tm_up_Nd_up)])
rownames(Streptococcus_pyogenes_expr_Tm_up_Nd_up) =Streptococcus_pyogenes_Tm_up_Nd_up[,1]

##################STAGA_Nd_down_Tm_up##########################
Streptococcus_agalactiae_Tm_up_Nd_down=read.csv("STAGA_Tm_1Tm_2Tm_3Vic_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Tm_up_Nd_down=as.data.frame(Streptococcus_agalactiae_Tm_up_Nd_down[,2:ncol(Streptococcus_agalactiae_Tm_up_Nd_down)])
rownames(Streptococcus_agalactiae_expr_Tm_up_Nd_down) =Streptococcus_agalactiae_Tm_up_Nd_down[,1]

##################STAGA_Nd_down_Vic_up##########################
Streptococcus_agalactiae_Vic_up_Nd_down=read.csv("STAGA_Vic_1Vic_2Vic_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Vic_up_Nd_down=as.data.frame(Streptococcus_agalactiae_Vic_up_Nd_down[,2:ncol(Streptococcus_agalactiae_Vic_up_Nd_down)])
rownames(Streptococcus_agalactiae_expr_Vic_up_Nd_down) =Streptococcus_agalactiae_Vic_up_Nd_down[,1]


##################Listeria_Nd1_down##########################
Listeria_monocytogenes_Nd1_down=read.csv("Listeria_Mig_2Mig_3Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Listeria_monocytogenes_expr_Nd1_down=as.data.frame(Listeria_monocytogenes_Nd1_down[,2:ncol(Listeria_monocytogenes_Nd1_down)])
rownames(Listeria_monocytogenes_expr_Nd1_down) =Listeria_monocytogenes_Nd1_down[,1]


###########Listeria_Nd2_down_createdbyCansu############
Listeria_monocytogenes_Nd2_down=read.csv("Listeria_Ctrl_1Ctrl_2Ctrl_3_UP_Nd_1Nd_2Nd_3_DOWN.csv",  sep=",")
Listeria_monocytogenes_expr_Nd2_down=as.data.frame(Listeria_monocytogenes_Nd2_down[,2:ncol(Listeria_monocytogenes_Nd2_down)])
rownames(Listeria_monocytogenes_expr_Nd2_down) =Listeria_monocytogenes_Nd2_down[,1]

#########Listeria_Nd3_down#########
Listeria_monocytogenes_Nd3_down = read.csv("Listeria_Mig_1Mig_2_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Listeria_monocytogenes_expr_Nd3_down =as.data.frame(Listeria_monocytogenes_Nd3_down[,2:ncol(Listeria_monocytogenes_Nd3_down)])
rownames(Listeria_monocytogenes_expr_Nd3_down) = Listeria_monocytogenes_Nd3_down[,1]




#########Listeria_Nd_down_As_Up#########
Listeria_monocytogenes_As_up_Nd_down = read.csv("Listeria_As_1As_2As_3Mig_1Mig_2_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Listeria_monocytogenes_expr_As_up_Nd_down =as.data.frame(Listeria_monocytogenes_As_up_Nd_down[,2:ncol(Listeria_monocytogenes_As_up_Nd_down)])
rownames(Listeria_monocytogenes_expr_As_up_Nd_down) = Listeria_monocytogenes_As_up_Nd_down[,1]



###########HPJ99_Nd1_up_createdbyCansu############
Helicobacter_pylori_J99_Nd_up=read.csv("HP_J99_Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Helicobacter_pylori_J99_expr_Nd_up=as.data.frame(Helicobacter_pylori_J99_Nd_up[,2:ncol(Helicobacter_pylori_J99_Nd_up)])
rownames(Helicobacter_pylori_J99_expr_Nd_up) =Helicobacter_pylori_J99_Nd_up[,1]

###########HPJ99_Nd2_up_createdbyCansu############
Helicobacter_pylori_J99_Nd2_up=read.csv("HP_J99_Nd_1Nd_2Nd_3Sp_2Sp_3_UP.csv",  sep=",")
Helicobacter_pylori_J99_expr_Nd2_up=as.data.frame(Helicobacter_pylori_J99_Nd2_up[,2:ncol(Helicobacter_pylori_J99_Nd2_up)])
rownames(Helicobacter_pylori_J99_expr_Nd2_up) =Helicobacter_pylori_J99_Nd2_up[,1]

##################BileStress##########################
##################EPEC_Bs##########################
Escherichia_coli_EPEC_0127_H6_E2348_69_Bs = read.csv("EPEC_Bs_1Bs_2Bs_3Ctrl_1Ctrl_3_UP.csv", sep=",")
Escherichia_coli_EPEC_0127_H6_E2348_69_expr_Bs =as.data.frame(Escherichia_coli_EPEC_0127_H6_E2348_69_Bs[,2:ncol(Escherichia_coli_EPEC_0127_H6_E2348_69_Bs)])
rownames(Escherichia_coli_EPEC_0127_H6_E2348_69_expr_Bs) = Escherichia_coli_EPEC_0127_H6_E2348_69_Bs[,1]

##################ETEC_Bs##########################
Escherichia_coli_ETEC_H10407_Bs = read.csv("ETEC_Bs_1Bs_2Bs_3Ctrl_1Ctrl_2Ctrl_3_UP.csv", sep=",")
Escherichia_coli_ETEC_H10407_expr_Bs=as.data.frame(Escherichia_coli_ETEC_H10407_Bs[,2:ncol(Escherichia_coli_ETEC_H10407_Bs)])
rownames(Escherichia_coli_ETEC_H10407_expr_Bs) = Escherichia_coli_ETEC_H10407_Bs[,1]

##################KLEBS_Bs_Sp##########################
Klebsiella_pneumoniae_Bs_Sp= read.csv("KLEBS_Bs_1Bs_2Bs_3Ctrl_1Ctrl_2Ctrl_3_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Klebsiella_pneumoniae_expr_Bs_Sp=as.data.frame(Klebsiella_pneumoniae_Bs_Sp[,2:ncol(Klebsiella_pneumoniae_Bs_Sp)])
rownames(Klebsiella_pneumoniae_expr_Bs_Sp) = Klebsiella_pneumoniae_Bs_Sp[,1]

##################KLEBS_Bs##########################
Klebsiella_pneumoniae_Bs=read.csv("KLEBS_Bs_1Bs_2Bs_3_UP.csv", sep=",") 
Klebsiella_pneumoniae_expr_Bs=as.data.frame(Klebsiella_pneumoniae_Bs[,2:ncol(Klebsiella_pneumoniae_Bs)])
rownames(Klebsiella_pneumoniae_expr_Bs) =Klebsiella_pneumoniae_Bs[,1]

##################SALMT_Bs##########################
Salmonella_enterica_Bs= read.csv("SALMT_Bs_1Bs_2Bs_3Ctrl_1Ctrl_2Ctrl_3_UP.csv", sep=",")
Salmonella_enterica_expr_Bs=as.data.frame(Salmonella_enterica_Bs[,2:ncol(Salmonella_enterica_Bs)])
rownames(Salmonella_enterica_expr_Bs) = Salmonella_enterica_Bs[,1]

##################UPEC_Bs##########################
Escherichia_coli_UPEC_Bs= read.csv("UPEC_Bs_1Bs_2Bs_3Ctrl_1Ctrl_2Ctrl_3_UP.csv", sep=",")
Escherichia_coli_UPEC_expr_Bs=as.data.frame(Escherichia_coli_UPEC_Bs[,2:ncol(Escherichia_coli_UPEC_Bs)])
rownames(Escherichia_coli_UPEC_expr_Bs) = Escherichia_coli_UPEC_Bs[,1]

##################FRAT_ Bs##########################
Francisella_tularensis_Bs=read.csv("FRAT_Bs_1Bs_2Bs_3_UP.csv", sep=",") 
Francisella_tularensis_expr_Bs= as.data.frame(Francisella_tularensis_Bs[,2:ncol(Francisella_tularensis_Bs)])
rownames(Francisella_tularensis_expr_Bs) = Francisella_tularensis_Bs[,1]

##################ENTFA_Bs##########################
Enterococcus_faecalis_Bs=read.csv("ENTFA_Bs_1Bs_2Bs_3_UP.csv", sep=",") 
Enterococcus_faecalis_expr_Bs=as.data.frame(Enterococcus_faecalis_Bs[,2:ncol(Enterococcus_faecalis_Bs)])
rownames(Enterococcus_faecalis_expr_Bs) =Enterococcus_faecalis_Bs[,1]


##################Listeria_Bs1##########################
Listeria_Bs1=read.csv("Listeria_Bs_1Bs_2Bs_3_UP.csv", sep=",") 
Listeria_expr_Bs1=as.data.frame(Listeria_Bs1[,2:ncol(Listeria_Bs1)])
rownames(Listeria_expr_Bs1) =Listeria_Bs1[,1]

##################Listeria_Bs2##########################
Listeria_Bs2=read.csv("Listeria_Bs_1Bs_2Bs_3_UPpink.csv", sep=",") 
Listeria_expr_Bs2=as.data.frame(Listeria_Bs2[,2:ncol(Listeria_Bs2)])
rownames(Listeria_expr_Bs2) =Listeria_Bs2[,1]


##################MRSA252_Bs##########################
Staphylococcus_aureus_MRSA252_Bs=read.csv("MRSA252_Bs_1Bs_2Bs_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Bs=as.data.frame(Staphylococcus_aureus_MRSA252_Bs[,2:ncol(Staphylococcus_aureus_MRSA252_Bs)])
rownames(Staphylococcus_aureus_MRSA252_expr_Bs) =Staphylococcus_aureus_MRSA252_Bs[,1]

##################NGON_Bs##########################
Neisseria_gonorrhoeae_Bs1=read.csv("NGON_Bs_1Bs_2Bs_3_UP.csv", sep=",")
Neisseria_gonorrhoeae_expr_Bs1=as.data.frame(Neisseria_gonorrhoeae_Bs1[,2:ncol(Neisseria_gonorrhoeae_Bs1)])
rownames(Neisseria_gonorrhoeae_expr_Bs1) =Neisseria_gonorrhoeae_Bs1[,1]

##################NGON_Bs_down##########################
Neisseria_gonorrhoeae_Bs2=read.csv("NGON_Bs_1Bs_2Bs_3_UP_Sp_3_DOWN.csv", sep=",")
Neisseria_gonorrhoeae_expr_Bs2=as.data.frame(Neisseria_gonorrhoeae_Bs2[,2:ncol(Neisseria_gonorrhoeae_Bs2)])
rownames(Neisseria_gonorrhoeae_expr_Bs2) =Neisseria_gonorrhoeae_Bs2[,1]

##################PSEUDO_Bs##########################
Pseudomonas_aeruginosa_Bs=read.csv("PSEUDO_Bs_1Bs_2Bs_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Bs=as.data.frame(Pseudomonas_aeruginosa_Bs[,2:ncol(Pseudomonas_aeruginosa_Bs)])
rownames(Pseudomonas_aeruginosa_expr_Bs) =Pseudomonas_aeruginosa_Bs[,1]

##################Staphylococcus_epidermidis_Bs##########################
Staphylococcus_epidermidis_Bs=read.csv("SEPI_Bs_1Bs_2Bs_3_UP.csv", sep=",")
Staphylococcus_epidermidis_expr_Bs=as.data.frame(Staphylococcus_epidermidis_Bs[,2:ncol(Staphylococcus_epidermidis_Bs)])
rownames(Staphylococcus_epidermidis_expr_Bs) =Staphylococcus_epidermidis_Bs[,1]

##################Achromobacter_xylosoxidans_Bs##########################
Achromobacter_xylosoxidans_Bs=read.csv("ACHX_Bs_1Bs_2Bs_3_UP.csv", sep=",")
Achromobacter_xylosoxidans_expr_Bs=as.data.frame(Achromobacter_xylosoxidans_Bs[,2:ncol(Achromobacter_xylosoxidans_Bs)])
rownames(Achromobacter_xylosoxidans_expr_Bs) =Achromobacter_xylosoxidans_Bs[,1]


##################Aggregatibacter_actinomycetemcomitans_Bs1##########################
Aggregatibacter_actinomycetemcomitans_Bs1=read.csv("AGGA_Bs_1Bs_2Bs_3_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Bs1=as.data.frame(Aggregatibacter_actinomycetemcomitans_Bs1[,2:ncol(Aggregatibacter_actinomycetemcomitans_Bs1)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Bs1) =Aggregatibacter_actinomycetemcomitans_Bs1[,1]

##################Aggregatibacter_actinomycetemcomitans_Bs2##########################
Aggregatibacter_actinomycetemcomitans_Bs2=read.csv("AGGA_Bs_1Bs_2Bs_3Nd_1_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Bs2=as.data.frame(Aggregatibacter_actinomycetemcomitans_Bs2[,2:ncol(Aggregatibacter_actinomycetemcomitans_Bs2)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Bs2) =Aggregatibacter_actinomycetemcomitans_Bs2[,1]

##################Aggregatibacter_actinomycetemcomitans_Bs3##########################
Aggregatibacter_actinomycetemcomitans_Bs3=read.csv("AGGA_Bs_1Bs_2Bs_3Oxs_3_UP_Nd_2_DOWN.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Bs3=as.data.frame(Aggregatibacter_actinomycetemcomitans_Bs3[,2:ncol(Aggregatibacter_actinomycetemcomitans_Bs3)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Bs3) =Aggregatibacter_actinomycetemcomitans_Bs3[,1]

##################Aggregatibacter_actinomycetemcomitans_Bs4##########################
Aggregatibacter_actinomycetemcomitans_Bs4=read.csv("AGGA_Bs_1Bs_2Bs_3_UP_Vic_1_DOWN.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Bs4=as.data.frame(Aggregatibacter_actinomycetemcomitans_Bs4[,2:ncol(Aggregatibacter_actinomycetemcomitans_Bs4)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Bs4) =Aggregatibacter_actinomycetemcomitans_Bs4[,1]

##################Campy_Bs##########################
Campylobacter_jejuni_Bs=read.csv("Campy_Bs_1Bs_2Bs_3_UP.csv", sep=",")
Campylobacter_jejuni_expr_Bs=as.data.frame(Campylobacter_jejuni_Bs[,2:ncol(Campylobacter_jejuni_Bs)])
rownames(Campylobacter_jejuni_expr_Bs) =Campylobacter_jejuni_Bs[,1]

##################Streptococcus_pneumoniae_Bs_Mig##########################
Streptococcus_pneumoniae_Bs_Mig=read.csv("STRPN_Bs_1Bs_2Bs_3Vic_1Vic_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Streptococcus_pneumoniae_expr_Bs_Mig=as.data.frame(Streptococcus_pneumoniae_Bs_Mig[,2:ncol(Streptococcus_pneumoniae_Bs_Mig)])
rownames(Streptococcus_pneumoniae_expr_Bs_Mig) =Streptococcus_pneumoniae_Bs_Mig[,1]

##################Vibrio_cholerae_Bs##########################
Vibrio_cholerae_Bs=read.csv("Vibrio_Bs_1Bs_2Bs_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Bs=as.data.frame(Vibrio_cholerae_Bs[,2:ncol(Vibrio_cholerae_Bs)])
rownames(Vibrio_cholerae_expr_Bs) =Vibrio_cholerae_Bs[,1]

##################Neisseria_meningitidis_Bs##########################
Neisseria_meningitidis_Bs=read.csv("NMEN_Bs_1Bs_2Bs_3_UP.csv", sep=",")
Neisseria_meningitidis_expr_Bs=as.data.frame(Neisseria_meningitidis_Bs[,2:ncol(Neisseria_meningitidis_Bs)])
rownames(Neisseria_meningitidis_expr_Bs) =Neisseria_meningitidis_Bs[,1]

##################LowIron##########################

##################ACHX_Li1_up##########################
Achromobacter_xylosoxidans_Li1_up=read.csv("ACHX_Li_1Li_2Li_3Oxs_2_UP.csv", sep=",")
Achromobacter_xylosoxidans_expr_Li1_up=as.data.frame(Achromobacter_xylosoxidans_Li1_up[,2:ncol(Achromobacter_xylosoxidans_Li1_up)])
rownames(Achromobacter_xylosoxidans_expr_Li1_up) =Achromobacter_xylosoxidans_Li1_up[,1]

##################ACHX_Li2_up##########################
Achromobacter_xylosoxidans_Li2_up=read.csv("ACHX_Li_1Li_2Li_3_UP.csv", sep=",")
Achromobacter_xylosoxidans_expr_Li2_up=as.data.frame(Achromobacter_xylosoxidans_Li2_up[,2:ncol(Achromobacter_xylosoxidans_Li2_up)])
rownames(Achromobacter_xylosoxidans_expr_Li2_up) =Achromobacter_xylosoxidans_Li2_up[,1]

##################ACIB_Li_up_Oxs_up##########################
Acinetobacter_baumannii_Li_up_Oxs_up=read.csv("ACIB_Li_1Li_2Li_3Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Acinetobacter_baumannii_expr_Li_up_Oxs_up=as.data.frame(Acinetobacter_baumannii_Li_up_Oxs_up[,2:ncol(Acinetobacter_baumannii_Li_up_Oxs_up)])
rownames(Acinetobacter_baumannii_expr_Li_up_Oxs_up) =Acinetobacter_baumannii_Li_up_Oxs_up[,1]

##################AGGA_Li1_up##########################
Aggregatibacter_actinomycetemcomitans_Li1_up=read.csv("AGGA_Li_1Li_2Li_3Nd_1Nd_3Oxs_1_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Li1_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Li1_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Li1_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Li1_up) =Aggregatibacter_actinomycetemcomitans_Li1_up[,1]

##################AGGA_Li2_up##########################
Aggregatibacter_actinomycetemcomitans_Li2_up=read.csv("AGGA_Li_1Li_2Li_3_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Li2_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Li2_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Li2_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Li2_up) =Aggregatibacter_actinomycetemcomitans_Li2_up[,1]

##################Burk_Li_up##########################
Burkholderia_pseudomallei_Li_up=read.csv("BURK_Li_1Li_2Li_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Li_up=as.data.frame(Burkholderia_pseudomallei_Li_up[,2:ncol(Burkholderia_pseudomallei_Li_up)])
rownames(Burkholderia_pseudomallei_expr_Li_up) =Burkholderia_pseudomallei_Li_up[,1]

##################Campy_Li_up##########################
Campylobacter_jejuni_Li_up=read.csv("Campy_Li_1Li_2Li_3_UP.csv", sep=",")
Campylobacter_jejuni_expr_Li_up=as.data.frame(Campylobacter_jejuni_Li_up[,2:ncol(Campylobacter_jejuni_Li_up)])
rownames(Campylobacter_jejuni_expr_Li_up) =Campylobacter_jejuni_Li_up[,1]

##################ENTFA_Li_up##########################
Enterococcus_faecalis_Li_up=read.csv("ENTFA_Li_1Li_2Li_3_UP.csv", sep=",")
Enterococcus_faecalis_expr_Li_up=as.data.frame(Enterococcus_faecalis_Li_up[,2:ncol(Enterococcus_faecalis_Li_up)])
rownames(Enterococcus_faecalis_expr_Li_up) =Enterococcus_faecalis_Li_up[,1]

##################EPEC_Li_up##########################
Escherichia_coli_EPEC_Li_up=read.csv("EPEC_Li_1Li_2Li_3_UP.csv", sep=",")
Escherichia_coli_EPEC_expr_Li_up=as.data.frame(Escherichia_coli_EPEC_Li_up[,2:ncol(Escherichia_coli_EPEC_Li_up)])
rownames(Escherichia_coli_EPEC_expr_Li_up) =Escherichia_coli_EPEC_Li_up[,1]

##################ETEC_Li_up##########################
Escherichia_coli_ETEC_Li_up=read.csv("ETEC_Li_1Li_2Li_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Li_up=as.data.frame(Escherichia_coli_ETEC_Li_up[,2:ncol(Escherichia_coli_ETEC_Li_up)])
rownames(Escherichia_coli_ETEC_expr_Li_up) =Escherichia_coli_ETEC_Li_up[,1]

##################FRAT_Li_up##########################
Francisella_tularensis_Li_up=read.csv("FRAT_Li_1Li_2Li_3_UP.csv", sep=",")
Francisella_tularensis_expr_Li_up=as.data.frame(Francisella_tularensis_Li_up[,2:ncol(Francisella_tularensis_Li_up)])
rownames(Francisella_tularensis_expr_Li_up) =Francisella_tularensis_Li_up[,1]


##################HP-G27_Li_up##########################
Helicobacter_pylori_G27_Li_up=read.csv("HP_G27_Li_1Li_2Li_3_UP.csv", sep=",")
Helicobacter_pylori_G27_expr_Li_up=as.data.frame(Helicobacter_pylori_G27_Li_up[,2:ncol(Helicobacter_pylori_G27_Li_up)])
rownames(Helicobacter_pylori_G27_expr_Li_up) =Helicobacter_pylori_G27_Li_up[,1]

##################HP-J99_Li_up##########################
Helicobacter_pylori_J99_Li_up=read.csv("HP_J99_Li_1Li_2Li_3_UP.csv", sep=",")
Helicobacter_pylori_J99_expr_Li_up=as.data.frame(Helicobacter_pylori_J99_Li_up[,2:ncol(Helicobacter_pylori_J99_Li_up)])
rownames(Helicobacter_pylori_J99_expr_Li_up) =Helicobacter_pylori_J99_Li_up[,1]

##################KLEBS_Li_up##########################
Klebsiella_pneumoniae_Li_up=read.csv("KLEBS_Li_1Li_2Li_3_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Li_up=as.data.frame(Klebsiella_pneumoniae_Li_up[,2:ncol(Klebsiella_pneumoniae_Li_up)])
rownames(Klebsiella_pneumoniae_expr_Li_up) =Klebsiella_pneumoniae_Li_up[,1]

##################Listeria_Li_up_Nd_down##########################
Listeria_monocytogenes_Li_up_Nd_down=read.csv("Listeria_Li_1Li_2Li_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Listeria_monocytogenes_expr_Li_up_Nd_down=as.data.frame(Listeria_monocytogenes_Li_up_Nd_down[,2:ncol(Listeria_monocytogenes_Li_up_Nd_down)])
rownames(Listeria_monocytogenes_expr_Li_up_Nd_down) =Listeria_monocytogenes_Li_up_Nd_down[,1]


##################MRSA252_Li_up##########################
Staphylococcus_aureus_MRSA252_Li_up=read.csv("MRSA252_Li_1Li_2Li_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Li_up=as.data.frame(Staphylococcus_aureus_MRSA252_Li_up[,2:ncol(Staphylococcus_aureus_MRSA252_Li_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Li_up) =Staphylococcus_aureus_MRSA252_Li_up[,1]

#################NGON_Li_up##########################
Neisseria_gonorrhoeae_Li_up=read.csv("NGON_Li_1Li_2Li_3_UP.csv", sep=",")
Neisseria_gonorrhoeae_expr_Li_up=as.data.frame(Neisseria_gonorrhoeae_Li_up[,2:ncol(Neisseria_gonorrhoeae_Li_up)])
rownames(Neisseria_gonorrhoeae_expr_Li_up) =Neisseria_gonorrhoeae_Li_up[,1]


#################PSEUDO_Li1_up##########################
Pseudomonas_aeruginosa_Li1_up=read.csv("PSEUDO_Li_1Li_2Li_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Li1_up=as.data.frame(Pseudomonas_aeruginosa_Li1_up[,2:ncol(Pseudomonas_aeruginosa_Li1_up)])
rownames(Pseudomonas_aeruginosa_expr_Li1_up) =Pseudomonas_aeruginosa_Li1_up[,1]

#################PSEUDO_Li2_up##########################
Pseudomonas_aeruginosa_Li2_up=read.csv("PSEUDO_Li_1Li_2Li_3_UPmagenta.csv", sep=",")
Pseudomonas_aeruginosa_expr_Li2_up=as.data.frame(Pseudomonas_aeruginosa_Li2_up[,2:ncol(Pseudomonas_aeruginosa_Li2_up)])
rownames(Pseudomonas_aeruginosa_expr_Li2_up) =Pseudomonas_aeruginosa_Li2_up[,1]

#################SALMT_Li_up##########################
Salmonella_enterica_Li_up=read.csv("SALMT_Li_1Li_2Li_3_UP.csv", sep=",")
Salmonella_enterica_expr_Li_up=as.data.frame(Salmonella_enterica_Li_up[,2:ncol(Salmonella_enterica_Li_up)])
rownames(Salmonella_enterica_expr_Li_up) =Salmonella_enterica_Li_up[,1]

#################SALMT_Li_up_Oxs_up##########################
Salmonella_enterica_Li_up_Oxs_up=read.csv("SALMT_Li_1Li_2Li_3Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Salmonella_enterica_expr_Li_up_Oxs_up=as.data.frame(Salmonella_enterica_Li_up_Oxs_up[,2:ncol(Salmonella_enterica_Li_up_Oxs_up)])
rownames(Salmonella_enterica_expr_Li_up_Oxs_up) =Salmonella_enterica_Li_up_Oxs_up[,1]

#################SEPI_Li_up##########################
Staphylococcus_epidermidis_Li_up=read.csv("SEPI_Li_1Li_2Li_3_UP.csv", sep=",")
Staphylococcus_epidermidis_expr_Li_up=as.data.frame(Staphylococcus_epidermidis_Li_up[,2:ncol(Staphylococcus_epidermidis_Li_up)])
rownames(Staphylococcus_epidermidis_expr_Li_up) =Staphylococcus_epidermidis_Li_up[,1]

#################SHIF_Li_up##########################
Shigella_flexneri_Li_up=read.csv("SHIF_Li_1Li_2Li_3Ns_1_UP.csv", sep=",")
Shigella_flexneri_expr_Li_up=as.data.frame(Shigella_flexneri_Li_up[,2:ncol(Shigella_flexneri_Li_up)])
rownames(Shigella_flexneri_expr_Li_up) =Shigella_flexneri_Li_up[,1]

#################SHIF_Li_up_Oxs_up##########################
Shigella_flexneri_Li_up_Oxs_up=read.csv("SHIF_Li_1Li_2Li_3Ns_1Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Shigella_flexneri_expr_Li_up_Oxs_up=as.data.frame(Shigella_flexneri_Li_up_Oxs_up[,2:ncol(Shigella_flexneri_Li_up_Oxs_up)])
rownames(Shigella_flexneri_expr_Li_up_Oxs_up) =Shigella_flexneri_Li_up_Oxs_up[,1]

#################STAGA_Li1_up##########################
Streptococcus_agalactiae_Li1_up=read.csv("STAGA_Li_1Li_2Li_3_UP.csv", sep=",")
Streptococcus_agalactiae_expr_Li1_up=as.data.frame(Streptococcus_agalactiae_Li1_up[,2:ncol(Streptococcus_agalactiae_Li1_up)])
rownames(Streptococcus_agalactiae_expr_Li1_up) =Streptococcus_agalactiae_Li1_up[,1]

#################STAGA_Li2_up##########################
Streptococcus_agalactiae_Li2_up=read.csv("STAGA_Li_1Li_2Li_3Ns_2Ns_3_UP.csv", sep=",")
Streptococcus_agalactiae_expr_Li2_up=as.data.frame(Streptococcus_agalactiae_Li2_up[,2:ncol(Streptococcus_agalactiae_Li2_up)])
rownames(Streptococcus_agalactiae_expr_Li2_up) =Streptococcus_agalactiae_Li2_up[,1]

#################STRPN_Li_up##########################
Streptococcus_pneumoniae_Li_up=read.csv("STRPN_Li_1Li_2Li_3_UP.csv", sep=",")
Streptococcus_pneumoniae_expr_Li_up=as.data.frame(Streptococcus_pneumoniae_Li_up[,2:ncol(Streptococcus_pneumoniae_Li_up)])
rownames(Streptococcus_pneumoniae_expr_Li_up) =Streptococcus_pneumoniae_Li_up[,1]

#################UPEC_Li_up##########################
Escherichia_coli_UPEC_Li_up=read.csv("UPEC_Li_1Li_2Li_3_UP.csv", sep=",")
Escherichia_coli_UPEC_expr_Li_up=as.data.frame(Escherichia_coli_UPEC_Li_up[,2:ncol(Escherichia_coli_UPEC_Li_up)])
rownames(Escherichia_coli_UPEC_expr_Li_up) =Escherichia_coli_UPEC_Li_up[,1]

#################Vibrio_Li1_up##########################
Vibrio_cholerae_Li1_up=read.csv("Vibrio_Li_1Li_2Li_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Li1_up=as.data.frame(Vibrio_cholerae_Li1_up[,2:ncol(Vibrio_cholerae_Li1_up)])
rownames(Vibrio_cholerae_expr_Li1_up) =Vibrio_cholerae_Li1_up[,1]

#################Vibrio_Li2_up##########################
Vibrio_cholerae_Li2_up=read.csv("Vibrio_Li_1Li_2Li_3Oss_2_UP.csv", sep=",")
Vibrio_cholerae_expr_Li2_up=as.data.frame(Vibrio_cholerae_Li2_up[,2:ncol(Vibrio_cholerae_Li2_up)])
rownames(Vibrio_cholerae_expr_Li2_up) =Vibrio_cholerae_Li2_up[,1]

#################YPSTB_Li_up##########################
Yersinia_pseudotuberculosis_Li_up=read.csv("YPSTB_Li_1Li_2Li_3_UP.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Li_up=as.data.frame(Yersinia_pseudotuberculosis_Li_up[,2:ncol(Yersinia_pseudotuberculosis_Li_up)])
rownames(Yersinia_pseudotuberculosis_expr_Li_up) =Yersinia_pseudotuberculosis_Li_up[,1]

################NMEN_Li_up_Sp_up_createdbyCansu###############
Neisseria_meningitidis_Sp_up_Li_up=read.csv("NMEN_Li_1Li_2Li_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Neisseria_meningitidis_expr_Sp_up_Li_up=as.data.frame(Neisseria_meningitidis_Sp_up_Li_up[,2:ncol(Neisseria_meningitidis_Sp_up_Li_up)])
rownames(Neisseria_meningitidis_expr_Sp_up_Li_up) =Neisseria_meningitidis_Sp_up_Li_up[,1]

################Vibrio_Li_up_Sp´_up_createdbyCansu###############
Vibrio_cholerae_Sp_up_Li_up=read.csv("Vibrio_Li_1Li_2Li_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Sp_up_Li_up=as.data.frame(Vibrio_cholerae_Sp_up_Li_up[,2:ncol(Vibrio_cholerae_Sp_up_Li_up)])
rownames(Vibrio_cholerae_expr_Sp_up_Li_up) =Vibrio_cholerae_Sp_up_Li_up[,1]

##################ETEC_Li_up_Tm_up##########################
Escherichia_coli_ETEC_Tm_up_Li_up=read.csv("ETEC_Li_1Li_2Li_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Tm_up_Li_up=as.data.frame(Escherichia_coli_ETEC_Tm_up_Li_up[,2:ncol(Escherichia_coli_ETEC_Tm_up_Li_up)])
rownames(Escherichia_coli_ETEC_expr_Tm_up_Li_up) =Escherichia_coli_ETEC_Tm_up_Li_up[,1]

#################STRPN_Li_up_Mig_down##########################
Streptococcus_pneumoniae_Li_up_Mig_down=read.csv("STRPN_Bs_2Bs_3Li_1Li_2Li_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Streptococcus_pneumoniae_expr_Li_up_Mig_down=as.data.frame(Streptococcus_pneumoniae_Li_up_Mig_down[,2:ncol(Streptococcus_pneumoniae_Li_up_Mig_down)])
rownames(Streptococcus_pneumoniae_expr_Li_up_Mig_down) =Streptococcus_pneumoniae_Li_up_Mig_down[,1]

########SPYO_Li_up_Sp_down_Mig_down###############
Streptococcus_pyogenes_Mig_Sp_Li = read.csv("SPYO_Li_1Li_2Li_3_UP_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Streptococcus_pyogenes_expr_Mig_Sp_Li =as.data.frame(Streptococcus_pyogenes_Mig_Sp_Li[,2:ncol(Streptococcus_pyogenes_Mig_Sp_Li)])
rownames(Streptococcus_pyogenes_expr_Mig_Sp_Li) = Streptococcus_pyogenes_Mig_Sp_Li[,1]


##################Microaerophilic(hypoxia)##########################

##################BURK_Mig_up_Sp_up##########################
Burkholderia_pseudomallei_Mig_up_Sp_up=read.csv("BURK_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Mig_up_Sp_up=as.data.frame(Burkholderia_pseudomallei_Mig_up_Sp_up[,2:ncol(Burkholderia_pseudomallei_Mig_up_Sp_up)])
rownames(Burkholderia_pseudomallei_expr_Mig_up_Sp_up) =Burkholderia_pseudomallei_Mig_up_Sp_up[,1]


##################ENTFA_Mig_up_Vic_up##########################
Enterococcus_faecalis_Mig_up_Vic_up=read.csv("ENTFA_Mig_1Mig_2Mig_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Enterococcus_faecalis_expr_Mig_up_Vic_up=as.data.frame(Enterococcus_faecalis_Mig_up_Vic_up[,2:ncol(Enterococcus_faecalis_Mig_up_Vic_up)])
rownames(Enterococcus_faecalis_expr_Mig_up_Vic_up) =Enterococcus_faecalis_Mig_up_Vic_up[,1]

##################SPYO_Mig_up_Sp_up##########################
Streptococcus_pyogenes_Mig_up_Sp_up=read.csv("SPYO_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Streptococcus_pyogenes_expr_Mig_up_Sp_up=as.data.frame(Streptococcus_pyogenes_Mig_up_Sp_up[,2:ncol(Streptococcus_pyogenes_Mig_up_Sp_up)])
rownames(Streptococcus_pyogenes_expr_Mig_up_Sp_up) =Streptococcus_pyogenes_Mig_up_Sp_up[,1]



##################Vibrio_Mig_up_Nd_up##########################
Vibrio_cholerae_Mig_up_Nd_up=read.csv("Vibrio_Mig_1Mig_2Mig_3Nd_1Nd_2Nd_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Mig_up_Nd_up=as.data.frame(Vibrio_cholerae_Mig_up_Nd_up[,2:ncol(Vibrio_cholerae_Mig_up_Nd_up)])
rownames(Vibrio_cholerae_expr_Mig_up_Nd_up) =Vibrio_cholerae_Mig_up_Nd_up[,1]

##################ACHX_Mig1_up##########################
Achromobacter_xylosoxidans_Mig1_up=read.csv("ACHX_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Achromobacter_xylosoxidans_expr_Mig1_up=as.data.frame(Achromobacter_xylosoxidans_Mig1_up[,2:ncol(Achromobacter_xylosoxidans_Mig1_up)])
rownames(Achromobacter_xylosoxidans_expr_Mig1_up) =Achromobacter_xylosoxidans_Mig1_up[,1]

##################ACHX_Mig2_up##########################
Achromobacter_xylosoxidans_Mig2_up=read.csv("ACHX_Mig_1Mig_2Mig_3Ns_1Ns_2_UP.csv", sep=",")
Achromobacter_xylosoxidans_expr_Mig2_up=as.data.frame(Achromobacter_xylosoxidans_Mig2_up[,2:ncol(Achromobacter_xylosoxidans_Mig2_up)])
rownames(Achromobacter_xylosoxidans_expr_Mig2_up) =Achromobacter_xylosoxidans_Mig2_up[,1]

##################AGGA_Mig_up_Nd_down##########################
Aggregatibacter_actinomycetemcomitans_Mig_up_Nd_down=read.csv("AGGA_Mig_1Mig_2Mig_3_UP_Nd_1Nd_2Nd_3Tm_2_DOWN.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Mig_up_Nd_down=as.data.frame(Aggregatibacter_actinomycetemcomitans_Mig_up_Nd_down[,2:ncol(Aggregatibacter_actinomycetemcomitans_Mig_up_Nd_down)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Mig_up_Nd_down) =Aggregatibacter_actinomycetemcomitans_Mig_up_Nd_down[,1]

##################BURK_Mig1_up##########################
Burkholderia_pseudomallei_Mig1_up=read.csv("BURK_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Mig1_up=as.data.frame(Burkholderia_pseudomallei_Mig1_up[,2:ncol(Burkholderia_pseudomallei_Mig1_up)])
rownames(Burkholderia_pseudomallei_expr_Mig1_up) =Burkholderia_pseudomallei_Mig1_up[,1]

##################BURK_Mig2_up##########################
Burkholderia_pseudomallei_Mig2_up=read.csv("BURK_Mig_1Mig_2Mig_3Oxs_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Mig2_up=as.data.frame(Burkholderia_pseudomallei_Mig2_up[,2:ncol(Burkholderia_pseudomallei_Mig2_up)])
rownames(Burkholderia_pseudomallei_expr_Mig2_up) =Burkholderia_pseudomallei_Mig2_up[,1]


##################ENTFA_Mig_up##########################
Enterococcus_faecalis_Mig_up=read.csv("ENTFA_Mig_1Mig_2Mig_3_UP_Nd_2_DOWN.csv", sep=",")
Enterococcus_faecalis_expr_Mig_up=as.data.frame(Enterococcus_faecalis_Mig_up[,2:ncol(Enterococcus_faecalis_Mig_up)])
rownames(Enterococcus_faecalis_expr_Mig_up) =Enterococcus_faecalis_Mig_up[,1]

##################EPEC_Mig_up##########################
Escherichia_coli_EPEC_Mig_up=read.csv("EPEC_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Escherichia_coli_EPEC_expr_Mig_up=as.data.frame(Escherichia_coli_EPEC_Mig_up[,2:ncol(Escherichia_coli_EPEC_Mig_up)])
rownames(Escherichia_coli_EPEC_expr_Mig_up) =Escherichia_coli_EPEC_Mig_up[,1]

##################ETEC_Mig_up##########################
Escherichia_coli_ETEC_Mig_up=read.csv("ETEC_Mig_1Mig_2Mig_3Sp_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Mig_up=as.data.frame(Escherichia_coli_ETEC_Mig_up[,2:ncol(Escherichia_coli_ETEC_Mig_up)])
rownames(Escherichia_coli_ETEC_expr_Mig_up) =Escherichia_coli_ETEC_Mig_up[,1]

##################ETEC_Mig_up_Oss_up##########################
Escherichia_coli_ETEC_Mig_up_Oss_up=read.csv("ETEC_Mig_1Mig_2Mig_3Oss_1Oss_2Oss_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Mig_up_Oss_up=as.data.frame(Escherichia_coli_ETEC_Mig_up_Oss_up[,2:ncol(Escherichia_coli_ETEC_Mig_up_Oss_up)])
rownames(Escherichia_coli_ETEC_expr_Mig_up_Oss_up) =Escherichia_coli_ETEC_Mig_up_Oss_up[,1]

##################HINF_Mig_up##########################
Haemophilus_influenzae_Mig_up=read.csv("HINF_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Haemophilus_influenzae_expr_Mig_up=as.data.frame(Haemophilus_influenzae_Mig_up[,2:ncol(Haemophilus_influenzae_Mig_up)])
rownames(Haemophilus_influenzae_expr_Mig_up) =Haemophilus_influenzae_Mig_up[,1]

##################HINF_Mig_up_Vic_up##########################
Haemophilus_influenzae_Mig_up_Vic_up=read.csv("HINF_Mig_1Mig_2Mig_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Haemophilus_influenzae_expr_Mig_up_Vic_up=as.data.frame(Haemophilus_influenzae_Mig_up_Vic_up[,2:ncol(Haemophilus_influenzae_Mig_up_Vic_up)])
rownames(Haemophilus_influenzae_expr_Mig_up_Vic_up) =Haemophilus_influenzae_Mig_up_Vic_up[,1]

##################KLEBS_Mig_up##########################
Klebsiella_pneumoniae_Mig_up=read.csv("KLEBS_Mig_1Mig_2Mig_3Tm_2_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Mig_up=as.data.frame(Klebsiella_pneumoniae_Mig_up[,2:ncol(Klebsiella_pneumoniae_Mig_up)])
rownames(Klebsiella_pneumoniae_expr_Mig_up) =Klebsiella_pneumoniae_Mig_up[,1]

##################Listeria_Mig_up##########################
Listeria_monocytogenes_Mig_up=read.csv("Listeria_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Listeria_monocytogenes_expr_Mig_up=as.data.frame(Listeria_monocytogenes_Mig_up[,2:ncol(Listeria_monocytogenes_Mig_up)])
rownames(Listeria_monocytogenes_expr_Mig_up) =Listeria_monocytogenes_Mig_up[,1]

##################MRSA252_Mig_up_Nd_up##########################
Staphylococcus_aureus_MRSA252_Mig_up_Nd_up=read.csv("MRSA252_Mig_1Mig_2Mig_3Nd_1Nd_2Nd_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Mig_up_Nd_up=as.data.frame(Staphylococcus_aureus_MRSA252_Mig_up_Nd_up[,2:ncol(Staphylococcus_aureus_MRSA252_Mig_up_Nd_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Mig_up_Nd_up) =Staphylococcus_aureus_MRSA252_Mig_up_Nd_up[,1]

##################MRSA252_Mig1_up##########################
Staphylococcus_aureus_MRSA252_Mig1_up=read.csv("MRSA252_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Mig1_up=as.data.frame(Staphylococcus_aureus_MRSA252_Mig1_up[,2:ncol(Staphylococcus_aureus_MRSA252_Mig1_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Mig1_up) =Staphylococcus_aureus_MRSA252_Mig1_up[,1]

##################MRSA252_Mig2_up##########################
Staphylococcus_aureus_MRSA252_Mig2_up=read.csv("MRSA252_Mig_1Mig_2Mig_3_UPsalmon.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Mig2_up=as.data.frame(Staphylococcus_aureus_MRSA252_Mig2_up[,2:ncol(Staphylococcus_aureus_MRSA252_Mig2_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Mig2_up) =Staphylococcus_aureus_MRSA252_Mig2_up[,1]

##################MSSA476_Mig_up##########################
Staphylococcus_aureus_MSSA476_Mig_up=read.csv("MSSA476_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Staphylococcus_aureus_MSSA476_expr_Mig_up=as.data.frame(Staphylococcus_aureus_MSSA476_Mig_up[,2:ncol(Staphylococcus_aureus_MSSA476_Mig_up)])
rownames(Staphylococcus_aureus_MSSA476_expr_Mig_up) =Staphylococcus_aureus_MSSA476_Mig_up[,1]

##################NGON_Mig1_up##########################
Neisseria_gonorrhoeae_Mig1_up=read.csv("NGON_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Neisseria_gonorrhoeae_expr_Mig1_up=as.data.frame(Neisseria_gonorrhoeae_Mig1_up[,2:ncol(Neisseria_gonorrhoeae_Mig1_up)])
rownames(Neisseria_gonorrhoeae_expr_Mig1_up) =Neisseria_gonorrhoeae_Mig1_up[,1]

##################NGON_Mig2_up##########################
Neisseria_gonorrhoeae_Mig2_up=read.csv("NGON_Mig_1Mig_2Mig_3Nd_1Nd_3_UP_Bs_1Bs_3_DOWN.csv", sep=",")
Neisseria_gonorrhoeae_expr_Mig2_up=as.data.frame(Neisseria_gonorrhoeae_Mig2_up[,2:ncol(Neisseria_gonorrhoeae_Mig2_up)])
rownames(Neisseria_gonorrhoeae_expr_Mig2_up) =Neisseria_gonorrhoeae_Mig2_up[,1]

##################NMEN_Mig_up##########################
Neisseria_meningitidis_Mig_up=read.csv("NMEN_Mig_1Mig_2Mig_3Vic_2_UP.csv", sep=",")
Neisseria_meningitidis_expr_Mig_up=as.data.frame(Neisseria_meningitidis_Mig_up[,2:ncol(Neisseria_meningitidis_Mig_up)])
rownames(Neisseria_meningitidis_expr_Mig_up) =Neisseria_meningitidis_Mig_up[,1]

##################PSEUDO_Mig_up_Oxs_up##########################
Pseudomonas_aeruginosa_Mig_up_Oxs_up=read.csv("PSEUDO_Mig_1Mig_2Mig_3Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Mig_up_Oxs_up=as.data.frame(Pseudomonas_aeruginosa_Mig_up_Oxs_up[,2:ncol(Pseudomonas_aeruginosa_Mig_up_Oxs_up)])
rownames(Pseudomonas_aeruginosa_expr_Mig_up_Oxs_up) =Pseudomonas_aeruginosa_Mig_up_Oxs_up[,1]

##################PSEUDO_Mig_up##########################
Pseudomonas_aeruginosa_Mig_up=read.csv("PSEUDO_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Mig_up=as.data.frame(Pseudomonas_aeruginosa_Mig_up[,2:ncol(Pseudomonas_aeruginosa_Mig_up)])
rownames(Pseudomonas_aeruginosa_expr_Mig_up) =Pseudomonas_aeruginosa_Mig_up[,1]

##################PSEUDO_Mig_up_Ns_down##########################
Pseudomonas_aeruginosa_Mig_up_Ns_down=read.csv("PSEUDO_Mig_1Mig_2Mig_3_UP_Nd_1Nd_2Ns_1Ns_2Ns_3_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_Mig_up_Ns_down=as.data.frame(Pseudomonas_aeruginosa_Mig_up_Ns_down[,2:ncol(Pseudomonas_aeruginosa_Mig_up_Ns_down)])
rownames(Pseudomonas_aeruginosa_expr_Mig_up_Ns_down) =Pseudomonas_aeruginosa_Mig_up_Ns_down[,1]

##################Salmt_Mig_up##########################
Salmonella_enterica_Mig_up=read.csv("SALMT_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Salmonella_enterica_expr_Mig_up=as.data.frame(Salmonella_enterica_Mig_up[,2:ncol(Salmonella_enterica_Mig_up)])
rownames(Salmonella_enterica_expr_Mig_up) =Salmonella_enterica_Mig_up[,1]

##################SEPI_Mig1_up##########################
Staphylococcus_epidermidis_Mig1_up=read.csv("SEPI_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Staphylococcus_epidermidis_expr_Mig1_up=as.data.frame(Staphylococcus_epidermidis_Mig1_up[,2:ncol(Staphylococcus_epidermidis_Mig1_up)])
rownames(Staphylococcus_epidermidis_expr_Mig1_up) =Staphylococcus_epidermidis_Mig1_up[,1]

##################SEPI_Mig2_up##########################
Staphylococcus_epidermidis_Mig2_up=read.csv("SEPI_Mig_1Mig_2Mig_3_UP_Vic_1Vic_2_DOWN.csv", sep=",")
Staphylococcus_epidermidis_expr_Mig2_up=as.data.frame(Staphylococcus_epidermidis_Mig2_up[,2:ncol(Staphylococcus_epidermidis_Mig2_up)])
rownames(Staphylococcus_epidermidis_expr_Mig2_up) =Staphylococcus_epidermidis_Mig2_up[,1]

##################SEPI_Mig_up_Tm_down##########################
Staphylococcus_epidermidis_Mig_up_Tm_down=read.csv("SEPI_Mig_1Mig_2Mig_3Nd_1_UP_Tm_1Tm_2Tm_3_DOWN.csv", sep=",")
Staphylococcus_epidermidis_expr_Mig_up_Tm_down=as.data.frame(Staphylococcus_epidermidis_Mig_up_Tm_down[,2:ncol(Staphylococcus_epidermidis_Mig_up_Tm_down)])
rownames(Staphylococcus_epidermidis_expr_Mig_up_Tm_down) =Staphylococcus_epidermidis_Mig_up_Tm_down[,1]

##################SEPI_Mig_up_As_up##########################
Staphylococcus_epidermidis_Mig_up_As_up=read.csv("SEPI_As_1As_2As_3Mig_1Mig_2Mig_3_UP.csv", sep=",")
Staphylococcus_epidermidis_expr_Mig_up_As_up=as.data.frame(Staphylococcus_epidermidis_Mig_up_As_up[,2:ncol(Staphylococcus_epidermidis_Mig_up_As_up)])
rownames(Staphylococcus_epidermidis_expr_Mig_up_As_up) =Staphylococcus_epidermidis_Mig_up_As_up[,1]

##################SHIF_Mig_up##########################
Shigella_flexneri_Mig_up=read.csv("SHIF_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Shigella_flexneri_expr_Mig_up=as.data.frame(Shigella_flexneri_Mig_up[,2:ncol(Shigella_flexneri_Mig_up)])
rownames(Shigella_flexneri_expr_Mig_up) =Shigella_flexneri_Mig_up[,1]

##################STAGA_Mig_up_Nd_down##########################
Streptococcus_agalactiae_Mig_up_Nd_down=read.csv("STAGA_Mig_1Mig_2Mig_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Mig_up_Nd_down=as.data.frame(Streptococcus_agalactiae_Mig_up_Nd_down[,2:ncol(Streptococcus_agalactiae_Mig_up_Nd_down)])
rownames(Streptococcus_agalactiae_expr_Mig_up_Nd_down) =Streptococcus_agalactiae_Mig_up_Nd_down[,1]

##################STRPN_Mig_up##########################
Streptococcus_pneumoniae_Mig_up=read.csv("STRPN_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Streptococcus_pneumoniae_expr_Mig_up=as.data.frame(Streptococcus_pneumoniae_Mig_up[,2:ncol(Streptococcus_pneumoniae_Mig_up)])
rownames(Streptococcus_pneumoniae_expr_Mig_up) =Streptococcus_pneumoniae_Mig_up[,1]

##################UPEC_Mig_up##########################
Escherichia_coli_UPEC_Mig_up=read.csv("UPEC_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Escherichia_coli_UPEC_expr_Mig_up=as.data.frame(Escherichia_coli_UPEC_Mig_up[,2:ncol(Escherichia_coli_UPEC_Mig_up)])
rownames(Escherichia_coli_UPEC_expr_Mig_up) =Escherichia_coli_UPEC_Mig_up[,1]

##################UPEC_Mig_up_Vic_up##########################
Escherichia_coli_UPEC_Mig_up_Vic_up=read.csv("UPEC_Mig_1Mig_2Mig_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Escherichia_coli_UPEC_expr_Mig_up_Vic_up=as.data.frame(Escherichia_coli_UPEC_Mig_up_Vic_up[,2:ncol(Escherichia_coli_UPEC_Mig_up_Vic_up)])
rownames(Escherichia_coli_UPEC_expr_Mig_up_Vic_up) =Escherichia_coli_UPEC_Mig_up_Vic_up[,1]

##################Vibrio_Mig1_up##########################
Vibrio_cholerae_Mig1_up=read.csv("Vibrio_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Mig1_up=as.data.frame(Vibrio_cholerae_Mig1_up[,2:ncol(Vibrio_cholerae_Mig1_up)])
rownames(Vibrio_cholerae_expr_Mig1_up) =Vibrio_cholerae_Mig1_up[,1]

##################Vibrio_Mig2_up##########################
Vibrio_cholerae_Mig2_up=read.csv("Vibrio_Ctrl_1Ctrl_2Ctrl_3Mig_1Mig_2Mig_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Mig2_up=as.data.frame(Vibrio_cholerae_Mig2_up[,2:ncol(Vibrio_cholerae_Mig2_up)])
rownames(Vibrio_cholerae_expr_Mig2_up) =Vibrio_cholerae_Mig2_up[,1]

##################YPSTB_Mig_up##########################
Yersinia_pseudotuberculosis_Mig_up=read.csv("YPSTB_Mig_1Mig_2Mig_3_UP.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Mig_up=as.data.frame(Yersinia_pseudotuberculosis_Mig_up[,2:ncol(Yersinia_pseudotuberculosis_Mig_up)])
rownames(Yersinia_pseudotuberculosis_expr_Mig_up) =Yersinia_pseudotuberculosis_Mig_up[,1]

########SPYO_Mig_Sp_down1###############
Streptococcus_pyogenes_Sp_Mg = read.csv("SPYO_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Streptococcus_pyogenes_expr_Sp_Mg =as.data.frame(Streptococcus_pyogenes_Sp_Mg[,2:ncol(Streptococcus_pyogenes_Sp_Mg)])
rownames(Streptococcus_pyogenes_expr_Sp_Mg) = Streptococcus_pyogenes_Sp_Mg[,1]

################SPYO_Mig_Sp_down2###############
Streptococcus_pyogenes_Mig2_down_Sp_down=read.csv("SPYO_Vic_1_UP_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Streptococcus_pyogenes_expr_Mig2_down_Sp_down =as.data.frame(Streptococcus_pyogenes_Mig2_down_Sp_down[,2:ncol(Streptococcus_pyogenes_Mig2_down_Sp_down)])
rownames(Streptococcus_pyogenes_expr_Mig2_down_Sp_down) =Streptococcus_pyogenes_Mig2_down_Sp_down[,1]


################Ypstb_Mig_down_Sp_up_createdbyCansu###############
Yersinia_pseudotuberculosis_Sp_up_Mig_down=read.csv("YPSTB_Sp_1Sp_2Sp_3Vic_1Vic_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Sp_up_Mig_down=as.data.frame(Yersinia_pseudotuberculosis_Sp_up_Mig_down[,2:ncol(Yersinia_pseudotuberculosis_Sp_up_Mig_down)])
rownames(Yersinia_pseudotuberculosis_expr_Sp_up_Mig_down) =Yersinia_pseudotuberculosis_Sp_up_Mig_down[,1]



########Entfa_Mig_up_Sp_up_createdbyCansu###############
Enterococcus_faecalis_Sp_up_Mig_up=read.csv("ENTFA_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Enterococcus_faecalis_expr_Sp_up_Mig_up =as.data.frame(Enterococcus_faecalis_Sp_up_Mig_up[,2:ncol(Enterococcus_faecalis_Sp_up_Mig_up)])
rownames(Enterococcus_faecalis_expr_Sp_up_Mig_up) =Enterococcus_faecalis_Sp_up_Mig_up[,1]




########STAGA_Mig_up_Sp_up_createdbyCansu###############
Streptococcus_agalactiae_Sp_up_Mig_up=read.csv("STAGA_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_UP.csv", sep=",")
Streptococcus_agalactiae_expr_Sp_up_Mig_up =as.data.frame(Streptococcus_agalactiae_Sp_up_Mig_up[,2:ncol(Streptococcus_agalactiae_Sp_up_Mig_up)])
rownames(Streptococcus_agalactiae_expr_Sp_up_Mig_up) =Streptococcus_agalactiae_Sp_up_Mig_up[,1]




###########MSSA476_Mig_Nd############
Staphylococcus_aureus_MSSA476_Nd_Mig =read.csv("MSSA476_Ns_3_UP_Mig_1Mig_2Mig_3Nd_1Nd_2Nd_3Sp_3_DOWN.csv",  sep=",")
Staphylococcus_aureus_MSSA476_expr_Nd_Mig =as.data.frame(Staphylococcus_aureus_MSSA476_Nd_Mig[,2:ncol(Staphylococcus_aureus_MSSA476_Nd_Mig)])
rownames(Staphylococcus_aureus_MSSA476_expr_Nd_Mig) = Staphylococcus_aureus_MSSA476_Nd_Mig[,1]

###########MRSA252_Mig_Nd############
Staphylococcus_aureus_MRSA252_Nd_Mig =read.csv("MRSA252_Mig_1Mig_2Mig_3Nd_1Nd_2Nd_3_DOWN.csv",  sep=",")
Staphylococcus_aureus_MRSA252_expr_Nd_Mig =as.data.frame(Staphylococcus_aureus_MRSA252_Nd_Mig[,2:ncol(Staphylococcus_aureus_MRSA252_Nd_Mig)])
rownames(Staphylococcus_aureus_MRSA252_expr_Nd_Mig) = Staphylococcus_aureus_MRSA252_Nd_Mig[,1]

###########ETEC_Mig_up_Nd_upcreatedbyCansu############
Escherichia_coli_ETEC_Nd_up_Mig_up=read.csv("ETEC_Mig_1Mig_2Mig_3Nd_1Nd_2Nd_3_UP.csv",  sep=",")
Escherichia_coli_ETEC_expr_Nd_up_Mig_up=as.data.frame(Escherichia_coli_ETEC_Nd_up_Mig_up[,2:ncol(Escherichia_coli_ETEC_Nd_up_Mig_up)])
rownames(Escherichia_coli_ETEC_expr_Nd_up_Mig_up) =Escherichia_coli_ETEC_Nd_up_Mig_up[,1]


##################STRPN_Mig_Bs##########################
Streptococcus_pneumoniae_Bs_Mig=read.csv("STRPN_Bs_1Bs_2Bs_3Vic_1Vic_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Streptococcus_pneumoniae_expr_Bs_Mig=as.data.frame(Streptococcus_pneumoniae_Bs_Mig[,2:ncol(Streptococcus_pneumoniae_Bs_Mig)])
rownames(Streptococcus_pneumoniae_expr_Bs_Mig) =Streptococcus_pneumoniae_Bs_Mig[,1]

##################HINF_Mig_down_Vic_up##########################
Haemophilus_influenzae_Vic_up_Mig_down=read.csv("HINF_Vic_1Vic_2Vic_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Haemophilus_influenzae_expr_Vic_up_Mig_down=as.data.frame(Haemophilus_influenzae_Vic_up_Mig_down[,2:ncol(Haemophilus_influenzae_Vic_up_Mig_down)])
rownames(Haemophilus_influenzae_expr_Vic_up_Mig_down) =Haemophilus_influenzae_Vic_up_Mig_down[,1]

##################STRPN_Mig_down_Vic_up_##########################
Streptococcus_pneumoniae_Vic_up_Mig_down=read.csv("STRPN_Vic_1Vic_2Vic_3_UP_Mig_1Mig_2Mig_3Nd_3_DOWN.csv", sep=",")
Streptococcus_pneumoniae_expr_Vic_up_Mig_down=as.data.frame(Streptococcus_pneumoniae_Vic_up_Mig_down[,2:ncol(Streptococcus_pneumoniae_Vic_up_Mig_down)])
rownames(Streptococcus_pneumoniae_expr_Vic_up_Mig_down) =Streptococcus_pneumoniae_Vic_up_Mig_down[,1]

################YPSTB_Mig_down_Vic_up_###############
Yersinia_pseudotuberculosis_Vic_up_Mig_down=read.csv("YPSTB_Vic_1Vic_2Vic_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Vic_up_Mig_down =as.data.frame(Yersinia_pseudotuberculosis_Vic_up_Mig_down[,2:ncol(Yersinia_pseudotuberculosis_Vic_up_Mig_down)])
rownames(Yersinia_pseudotuberculosis_expr_Vic_up_Mig_down) =Yersinia_pseudotuberculosis_Vic_up_Mig_down[,1]

#################STRPN_Mig_down_Li_up##########################
Streptococcus_pneumoniae_Li_up_Mig_down=read.csv("STRPN_Bs_2Bs_3Li_1Li_2Li_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Streptococcus_pneumoniae_expr_Li_up_Mig_down=as.data.frame(Streptococcus_pneumoniae_Li_up_Mig_down[,2:ncol(Streptococcus_pneumoniae_Li_up_Mig_down)])
rownames(Streptococcus_pneumoniae_expr_Li_up_Mig_down) =Streptococcus_pneumoniae_Li_up_Mig_down[,1]




##################HINF_Mig1_down##########################
Haemophilus_influenzae_Mig_down=read.csv("HINF_Mig_1Mig_2Mig_3Ns_3_DOWN.csv", sep=",")
Haemophilus_influenzae_expr_Mig_down=as.data.frame(Haemophilus_influenzae_Mig_down[,2:ncol(Haemophilus_influenzae_Mig_down)])
rownames(Haemophilus_influenzae_expr_Mig_down) =Haemophilus_influenzae_Mig_down[,1]

##################HINF_Mig2_down##########################
Haemophilus_influenzae_Mig2_down=read.csv("HINF_Nd_1Nd_2_UP_Mig_1Mig_2Mig_3Oss_2_DOWN.csv", sep=",")
Haemophilus_influenzae_expr_Mig2_down=as.data.frame(Haemophilus_influenzae_Mig2_down[,2:ncol(Haemophilus_influenzae_Mig2_down)])
rownames(Haemophilus_influenzae_expr_Mig2_down) =Haemophilus_influenzae_Mig2_down[,1]


########SPYO_Mig_down_Sp_down_Li_up###############
Streptococcus_pyogenes_Mig_Sp_Li = read.csv("SPYO_Li_1Li_2Li_3_UP_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Streptococcus_pyogenes_expr_Mig_Sp_Li =as.data.frame(Streptococcus_pyogenes_Mig_Sp_Li[,2:ncol(Streptococcus_pyogenes_Mig_Sp_Li)])
rownames(Streptococcus_pyogenes_expr_Mig_Sp_Li) = Streptococcus_pyogenes_Mig_Sp_Li[,1]

#################UPEC_Mig_Sp_down##########################
Escherichia_coli_UPEC_Mig_Sp_down=read.csv("UPEC_Li_1Li_2Ns_2_UP_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Escherichia_coli_UPEC_expr_Mig_Sp_down=as.data.frame(Escherichia_coli_UPEC_Mig_Sp_down[,2:ncol(Escherichia_coli_UPEC_Mig_Sp_down)])
rownames(Escherichia_coli_UPEC_expr_Mig_Sp_down) =Escherichia_coli_UPEC_Mig_Sp_down[,1]

##################SHIF_Mig_Sp_down##########################
Shigella_flexneri_Mig_Sp_down=read.csv("SHIF_Mig_1Mig_2Mig_3Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Shigella_flexneri_expr_Mig_Sp_down=as.data.frame(Shigella_flexneri_Mig_Sp_down[,2:ncol(Shigella_flexneri_Mig_Sp_down)])
rownames(Shigella_flexneri_expr_Mig_Sp_down) =Shigella_flexneri_Mig_Sp_down[,1]


##################Listeria_Mig_down##########################
Listeria_monocytogenes_Mig_down=read.csv("Listeria_Oxs_1Sp_2Sp_3Vic_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Listeria_monocytogenes_expr_Mig_down=as.data.frame(Listeria_monocytogenes_Mig_down[,2:ncol(Listeria_monocytogenes_Mig_down)])
rownames(Listeria_monocytogenes_expr_Mig_down) =Listeria_monocytogenes_Mig_down[,1]

##################Nitrosative_stress##########################

##################ACHX_Ns_up##########################
Achromobacter_xylosoxidans_Ns_up=read.csv("ACHX_Ns_1Ns_2Ns_3_UP.csv", sep=",")
Achromobacter_xylosoxidans_expr_Ns_up=as.data.frame(Achromobacter_xylosoxidans_Ns_up[,2:ncol(Achromobacter_xylosoxidans_Ns_up)])
rownames(Achromobacter_xylosoxidans_expr_Ns_up) =Achromobacter_xylosoxidans_Ns_up[,1]

##################AGGA_Ns_up##########################
Aggregatibacter_actinomycetemcomitans_Ns_up=read.csv("AGGA_Li_1Ns_1Ns_2Ns_3_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Ns_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Ns_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Ns_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Ns_up) =Aggregatibacter_actinomycetemcomitans_Ns_up[,1]

##################ENTFA_Ns_up##########################
Enterococcus_faecalis_Ns_up=read.csv("ENTFA_Ns_1Ns_2Ns_3_UP.csv", sep=",")
Enterococcus_faecalis_expr_Ns_up=as.data.frame(Enterococcus_faecalis_Ns_up[,2:ncol(Enterococcus_faecalis_Ns_up)])
rownames(Enterococcus_faecalis_expr_Ns_up) =Enterococcus_faecalis_Ns_up[,1]

##################EPEC_Ns_up##########################
Escherichia_coli_EPEC_Ns_up=read.csv("EPEC_Ns_1Ns_2Ns_3_UP.csv", sep=",")
Escherichia_coli_EPEC_expr_Ns_up=as.data.frame(Escherichia_coli_EPEC_Ns_up[,2:ncol(Escherichia_coli_EPEC_Ns_up)])
rownames(Escherichia_coli_EPEC_expr_Ns_up) =Escherichia_coli_EPEC_Ns_up[,1]

##################ETEC_Ns_up##########################
Escherichia_coli_ETEC_Ns_up=read.csv("ETEC_Ns_1Ns_2Ns_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Ns_up=as.data.frame(Escherichia_coli_ETEC_Ns_up[,2:ncol(Escherichia_coli_ETEC_Ns_up)])
rownames(Escherichia_coli_ETEC_expr_Ns_up) =Escherichia_coli_ETEC_Ns_up[,1]

##################HINF_Ns_up##########################
Haemophilus_influenzae_Ns_up=read.csv("HINF_Ns_1Ns_2Ns_3_UP.csv", sep=",")
Haemophilus_influenzae_expr_Ns_up=as.data.frame(Haemophilus_influenzae_Ns_up[,2:ncol(Haemophilus_influenzae_Ns_up)])
rownames(Haemophilus_influenzae_expr_Ns_up) =Haemophilus_influenzae_Ns_up[,1]

##################KLEBS_Ns_up##########################
Klebsiella_pneumoniae_Ns_up=read.csv("KLEBS_Ns_1Ns_2Ns_3_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Ns_up=as.data.frame(Klebsiella_pneumoniae_Ns_up[,2:ncol(Klebsiella_pneumoniae_Ns_up)])
rownames(Klebsiella_pneumoniae_expr_Ns_up) =Klebsiella_pneumoniae_Ns_up[,1]

##################Listeria_Ns1_up##########################
Listeria_monocytogenes_Ns_up=read.csv("Listeria_Ctrl_1Ns_1Ns_2Ns_3_UP.csv", sep=",")
Listeria_monocytogenes_expr_Ns_up=as.data.frame(Listeria_monocytogenes_Ns_up[,2:ncol(Listeria_monocytogenes_Ns_up)])
rownames(Listeria_monocytogenes_expr_Ns_up) =Listeria_monocytogenes_Ns_up[,1]


################Listeria_Ns2_up###############
Listeria_monocytogenes_Ns2_up=read.csv("Listeria_Ns_1Ns_2Ns_3Oxs_1Oxs_2_UP.csv", sep=",")
Listeria_monocytogenes_expr_Ns2_up=as.data.frame(Listeria_monocytogenes_Ns2_up[,2:ncol(Listeria_monocytogenes_Ns2_up)])
rownames(Listeria_monocytogenes_expr_Ns2_up) =Listeria_monocytogenes_Ns2_up[,1]

##################Pseudo_Ns_up##########################
Pseudomonas_aeruginosa_Ns_up=read.csv("PSEUDO_Ns_1Ns_2Ns_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Ns_up=as.data.frame(Pseudomonas_aeruginosa_Ns_up[,2:ncol(Pseudomonas_aeruginosa_Ns_up)])
rownames(Pseudomonas_aeruginosa_expr_Ns_up) =Pseudomonas_aeruginosa_Ns_up[,1]

##################Pseudo_Ns1_up_Sp_down##########################
Pseudomonas_aeruginosa_Ns1_up_Sp_down=read.csv("PSEUDO_Li_2Ns_1Ns_2Ns_3_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_Ns1_up_Sp_down=as.data.frame(Pseudomonas_aeruginosa_Ns1_up_Sp_down[,2:ncol(Pseudomonas_aeruginosa_Ns1_up_Sp_down)])
rownames(Pseudomonas_aeruginosa_expr_Ns1_up_Sp_down) =Pseudomonas_aeruginosa_Ns1_up_Sp_down[,1]

##################Pseudo_Ns2_up_Sp_down##########################
Pseudomonas_aeruginosa_Ns2_up_Sp_down=read.csv("PSEUDO_Ns_1Ns_2Ns_3_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_Ns2_up_Sp_down=as.data.frame(Pseudomonas_aeruginosa_Ns2_up_Sp_down[,2:ncol(Pseudomonas_aeruginosa_Ns2_up_Sp_down)])
rownames(Pseudomonas_aeruginosa_expr_Ns2_up_Sp_down) =Pseudomonas_aeruginosa_Ns2_up_Sp_down[,1]

##################SALMT_Ns_up##########################
Salmonella_enterica_Ns_up=read.csv("SALMT_Ns_1Ns_2Ns_3_UP.csv", sep=",")
Salmonella_enterica_expr_Ns_up=as.data.frame(Salmonella_enterica_Ns_up[,2:ncol(Salmonella_enterica_Ns_up)])
rownames(Salmonella_enterica_expr_Ns_up) =Salmonella_enterica_Ns_up[,1]

##################SEPI_Ns_up##########################
Staphylococcus_epidermidis_Ns_up=read.csv("SEPI_Li_1Li_3Ns_1Ns_2Ns_3_UP.csv", sep=",")
Staphylococcus_epidermidis_expr_Ns_up=as.data.frame(Staphylococcus_epidermidis_Ns_up[,2:ncol(Staphylococcus_epidermidis_Ns_up)])
rownames(Staphylococcus_epidermidis_expr_Ns_up) =Staphylococcus_epidermidis_Ns_up[,1]

##################SEPI_Ns_up_Tm_down##########################
Staphylococcus_epidermidis_Ns_up_Tm_down=read.csv("SEPI_Ns_1Ns_2Ns_3_UP_Tm_1Tm_2Tm_3_DOWN.csv", sep=",")
Staphylococcus_epidermidis_expr_Ns_up_Tm_down=as.data.frame(Staphylococcus_epidermidis_Ns_up_Tm_down[,2:ncol(Staphylococcus_epidermidis_Ns_up_Tm_down)])
rownames(Staphylococcus_epidermidis_expr_Ns_up_Tm_down) =Staphylococcus_epidermidis_Ns_up_Tm_down[,1]

##################SHIF_Ns_up##########################
Shigella_flexneri_Ns_up=read.csv("SHIF_Ns_1Ns_2Ns_3_UP.csv", sep=",")
Shigella_flexneri_expr_Ns_up=as.data.frame(Shigella_flexneri_Ns_up[,2:ncol(Shigella_flexneri_Ns_up)])
rownames(Shigella_flexneri_expr_Ns_up) =Shigella_flexneri_Ns_up[,1]

##################STAGA_Ns_up##########################
Streptococcus_agalactiae_Ns_up=read.csv("STAGA_Ctrl_2Ctrl_3Ns_1Ns_2Ns_3_UP.csv", sep=",")
Streptococcus_agalactiae_expr_Ns_up=as.data.frame(Streptococcus_agalactiae_Ns_up[,2:ncol(Streptococcus_agalactiae_Ns_up)])
rownames(Streptococcus_agalactiae_expr_Ns_up) =Streptococcus_agalactiae_Ns_up[,1]

##################STRPN_Ns_up##########################
Streptococcus_pneumoniae_Ns_up=read.csv("STRPN_Ns_1Ns_2Ns_3Oxs_1_UP.csv", sep=",")
Streptococcus_pneumoniae_expr_Ns_up=as.data.frame(Streptococcus_pneumoniae_Ns_up[,2:ncol(Streptococcus_pneumoniae_Ns_up)])
rownames(Streptococcus_pneumoniae_expr_Ns_up) =Streptococcus_pneumoniae_Ns_up[,1]

###########_MSSA476_Ns_As############
Staphylococcus_aureus_MSSA476_As_Ns =read.csv("MSSA476_As_1As_2As_3_UP_Ns_1Ns_2Ns_3_DOWN.csv",  sep=",")
Staphylococcus_aureus_MSSA476_expr_As_Ns =as.data.frame(Staphylococcus_aureus_MSSA476_As_Ns[,2:ncol(Staphylococcus_aureus_MSSA476_As_Ns)])
rownames(Staphylococcus_aureus_MSSA476_expr_As_Ns) = Staphylococcus_aureus_MSSA476_As_Ns[,1]

###########MRSA252_Ns_down_Nd_up_createdbyCansu############
Staphylococcus_aureus_MRSA252_Nd_up_Ns_down=read.csv("MRSA252_As_2As_3Nd_1Nd_2Nd_3_UP_Ns_1Ns_2Ns_3_DOWN.csv",  sep=",")
Staphylococcus_aureus_MRSA252_expr_Nd_up_Ns_down=as.data.frame(Staphylococcus_aureus_MRSA252_Nd_up_Ns_down[,2:ncol(Staphylococcus_aureus_MRSA252_Nd_up_Ns_down)])
rownames(Staphylococcus_aureus_MRSA252_expr_Nd_up_Ns_down) =Staphylococcus_aureus_MRSA252_Nd_up_Ns_down[,1]

##################PSEUDO_Ns_down_Mig_up##########################
Pseudomonas_aeruginosa_Mig_up_Ns_down=read.csv("PSEUDO_Mig_1Mig_2Mig_3_UP_Nd_1Nd_2Ns_1Ns_2Ns_3_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_Mig_up_Ns_down=as.data.frame(Pseudomonas_aeruginosa_Mig_up_Ns_down[,2:ncol(Pseudomonas_aeruginosa_Mig_up_Ns_down)])
rownames(Pseudomonas_aeruginosa_expr_Mig_up_Ns_down) =Pseudomonas_aeruginosa_Mig_up_Ns_down[,1]

################SPYO_Ns_up_Oss_up###############
Streptococcus_pyogenes_Oss_up_Ns_up=read.csv("SPYO_Ns_1Ns_2Ns_3Oss_1Oss_2Oss_3_UP.csv", sep=",")
Streptococcus_pyogenes_expr_Oss_up_Ns_up =as.data.frame(Streptococcus_pyogenes_Oss_up_Ns_up[,2:ncol(Streptococcus_pyogenes_Oss_up_Ns_up)])
rownames(Streptococcus_pyogenes_expr_Oss_up_Ns_up) =Streptococcus_pyogenes_Oss_up_Ns_up[,1]


###########MSSA476_Ns_down############
Staphylococcus_aureus_MSSA476_Ns_down =read.csv("MSSA476_Vic_1_UP_Li_2Nd_1Nd_3Ns_1Ns_2Ns_3_DOWN.csv",  sep=",")
Staphylococcus_aureus_MSSA476_expr_Ns_down =as.data.frame(Staphylococcus_aureus_MSSA476_Ns_down[,2:ncol(Staphylococcus_aureus_MSSA476_Ns_down)])
rownames(Staphylococcus_aureus_MSSA476_expr_Ns_down) = Staphylococcus_aureus_MSSA476_Ns_down[,1]

##################Pseudo_Ns_down##########################
Pseudomonas_aeruginosa_Ns_down=read.csv("PSEUDO_Mig_1Mig_3Sp_3_UP_Nd_1Nd_2Ns_1Ns_2Ns_3_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_Ns_down=as.data.frame(Pseudomonas_aeruginosa_Ns_down[,2:ncol(Pseudomonas_aeruginosa_Ns_down)])
rownames(Pseudomonas_aeruginosa_expr_Ns_down) =Pseudomonas_aeruginosa_Ns_down[,1]

###########Pseudo_Ns_up_Nd_up_createdbyCansu############
Pseudomonas_aeruginosa_Nd_up_Ns_up=read.csv("PSEUDO_Nd_1Nd_2Nd_3Ns_1Ns_2Ns_3_UP.csv",  sep=",")
Pseudomonas_aeruginosa_expr_Nd_up_Ns_up=as.data.frame(Pseudomonas_aeruginosa_Nd_up_Ns_up[,2:ncol(Pseudomonas_aeruginosa_Nd_up_Ns_up)])
rownames(Pseudomonas_aeruginosa_expr_Nd_up_Ns_up) =Pseudomonas_aeruginosa_Nd_up_Ns_up[,1]

##################Osmotic_stress##########################

##################ACHX_Oss_up##########################
Achromobacter_xylosoxidans_Oss_up=read.csv("ACHX_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Achromobacter_xylosoxidans_expr_Oss_up=as.data.frame(Achromobacter_xylosoxidans_Oss_up[,2:ncol(Achromobacter_xylosoxidans_Oss_up)])
rownames(Achromobacter_xylosoxidans_expr_Oss_up) =Achromobacter_xylosoxidans_Oss_up[,1]

##################ACIB_Oss_up##########################
Acinetobacter_baumannii_Oss_up=read.csv("ACIB_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Acinetobacter_baumannii_expr_Oss_up=as.data.frame(Acinetobacter_baumannii_Oss_up[,2:ncol(Acinetobacter_baumannii_Oss_up)])
rownames(Acinetobacter_baumannii_expr_Oss_up) =Acinetobacter_baumannii_Oss_up[,1]

##################AGGA_Oss_up##########################
Aggregatibacter_actinomycetemcomitans_Oss_up=read.csv("AGGA_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Oss_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Oss_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Oss_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Oss_up) =Aggregatibacter_actinomycetemcomitans_Oss_up[,1]

##################BBURG_Oss_up##########################
Borrelia_burgdorferi_Oss_up=read.csv("BBURG_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Borrelia_burgdorferi_expr_Oss_up=as.data.frame(Borrelia_burgdorferi_Oss_up[,2:ncol(Borrelia_burgdorferi_Oss_up)])
rownames(Borrelia_burgdorferi_expr_Oss_up) =Borrelia_burgdorferi_Oss_up[,1]

##################BURK_Oss_up##########################
Burkholderia_pseudomallei_Oss_up=read.csv("BURK_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Oss_up=as.data.frame(Burkholderia_pseudomallei_Oss_up[,2:ncol(Burkholderia_pseudomallei_Oss_up)])
rownames(Burkholderia_pseudomallei_expr_Oss_up) =Burkholderia_pseudomallei_Oss_up[,1]

##################ENTFA_Oss_up##########################
Enterococcus_faecalis_Oss_up=read.csv("ENTFA_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Enterococcus_faecalis_expr_Oss_up=as.data.frame(Enterococcus_faecalis_Oss_up[,2:ncol(Enterococcus_faecalis_Oss_up)])
rownames(Enterococcus_faecalis_expr_Oss_up) =Enterococcus_faecalis_Oss_up[,1]

##################EPEC_Oss_up##########################
Escherichia_coli_EPEC_Oss_up=read.csv("EPEC_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Escherichia_coli_EPEC_expr_Oss_up=as.data.frame(Escherichia_coli_EPEC_Oss_up[,2:ncol(Escherichia_coli_EPEC_Oss_up)])
rownames(Escherichia_coli_EPEC_expr_Oss_up) =Escherichia_coli_EPEC_Oss_up[,1]

##################ETEC_Oss_up##########################
Escherichia_coli_ETEC_Oss_up=read.csv("ETEC_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Oss_up=as.data.frame(Escherichia_coli_ETEC_Oss_up[,2:ncol(Escherichia_coli_ETEC_Oss_up)])
rownames(Escherichia_coli_ETEC_expr_Oss_up) =Escherichia_coli_ETEC_Oss_up[,1]


##################FRAT_Oss_up##########################
Francisella_tularensis_Oss_up=read.csv("FRAT_Mig_2Mig_3Oss_1Oss_2Oss_3_UP.csv", sep=",")
Francisella_tularensis_expr_Oss_up=as.data.frame(Francisella_tularensis_Oss_up[,2:ncol(Francisella_tularensis_Oss_up)])
rownames(Francisella_tularensis_expr_Oss_up) =Francisella_tularensis_Oss_up[,1]

##################HINF_Oss_up##########################
Haemophilus_influenzae_Oss_up=read.csv("HINF_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Haemophilus_influenzae_expr_Oss_up=as.data.frame(Haemophilus_influenzae_Oss_up[,2:ncol(Haemophilus_influenzae_Oss_up)])
rownames(Haemophilus_influenzae_expr_Oss_up) =Haemophilus_influenzae_Oss_up[,1]

################HP_G27_Oss_up###############
Helicobacter_pylori_G27_Oss_up=read.csv("HP_G27_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Helicobacter_pylori_G27_expr_Oss_up =as.data.frame(Helicobacter_pylori_G27_Oss_up[,2:ncol(Helicobacter_pylori_G27_Oss_up)])
rownames(Helicobacter_pylori_G27_expr_Oss_up) =Helicobacter_pylori_G27_Oss_up[,1]

################HP_J99_Oss_up###############
Helicobacter_pylori_J99_Oss_up=read.csv("HP_J99_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Helicobacter_pylori_J99_expr_Oss_up =as.data.frame(Helicobacter_pylori_J99_Oss_up[,2:ncol(Helicobacter_pylori_J99_Oss_up)])
rownames(Helicobacter_pylori_J99_expr_Oss_up) =Helicobacter_pylori_J99_Oss_up[,1]

################KLEBS_Oss1_up###############
Klebsiella_pneumoniae_Oss1_up=read.csv("KLEBS_Mig_1Mig_3Oss_1Oss_2Oss_3_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Oss1_up =as.data.frame(Klebsiella_pneumoniae_Oss1_up[,2:ncol(Klebsiella_pneumoniae_Oss1_up)])
rownames(Klebsiella_pneumoniae_expr_Oss1_up) =Klebsiella_pneumoniae_Oss1_up[,1]

################KLEBS_Oss2_up###############
Klebsiella_pneumoniae_Oss2_up=read.csv("KLEBS_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Oss2_up =as.data.frame(Klebsiella_pneumoniae_Oss2_up[,2:ncol(Klebsiella_pneumoniae_Oss2_up)])
rownames(Klebsiella_pneumoniae_expr_Oss2_up) =Klebsiella_pneumoniae_Oss2_up[,1]

################KLEBS_Oss_up_Tm_up###############
Klebsiella_pneumoniae_Oss_up_Tm_up=read.csv("KLEBS_Oss_1Oss_2Oss_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Oss_up_Tm_up=as.data.frame(Klebsiella_pneumoniae_Oss_up_Tm_up[,2:ncol(Klebsiella_pneumoniae_Oss_up_Tm_up)])
rownames(Klebsiella_pneumoniae_expr_Oss_up_Tm_up) =Klebsiella_pneumoniae_Oss_up_Tm_up[,1]

################Listeria_Oss_up_Nd_down###############
Listeria_monocytogenes_Oss_up_Nd_down=read.csv("Listeria_Oss_1Oss_2Oss_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Listeria_monocytogenes_expr_Oss_up_Nd_down=as.data.frame(Listeria_monocytogenes_Oss_up_Nd_down[,2:ncol(Listeria_monocytogenes_Oss_up_Nd_down)])
rownames(Listeria_monocytogenes_expr_Oss_up_Nd_down) =Listeria_monocytogenes_Oss_up_Nd_down[,1]

################MRSA252_Oss_up###############
Staphylococcus_aureus_MRSA252_Oss_up=read.csv("MRSA252_Ctrl_1Ctrl_2Ctrl_3Oss_1Oss_2Oss_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Oss_up =as.data.frame(Staphylococcus_aureus_MRSA252_Oss_up[,2:ncol(Staphylococcus_aureus_MRSA252_Oss_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Oss_up) =Staphylococcus_aureus_MRSA252_Oss_up[,1]

################MSSA476_Oss1_up###############
Staphylococcus_aureus_MSSA476_Oss1_up=read.csv("MSSA476_Ctrl_2Ctrl_3Oss_1Oss_2Oss_3_UP.csv", sep=",")
Staphylococcus_aureus_MSSA476_expr_Oss1_up =as.data.frame(Staphylococcus_aureus_MSSA476_Oss1_up[,2:ncol(Staphylococcus_aureus_MSSA476_Oss1_up)])
rownames(Staphylococcus_aureus_MSSA476_expr_Oss1_up) =Staphylococcus_aureus_MSSA476_Oss1_up[,1]

################MSSA476_Oss2_up###############
Staphylococcus_aureus_MSSA476_Oss2_up=read.csv("MSSA476_Li_2Oss_1Oss_2Oss_3_UP.csv", sep=",")
Staphylococcus_aureus_MSSA476_expr_Oss2_up =as.data.frame(Staphylococcus_aureus_MSSA476_Oss2_up[,2:ncol(Staphylococcus_aureus_MSSA476_Oss2_up)])
rownames(Staphylococcus_aureus_MSSA476_expr_Oss2_up) =Staphylococcus_aureus_MSSA476_Oss2_up[,1]

################NMEN_Oss_up###############
Neisseria_meningitidis_Oss_up=read.csv("NMEN_Oss_1Oss_2Oss_3_UP_Oxs_1_DOWN.csv", sep=",")
Neisseria_meningitidis_expr_Oss_up =as.data.frame(Neisseria_meningitidis_Oss_up[,2:ncol(Neisseria_meningitidis_Oss_up)])
rownames(Neisseria_meningitidis_expr_Oss_up) =Neisseria_meningitidis_Oss_up[,1]

################PSEUDO_Oss_up###############
Pseudomonas_aeruginosa_Oss_up=read.csv("PSEUDO_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Oss_up =as.data.frame(Pseudomonas_aeruginosa_Oss_up[,2:ncol(Pseudomonas_aeruginosa_Oss_up)])
rownames(Pseudomonas_aeruginosa_expr_Oss_up) =Pseudomonas_aeruginosa_Oss_up[,1]

################SALMT_Oss_up###############
Salmonella_enterica_Oss_up=read.csv("SALMT_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Salmonella_enterica_expr_Oss_up =as.data.frame(Salmonella_enterica_Oss_up[,2:ncol(Salmonella_enterica_Oss_up)])
rownames(Salmonella_enterica_expr_Oss_up) =Salmonella_enterica_Oss_up[,1]

################SHIF_Oss_up1###############
Shigella_flexneri_Oss_up=read.csv("SHIF_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Shigella_flexneri_expr_Oss_up =as.data.frame(Shigella_flexneri_Oss_up[,2:ncol(Shigella_flexneri_Oss_up)])
rownames(Shigella_flexneri_expr_Oss_up) =Shigella_flexneri_Oss_up[,1]

################SHIF_Oss_up2###############
Shigella_flexneri_Oss2_up=read.csv("SHIF_Ns_1Ns_2Oss_1Oss_2Oss_3_UP_Sp_1Sp_2_DOWN.csv", sep=",")
Shigella_flexneri_expr_Oss2_up =as.data.frame(Shigella_flexneri_Oss2_up[,2:ncol(Shigella_flexneri_Oss2_up)])
rownames(Shigella_flexneri_expr_Oss2_up) =Shigella_flexneri_Oss2_up[,1]

################SHIF_Oss_up_Tm_up###############
Shigella_flexneri_Oss_up_Tm_up=read.csv("SHIF_Oss_1Oss_2Oss_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Shigella_flexneri_expr_Oss_up_Tm_up =as.data.frame(Shigella_flexneri_Oss_up_Tm_up[,2:ncol(Shigella_flexneri_Oss_up_Tm_up)])
rownames(Shigella_flexneri_expr_Oss_up_Tm_up) =Shigella_flexneri_Oss_up_Tm_up[,1]

################SPYO_Oss_up###############
Streptococcus_pyogenes_Oss_up=read.csv("SPYO_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Streptococcus_pyogenes_expr_Oss_up =as.data.frame(Streptococcus_pyogenes_Oss_up[,2:ncol(Streptococcus_pyogenes_Oss_up)])
rownames(Streptococcus_pyogenes_expr_Oss_up) =Streptococcus_pyogenes_Oss_up[,1]

################SPYO_Oss_up_Ns_up###############
Streptococcus_pyogenes_Oss_up_Ns_up=read.csv("SPYO_Ns_1Ns_2Ns_3Oss_1Oss_2Oss_3_UP.csv", sep=",")
Streptococcus_pyogenes_expr_Oss_up_Ns_up =as.data.frame(Streptococcus_pyogenes_Oss_up_Ns_up[,2:ncol(Streptococcus_pyogenes_Oss_up_Ns_up)])
rownames(Streptococcus_pyogenes_expr_Oss_up_Ns_up) =Streptococcus_pyogenes_Oss_up_Ns_up[,1]

################STAGA_Oss_up###############
Streptococcus_agalactiae_Oss_up=read.csv("STAGA_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Streptococcus_agalactiae_expr_Oss_up =as.data.frame(Streptococcus_agalactiae_Oss_up[,2:ncol(Streptococcus_agalactiae_Oss_up)])
rownames(Streptococcus_agalactiae_expr_Oss_up) =Streptococcus_agalactiae_Oss_up[,1]

################STAGA_Oss_up_Nd_down###############
Streptococcus_agalactiae_Oss_up_Nd_down=read.csv("STAGA_Oss_1Oss_2Oss_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Oss_up_Nd_down =as.data.frame(Streptococcus_agalactiae_Oss_up_Nd_down[,2:ncol(Streptococcus_agalactiae_Oss_up_Nd_down)])
rownames(Streptococcus_agalactiae_expr_Oss_up_Nd_down) =Streptococcus_agalactiae_Oss_up_Nd_down[,1]

################STAGA_Oss_up_Vic_down###############
Streptococcus_agalactiae_Oss_up_Vic_down=read.csv("STAGA_Oss_1Oss_2Oss_3_UP_Vic_1Vic_2Vic_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Oss_up_Vic_down =as.data.frame(Streptococcus_agalactiae_Oss_up_Vic_down[,2:ncol(Streptococcus_agalactiae_Oss_up_Vic_down)])
rownames(Streptococcus_agalactiae_expr_Oss_up_Vic_down) =Streptococcus_agalactiae_Oss_up_Vic_down[,1]

################UPEC_Oss_up###############
Escherichia_coli_UPEC_Oss_up=read.csv("UPEC_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Escherichia_coli_UPEC_expr_Oss_up =as.data.frame(Escherichia_coli_UPEC_Oss_up[,2:ncol(Escherichia_coli_UPEC_Oss_up)])
rownames(Escherichia_coli_UPEC_expr_Oss_up) =Escherichia_coli_UPEC_Oss_up[,1]

################Vibrio_Oss_up###############
Vibrio_cholerae_Oss_up=read.csv("Vibrio_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Oss_up =as.data.frame(Vibrio_cholerae_Oss_up[,2:ncol(Vibrio_cholerae_Oss_up)])
rownames(Vibrio_cholerae_expr_Oss_up) =Vibrio_cholerae_Oss_up[,1]

################YPSTB_Oss_up###############
Yersinia_pseudotuberculosis_Oss_up=read.csv("YPSTB_Oss_1Oss_2Oss_3_UP.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Oss_up =as.data.frame(Yersinia_pseudotuberculosis_Oss_up[,2:ncol(Yersinia_pseudotuberculosis_Oss_up)])
rownames(Yersinia_pseudotuberculosis_expr_Oss_up) =Yersinia_pseudotuberculosis_Oss_up[,1]

###########EPEC_Oss_up_Nd_up_createdbyCansu############
Escherichia_coli_EPEC_Nd_up_Oss_up=read.csv("EPEC_Nd_1Nd_2Nd_3Oss_1Oss_2Oss_3_UP.csv",  sep=",")
Escherichia_coli_EPEC_expr_Nd_up_Oss_up=as.data.frame(Escherichia_coli_EPEC_Nd_up_Oss_up[,2:ncol(Escherichia_coli_EPEC_Nd_up_Oss_up)])
rownames(Escherichia_coli_EPEC_expr_Nd_up_Oss_up) =Escherichia_coli_EPEC_Nd_up_Oss_up[,1]

###########Pseudo_Oss_up_Nd_up_createdbyCansu############
Pseudomonas_aeruginosa_Nd_up_Oss_up=read.csv("PSEUDO_Nd_1Nd_2Nd_3Oss_1Oss_2Oss_3_UP.csv",  sep=",")
Pseudomonas_aeruginosa_expr_Nd_up_Oss_up=as.data.frame(Pseudomonas_aeruginosa_Nd_up_Oss_up[,2:ncol(Pseudomonas_aeruginosa_Nd_up_Oss_up)])
rownames(Pseudomonas_aeruginosa_expr_Nd_up_Oss_up) =Pseudomonas_aeruginosa_Nd_up_Oss_up[,1]

##################ETEC_Oss_up_Mig_up_##########################
Escherichia_coli_ETEC_Mig_up_Oss_up=read.csv("ETEC_Mig_1Mig_2Mig_3Oss_1Oss_2Oss_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Mig_up_Oss_up=as.data.frame(Escherichia_coli_ETEC_Mig_up_Oss_up[,2:ncol(Escherichia_coli_ETEC_Mig_up_Oss_up)])
rownames(Escherichia_coli_ETEC_expr_Mig_up_Oss_up) =Escherichia_coli_ETEC_Mig_up_Oss_up[,1]

##################BURK_Oss_down##########################
Burkholderia_pseudomallei_Oss_down=read.csv("BURK_Li_1Li_2Nd_1_UP_Oss_1Oss_2Oss_3_DOWN.csv", sep=",")
Burkholderia_pseudomallei_expr_Oss_down=as.data.frame(Burkholderia_pseudomallei_Oss_down[,2:ncol(Burkholderia_pseudomallei_Oss_down)])
rownames(Burkholderia_pseudomallei_expr_Oss_down) =Burkholderia_pseudomallei_Oss_down[,1]

################PSEUDO_Oss_down###############
Pseudomonas_aeruginosa_Oss_down=read.csv("PSEUDO_Nd_1Nd_3_UP_Oss_1Oss_2Oss_3_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_Oss_down =as.data.frame(Pseudomonas_aeruginosa_Oss_down[,2:ncol(Pseudomonas_aeruginosa_Oss_down)])
rownames(Pseudomonas_aeruginosa_expr_Oss_down) =Pseudomonas_aeruginosa_Oss_down[,1]

##################Oxidative_stress##########################

##################ACIB_Oxs_up##########################
Acinetobacter_baumannii_Oxs_up=read.csv("ACIB_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Acinetobacter_baumannii_expr_Oxs_up=as.data.frame(Acinetobacter_baumannii_Oxs_up[,2:ncol(Acinetobacter_baumannii_Oxs_up)])
rownames(Acinetobacter_baumannii_expr_Oxs_up) =Acinetobacter_baumannii_Oxs_up[,1]

##################AGGA_Oxs_up##########################
Aggregatibacter_actinomycetemcomitans_Oxs_up=read.csv("AGGA_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Oxs_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Oxs_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Oxs_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Oxs_up) =Aggregatibacter_actinomycetemcomitans_Oxs_up[,1]

##################BURK_Oxs_up##########################
Burkholderia_pseudomallei_Oxs_up=read.csv("BURK_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Oxs_up=as.data.frame(Burkholderia_pseudomallei_Oxs_up[,2:ncol(Burkholderia_pseudomallei_Oxs_up)])
rownames(Burkholderia_pseudomallei_expr_Oxs_up) =Burkholderia_pseudomallei_Oxs_up[,1]

##################CAMPY_Oxs1_up##########################
Campylobacter_jejuni_Oxs1_up=read.csv("Campy_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Campylobacter_jejuni_expr_Oxs1_up=as.data.frame(Campylobacter_jejuni_Oxs1_up[,2:ncol(Campylobacter_jejuni_Oxs1_up)])
rownames(Campylobacter_jejuni_expr_Oxs1_up) =Campylobacter_jejuni_Oxs1_up[,1]

##################CAMPY_Oxs2_up##########################
Campylobacter_jejuni_Oxs2_up=read.csv("Campy_Li_1Li_2Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Campylobacter_jejuni_expr_Oxs2_up=as.data.frame(Campylobacter_jejuni_Oxs2_up[,2:ncol(Campylobacter_jejuni_Oxs2_up)])
rownames(Campylobacter_jejuni_expr_Oxs2_up) =Campylobacter_jejuni_Oxs2_up[,1]

##################ENTFA_Oxs_up_Vic_down##########################
Enterococcus_faecalis_Oxs_up_Vic_down=read.csv("ENTFA_Oxs_1Oxs_2Oxs_3_UP_Vic_1Vic_2Vic_3_DOWN.csv", sep=",")
Enterococcus_faecalis_expr_Oxs_up_Vic_down=as.data.frame(Enterococcus_faecalis_Oxs_up_Vic_down[,2:ncol(Enterococcus_faecalis_Oxs_up_Vic_down)])
rownames(Enterococcus_faecalis_expr_Oxs_up_Vic_down) =Enterococcus_faecalis_Oxs_up_Vic_down[,1]

##################EPEC_Oxs1_up##########################
Escherichia_coli_EPEC_Oxs1_up=read.csv("EPEC_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Escherichia_coli_EPEC_expr_Oxs1_up=as.data.frame(Escherichia_coli_EPEC_Oxs1_up[,2:ncol(Escherichia_coli_EPEC_Oxs1_up)])
rownames(Escherichia_coli_EPEC_expr_Oxs1_up) =Escherichia_coli_EPEC_Oxs1_up[,1]

##################EPEC_Oxs2_up##########################
Escherichia_coli_EPEC_Oxs2_up=read.csv("EPEC_As_1Oxs_1Oxs_2Oxs_3_UP_Sp_3_DOWN.csv", sep=",")
Escherichia_coli_EPEC_expr_Oxs2_up=as.data.frame(Escherichia_coli_EPEC_Oxs2_up[,2:ncol(Escherichia_coli_EPEC_Oxs2_up)])
rownames(Escherichia_coli_EPEC_expr_Oxs2_up) =Escherichia_coli_EPEC_Oxs2_up[,1]

##################ETEC_Oxs_up##########################
Escherichia_coli_ETEC_Oxs_up=read.csv("ETEC_Oxs_1Oxs_2Oxs_3_UPpink.csv", sep=",")
Escherichia_coli_ETEC_expr_Oxs_up=as.data.frame(Escherichia_coli_ETEC_Oxs_up[,2:ncol(Escherichia_coli_ETEC_Oxs_up)])
rownames(Escherichia_coli_ETEC_expr_Oxs_up) =Escherichia_coli_ETEC_Oxs_up[,1]

##################ETEC_Oxs2_up##########################
Escherichia_coli_ETEC_Oxs2_up=read.csv("ETEC_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Oxs2_up=as.data.frame(Escherichia_coli_ETEC_Oxs2_up[,2:ncol(Escherichia_coli_ETEC_Oxs2_up)])
rownames(Escherichia_coli_ETEC_expr_Oxs2_up) =Escherichia_coli_ETEC_Oxs2_up[,1]


################HP_G27_Oxs_up###############
Helicobacter_pylori_G27_Oxs_up=read.csv("HP_G27_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Helicobacter_pylori_G27_expr_Oxs_up =as.data.frame(Helicobacter_pylori_G27_Oxs_up[,2:ncol(Helicobacter_pylori_G27_Oxs_up)])
rownames(Helicobacter_pylori_G27_expr_Oxs_up) =Helicobacter_pylori_G27_Oxs_up[,1]

################HP_J99_Oxs_up###############
Helicobacter_pylori_J99_Oxs_up=read.csv("HP_J99_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Helicobacter_pylori_J99_expr_Oxs_up =as.data.frame(Helicobacter_pylori_J99_Oxs_up[,2:ncol(Helicobacter_pylori_J99_Oxs_up)])
rownames(Helicobacter_pylori_J99_expr_Oxs_up) =Helicobacter_pylori_J99_Oxs_up[,1]

##################KLEBS_Oxs_up##########################
Klebsiella_pneumoniae_Oxs_up=read.csv("KLEBS_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Oxs_up=as.data.frame(Klebsiella_pneumoniae_Oxs_up[,2:ncol(Klebsiella_pneumoniae_Oxs_up)])
rownames(Klebsiella_pneumoniae_expr_Oxs_up) =Klebsiella_pneumoniae_Oxs_up[,1]

################Listeria_Oxs_up_As_up###############
Listeria_monocytogenes_Oxs_up_As_up=read.csv("Listeria_As_1As_2As_3Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Listeria_monocytogenes_expr_Oxs_up_As_up=as.data.frame(Listeria_monocytogenes_Oxs_up_As_up[,2:ncol(Listeria_monocytogenes_Oxs_up_As_up)])
rownames(Listeria_monocytogenes_expr_Oxs_up_As_up) =Listeria_monocytogenes_Oxs_up_As_up[,1]

################MRSA252_Oxs_up###############
Staphylococcus_aureus_MRSA252_Oxs_up=read.csv("MRSA252_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Oxs_up =as.data.frame(Staphylococcus_aureus_MRSA252_Oxs_up[,2:ncol(Staphylococcus_aureus_MRSA252_Oxs_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Oxs_up) =Staphylococcus_aureus_MRSA252_Oxs_up[,1]

################MRSA252_Oxs_up_Vic_up###############
Staphylococcus_aureus_MRSA252_Oxs_up_Vic_up=read.csv("MRSA252_Oxs_1Oxs_2Oxs_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Oxs_up_Vic_up =as.data.frame(Staphylococcus_aureus_MRSA252_Oxs_up_Vic_up[,2:ncol(Staphylococcus_aureus_MRSA252_Oxs_up_Vic_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Oxs_up_Vic_up) =Staphylococcus_aureus_MRSA252_Oxs_up_Vic_up[,1]

################MSSA476_Oxs_up###############
Staphylococcus_aureus_MSSA476_Oxs_up=read.csv("MSSA476_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Staphylococcus_aureus_MSSA476_expr_Oxs_up =as.data.frame(Staphylococcus_aureus_MSSA476_Oxs_up[,2:ncol(Staphylococcus_aureus_MSSA476_Oxs_up)])
rownames(Staphylococcus_aureus_MSSA476_expr_Oxs_up) =Staphylococcus_aureus_MSSA476_Oxs_up[,1]

################PSEUDO_Oxs_up###############
Pseudomonas_aeruginosa_Oxs_up=read.csv("PSEUDO_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Oxs_up =as.data.frame(Pseudomonas_aeruginosa_Oxs_up[,2:ncol(Pseudomonas_aeruginosa_Oxs_up)])
rownames(Pseudomonas_aeruginosa_expr_Oxs_up) =Pseudomonas_aeruginosa_Oxs_up[,1]

################SALMT_Oxs_up###############
Salmonella_enterica_Oxs_up=read.csv("SALMT_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Salmonella_enterica_expr_Oxs_up =as.data.frame(Salmonella_enterica_Oxs_up[,2:ncol(Salmonella_enterica_Oxs_up)])
rownames(Salmonella_enterica_expr_Oxs_up) =Salmonella_enterica_Oxs_up[,1]

################SHIF_Oxs_up###############
Shigella_flexneri_Oxs_up=read.csv("SHIF_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Shigella_flexneri_expr_Oxs_up =as.data.frame(Shigella_flexneri_Oxs_up[,2:ncol(Shigella_flexneri_Oxs_up)])
rownames(Shigella_flexneri_expr_Oxs_up) =Shigella_flexneri_Oxs_up[,1]

################SPYO_Oxs1_up###############
Streptococcus_pyogenes_Oxs_up=read.csv("SPYO_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Streptococcus_pyogenes_expr_Oxs_up =as.data.frame(Streptococcus_pyogenes_Oxs_up[,2:ncol(Streptococcus_pyogenes_Oxs_up)])
rownames(Streptococcus_pyogenes_expr_Oxs_up) =Streptococcus_pyogenes_Oxs_up[,1]

################SPYO_Oxs2_up###############
Streptococcus_pyogenes_Oxs2_up=read.csv("SPYO_Oxs_1Oxs_2Oxs_3_UPsalmon.csv", sep=",")
Streptococcus_pyogenes_expr_Oxs2_up =as.data.frame(Streptococcus_pyogenes_Oxs2_up[,2:ncol(Streptococcus_pyogenes_Oxs2_up)])
rownames(Streptococcus_pyogenes_expr_Oxs2_up) =Streptococcus_pyogenes_Oxs2_up[,1]

################UPEC_Oxs_up###############
Escherichia_coli_UPEC_Oxs_up=read.csv("UPEC_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Escherichia_coli_UPEC_expr_Oxs_up =as.data.frame(Escherichia_coli_UPEC_Oxs_up[,2:ncol(Escherichia_coli_UPEC_Oxs_up)])
rownames(Escherichia_coli_UPEC_expr_Oxs_up) =Escherichia_coli_UPEC_Oxs_up[,1]

################YPSTB_Oxs_up###############
Yersinia_pseudotuberculosis_Oxs_up=read.csv("YPSTB_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Oxs_up =as.data.frame(Yersinia_pseudotuberculosis_Oxs_up[,2:ncol(Yersinia_pseudotuberculosis_Oxs_up)])
rownames(Yersinia_pseudotuberculosis_expr_Oxs_up) =Yersinia_pseudotuberculosis_Oxs_up[,1]

###########ENTFA_Oxs_up_Nd_up_Vic_down_createdbyCansu############
Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down=read.csv("ENTFA_Nd_1Nd_2Nd_3Oxs_1Oxs_2Oxs_3_UP_Vic_1Vic_2Vic_3_DOWN.csv",  sep=",")
Enterococcus_faecalis_expr_Nd_up_Oxs_up_Vic_down=as.data.frame(Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down[,2:ncol(Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down)])
rownames(Enterococcus_faecalis_expr_Nd_up_Oxs_up_Vic_down) =Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down[,1]

###########KLEBS_Oxs_up_Nd_up_createdbyCansu############
Klebsiella_pneumoniae_Nd_up_Oxs_up=read.csv("KLEBS_Nd_1Nd_2Nd_3Oxs_1Oxs_2Oxs_3_UP.csv",  sep=",")
Klebsiella_pneumoniae_expr_Nd_up_Oxs_up=as.data.frame(Klebsiella_pneumoniae_Nd_up_Oxs_up[,2:ncol(Klebsiella_pneumoniae_Nd_up_Oxs_up)])
rownames(Klebsiella_pneumoniae_expr_Nd_up_Oxs_up) =Klebsiella_pneumoniae_Nd_up_Oxs_up[,1]

###########SALMT_Oxs_up_Nd_up_createdbyCansu############
Salmonella_enterica_Nd_up_Oxs_up=read.csv("SALMT_Nd_1Nd_2Nd_3Oxs_1Oxs_2Oxs_3_UP.csv",  sep=",")
Salmonella_enterica_expr_Nd_up_Oxs_up=as.data.frame(Salmonella_enterica_Nd_up_Oxs_up[,2:ncol(Salmonella_enterica_Nd_up_Oxs_up)])
rownames(Salmonella_enterica_expr_Nd_up_Oxs_up) =Salmonella_enterica_Nd_up_Oxs_up[,1]

###########SHIF_Oxs_up_Nd_up_createdbyCansu############
Shigella_flexneri_Nd_up_Oxs_up=read.csv("SHIF_Nd_1Nd_2Nd_3Oxs_1Oxs_2Oxs_3_UP.csv",  sep=",")
Shigella_flexneri_expr_Nd_up_Oxs_up=as.data.frame(Shigella_flexneri_Nd_up_Oxs_up[,2:ncol(Shigella_flexneri_Nd_up_Oxs_up)])
rownames(Shigella_flexneri_expr_Nd_up_Oxs_up) =Shigella_flexneri_Nd_up_Oxs_up[,1]

##################ACIB_Oxs_up_Li_up##########################
Acinetobacter_baumannii_Li_up_Oxs_up=read.csv("ACIB_Li_1Li_2Li_3Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Acinetobacter_baumannii_expr_Li_up_Oxs_up=as.data.frame(Acinetobacter_baumannii_Li_up_Oxs_up[,2:ncol(Acinetobacter_baumannii_Li_up_Oxs_up)])
rownames(Acinetobacter_baumannii_expr_Li_up_Oxs_up) =Acinetobacter_baumannii_Li_up_Oxs_up[,1]


#################SALMT_Oxs_up_Li_up##########################
Salmonella_enterica_Li_up_Oxs_up=read.csv("SALMT_Li_1Li_2Li_3Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Salmonella_enterica_expr_Li_up_Oxs_up=as.data.frame(Salmonella_enterica_Li_up_Oxs_up[,2:ncol(Salmonella_enterica_Li_up_Oxs_up)])
rownames(Salmonella_enterica_expr_Li_up_Oxs_up) =Salmonella_enterica_Li_up_Oxs_up[,1]


#################SHIF_Oxs_up_Li_up##########################
Shigella_flexneri_Li_up_Oxs_up=read.csv("SHIF_Li_1Li_2Li_3Ns_1Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Shigella_flexneri_expr_Li_up_Oxs_up=as.data.frame(Shigella_flexneri_Li_up_Oxs_up[,2:ncol(Shigella_flexneri_Li_up_Oxs_up)])
rownames(Shigella_flexneri_expr_Li_up_Oxs_up) =Shigella_flexneri_Li_up_Oxs_up[,1]

##################PSEUDO_Oxs_up_Mig_up##########################
Pseudomonas_aeruginosa_Mig_up_Oxs_up=read.csv("PSEUDO_Mig_1Mig_2Mig_3Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Mig_up_Oxs_up=as.data.frame(Pseudomonas_aeruginosa_Mig_up_Oxs_up[,2:ncol(Pseudomonas_aeruginosa_Mig_up_Oxs_up)])
rownames(Pseudomonas_aeruginosa_expr_Mig_up_Oxs_up) =Pseudomonas_aeruginosa_Mig_up_Oxs_up[,1]

################NMEN_Oxs_down_Tm_up_###############
Neisseria_meningitidis_Tm_up_Oxs_down=read.csv("NMEN_Tm_1Tm_2Tm_3_UP_Oxs_1Oxs_2Oxs_3_DOWN.csv", sep=",")
Neisseria_meningitidis_expr_Tm_up_Oxs_down=as.data.frame(Neisseria_meningitidis_Tm_up_Oxs_down[,2:ncol(Neisseria_meningitidis_Tm_up_Oxs_down)])
rownames(Neisseria_meningitidis_expr_Tm_up_Oxs_down) =Neisseria_meningitidis_Tm_up_Oxs_down[,1]

##################ACHX_Oxs_down##########################
Achromobacter_xylosoxidans_Oxs_down=read.csv("ACHX_Oxs_1Oxs_2Oxs_3_DOWN.csv", sep=",")
Achromobacter_xylosoxidans_expr_Oxs_down=as.data.frame(Achromobacter_xylosoxidans_Oxs_down[,2:ncol(Achromobacter_xylosoxidans_Oxs_down)])
rownames(Achromobacter_xylosoxidans_expr_Oxs_down) =Achromobacter_xylosoxidans_Oxs_down[,1]

################STAGA_Oxs_up###############
Streptococcus_agalactiae_Oxs_up=read.csv("STAGA_Nd_1Nd_3Oxs_1Oxs_2Oxs_3Tm_3_UP.csv", sep=",")
Streptococcus_agalactiae_expr_Oxs_up =as.data.frame(Streptococcus_agalactiae_Oxs_up[,2:ncol(Streptococcus_agalactiae_Oxs_up)])
rownames(Streptococcus_agalactiae_expr_Oxs_up)=Streptococcus_agalactiae_Oxs_up[,1]

################STAGA_Oxs_up2###############
Streptococcus_agalactiae_Oxs2_up=read.csv("STAGA_Nd_3Oxs_1Oxs_2Oxs_3_UP_Vic_1Vic_2_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Oxs2_up =as.data.frame(Streptococcus_agalactiae_Oxs2_up[,2:ncol(Streptococcus_agalactiae_Oxs2_up)])
rownames(Streptococcus_agalactiae_expr_Oxs2_up) =Streptococcus_agalactiae_Oxs2_up[,1]


##################Temperature_stress##########################

##################ACHX_Tm1_up##########################
Achromobacter_xylosoxidans_Tm1_up=read.csv("ACHX_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Achromobacter_xylosoxidans_expr_Tm1_up=as.data.frame(Achromobacter_xylosoxidans_Tm1_up[,2:ncol(Achromobacter_xylosoxidans_Tm1_up)])
rownames(Achromobacter_xylosoxidans_expr_Tm1_up) =Achromobacter_xylosoxidans_Tm1_up[,1]

##################ACHX_Tm2_up##########################
Achromobacter_xylosoxidans_Tm2_up=read.csv("ACHX_Tm_1Tm_2Tm_3_UP_Ctrl_1_DOWN.csv", sep=",")
Achromobacter_xylosoxidans_expr_Tm2_up=as.data.frame(Achromobacter_xylosoxidans_Tm2_up[,2:ncol(Achromobacter_xylosoxidans_Tm2_up)])
rownames(Achromobacter_xylosoxidans_expr_Tm2_up) =Achromobacter_xylosoxidans_Tm2_up[,1]

##################ACHX_Tm3_up##########################
Achromobacter_xylosoxidans_Tm3_up=read.csv("ACHX_Oss_2Tm_1Tm_2Tm_3_UP.csv", sep=",")
Achromobacter_xylosoxidans_expr_Tm3_up=as.data.frame(Achromobacter_xylosoxidans_Tm3_up[,2:ncol(Achromobacter_xylosoxidans_Tm3_up)])
rownames(Achromobacter_xylosoxidans_expr_Tm3_up) =Achromobacter_xylosoxidans_Tm3_up[,1]
##################AGGA_Tm_up##########################
Aggregatibacter_actinomycetemcomitans_Tm_up=read.csv("AGGA_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Tm_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Tm_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Tm_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Tm_up) =Aggregatibacter_actinomycetemcomitans_Tm_up[,1]

################BBURG_Tm_up###############
Borrelia_burgdorferi_Tm_up=read.csv("BBURG_Tm_1Tm_2Tm_3_UP_Li_1Li_2_DOWN.csv", sep=",")
Borrelia_burgdorferi_expr_Tm_up =as.data.frame(Borrelia_burgdorferi_Tm_up[,2:ncol(Borrelia_burgdorferi_Tm_up)])
rownames(Borrelia_burgdorferi_expr_Tm_up) =Borrelia_burgdorferi_Tm_up[,1]

##################BURK_Tm_up##########################
Burkholderia_pseudomallei_Tm_up=read.csv("BURK_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Tm_up=as.data.frame(Burkholderia_pseudomallei_Tm_up[,2:ncol(Burkholderia_pseudomallei_Tm_up)])
rownames(Burkholderia_pseudomallei_expr_Tm_up) =Burkholderia_pseudomallei_Tm_up[,1]

##################ENTFA_Tm_up##########################
Enterococcus_faecalis_Tm_up=read.csv("ENTFA_As_1As_2Tm_1Tm_2Tm_3_UP.csv", sep=",")
Enterococcus_faecalis_expr_Tm_up=as.data.frame(Enterococcus_faecalis_Tm_up[,2:ncol(Enterococcus_faecalis_Tm_up)])
rownames(Enterococcus_faecalis_expr_Tm_up) =Enterococcus_faecalis_Tm_up[,1]

##################EPEC_Tm_up##########################
Escherichia_coli_EPEC_Tm_up=read.csv("EPEC_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Escherichia_coli_EPEC_expr_Tm_up=as.data.frame(Escherichia_coli_EPEC_Tm_up[,2:ncol(Escherichia_coli_EPEC_Tm_up)])
rownames(Escherichia_coli_EPEC_expr_Tm_up) =Escherichia_coli_EPEC_Tm_up[,1]

##################ETEC_Tm_up##########################
Escherichia_coli_ETEC_Tm_up=read.csv("ETEC_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Tm_up=as.data.frame(Escherichia_coli_ETEC_Tm_up[,2:ncol(Escherichia_coli_ETEC_Tm_up)])
rownames(Escherichia_coli_ETEC_expr_Tm_up) =Escherichia_coli_ETEC_Tm_up[,1]

##################ETEC_Tm_up_Li_up##########################
Escherichia_coli_ETEC_Tm_up_Li_up=read.csv("ETEC_Li_1Li_2Li_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Tm_up_Li_up=as.data.frame(Escherichia_coli_ETEC_Tm_up_Li_up[,2:ncol(Escherichia_coli_ETEC_Tm_up_Li_up)])
rownames(Escherichia_coli_ETEC_expr_Tm_up_Li_up) =Escherichia_coli_ETEC_Tm_up_Li_up[,1]

##################FRAT_Tm_up##########################
Francisella_tularensis_Tm_up=read.csv("FRAT_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Francisella_tularensis_expr_Tm_up=as.data.frame(Francisella_tularensis_Tm_up[,2:ncol(Francisella_tularensis_Tm_up)])
rownames(Francisella_tularensis_expr_Tm_up) =Francisella_tularensis_Tm_up[,1]

##################HINF_Tm1_up##########################
Haemophilus_influenzae_Tm1_up=read.csv("HINF_Sp_2Tm_1Tm_2Tm_3_UP.csv", sep=",")
Haemophilus_influenzae_expr_Tm1_up=as.data.frame(Haemophilus_influenzae_Tm1_up[,2:ncol(Haemophilus_influenzae_Tm1_up)])
rownames(Haemophilus_influenzae_expr_Tm1_up) =Haemophilus_influenzae_Tm1_up[,1]

##################HINF_Tm2_up##########################
Haemophilus_influenzae_Tm2_up=read.csv("HINF_Tm_1Tm_2Tm_3_UP_Mig_3_DOWN.csv", sep=",")
Haemophilus_influenzae_expr_Tm2_up=as.data.frame(Haemophilus_influenzae_Tm2_up[,2:ncol(Haemophilus_influenzae_Tm2_up)])
rownames(Haemophilus_influenzae_expr_Tm2_up) =Haemophilus_influenzae_Tm2_up[,1]

##################KLEBS_Tm_up##########################
Klebsiella_pneumoniae_Tm_up=read.csv("KLEBS_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Tm_up=as.data.frame(Klebsiella_pneumoniae_Tm_up[,2:ncol(Klebsiella_pneumoniae_Tm_up)])
rownames(Klebsiella_pneumoniae_expr_Tm_up) =Klebsiella_pneumoniae_Tm_up[,1]

##################Listeria_Tm1_up##########################
Listeria_monocytogenes_Tm1_up=read.csv("Listeria_Ns_2Ns_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Listeria_monocytogenes_expr_Tm1_up=as.data.frame(Listeria_monocytogenes_Tm1_up[,2:ncol(Listeria_monocytogenes_Tm1_up)])
rownames(Listeria_monocytogenes_expr_Tm1_up) =Listeria_monocytogenes_Tm1_up[,1]

##################Listeria_Tm2_up##########################
Listeria_monocytogenes_Tm2_up=read.csv("Listeria_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Listeria_monocytogenes_expr_Tm2_up=as.data.frame(Listeria_monocytogenes_Tm2_up[,2:ncol(Listeria_monocytogenes_Tm2_up)])
rownames(Listeria_monocytogenes_expr_Tm2_up) =Listeria_monocytogenes_Tm2_up[,1]

################MRSA252_Tm_up###############
Staphylococcus_aureus_MRSA252_Tm_up=read.csv("MRSA252_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Tm_up =as.data.frame(Staphylococcus_aureus_MRSA252_Tm_up[,2:ncol(Staphylococcus_aureus_MRSA252_Tm_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Tm_up) =Staphylococcus_aureus_MRSA252_Tm_up[,1]

##################MSSA476_Tm_up##########################
Staphylococcus_aureus_MSSA476_Tm_up=read.csv("MSSA476_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Staphylococcus_aureus_MSSA476_expr_Tm_up=as.data.frame(Staphylococcus_aureus_MSSA476_Tm_up[,2:ncol(Staphylococcus_aureus_MSSA476_Tm_up)])
rownames(Staphylococcus_aureus_MSSA476_expr_Tm_up) =Staphylococcus_aureus_MSSA476_Tm_up[,1]

################NMEN_Tm_up###############
Neisseria_meningitidis_Tm_up=read.csv("NMEN_Mig_2Nd_2Nd_3Tm_1Tm_2Tm_3_UP_Oss_3_DOWN.csv", sep=",")
Neisseria_meningitidis_expr_Tm_up =as.data.frame(Neisseria_meningitidis_Tm_up[,2:ncol(Neisseria_meningitidis_Tm_up)])
rownames(Neisseria_meningitidis_expr_Tm_up) =Neisseria_meningitidis_Tm_up[,1]

################NMEN_Tm_up_Oxs_down###############
Neisseria_meningitidis_Tm_up_Oxs_down=read.csv("NMEN_Tm_1Tm_2Tm_3_UP_Oxs_1Oxs_2Oxs_3_DOWN.csv", sep=",")
Neisseria_meningitidis_expr_Tm_up_Oxs_down=as.data.frame(Neisseria_meningitidis_Tm_up_Oxs_down[,2:ncol(Neisseria_meningitidis_Tm_up_Oxs_down)])
rownames(Neisseria_meningitidis_expr_Tm_up_Oxs_down) =Neisseria_meningitidis_Tm_up_Oxs_down[,1]

################SALMT_Tm_up###############
Salmonella_enterica_Tm_up=read.csv("SALMT_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Salmonella_enterica_expr_Tm_up =as.data.frame(Salmonella_enterica_Tm_up[,2:ncol(Salmonella_enterica_Tm_up)])
rownames(Salmonella_enterica_expr_Tm_up) =Salmonella_enterica_Tm_up[,1]

##################SEPI_Tm_up##########################
Staphylococcus_epidermidis_Tm_up=read.csv("SEPI_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Staphylococcus_epidermidis_expr_Tm_up=as.data.frame(Staphylococcus_epidermidis_Tm_up[,2:ncol(Staphylococcus_epidermidis_Tm_up)])
rownames(Staphylococcus_epidermidis_expr_Tm_up) =Staphylococcus_epidermidis_Tm_up[,1]

##################SEPI_Tm_down##########################
Staphylococcus_epidermidis_Tm_down=read.csv("SEPI_Bs_2_UP_Tm_1Tm_2Tm_3_DOWN.csv", sep=",")
Staphylococcus_epidermidis_expr_Tm_down=as.data.frame(Staphylococcus_epidermidis_Tm_down[,2:ncol(Staphylococcus_epidermidis_Tm_down)])
rownames(Staphylococcus_epidermidis_expr_Tm_down) =Staphylococcus_epidermidis_Tm_down[,1]

################SHIF_Tm_up###############
Shigella_flexneri_Tm_up=read.csv("SHIF_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Shigella_flexneri_expr_Tm_up =as.data.frame(Shigella_flexneri_Tm_up[,2:ncol(Shigella_flexneri_Tm_up)])
rownames(Shigella_flexneri_expr_Tm_up) =Shigella_flexneri_Tm_up[,1]

################SPYO_Tm_up_Nd_up###############
Streptococcus_pyogenes_Tm_up_Nd_up=read.csv("SPYO_Nd_1Nd_2Nd_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Streptococcus_pyogenes_expr_Tm_up_Nd_up =as.data.frame(Streptococcus_pyogenes_Tm_up_Nd_up[,2:ncol(Streptococcus_pyogenes_Tm_up_Nd_up)])
rownames(Streptococcus_pyogenes_expr_Tm_up_Nd_up) =Streptococcus_pyogenes_Tm_up_Nd_up[,1]

##################STAGA_Tm_up##########################
Streptococcus_agalactiae_Tm_up=read.csv("STAGA_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Streptococcus_agalactiae_expr_Tm_up=as.data.frame(Streptococcus_agalactiae_Tm_up[,2:ncol(Streptococcus_agalactiae_Tm_up)])
rownames(Streptococcus_agalactiae_expr_Tm_up) =Streptococcus_agalactiae_Tm_up[,1]

##################STAGA_Tm_up_Nd_down##########################
Streptococcus_agalactiae_Tm_up_Nd_down=read.csv("STAGA_Tm_1Tm_2Tm_3Vic_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Tm_up_Nd_down=as.data.frame(Streptococcus_agalactiae_Tm_up_Nd_down[,2:ncol(Streptococcus_agalactiae_Tm_up_Nd_down)])
rownames(Streptococcus_agalactiae_expr_Tm_up_Nd_down) =Streptococcus_agalactiae_Tm_up_Nd_down[,1]

##################STRPN_Tm_up_Vic_up##########################
Streptococcus_pneumoniae_Tm_up_Vic_up=read.csv("STRPN_Tm_1Tm_2Tm_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Streptococcus_pneumoniae_expr_Tm_up_Vic_up=as.data.frame(Streptococcus_pneumoniae_Tm_up_Vic_up[,2:ncol(Streptococcus_pneumoniae_Tm_up_Vic_up)])
rownames(Streptococcus_pneumoniae_expr_Tm_up_Vic_up) =Streptococcus_pneumoniae_Tm_up_Vic_up[,1]

################UPEC_Tm_up###############
Escherichia_coli_UPEC_Tm_up=read.csv("UPEC_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Escherichia_coli_UPEC_expr_Tm_up =as.data.frame(Escherichia_coli_UPEC_Tm_up[,2:ncol(Escherichia_coli_UPEC_Tm_up)])
rownames(Escherichia_coli_UPEC_expr_Tm_up) =Escherichia_coli_UPEC_Tm_up[,1]

################Vibrio_Tm_up###############
Vibrio_cholerae_Tm_up=read.csv("Vibrio_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Tm_up =as.data.frame(Vibrio_cholerae_Tm_up[,2:ncol(Vibrio_cholerae_Tm_up)])
rownames(Vibrio_cholerae_expr_Tm_up) =Vibrio_cholerae_Tm_up[,1]

################YPSTB_Tm_up###############
Yersinia_pseudotuberculosis_Tm_up=read.csv("YPSTB_Tm_1Tm_2Tm_3_UP.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Tm_up =as.data.frame(Yersinia_pseudotuberculosis_Tm_up[,2:ncol(Yersinia_pseudotuberculosis_Tm_up)])
rownames(Yersinia_pseudotuberculosis_expr_Tm_up) =Yersinia_pseudotuberculosis_Tm_up[,1]

#########PSEUDO_Tm_As#########
Pseudomonas_aeruginosa_As_Tm = read.csv("PSEUDO_As_1As_2As_3Tm_1Tm_2Tm_3_UP_Nd_1Nd_2_DOWN.csv", sep=",")
Pseudomonas_aeruginosa_expr_As_Tm =as.data.frame(Pseudomonas_aeruginosa_As_Tm[,2:ncol(Pseudomonas_aeruginosa_As_Tm)])
rownames(Pseudomonas_aeruginosa_expr_As_Tm) = Pseudomonas_aeruginosa_As_Tm[,1]

################BURK_Tm_Sp_createdbyCansu###############
Burkholderia_pseudomallei_Sp_Tm=read.csv("BURK_Sp_1Sp_2Sp_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Sp_Tm =as.data.frame(Burkholderia_pseudomallei_Sp_Tm[,2:ncol(Burkholderia_pseudomallei_Sp_Tm)])
rownames(Burkholderia_pseudomallei_expr_Sp_Tm) =Burkholderia_pseudomallei_Sp_Tm[,1]

################SEPI_Tm_up_Sp_up_createdbyCansu###############
Staphylococcus_epidermidis_Sp_up_Tm_up=read.csv("SEPI_Sp_1Sp_2Sp_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Staphylococcus_epidermidis_expr_Sp_up_Tm_up=as.data.frame(Staphylococcus_epidermidis_Sp_up_Tm_up[,2:ncol(Staphylococcus_epidermidis_Sp_up_Tm_up)])
rownames(Staphylococcus_epidermidis_expr_Sp_up_Tm_up) =Staphylococcus_epidermidis_Sp_up_Tm_up[,1]

###########SEPI_Tm_down_Nd_up_createdbyCansu############
Staphylococcus_epidermidis_Nd_up_Tm_down=read.csv("SEPI_Nd_1Nd_2Nd_3_UP_Tm_1Tm_2Tm_3_DOWN.csv",  sep=",")
Staphylococcus_epidermidis_expr_Nd_up_Tm_down=as.data.frame(Staphylococcus_epidermidis_Nd_up_Tm_down[,2:ncol(Staphylococcus_epidermidis_Nd_up_Tm_down)])
rownames(Staphylococcus_epidermidis_expr_Nd_up_Tm_down) =Staphylococcus_epidermidis_Nd_up_Tm_down[,1]

##################SEPI_Tm_down_Mig_up##########################
Staphylococcus_epidermidis_Mig_up_Tm_down=read.csv("SEPI_Mig_1Mig_2Mig_3Nd_1_UP_Tm_1Tm_2Tm_3_DOWN.csv", sep=",")
Staphylococcus_epidermidis_expr_Mig_up_Tm_down=as.data.frame(Staphylococcus_epidermidis_Mig_up_Tm_down[,2:ncol(Staphylococcus_epidermidis_Mig_up_Tm_down)])
rownames(Staphylococcus_epidermidis_expr_Mig_up_Tm_down) =Staphylococcus_epidermidis_Mig_up_Tm_down[,1]

##################SEPI_Tm_down_Ns_up##########################
Staphylococcus_epidermidis_Ns_up_Tm_down=read.csv("SEPI_Ns_1Ns_2Ns_3_UP_Tm_1Tm_2Tm_3_DOWN.csv", sep=",")
Staphylococcus_epidermidis_expr_Ns_up_Tm_down=as.data.frame(Staphylococcus_epidermidis_Ns_up_Tm_down[,2:ncol(Staphylococcus_epidermidis_Ns_up_Tm_down)])
rownames(Staphylococcus_epidermidis_expr_Ns_up_Tm_down) =Staphylococcus_epidermidis_Ns_up_Tm_down[,1]

##################SEPI_Tm_down_As_up##########################
Staphylococcus_epidermidis_As_up_Tm_down=read.csv("SEPI_As_1As_2As_3Li_1Li_3_UP_Sp_2Tm_1Tm_2Tm_3_DOWN.csv", sep=",")
Staphylococcus_epidermidis_expr_As_up_Tm_down=as.data.frame(Staphylococcus_epidermidis_As_up_Tm_down[,2:ncol(Staphylococcus_epidermidis_As_up_Tm_down)])
rownames(Staphylococcus_epidermidis_expr_As_up_Tm_down) =Staphylococcus_epidermidis_As_up_Tm_down[,1]

################KLEBS_Tm_up_Oss_up###############
Klebsiella_pneumoniae_Oss_up_Tm_up=read.csv("KLEBS_Oss_1Oss_2Oss_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Oss_up_Tm_up=as.data.frame(Klebsiella_pneumoniae_Oss_up_Tm_up[,2:ncol(Klebsiella_pneumoniae_Oss_up_Tm_up)])
rownames(Klebsiella_pneumoniae_expr_Oss_up_Tm_up) =Klebsiella_pneumoniae_Oss_up_Tm_up[,1]

################SHIF_Tm_up_Oss_up###############
Shigella_flexneri_Oss_up_Tm_up=read.csv("SHIF_Oss_1Oss_2Oss_3Tm_1Tm_2Tm_3_UP.csv", sep=",")
Shigella_flexneri_expr_Oss_up_Tm_up =as.data.frame(Shigella_flexneri_Oss_up_Tm_up[,2:ncol(Shigella_flexneri_Oss_up_Tm_up)])
rownames(Shigella_flexneri_expr_Oss_up_Tm_up) =Shigella_flexneri_Oss_up_Tm_up[,1]

################SHIF_Tm_up_Vic_up###############
Shigella_flexneri_Vic_up_Tm_up=read.csv("SHIF_Tm_1Tm_2Tm_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Shigella_flexneri_expr_Vic_up_Tm_up =as.data.frame(Shigella_flexneri_Vic_up_Tm_up[,2:ncol(Shigella_flexneri_Vic_up_Tm_up)])
rownames(Shigella_flexneri_expr_Vic_up_Tm_up) =Shigella_flexneri_Vic_up_Tm_up[,1]




##################Virulence_inducing_condition##########################

##################AGGA_Vic1_up##########################
Aggregatibacter_actinomycetemcomitans_Vic1_up=read.csv("AGGA_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Vic1_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Vic1_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Vic1_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Vic1_up) =Aggregatibacter_actinomycetemcomitans_Vic1_up[,1]

##################AGGA_Vic2_up##########################
Aggregatibacter_actinomycetemcomitans_Vic2_up=read.csv("AGGA_Bs_1Vic_1Vic_2Vic_3_UP.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Vic2_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Vic2_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Vic2_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Vic2_up) =Aggregatibacter_actinomycetemcomitans_Vic2_up[,1]

##################AGGA_Vic3_up##########################
Aggregatibacter_actinomycetemcomitans_Vic3_up=read.csv("AGGA_Vic_1Vic_2Vic_3_UP_Ns_2_DOWN.csv", sep=",")
Aggregatibacter_actinomycetemcomitans_expr_Vic3_up=as.data.frame(Aggregatibacter_actinomycetemcomitans_Vic3_up[,2:ncol(Aggregatibacter_actinomycetemcomitans_Vic3_up)])
rownames(Aggregatibacter_actinomycetemcomitans_expr_Vic3_up) =Aggregatibacter_actinomycetemcomitans_Vic3_up[,1]

##################BURK_Vic_up##########################
Burkholderia_pseudomallei_Vic_up=read.csv("BURK_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Burkholderia_pseudomallei_expr_Vic_up=as.data.frame(Burkholderia_pseudomallei_Vic_up[,2:ncol(Burkholderia_pseudomallei_Vic_up)])
rownames(Burkholderia_pseudomallei_expr_Vic_up) =Burkholderia_pseudomallei_Vic_up[,1]

##################EPEC_Vic_up##########################
Escherichia_coli_EPEC_Vic_up=read.csv("EPEC_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Escherichia_coli_EPEC_expr_Vic_up=as.data.frame(Escherichia_coli_EPEC_Vic_up[,2:ncol(Escherichia_coli_EPEC_Vic_up)])
rownames(Escherichia_coli_EPEC_expr_Vic_up) =Escherichia_coli_EPEC_Vic_up[,1]

##################ETEC_Vic1_up##########################
Escherichia_coli_ETEC_Vic1_up=read.csv("ETEC_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Vic1_up=as.data.frame(Escherichia_coli_ETEC_Vic1_up[,2:ncol(Escherichia_coli_ETEC_Vic1_up)])
rownames(Escherichia_coli_ETEC_expr_Vic1_up) =Escherichia_coli_ETEC_Vic1_up[,1]

##################ETEC_Vic2_up##########################
Escherichia_coli_ETEC_Vic2_up=read.csv("ETEC_As_1Vic_1Vic_2Vic_3_UP.csv", sep=",")
Escherichia_coli_ETEC_expr_Vic2_up=as.data.frame(Escherichia_coli_ETEC_Vic2_up[,2:ncol(Escherichia_coli_ETEC_Vic2_up)])
rownames(Escherichia_coli_ETEC_expr_Vic2_up) =Escherichia_coli_ETEC_Vic2_up[,1]

##################HINF_Vic_up##########################
Haemophilus_influenzae_Vic_up=read.csv("HINF_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Haemophilus_influenzae_expr_Vic_up=as.data.frame(Haemophilus_influenzae_Vic_up[,2:ncol(Haemophilus_influenzae_Vic_up)])
rownames(Haemophilus_influenzae_expr_Vic_up) =Haemophilus_influenzae_Vic_up[,1]

##################HINF_Vic_up_Mig_down##########################
Haemophilus_influenzae_Vic_up_Mig_down=read.csv("HINF_Vic_1Vic_2Vic_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Haemophilus_influenzae_expr_Vic_up_Mig_down=as.data.frame(Haemophilus_influenzae_Vic_up_Mig_down[,2:ncol(Haemophilus_influenzae_Vic_up_Mig_down)])
rownames(Haemophilus_influenzae_expr_Vic_up_Mig_down) =Haemophilus_influenzae_Vic_up_Mig_down[,1]

##################HINF_Vic_up_Mig_up##########################
Haemophilus_influenzae_Mig_up_Vic_up=read.csv("HINF_Mig_1Mig_2Mig_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Haemophilus_influenzae_expr_Mig_up_Vic_up=as.data.frame(Haemophilus_influenzae_Mig_up_Vic_up[,2:ncol(Haemophilus_influenzae_Mig_up_Vic_up)])
rownames(Haemophilus_influenzae_expr_Mig_up_Vic_up) =Haemophilus_influenzae_Mig_up_Vic_up[,1]

################HP_J99_Vic_up###############
Helicobacter_pylori_J99_Vic_up=read.csv("HP_J99_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Helicobacter_pylori_J99_expr_Vic_up =as.data.frame(Helicobacter_pylori_J99_Vic_up[,2:ncol(Helicobacter_pylori_J99_Vic_up)])
rownames(Helicobacter_pylori_J99_expr_Vic_up) =Helicobacter_pylori_J99_Vic_up[,1]

##################KLEBS_Vic_up##########################
Klebsiella_pneumoniae_Vic_up=read.csv("KLEBS_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Klebsiella_pneumoniae_expr_Vic_up=as.data.frame(Klebsiella_pneumoniae_Vic_up[,2:ncol(Klebsiella_pneumoniae_Vic_up)])
rownames(Klebsiella_pneumoniae_expr_Vic_up) =Klebsiella_pneumoniae_Vic_up[,1]

################MRSA252_Vic_up###############
Staphylococcus_aureus_MRSA252_Vic_up=read.csv("MRSA252_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Vic_up =as.data.frame(Staphylococcus_aureus_MRSA252_Vic_up[,2:ncol(Staphylococcus_aureus_MRSA252_Vic_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Vic_up) =Staphylococcus_aureus_MRSA252_Vic_up[,1]

##################MSSA476_Vic_up##########################
Staphylococcus_aureus_MSSA476_Vic_up=read.csv("MSSA476_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Staphylococcus_aureus_MSSA476_expr_Vic_up=as.data.frame(Staphylococcus_aureus_MSSA476_Vic_up[,2:ncol(Staphylococcus_aureus_MSSA476_Vic_up)])
rownames(Staphylococcus_aureus_MSSA476_expr_Vic_up) =Staphylococcus_aureus_MSSA476_Vic_up[,1]

################PSEUDO_Vic_up###############
Pseudomonas_aeruginosa_Vic_up=read.csv("PSEUDO_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Pseudomonas_aeruginosa_expr_Vic_up =as.data.frame(Pseudomonas_aeruginosa_Vic_up[,2:ncol(Pseudomonas_aeruginosa_Vic_up)])
rownames(Pseudomonas_aeruginosa_expr_Vic_up) =Pseudomonas_aeruginosa_Vic_up[,1]

################SALMT_Vic_up###############
Salmonella_enterica_Vic_up=read.csv("SALMT_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Salmonella_enterica_expr_Vic_up =as.data.frame(Salmonella_enterica_Vic_up[,2:ncol(Salmonella_enterica_Vic_up)])
rownames(Salmonella_enterica_expr_Vic_up) =Salmonella_enterica_Vic_up[,1]

################SHIF_Vic1_up###############
Shigella_flexneri_Vic1_up=read.csv("SHIF_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Shigella_flexneri_expr_Vic1_up =as.data.frame(Shigella_flexneri_Vic1_up[,2:ncol(Shigella_flexneri_Vic1_up)])
rownames(Shigella_flexneri_expr_Vic1_up) =Shigella_flexneri_Vic1_up[,1]

################SHIF_Vic2_up###############
Shigella_flexneri_Vic2_up=read.csv("SHIF_Vic_1Vic_2Vic_3_UPlightyellow.csv", sep=",")
Shigella_flexneri_expr_Vic2_up =as.data.frame(Shigella_flexneri_Vic2_up[,2:ncol(Shigella_flexneri_Vic2_up)])
rownames(Shigella_flexneri_expr_Vic2_up) =Shigella_flexneri_Vic2_up[,1]

################SHIF_Vic_up_Tm_up###############
Shigella_flexneri_Vic_up_Tm_up=read.csv("SHIF_Tm_1Tm_2Tm_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Shigella_flexneri_expr_Vic_up_Tm_up =as.data.frame(Shigella_flexneri_Vic_up_Tm_up[,2:ncol(Shigella_flexneri_Vic_up_Tm_up)])
rownames(Shigella_flexneri_expr_Vic_up_Tm_up) =Shigella_flexneri_Vic_up_Tm_up[,1]

##################STAGA_Vic_up_Nd_down##########################
Streptococcus_agalactiae_Vic_up_Nd_down=read.csv("STAGA_Vic_1Vic_2Vic_3_UP_Nd_1Nd_2Nd_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Vic_up_Nd_down=as.data.frame(Streptococcus_agalactiae_Vic_up_Nd_down[,2:ncol(Streptococcus_agalactiae_Vic_up_Nd_down)])
rownames(Streptococcus_agalactiae_expr_Vic_up_Nd_down) =Streptococcus_agalactiae_Vic_up_Nd_down[,1]

##################STRPN_Vic_up_Mig_down##########################
Streptococcus_pneumoniae_Vic_up_Mig_down=read.csv("STRPN_Vic_1Vic_2Vic_3_UP_Mig_1Mig_2Mig_3Nd_3_DOWN.csv", sep=",")
Streptococcus_pneumoniae_expr_Vic_up_Mig_down=as.data.frame(Streptococcus_pneumoniae_Vic_up_Mig_down[,2:ncol(Streptococcus_pneumoniae_Vic_up_Mig_down)])
rownames(Streptococcus_pneumoniae_expr_Vic_up_Mig_down) =Streptococcus_pneumoniae_Vic_up_Mig_down[,1]

################YPSTB_Vic1_up###############
Yersinia_pseudotuberculosis_Vic1_up=read.csv("YPSTB_Vic_1Vic_2Vic_3_UP.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Vic1_up =as.data.frame(Yersinia_pseudotuberculosis_Vic1_up[,2:ncol(Yersinia_pseudotuberculosis_Vic1_up)])
rownames(Yersinia_pseudotuberculosis_expr_Vic1_up) =Yersinia_pseudotuberculosis_Vic1_up[,1]

################YPSTB_Vic2_up###############
Yersinia_pseudotuberculosis_Vic2_up=read.csv("YPSTB_Vic_1Vic_2Vic_3_UPcyan.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Vic2_up =as.data.frame(Yersinia_pseudotuberculosis_Vic2_up[,2:ncol(Yersinia_pseudotuberculosis_Vic2_up)])
rownames(Yersinia_pseudotuberculosis_expr_Vic2_up) =Yersinia_pseudotuberculosis_Vic2_up[,1]

################YPSTB_Vic_up_Mig_down###############
Yersinia_pseudotuberculosis_Vic_up_Mig_down=read.csv("YPSTB_Vic_1Vic_2Vic_3_UP_Mig_1Mig_2Mig_3_DOWN.csv", sep=",")
Yersinia_pseudotuberculosis_expr_Vic_up_Mig_down =as.data.frame(Yersinia_pseudotuberculosis_Vic_up_Mig_down[,2:ncol(Yersinia_pseudotuberculosis_Vic_up_Mig_down)])
rownames(Yersinia_pseudotuberculosis_expr_Vic_up_Mig_down) =Yersinia_pseudotuberculosis_Vic_up_Mig_down[,1]

################Listeria_Vic_Sp_up_createdbyCansu###############
Listeria_monocytogenes_Sp_Vic_up=read.csv("Listeria_Sp_1Sp_2Sp_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Listeria_monocytogenes_expr_Sp_Vic_up=as.data.frame(Listeria_monocytogenes_Sp_Vic_up[,2:ncol(Listeria_monocytogenes_Sp_Vic_up)])
rownames(Listeria_monocytogenes_expr_Sp_Vic_up) =Listeria_monocytogenes_Sp_Vic_up[,1]

###########YPSTB_Vic_down_Nd_up_createdbyCansu############
Yersinia_pseudotuberculosis_Nd_up_Vic_down=read.csv("YPSTB_Nd_1Nd_2Nd_3_UP_Mig_1Mig_3Vic_1Vic_2Vic_3_DOWN.csv",  sep=",")
Yersinia_pseudotuberculosis_expr_Nd_up_Vic_down=as.data.frame(Yersinia_pseudotuberculosis_Nd_up_Vic_down[,2:ncol(Yersinia_pseudotuberculosis_Nd_up_Vic_down)])
rownames(Yersinia_pseudotuberculosis_expr_Nd_up_Vic_down) =Yersinia_pseudotuberculosis_Nd_up_Vic_down[,1]

###########ENTFA_Vic_down_Nd_up_Oxs_up_createdbyCansu############
Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down=read.csv("ENTFA_Nd_1Nd_2Nd_3Oxs_1Oxs_2Oxs_3_UP_Vic_1Vic_2Vic_3_DOWN.csv",  sep=",")
Enterococcus_faecalis_expr_Nd_up_Oxs_up_Vic_down=as.data.frame(Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down[,2:ncol(Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down)])
rownames(Enterococcus_faecalis_expr_Nd_up_Oxs_up_Vic_down) =Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down[,1]

##################ENTFA_Vic_up_Mig_up##########################
Enterococcus_faecalis_Mig_up_Vic_up=read.csv("ENTFA_Mig_1Mig_2Mig_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Enterococcus_faecalis_expr_Mig_up_Vic_up=as.data.frame(Enterococcus_faecalis_Mig_up_Vic_up[,2:ncol(Enterococcus_faecalis_Mig_up_Vic_up)])
rownames(Enterococcus_faecalis_expr_Mig_up_Vic_up) =Enterococcus_faecalis_Mig_up_Vic_up[,1]





##################UPEC_Vic_up_Mig_up##########################
Escherichia_coli_UPEC_Mig_up_Vic_up=read.csv("UPEC_Mig_1Mig_2Mig_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Escherichia_coli_UPEC_expr_Mig_up_Vic_up=as.data.frame(Escherichia_coli_UPEC_Mig_up_Vic_up[,2:ncol(Escherichia_coli_UPEC_Mig_up_Vic_up)])
rownames(Escherichia_coli_UPEC_expr_Mig_up_Vic_up) =Escherichia_coli_UPEC_Mig_up_Vic_up[,1]

################STAGA__Vic_down_Oss_up###############
Streptococcus_agalactiae_Oss_up_Vic_down=read.csv("STAGA_Oss_1Oss_2Oss_3_UP_Vic_1Vic_2Vic_3_DOWN.csv", sep=",")
Streptococcus_agalactiae_expr_Oss_up_Vic_down =as.data.frame(Streptococcus_agalactiae_Oss_up_Vic_down[,2:ncol(Streptococcus_agalactiae_Oss_up_Vic_down)])
rownames(Streptococcus_agalactiae_expr_Oss_up_Vic_down) =Streptococcus_agalactiae_Oss_up_Vic_down[,1]

##################ENTFA_Vic_down_Oxs_up##########################
Enterococcus_faecalis_Oxs_up_Vic_down=read.csv("ENTFA_Oxs_1Oxs_2Oxs_3_UP_Vic_1Vic_2Vic_3_DOWN.csv", sep=",")
Enterococcus_faecalis_expr_Oxs_up_Vic_down=as.data.frame(Enterococcus_faecalis_Oxs_up_Vic_down[,2:ncol(Enterococcus_faecalis_Oxs_up_Vic_down)])
rownames(Enterococcus_faecalis_expr_Oxs_up_Vic_down) =Enterococcus_faecalis_Oxs_up_Vic_down[,1]

################MRSA252_Vic_up_Oxs_up###############
Staphylococcus_aureus_MRSA252_Oxs_up_Vic_up=read.csv("MRSA252_Oxs_1Oxs_2Oxs_3_UP.csv", sep=",")
Staphylococcus_aureus_MRSA252_expr_Oxs_up_Vic_up =as.data.frame(Staphylococcus_aureus_MRSA252_Oxs_up_Vic_up[,2:ncol(Staphylococcus_aureus_MRSA252_Oxs_up_Vic_up)])
rownames(Staphylococcus_aureus_MRSA252_expr_Oxs_up_Vic_up) =Staphylococcus_aureus_MRSA252_Oxs_up_Vic_up[,1]

##################STRPN_Vic_up_Tm_up##########################
Streptococcus_pneumoniae_Tm_up_Vic_up=read.csv("STRPN_Tm_1Tm_2Tm_3Vic_1Vic_2Vic_3_UP.csv", sep=",")
Streptococcus_pneumoniae_expr_Tm_up_Vic_up=as.data.frame(Streptococcus_pneumoniae_Tm_up_Vic_up[,2:ncol(Streptococcus_pneumoniae_Tm_up_Vic_up)])
rownames(Streptococcus_pneumoniae_expr_Tm_up_Vic_up) =Streptococcus_pneumoniae_Tm_up_Vic_up[,1]

################SPYO_Vic_down###############
Streptococcus_pyogenes_Vic_down=read.csv("SPYO_As_3_UP_Bs_1Bs_3Vic_1Vic_2Vic_3_DOWN.csv", sep=",")
Streptococcus_pyogenes_expr_Vic_down =as.data.frame(Streptococcus_pyogenes_Vic_down[,2:ncol(Streptococcus_pyogenes_Vic_down)])
rownames(Streptococcus_pyogenes_expr_Vic_down) =Streptococcus_pyogenes_Vic_down[,1]

################Vibrio_Vic_up###############
Vibrio_cholerae_Vic_up=read.csv("Vibrio_Nd_1Vic_1Vic_2Vic_3_UP.csv", sep=",")
Vibrio_cholerae_expr_Vic_up =as.data.frame(Vibrio_cholerae_Vic_up[,2:ncol(Vibrio_cholerae_Vic_up)])
rownames(Vibrio_cholerae_expr_Vic_up) =Vibrio_cholerae_Vic_up[,1]



##################STRPN_Vic_down##########################
Streptococcus_pneumoniae_Vic_down=read.csv("STRPN_Mig_3Nd_2_UP_Vic_1Vic_2Vic_3_DOWN.csv", sep=",")
Streptococcus_pneumoniae_expr_Vic_down=as.data.frame(Streptococcus_pneumoniae_Vic_down[,2:ncol(Streptococcus_pneumoniae_Vic_down)])
rownames(Streptococcus_pneumoniae_expr_Vic_down) =Streptococcus_pneumoniae_Vic_down[,1]




####################################UI#########################

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("cosmo"),

    # Application title
    titlePanel(windowTitle = "COPATHOgenex",title=div(img(height = 50, width = 800, src = "Co-pathogenex logo V2.png") 
                         )),

    # Sidebar with a slider input for number of bins 
    # Copy the chunk below to make a group of checkboxes
    
    navbarPage("",
               tabPanel("Home",
               #          h2("Welcome to COPATHOgenex!", align="center", color="darkblue"),
               #           h4("General workflow explained below:"),
               # tags$video(id="video2", type = "video/mp4",src = "General_workflow_explained.mp4", autoplay=TRUE,
               #            muted=TRUE,
               #            playsinline=TRUE,
               #            loop=TRUE,
               #            width = 1300,
               #            top="2px",
               #            controls = "controls")
               img(height = 800, width = 900,src='Co-pathogenex homepage.png', align = "center"),
               ),
               navbarMenu("Gene co-expression network construction (GCN)",
    
    tabPanel("Pathogen-specific GCN",
             
             titlePanel("Pathogen specific gene co-expression network (GCN)"),
                 #sliderInput("myslider", "Number of rows to display", min=1, max = 32, value = 5),
                 #DT::dataTableOutput("mytable")
                 
                 
             
            
             mainPanel(
                 tags$style(type="text/css",
                            ".shiny-output-error { visibility: hidden; }",
                            ".shiny-output-error:before { visibility: hidden; }"
                            ),
                 
                 tabsetPanel(
                   tabPanel("Data Overview",
                            sidebarLayout(position = "left",
                                          sidebarPanel(      
                            selectInput(inputId = "dataset",
                                        label = "Choose a dataset:",
                                        choices = c("Aggregatibacter actinomycetemcomitans","Acinetobacter baumannii","Achromobacter xylosoxidans","Borrelia burgdorferi","Burkholderia pseudomallei",
                                                    "Campylobacter jejuni", "Enterococcus faecalis",
                                                    "Escherichia coli EPEC", "Escherichia coli ETEC", "Escherichia coli UPEC", "Francisella tularensis", "Haemophilus influenzae",
                                                    "Helicobacter pylori G27",  "Helicobacter pylori J99",  "Klebsiella pneumoniae", "Legionella pneumophila", "Listeria monocytogenes",
                                                    "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Neisseria meningitidis","Pseudomonas aeruginosa",
                                                    "Salmonella enterica", "Shigella flexneri" , "Staphylococcus aureus MRSA252", "Staphylococcus aureus MSSA476", "Staphylococcus epidermidis",
                                                    "Streptococcus agalactiae", "Streptococcus pneumoniae",   "Streptococcus pyogenes", "Streptococcus suis","Yersinia pseudotuberculosis"
                                                    , "Vibrio cholerae"
                                                    
                                                    
                                                     
                                                     
                                                    
                                        )),
                           
                            column(5,
                                   br(), 
                                   br(),
                                   submitButton("Update View", icon("refresh")),
                                   br(), 
                                   br()
                            ),
                            
                            ),
                                        
                            mainPanel("Display species information:", 
                            
                                      textOutput(outputId ="summary_datasetin1" )%>% withSpinner(type=5),
                            fluidRow(
                              #splitLayout(cellWidths = c("50%", "50%"), plotOutput("clusterone"), plotOutput("powerplot1"))
                              column(2, align="right",
                                     plotOutput(outputId = "clusterone", width  = "600px",height = "400px"),  
                                     plotOutput(outputId = "powerplot1", width  = "600px",height = "400px")  
                                   
                                     #plotOutput(outputId = "clusterone", width  = "500px",height = "400px")
                              ),
                              
                          
                             
                              
                            ))) 
                            
                            
                            
                            
                            ), 
                   
                   
                     tabPanel("Network construction",
                              # actionButton("update" ,"Update View", icon("arrows-rotate"),
                              #              class = "btn btn-primary"),
                              # helpText("When you click the button above, you should see",
                              #          "the output below update to reflect the values you",
                              #          "entered at the Construction tab"),
                              sidebarLayout(
                                  
                                  # Sidebar panel for inputs ----
                                  sidebarPanel(
                                      
                                      
                                      # Input: Selector for choosing dataset ----
                                      
                                      # selectInput(inputId = "dataset",
                                      #             label = "Choose a dataset:",
                                      #             choices = c("Escherichia coli EPEC", "Escherichia coli ETEC", "Escherichia coli UPEC","Helicobacter pylori G27",  "Helicobacter pylori J99", "Salmonella enterica", "Acinetobacter baumannii",
                                      #                         "Aggregatibacter actinomycetemcomitans", "Achromobacter xylosoxidans", "Burkholderia pseudomallei", "Borrelia burgdorferi", 
                                      #                         "Campylobacter jejuni", "Enterococcus faecalis", "Francisella tularensis", "Haemophilus influenzae", "Klebsiella pneumoniae", 
                                      #                         "Legionella pneumophila", "Listeria monocytogenes", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", 
                                      #                         "Neisseria meningitidis", "Pseudomonas aeruginosa", "Staphylococcus aureus MRSA252", 
                                      #                         "Staphylococcus aureus MSSA476", "Staphylococcus epidermidis", "Streptococcus agalactiae", "Streptococcus pneumoniae",
                                      #                         "Streptococcus pyogenes", "Streptococcus suis","Shigella flexneri" ,"Yersinia pseudotuberculosis", "Vibrio cholerae"
                                      #                         )),
                                      
                                      #Input: Selector for choosing method ----
                                      #checkboxInput( "graph",
                                       #            "Network construction"),
                                                   
                                      
                                      #          sidebarPanel(       
                                      # Input: Slider for the number of bins ----
                                      
                                      strong("Tuning parameters for network construction:"),
                                      sliderInput(inputId = "powerin",
                                                  label = "1. Power",
                                                  min = 2,
                                                  max = 20,
                                                  value = 6),
                                      
                                      sliderInput(inputId = "mingene",
                                                  label = "2. Genes per module (minimum)",
                                                  min = 2,
                                                  max = 30,
                                                  value = 15),
                                      
                                      column(5,
                                             br(),
                                             br(),
                                             submitButton("Update View", icon("refresh")),
                                             br(),
                                             br()
                                      )
                                      
                                      ),
                                  
                                  # Main panel for displaying outputs ----
                                  
                                  
                                  mainPanel("Clustering dendrogram of genes based on topological overlap",
                                            textOutput(outputId ="summary_datasetin3")%>% withSpinner(type=5)     ,
                                                   fluidRow(
                                                     #splitLayout(cellWidths = c("50%", "50%"), plotOutput("clusterone"), plotOutput("powerplot1"))
                                                     column(2, align="right",
                                                            
                                                            plotOutput(outputId = "selected_graph", width  = "600px",height = "400px" )
                                                              
                                                            
                                                            #plotOutput(outputId = "clusterone", width  = "500px",height = "400px")
                                                     ),
                                                     
                                                         
                                                    
                                                     
                                                   ),   
                                                   
                                            verbatimTextOutput(outputId ="summary_datasetin2" )%>% withSpinner(type=5)      
                                                   
                                                          
                                                   #plotOutput("selected_graph")), 
                                        
                                          #tabPanel("Summary",verbatimTextOutput("summary_datasetin1")), 
                                      
                                      
                                      
                                  
                                  
                              ))),       
                     
                     
                     tabPanel("Modules expression",
                              h4("Each module is represented by an expression value which belongs to the module 'eigengene'. This value is identified from the PCA"),
                              # submitButton("Update View", icon("refresh")),
                              # helpText("When you click the button above, you should see",
                              #          "the output below update to reflect the values you",
                              #          "entered at the Data Overview and Construction tab"),
                              textOutput(outputId ="summary_datasetin4" ),
                              dataTableOutput(outputId = "mytable") %>% withSpinner(type=5),   
                   br(),
                   hr(),
                   p(" The stress abbreviations used in the table are listed as follows:"),
                   tags$ul(
                     tags$li(tags$b("As"), " - Acidic stress"),
                     tags$li(tags$b("Bs"), " - Bile stress"),
                     tags$li(tags$b("Ctrl"), " - Control"),
                     tags$li(tags$b("Li"), " - Low iron"),
                     tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                     tags$li(tags$b("Nd"), " - Nutritional downshift"),
                     tags$li(tags$b("Ns"), " - Nitrosative stress"),
                     tags$li(tags$b("Oss"), " - Osmotic stress"),
                     tags$li(tags$b("Oxs"), " - Oxidative stress"),
                     tags$li(tags$b("Sp"), " - Stationary phase"),
                     tags$li(tags$b("Tm"), " - Temperature"),
                     tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                   )
                   
                     ),  
                     tabPanel("Gene expression",
                              # actionButton("update" ,"Update View", icon("arrows-rotate"),
                              #              class = "refresh"),
                              # helpText("When you click the button above, you should see",
                              #          "the output below update to reflect the values you",
                              #          "entered at the  Data Overview and Construction tab"),
                              # submitButton("Update View", icon("refresh")),
                              # helpText("When you click the button above, you should see",
                              #          "the output below update to reflect the values you",
                              #          "entered at the Data Overview and Construction tab"),
                              h4("Gene expression matrix represented as standarized (Z-Score) TPM values"),
                              textOutput(outputId ="summary_datasetin5" ),
                              dataTableOutput("mytable2")%>% withSpinner(type=5) ),
                     tabPanel("Module visualization",
                              
                              sidebarLayout(
                                
                                # Sidebar panel for inputs ----
                                sidebarPanel(
                                  
                                            textInput("textforcolor", "Choose a module for analysis"),
                                             
                                             radioButtons(inputId = "graph2",
                                                          label = "Plots for module overview",
                                                          choices = c("Gene expression heatmap", 
                                                                      #"Gene expression box plots", 
                                                                      "Module eigengene barplot")),
                                             
                                             column(5,
                                                    br(),
                                                    submitButton("Submit", icon("refresh")),
                                                    
                                                    
                                                   
                                                    
                                                    
                                                    
                                                    
                                                    
                                                    
                                             ),                      
                                            br(),
                                            br(),
                                            br(),
                                            helpText("Please enter a module color and hit Submit"),     
                                             
                                ),                 
                              
                              mainPanel(

                                
                                plotOutput("myresultscolor"),
                                br(),
                                br(), 
                                br(),
                                br(),
                                br(),
                                br(),
                                hr(), 
                                downloadButton(
                                  outputId = "download_heatmap",
                                  label = "Download heatmap"
                                ),
                                
                                downloadButton(
                                  outputId = "download_barplot",
                                  label = "Download barplot"
                                )           
                                  
                              )
                              
                     )
                                  
                                  
                     
                              
            
             )
                 )        
             )
             ),
    #tabPanel("GO overview",parcatsOutput("selected_graph2",height = "1000px",width = "1000px", inline = FALSE))
  tabPanel("Core GCN",
  titlePanel("Core transcriptome gene co-expression network (GCN)"),
  #sliderInput("myslider", "Number of rows to display", min=1, max = 32, value = 5),
  #DT::dataTableOutput("mytable")
  
  
  
  
  mainPanel(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),
      tabsetPanel(
        
        
        tabPanel("Data Overview",
                 sidebarLayout(position = "left",
                               sidebarPanel(      
                                 selectInput(inputId = "dataset_core",
                                             label = "Choose a dataset:",
                                                       choices = c("HP_J99 vs HP_G27", "UPEC vs EPEC", "UPEC vs ETEC", "EPEC vs ETEC",  "MSSA476 vs MRSA252"
                                                                   #,"EPEC vs ETEC vs UPEC"
                                                       )),
                                 
                                 column(5,
                                        br(), 
                                        br(),
                                        submitButton("Update View", icon("refresh")),
                                        br(), 
                                        br()
                                 ),
                                 
                               ),
                               
                               mainPanel("Display species information:", 
                                         
                                         textOutput(outputId ="summary_datasetin_core1" )%>% withSpinner(type=5),
                                         fluidRow(
                                           #splitLayout(cellWidths = c("50%", "50%"), plotOutput("clusterone"), plotOutput("powerplot1"))
                                           column(2, align="right",
                                                  plotOutput(outputId = "clustertwo", width  = "600px",height = "400px"),  
                                                  plotOutput(outputId = "powerplot2", width  = "600px",height = "400px")  
                                                  
                                                  #plotOutput(outputId = "clusterone", width  = "500px",height = "400px")
                                           ),
                                           
                                           
                                           
                                           
                                         ))) 
                 
                 
                 
                 
        ), 
        
        
        
        
          tabPanel("Network construction",
                   sidebarLayout(
                       
                       # Sidebar panel for inputs ----
                       sidebarPanel(
                           
                           
                         
                           
                           #Input: Selector for choosing method ----
                           # radioButtons(inputId = "graph_core",
                           #              label = "Inspection and modelling:",
                           #              choices = c("Sample cluster", "Power analysis", 
                           #                          "Network construction")),
                           
                           
                           #          sidebarPanel(       
                           # Input: Slider for the number of bins ----
                           
                           strong("Tuning parameters for network construction:"),
                           sliderInput(inputId = "powerin_core",
                                       label = "1. Power",
                                       min = 2,
                                       max = 20,
                                       value = 6),
                           
                           sliderInput(inputId = "mingene_core",
                                       label = "2. Genes per module (minimum)",
                                       min = 2,
                                       max = 30,
                                       value = 15),
                           
                           column(5,
                                  br(), 
                                  br(),
                                  submitButton("Submit"),
                                  br(), 
                                  br()
                           )),
                       
                       # Main panel for displaying outputs ----
                       
                   #     
                   #     mainPanel(
                   #         
                   #         tabsetPanel(
                   #             tabPanel("Plot",plotOutput("selected_graph_core")), 
                   #             tabPanel("Summary",verbatimTextOutput("summary_datasetin_core")), 
                   #         ),
                   #         
                   #         
                   #     )
                   #     
                   # )),       
                   mainPanel("Clustering dendrogram of genes based on topological overlap",
                             textOutput(outputId ="summary_datasetin_core3")%>% withSpinner(type=5)     ,
                             fluidRow(
                               #splitLayout(cellWidths = c("50%", "50%"), plotOutput("clusterone"), plotOutput("powerplot1"))
                               column(2, align="right",
                                      
                                      plotOutput(outputId = "selected_graph_core", width  = "600px",height = "400px" )
                                      
                                      
                                      #plotOutput(outputId = "clusterone", width  = "500px",height = "400px")
                               ),
                               
                               
                               
                               
                             ),   
                             
                             verbatimTextOutput(outputId ="summary_datasetin_core2" )%>% withSpinner(type=5)      
                             
                             
                             #plotOutput("selected_graph")), 
                             
                             #tabPanel("Summary",verbatimTextOutput("summary_datasetin1")), 
                             
                             
                             
                             
                             
                   ))),       
        
        tabPanel("Modules expression",
                 h4("Each module is represented by an expression value which belongs to the module 'eigengene'. This value is identified from the PCA"),
                 # submitButton("Update View", icon("refresh")),
                 # helpText("When you click the button above, you should see",
                 #          "the output below update to reflect the values you",
                 #          "entered at the Data Overview and Construction tab"),
                 textOutput(outputId ="summary_datasetin_core4" ),
                 dataTableOutput(outputId = "mytable_core") %>% withSpinner(type = 5)   ),  
          #tabPanel("Modules eigengene",dataTableOutput("mytable_core")), 
          tabPanel("Consensus between strains",
                   h4("Relating the consensus modules of two selected strains:"),
                   # actionButton("update" ,"Update View", icon("arrows-rotate"),
                   #              class = "btn btn-primary"),
                   # helpText("When you click the button above, you should see",
                   #          "the output below update to reflect the values you",
                   #          "entered at the Construction tab"),
                   #verbatimTextOutput("powerin_core"),
                   plotOutput("consensus_graph_core") %>% withSpinner(type = 5)),
        tabPanel("Gene expression",
                 # actionButton("update" ,"Update View", icon("arrows-rotate"),
                 #              class = "refresh"),
                 # helpText("When you click the button above, you should see",
                 #          "the output below update to reflect the values you",
                 #          "entered at the  Data Overview and Construction tab"),
                 # submitButton("Update View", icon("refresh")),
                 # helpText("When you click the button above, you should see",
                 #          "the output below update to reflect the values you",
                 #          "entered at the Data Overview and Construction tab"),
                 h4("Gene expression matrix represented as standarized (Z-Score) TPM values"),
                 textOutput(outputId ="summary_datasetin_core5" ),
                 dataTableOutput("mytable2_core")%>% withSpinner(type=5) ),
        
          #tabPanel("Gene expression",dataTableOutput("mytable2_core")),
        tabPanel("Module visualization",
                 
                 sidebarLayout(
                   
                   # Sidebar panel for inputs ----
                   sidebarPanel(
                     
                     textInput("textforcolor_core", "Choose a module for analysis"),
                     
                     radioButtons(inputId = "graph2_core",
                                  label = "Plots for module overview",
                                  choices = c("Gene expression heatmap", 
                                              #"Gene expression box plots", 
                                              "Module eigengene barplot")),
                     
                     column(5,
                            br(),
                            submitButton("Submit", icon("refresh")),
                            
                            
                            
                            
                            
                            
                            
                            
                            
                     ),                      
                     br(),
                     br(),
                     br(),
                     helpText("Please enter a module color and hit Submit"),     
                     
                   ),                 
                   
                   mainPanel(plotOutput("myresultscolor_core")),
                     
                     
                     
                     # br(),
                     # br(), 
                     # br(),
                     # br(),
                     # br(),
                     # br(),
                     # hr(), 
                     # downloadButton(
                     #   outputId = "download_heatmap",
                     #   label = "Download heatmap"
                     # ),
                     # 
                     # downloadButton(
                     #   outputId = "download_barplot",
                     #   label = "Download barplot"
                     # )           
                     
                   
                   
                 )
                 
                 
                 
                 
                 
        )
      )        
  )
  )),
        
  #         tabPanel("Module visualization",fluidRow(column(5,textInput("textforcolor_core", "Choose a module for analysis")),
  #                                             column(5,  
  #                                                    radioButtons(inputId = "graph2_core",
  #                                                                 label = "Plots for module overview",
  #                                                                 choices = c("Gene expression heatmap", 
  #                                                                             #"Gene expression box plots", 
  #                                                                             "Module eigengene barplot"))),
  #                                             
  #                                             
  #                                             column(5,
  #                                                    br(), 
  #                                                    br(),
  #                                                    submitButton("Submit"),
  #                                                    br(), 
  #                                                    br()
  #                                             )),
  #                  
  #                  plotOutput("myresultscolor_core")
  #                  
  #                  
  #         )
  #     )        
  # )
  # 
  #   )),
  navbarMenu("Stress-specific stimulons",
             
             ####Data_entry_acid_stress#########             
             
             tabPanel("Acid stress",
                      
                      titlePanel("Acid stress stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_acid_stress",
                            
                            label = "Select Species",
                            
                            choices = c(
                              
                              "Burkholderia pseudomallei As (+)",            
                              "Enterococcus faecalis As (+)",               
                              "Escherichia coli EPEC As (+) (I)",            
                              "Escherichia coli EPEC As (+) (II)",          
                              "Escherichia coli ETEC As (+)",               
                              "Escherichia coli UPEC As (+) Sp (-)",
                              "Klebsiella pneumoniae As (+) Sp (-)",         
                              "Listeria monocytogenes As (+)",              
                              "Listeria monocytogenes As (+) Nd (-)",        
                              "Listeria monocytogenes As (+) Oxs (+)",      
                              "Pseudomonas aeruginosa As (+)",               
                              "Pseudomonas aeruginosa As (+) Tm (+)",       
                              "Salmonella enterica As (+)",                  
                              "Shigella flexneri As (+)",                   
                              "Staphylococcus aureus MRSA252 As (+)",        
                              "Staphylococcus aureus MSSA476 As (+)",       
                              "Staphylococcus aureus MSSA476 As (+) Ns (-)", 
                              
                              "Staphylococcus epidermidis As (+)",          
                              "Staphylococcus epidermidis As (+) Mig (+) ", 
                              "Staphylococcus epidermidis As (+) Tm (-)"
                              
                              
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,

                                 br(),

                                 br(),

                                 submitButton("Submit"),

                                 br(),

                                 br()

                          )



                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Acidstress_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                              
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableAs")
                              
                            )
                            
                            
                            
                          ),
                          
                        ),
                        
                      )),
             
             ####Data_entry_sp#########
             
             tabPanel("Stationary phase",
                      
                      titlePanel("Stationary phase stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_sp",
                            
                            label = "Select Species",
                            
                            choices = c(
                              
                              "Achromobacter xylosoxidans Sp (+)",            
                              "Acinetobacter baumannii Sp (-) (I)",          
                              "Acinetobacter baumannii Sp (-) (II)" ,         
                              "Acinetobacter baumannii Sp (-) (III)",        
                              "Acinetobacter baumannii Sp (-) (IV)" ,         
                              "Aggregatibacter actinomycetemcomitans Sp (+)",
                              "Borrelia burgdorferi Sp (+) (I)",   
                              "Borrelia burgdorferi Sp (+) (II)", 
                              "Burkholderia pseudomallei Sp (-)",            
                              "Burkholderia pseudomallei Sp (+)",           
                              "Burkholderia pseudomallei Sp (+) Mig (+)",    
                              "Burkholderia pseudomallei Sp (+) Tm (+)",    
                              "Campylobacter jejuni Sp (+)",                 
                              "Enterococcus faecalis Sp (+)" ,                
                              "Enterococcus faecalis Sp (+) Mig (+)", 
                              "Escherichia coli EPEC Sp (-)",                
                              "Escherichia coli EPEC Sp (+)",              
                              "Escherichia coli ETEC Sp (+)",                 
                              "Escherichia coli UPEC Sp (-) As (+)",         
                              "Escherichia coli UPEC Sp (-) Mig (-)",        
                              "Escherichia coli UPEC Sp (+)",            
                              "Francisella tularensis Sp (+)",               
                              "Haemophilus influenzae Sp (+) (I)",    
                              "Haemophilus influenzae Sp (+) (II)",  
                              "Helicobacter pylori G27 Sp (+)",               
                              "Klebsiella pneumoniae Sp (-)",                
                              "Klebsiella pneumoniae Sp (-) As (+)",          
                              "Klebsiella pneumoniae Sp (-) Bs (+)",         
                              "Klebsiella pneumoniae Sp (+)",                 
                              "Listeria monocytogenes Sp (+) Vic (+)",
                              "Neisseria gonorrhoeae Sp (-)",                
                              "Neisseria gonorrhoeae Sp (+) (I)",            
                              "Neisseria gonorrhoeae Sp (+) (II)",            
                              "Neisseria meningitidis Sp (+)",               
                              "Neisseria meningitidis Sp (+) Li (+)",         
                              "Pseudomonas aeruginosa Sp (-) (I)",           
                              "Pseudomonas aeruginosa Sp (-) (II)",           
                              "Pseudomonas aeruginosa Sp (-) Ns (+) (I)",   
                              "Pseudomonas aeruginosa Sp (-) Ns (+) (II)",    
                              "Pseudomonas aeruginosa Sp (+)",               
                              "Pseudomonas aeruginosa Sp (+) Nd (+)",         
                              "Salmonella enterica Sp (+) (I)",              
                              "Salmonella enterica Sp (+) (II)",              
                              "Shigella flexneri Sp (-) Mig (-)",           
                              "Shigella flexneri Sp (+)",    
                              
                              "Staphylococcus aureus MRSA252 Sp (+)",        
                              "Staphylococcus aureus MSSA476 Sp (+)",         
                              "Staphylococcus aureus MSSA476 Sp (+) Nd (+)", 
                              "Staphylococcus epidermidis Sp (+)",            
                              "Staphylococcus epidermidis Sp (+) Tm (+)",  
                              "Streptococcus agalactiae Sp (+)", 
                              "Streptococcus agalactiae Sp (+) Mig (+)",
                              
                              "Streptococcus pneumoniae Sp (+) (I)",         
                              "Streptococcus pneumoniae Sp (+) (II)",  
                              
                              "Streptococcus pyogenes Sp (-) Mig (-) (I)",   
                              "Streptococcus pyogenes Sp (-) Mig (-) (II)",   
                              "Streptococcus pyogenes Sp (-) Mig (-) Li (+)",
                              "Streptococcus pyogenes Sp (+) Mig (+)",        
                              "Vibrio cholerae Sp (+)",                      
                              "Vibrio cholerae Sp (+) Li (+)",                
                              "Yersinia pseudotuberculosis Sp (+) (I)",      
                              "Yersinia pseudotuberculosis Sp (+) (II)",      
                              "Yersinia pseudotuberculosis Sp (+) Mig (-)"  
                              
                              
                              
                              
                              
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,
                                 
                                 br(), 
                                 
                                 br(),
                                 
                                 submitButton("Submit"),
                                 
                                 br(), 
                                 
                                 br()
                                 
                          )
                          
                          
                          
                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Sp_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableSp")
                              
                            )
                            
                            
                            
                          ),
                          
                        ),
                        
                      )),
             
             
             ####Data_entry_nd#########
             
             tabPanel("Nutritional downwshift",
                      
                      titlePanel("Nutritional downwshift stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_nd",
                            
                            label = "Select Species",
                            
                            choices = c(
                              
                              "Achromobacter xylosoxidans Nd (+)",                     
                              "Acinetobacter baumannii Nd (+)",                       
                              "Aggregatibacter actinomycetemcomitans  Nd (-) Mig (+)", 
                              "Aggregatibacter actinomycetemcomitans Nd (+) (I)",         
                              "Aggregatibacter actinomycetemcomitans Nd (+) (II)" ,         
                              "Borrelia burgdorferi Nd (+)",                          
                              "Burkholderia pseudomallei Nd (+) (I)",                  
                              "Burkholderia pseudomallei Nd (+) (II)",                
                              "Burkholderia pseudomallei Nd (+) (III)",                
                              "Campylobacter jejuni Nd (+)",                          
                              "Enterococcus faecalis Nd (-)",                         
                              "Enterococcus faecalis Nd (+) (I)",                     
                              "Enterococcus faecalis Nd (+) (II)",                     
                              "Enterococcus faecalis Nd (+) Oxs (+) Vic(-)",          
                              "Escherichia coli EPEC Nd (+)",                          
                              "Escherichia coli EPEC Nd (+) Oss (+)",                 
                              "Escherichia coli ETEC Nd (+) (I)",                      
                              "Escherichia coli ETEC Nd (+) (II)",                    
                              "Escherichia coli ETEC Nd (+) Mig (+)", 
                              
                              "Escherichia coli UPEC Nd (+) (I)",                     
                              "Escherichia coli UPEC Nd (+) (II)",                     
                              "Escherichia coli UPEC Nd (+) (III)",   
                              
                              "Francisella tularensis Nd (-) (I)",                     
                              "Francisella tularensis Nd (-) (II)", 
                              
                              "Francisella tularensis Nd (+) (I)",                     
                              "Francisella tularensis Nd (+) (II)",
                              
                              "Haemophilus influenzae Nd (+)",                        
                              "Helicobacter pylori G27 Nd (-) (I)",                       
                              "Helicobacter pylori G27 Nd (-) (II)", 
                              
                              "Helicobacter pylori G27 Nd (+)",                       
                              "Helicobacter pylori J99 Nd (-)",   
                              
                              "Helicobacter pylori J99 Nd (+) (I)",                   
                              "Helicobacter pylori J99 Nd (+) (II)",                  
                              "Klebsiella pneumoniae Nd (+) (I)",
                              "Klebsiella pneumoniae Nd (+) (II)",
                              "Klebsiella pneumoniae Nd (+) Oxs (+)",                 
                              "Listeria monocytogenes  Nd (-) Li (+)", 
                              
                              "Listeria monocytogenes Nd (-) (I)",                     
                              "Listeria monocytogenes Nd (-) (II)",                   
                              "Listeria monocytogenes Nd (-) (III)",   
                              
                              "Listeria monocytogenes Nd (-) As (+)",                 
                              "Listeria monocytogenes Nd (-) Oss (+) ",                
                              "Listeria monocytogenes Nd (+) (I)",                    
                              "Listeria monocytogenes Nd (+) (II)",                   
                              "Listeria monocytogenes Nd (+) (III)",                  
                              "Neisseria meningitidis Nd (+)",                         
                              "Pseudomonas aeruginosa Nd (+)",                        
                              "Pseudomonas aeruginosa Nd (+) Ns (+)",                  
                              "Pseudomonas aeruginosa Nd (+) Oss (+)",                
                              "Pseudomonas aeruginosa Nd (+) Sp (+)",                 
                              "Salmonella enterica Nd (+)",                           
                              "Salmonella enterica Nd (+) Oxs (+)",                    
                              "Shigella flexneri Nd (+)",                             
                              "Shigella flexneri Nd (+) Oxs (+)",                      
                              "Staphylococcus aureus MRSA252 Nd (-) Mig (-)",         
                              "Staphylococcus aureus MRSA252 Nd (+)",                  
                              "Staphylococcus aureus MRSA252 Nd (+) Mig (+) ",        
                              "Staphylococcus aureus MRSA252 Nd (+) Ns (-)",           
                              "Staphylococcus aureus MSSA476 Nd (-) Mig (-)",         
                              "Staphylococcus aureus MSSA476 Nd (+)",                  
                              "Staphylococcus aureus MSSA476 Nd (+) Sp (+)",          
                              "Staphylococcus epidermidis Nd (+) (I)",                
                              "Staphylococcus epidermidis Nd (+) (II)",               
                              "Staphylococcus epidermidis Nd (+) Tm (-)",              
                              "Streptococcus agalactiae Nd (-) (I)",                  
                              "Streptococcus agalactiae Nd (-) (II)",                 
                              "Streptococcus agalactiae Nd (-) Mig (+) ",             
                              "Streptococcus agalactiae Nd (-) Oss (+) ",              
                              "Streptococcus agalactiae Nd (-) Tm (+) ",              
                              "Streptococcus agalactiae Nd (-) Vic (+) ",              
                              "Streptococcus agalactiae Nd (+)",                      
                              "Streptococcus pyogenes Nd (+)",                         
                              "Streptococcus pyogenes Nd (+) Tm (+) ",                
                              "Vibrio cholerae Nd (+)",                               
                              "Vibrio cholerae Nd (+) Mig (+) ",                      
                              "Yersinia pseudotuberculosis Nd (+)",                   
                              "Yersinia pseudotuberculosis Nd (+) Vic (-)"   
                              
                              
                              
                              
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,
                                 
                                 br(), 
                                 
                                 br(),
                                 
                                 submitButton("Submit"),
                                 
                                 br(), 
                                 
                                 br()
                                 
                          )
                          
                          
                          
                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Nd_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableNd")
                              
                            )
                          ),
                        ),
                      )),
             
             ####Data_entry_bs#########
             
             tabPanel("Bile stress",
                      
                      titlePanel("Bile stress stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_bs",
                            
                            label = "Select Species",
                            
                            choices = c(
                              "Achromobacter xylosoxidans Bs (+)",                  
                              "Aggregatibacter actinomycetemcomitans Bs (+) (I)",  
                              "Aggregatibacter actinomycetemcomitans Bs (+) (II)",  
                              "Aggregatibacter actinomycetemcomitans Bs (+) (III)",
                              "Aggregatibacter actinomycetemcomitans Bs (+) (IV)",  
                              "Campylobacter jejuni Bs (+)",                       
                              "Enterococcus faecalis Bs (+)",                       
                              "Escherichia coli EPEC Bs (+)",                      
                              "Escherichia coli ETEC Bs (+)",                       
                              "Escherichia coli UPEC Bs (+)",                      
                              "Francisella tularensis Bs (+)",                      
                              "Klebsiella pneumoniae Bs (+)",                      
                              "Klebsiella pneumoniae Bs (+) Sp (-)",                
                              "Listeria monocytogenes Bs (+) (I) ",                
                              "Listeria monocytogenes Bs (+) (II)",                 
                              "Neisseria gonorrhoeae Bs (+) (I)",                  
                              "Neisseria gonorrhoeae Bs (+) (II)",                  
                              "Neisseria meningitidis Bs (+)",                     
                              "Pseudomonas aeruginosa Bs (+)",                      
                              "Salmonella enterica Bs (+)",                        
                              "Staphylococcus aureus MRSA252 Bs (+)",               
                              "Staphylococcus epidermidis Bs (+)",                 
                              "Streptococcus pneumoniae Bs (+) Mig (-)",            
                              "Vibrio cholerae Bs (+)"              
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,
                                 
                                 br(), 
                                 
                                 br(),
                                 
                                 submitButton("Submit"),
                                 
                                 br(), 
                                 
                                 br()
                                 
                          )
                          
                          
                          
                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Bs_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableBs")
                              
                            )
                            
                            
                          ),
                          
                        ),
                        
                      )),
             
             
             ####Data_entry_Li#########
             
             tabPanel("Low iron",
                      
                      titlePanel("Low Iron stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_li",
                            
                            label = "Select Species",
                            
                            choices = c(
                              "Achromobacter xylosoxidans Li (+) (I)",             
                              "Achromobacter xylosoxidans Li (+) (II)",           
                              "Acinetobacter baumannii Li (+) Oxs (+)",           
                              "Aggregatibacter actinomycetemcomitans Li (+) (I)", 
                              "Aggregatibacter actinomycetemcomitans Li (+) (II)", 
                              "Burkholderia pseudomallei Li (+)",                 
                              "Campylobacter jejuni Li (+)",                       
                              "Enterococcus faecalis Li (+)",                     
                              "Escherichia coli EPEC Li (+)",                      
                              "Escherichia coli ETEC Li (+)",                     
                              "Escherichia coli ETEC Li (+) Tm (+) ",              
                              "Escherichia coli UPEC Li (+)",                     
                              "Francisella tularensis Li (+)",                     
                              "Helicobacter pylori G27 Li (+)",                   
                              "Helicobacter pylori J99 Li (+)",                    
                              "Klebsiella pneumoniae Li (+)",                     
                              "Listeria monocytogenes Li (+) Nd (-)",              
                              "Neisseria gonorrhoeae Li (+)",                     
                              "Neisseria meningitidis Li (+) Sp (+)",              
                              "Pseudomonas aeruginosa Li (+) (I)",                
                              "Pseudomonas aeruginosa Li (+) (II)",                
                              "Salmonella enterica Li (+)",                       
                              "Salmonella enterica Li (+) Oxs (+)",                
                              "Shigella flexneri Li (+)",                         
                              "Shigella flexneri Li (+) Oxs (+)",                  
                              "Staphylococcus aureus MRSA252 Li (+)",             
                              "Staphylococcus epidermidis Li (+)",                 
                              "Streptococcus agalactiae Li (+) (I)",              
                              "Streptococcus agalactiae Li (+) (II)",              
                              "Streptococcus pneumoniae Li (+)",                  
                              "Streptococcus pneumoniae Li (+) Mig (-)",           
                              "Streptococcus pyogenes Li (+) Mig (-) Sp (-) ",    
                              "Vibrio cholerae Li (+) (I)",                        
                              "Vibrio cholerae Li (+) (II)",                      
                              "Vibrio cholerae Li (+) Sp (+)",                     
                              "Yersinia pseudotuberculosis Li (+)"   
                              
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,
                                 
                                 br(), 
                                 
                                 br(),
                                 
                                 submitButton("Submit"),
                                 
                                 br(), 
                                 
                                 br()
                                 
                          )
                          
                          
                          
                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Li_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableLi")
                              
                            )
                            
                            
                          ),
                          
                        ),
                        
                      )),
             
             
             ####Data_entry_Mig#########
             
             tabPanel("Microaerophilic growth",
                      
                      titlePanel("Microaerophilic growth stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_mig",
                            
                            label = "Select Species",
                            
                            choices = c(
                              "Achromobacter xylosoxidans Mig (+) (I)",               
                              "Achromobacter xylosoxidans Mig (+) (II)",             
                              "Aggregatibacter actinomycetemcomitans Mig (+) Nd (-)", 
                              "Burkholderia pseudomallei Mig (+) (I)",               
                              "Burkholderia pseudomallei Mig (+) (II)",               
                              "Burkholderia pseudomallei Mig (+) Sp (+)",        
                              
                              "Enterococcus faecalis Mig (+)",                       
                              "Enterococcus faecalis Mig (+) Sp (+)",                
                              "Enterococcus faecalis Mig (+) Vic (+)",               
                              "Escherichia coli EPEC Mig (+)",                        
                              "Escherichia coli ETEC Mig (+)",                       
                              "Escherichia coli ETEC Mig (+) Nd (+)",                 
                              "Escherichia coli ETEC Mig (+) Oss (+)",               
                              "Escherichia coli UPEC Mig (-) Sp (-)",                 
                              "Escherichia coli UPEC Mig (+)",                       
                              "Escherichia coli UPEC Mig (+) Vic (+)",                
                              "Haemophilus influenzae Mig (-) (I)",                  
                              "Haemophilus influenzae Mig (-) (II)",                  
                              "Haemophilus influenzae Mig (-) Vic (+)",              
                              "Haemophilus influenzae Mig (+)",                       
                              "Haemophilus influenzae Mig (+) Vic (+)",              
                              "Klebsiella pneumoniae Mig (+) ",                       
                              "Listeria monocytogenes Mig (-)",                      
                              "Listeria monocytogenes Mig (+)",                       
                              "Neisseria gonorrhoeae Mig (+) (I)",                   
                              "Neisseria gonorrhoeae Mig (+) (II)",                   
                              "Neisseria meningitidis Mig (+)",                      
                              "Pseudomonas aeruginosa Mig (+)",                       
                              "Pseudomonas aeruginosa Mig (+) Ns (-)",               
                              "Pseudomonas aeruginosa Mig (+) Oxs (+)",               
                              "Salmonella enterica Mig (+)",                         
                              "Shigella flexneri Mig (-) Sp (-)",                     
                              "Shigella flexneri Mig (+)",                           
                              "Staphylococcus aureus MRSA252 Mig (-) Nd (-)",         
                              "Staphylococcus aureus MRSA252 Mig (+) (I)",           
                              "Staphylococcus aureus MRSA252 Mig (+) (II)",           
                              "Staphylococcus aureus MRSA252 Mig (+) Nd (+)",        
                              "Staphylococcus aureus MSSA476 Mig (-) Nd (-)",         
                              "Staphylococcus aureus MSSA476 Mig (+)",               
                              "Staphylococcus epidermidis Mig (+) (I)",               
                              "Staphylococcus epidermidis Mig (+) (II)",             
                              "Staphylococcus epidermidis Mig (+) As (+)",            
                              "Staphylococcus epidermidis Mig (+) Tm (-)",    
                              
                              "Streptococcus agalactiae Mig (+) Nd (-)",              
                              "Streptococcus agalactiae Mig (+) Sp (+)",         
                              
                              "Streptococcus pneumoniae Mig (-) Bs (+)",             
                              "Streptococcus pneumoniae Mig (-) Li (+)",              
                              "Streptococcus pneumoniae Mig (-) Vic (+) ",           
                              "Streptococcus pneumoniae Mig (+)",                     
                              "Streptococcus pyogenes Mig (-) Sp (-) (I)",               
                              "Streptococcus pyogenes Mig (-) Sp (-) (II)",                
                              "Streptococcus pyogenes Mig (-) Sp (-) Li (+)",      
                              
                              "Streptococcus pyogenes Mig (+) Sp (+)",            
                              
                              "Vibrio cholerae Mig (+) (I)",                          
                              "Vibrio cholerae Mig (+) (II)",                        
                              "Vibrio cholerae Mig (+) Nd (+)",                       
                              "Yersinia pseudotuberculosis Mig (-) Sp (+)",          
                              "Yersinia pseudotuberculosis Mig (-) Vic (+) ",         
                              "Yersinia pseudotuberculosis Mig (+)"    
                              
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,
                                 
                                 br(), 
                                 
                                 br(),
                                 
                                 submitButton("Submit"),
                                 
                                 br(), 
                                 
                                 br()
                                 
                          )
                          
                          
                          
                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Mig_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableMig")
                              
                            )
                            
                            
                          ),
                          
                        ),
                        
                      )),
             
             ####Data_entry_Ns#########
             
             tabPanel("Nitrosative stress",
                      
                      titlePanel("Nitrosative stress stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_ns",
                            
                            label = "Select Species",
                            
                            choices = c(
                              "Achromobacter xylosoxidans Ns (+)",            
                              "Aggregatibacter actinomycetemcomitans Ns (+)",
                              "Enterococcus faecalis Ns (+)",                 
                              "Escherichia coli EPEC Ns (+)",                
                              "Escherichia coli ETEC Ns (+)",                 
                              "Haemophilus influenzae Ns (+)",               
                              "Klebsiella pneumoniae Ns (+)",                 
                              "Listeria monocytogenes Ns (+) (I)",               
                              "Listeria monocytogenes Ns (+) (II)",                
                              "Pseudomonas aeruginosa Ns (-)",               
                              "Pseudomonas aeruginosa Ns (-) Mig (+) ",       
                              "Pseudomonas aeruginosa Ns (+)",  
                              "Pseudomonas aeruginosa Ns (+) Nd (+)",
                              "Pseudomonas aeruginosa Ns (+) Sp (-) (I)",     
                              "Pseudomonas aeruginosa Ns (+) Sp (-) (II)",   
                              "Salmonella enterica Ns (+)",                   
                              "Shigella flexneri Ns (+)",                    
                              "Staphylococcus aureus MRSA252 Ns (-) Nd (+)",  
                              "Staphylococcus aureus MSSA476 Ns (-)",        
                              "Staphylococcus aureus MSSA476 Ns (-) As (+)",  
                              "Staphylococcus epidermidis Ns (+)",           
                              "Staphylococcus epidermidis Ns (+) Tm (-)",     
                              "Streptococcus agalactiae Ns (+)",             
                              "Streptococcus pneumoniae Ns (+)",              
                              "Streptococcus pyogenes Ns (+) Oss (+) "      
                              
                              
                              
                              
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,
                                 
                                 br(), 
                                 
                                 br(),
                                 
                                 submitButton("Submit"),
                                 
                                 br(), 
                                 
                                 br()
                                 
                          )
                          
                          
                          
                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Ns_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableNs")
                              
                            )
                            
                            
                          ),
                          
                        ),
                        
                      )),
             
             ####Data_entry_Oss#########
             
             tabPanel("Osmotic stress",
                      
                      titlePanel("Osmotic stress stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_oss",
                            
                            label = "Select Species",
                            
                            choices = c(
                              "Achromobacter xylosoxidans Oss (+)",            
                              "Acinetobacter baumannii Oss (+)",              
                              "Aggregatibacter actinomycetemcomitans Oss (+)", 
                              "Borrelia burgdorferi Oss (+)",                 
                              "Burkholderia pseudomallei Oss (-)",             
                              "Burkholderia pseudomallei Oss (+)",            
                              "Enterococcus faecalis Oss (+)",                 
                              "Escherichia coli EPEC Oss (+)",                
                              "Escherichia coli EPEC Oss (+) Nd (+)",          
                              "Escherichia coli ETEC Oss (+)",                
                              "Escherichia coli ETEC Oss (+) Mig (+) ",        
                              "Escherichia coli UPEC Oss (+)",                
                              "Francisella tularensis Oss (+)",                
                              "Haemophilus influenzae Oss (+)",              
                              "Helicobacter pylori G27 Oss (+)",               
                              "Helicobacter pylori J99 Oss (+)",              
                              "Klebsiella pneumoniae Oss (+) (I)",             
                              "Klebsiella pneumoniae Oss (+) (II)",           
                              "Klebsiella pneumoniae Oss (+) Tm (+)",          
                              "Listeria monocytogenes Oss (+) Nd (-)",        
                              "Neisseria meningitidis Oss (+)",               
                              "Pseudomonas aeruginosa Oss (-)",               
                              "Pseudomonas aeruginosa Oss (+)",                
                              "Pseudomonas aeruginosa Oss (+) Nd (+)",        
                              "Salmonella enterica Oss (+)",                   
                              "Shigella flexneri Oss (+) (I)",                
                              "Shigella flexneri Oss (+) (II)",               
                              "Shigella flexneri Oss (+) Tm (+)",             
                              "Staphylococcus aureus MRSA252 Oss (+)",         
                              "Staphylococcus aureus MSSA476 Oss (+) (I)",    
                              "Staphylococcus aureus MSSA476 Oss (+) (II)",    
                              "Streptococcus agalactiae Oss (+)",             
                              "Streptococcus agalactiae Oss (+) Nd (-)",       
                              "Streptococcus agalactiae Oss (+) Vic (-)",     
                              "Streptococcus pyogenes Oss (+)",                
                              "Streptococcus pyogenes Oss (+) Ns (+)",        
                              "Vibrio cholerae Oss (+)",                       
                              "Yersinia pseudotuberculosis Oss (+)"   
                              
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,
                                 
                                 br(), 
                                 
                                 br(),
                                 
                                 submitButton("Submit"),
                                 
                                 br(), 
                                 
                                 br()
                                 
                          )
                          
                          
                          
                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Oss_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableOss")
                              
                            )
                            
                            
                          ),
                          
                        ),
                        
                      )),
             
             ####Data_entry_Oxs#########
             
             tabPanel("Oxidative stress",
                      
                      titlePanel("Oxidative stress stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_oxs",
                            
                            label = "Select Species",
                            
                            choices = c(
                              
                              "Achromobacter xylosoxidans Oxs (-)",            
                              "Acinetobacter baumannii  Oxs (+) Li (+)",      
                              "Acinetobacter baumannii Oxs (+)",               
                              "Aggregatibacter actinomycetemcomitans Oxs (+)",
                              "Burkholderia pseudomallei Oxs (+)",             
                              "Campylobacter jejuni Oxs (+) (I)",             
                              "Campylobacter jejuni Oxs (+) (II)",             
                              "Enterococcus faecalis  Oxs (+) Nd (+) Vic(-)", 
                              "Enterococcus faecalis Oxs (+) Vic (-)",         
                              "Escherichia coli EPEC Oxs (+) (I)",            
                              "Escherichia coli EPEC Oxs (+) (II)",            
                              "Escherichia coli ETEC Oxs (+) (I)",            
                              "Escherichia coli ETEC Oxs (+) (II)",            
                              "Escherichia coli UPEC Oxs (+)",                
                              "Helicobacter pylori G27 Oxs (+)",               
                              "Helicobacter pylori J99 Oxs (+)",              
                              "Klebsiella pneumoniae  Oxs (+) Nd (+)",         
                              "Klebsiella pneumoniae Oxs (+)",                
                              "Listeria monocytogenes Oxs (+) As (+)",         
                              "Neisseria meningitidis Oxs (-) Tm (+) ",       
                              "Pseudomonas aeruginosa Oxs (+)",                
                              "Pseudomonas aeruginosa Oxs (+) Mig (+) ",      
                              "Salmonella enterica  Oxs (+) Li (+)",           
                              "Salmonella enterica  Oxs (+) Nd (+)",          
                              "Salmonella enterica Oxs (+)",                   
                              "Shigella flexneri  Oxs (+) Li (+)",            
                              "Shigella flexneri Oxs (+)",                     
                              "Shigella flexneri Oxs (+) Nd (+) ",            
                              "Staphylococcus aureus MRSA252 Oxs (+)",         
                              "Staphylococcus aureus MRSA252 Oxs (+) Vic (+)",
                              "Staphylococcus aureus MSSA476 Oxs (+)",         
                              "Streptococcus agalactiae Oxs (+) (I)",             
                              "Streptococcus agalactiae Oxs (+) (II)",      
                              "Streptococcus pyogenes Oxs (+) (I)",           
                              "Streptococcus pyogenes Oxs (+) (II)",           
                              "Yersinia pseudotuberculosis Oxs (+)"         
                              
                              
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,
                                 
                                 br(), 
                                 
                                 br(),
                                 
                                 submitButton("Submit"),
                                 
                                 br(), 
                                 
                                 br()
                                 
                          )
                          
                          
                          
                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Oxs_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableOxs")
                              
                            )
                            
                            
                          ),
                          
                        ),
                        
                      )),
             ####Data_entry_Tm#########
             
             tabPanel("Temperature stress",
                      
                      titlePanel("Temperature stress stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_tm",
                            
                            label = "Select Species",
                            
                            choices = c(
                              
                              "Achromobacter xylosoxidans Tm (+) (I)",        
                              "Achromobacter xylosoxidans Tm (+) (II)",      
                              "Achromobacter xylosoxidans Tm (+) (III)",      
                              "Aggregatibacter actinomycetemcomitans Tm (+)",
                              "Borrelia burgdorferi Tm (+)",                  
                              "Burkholderia pseudomallei Tm (+)",            
                              "Burkholderia pseudomallei Tm (+) Sp (+)",      
                              "Enterococcus faecalis Tm (+)",                
                              "Escherichia coli EPEC Tm (+)",                 
                              "Escherichia coli ETEC Tm (+)",                
                              "Escherichia coli ETEC Tm (+) Li (+)",          
                              "Escherichia coli UPEC Tm (+)",                
                              "Francisella tularensis Tm (+)",                
                              "Haemophilus influenzae Tm (+) (I)",           
                              "Haemophilus influenzae Tm (+) (II)",           
                              "Klebsiella pneumoniae  Tm (+) Oss (+)",       
                              "Klebsiella pneumoniae Tm (+)",                 
                              "Listeria monocytogenes Tm (+) (I)",           
                              "Listeria monocytogenes Tm (+) (II)",           
                              "Neisseria meningitidis Tm (+)",               
                              "Neisseria meningitidis Tm (+) Oxs (-)",        
                              "Pseudomonas aeruginosa Tm (+) As (+)",        
                              "Salmonella enterica Tm (+)",                   
                              "Shigella flexneri Tm (+)",                    
                              "Shigella flexneri Tm (+) Oss (+) ",            
                              "Shigella flexneri Tm (+) Vic (+) ",           
                              "Staphylococcus aureus MRSA252 Tm (+)",         
                              "Staphylococcus aureus MSSA476 Tm (+)",        
                              "Staphylococcus epidermidis Tm (-)", 
                              "Staphylococcus epidermidis Tm (-) As (+)",           
                              "Staphylococcus epidermidis Tm (-) Mig (+) ",  
                              "Staphylococcus epidermidis Tm (-) Nd (+)",     
                              "Staphylococcus epidermidis Tm (-) Ns (+) ",   
                              "Staphylococcus epidermidis Tm (+)", 
                              
                              "Staphylococcus epidermidis Tm (+) Sp (+)",    
                              "Streptococcus agalactiae Tm (+)",              
                              "Streptococcus agalactiae Tm (+) Nd (-)",      
                              "Streptococcus pneumoniae Tm (+) Vic (+)",      
                              "Streptococcus pyogenes Tm (+) Nd (+)",   
                              
                              "Vibrio cholerae Tm (+)",                       
                              "Yersinia pseudotuberculosis Tm (+)"   
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,
                                 
                                 br(), 
                                 
                                 br(),
                                 
                                 submitButton("Submit"),
                                 
                                 br(), 
                                 
                                 br()
                                 
                          )
                          
                          
                          
                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Tm_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableTm")
                              
                            )
                            
                            
                          ),
                          
                        ),
                        
                      )),
             ####Data_entry_Vic#########
             
             tabPanel("Virulence inducing condition",
                      
                      titlePanel("Virulence inducing condition stimulons per species"),
                      
                      
                      
                      sidebarLayout(
                        
                        
                        
                        sidebarPanel(
                          
                          # selector for district
                          
                          selectInput(
                            
                            inputId = "select_vic",
                            
                            label = "Select Species",
                            
                            choices = c(
                              "Aggregatibacter actinomycetemcomitans Vic (+) (I)",   
                              "Aggregatibacter actinomycetemcomitans Vic (+) (II)", 
                              "Aggregatibacter actinomycetemcomitans Vic (+) (III)", 
                              "Burkholderia pseudomallei Vic (+)",                  
                              "Enterococcus faecalis  Vic (+) Mig (+)",              
                              "Enterococcus faecalis Vic (-) Oxs (+) ",             
                              "Enterococcus faecalis Vic(-) Nd (+) Oxs (+)",         
                              "Escherichia coli EPEC Vic (+)",                      
                              "Escherichia coli ETEC Vic (+) (I)",                   
                              "Escherichia coli ETEC Vic (+) (II)",                 
                              "Escherichia coli UPEC  Vic (+) Mig (+)",             
                              
                              "Haemophilus influenzae Vic (+)",                      
                              "Haemophilus influenzae Vic (+) Mig (-)",             
                              "Haemophilus influenzae  Vic (+) Mig (+)",
                              
                              "Helicobacter pylori J99 Vic (+)",                     
                              "Klebsiella pneumoniae Vic (+)",                      
                              "Listeria monocytogenes Vic (+) Sp (+)",               
                              "Pseudomonas aeruginosa Vic (+)",                     
                              "Salmonella enterica Vic (+)",                         
                              "Shigella flexneri Vic (+) (I)",                      
                              "Shigella flexneri Vic (+) (II)",                      
                              "Shigella flexneri Vic (+) Tm (+)",                   
                              "Staphylococcus aureus MRSA252 Vic (+)",               
                              "Staphylococcus aureus MRSA252 Vic (+) Oxs (+) ",     
                              "Staphylococcus aureus MSSA476 Vic (+)",               
                              "Streptococcus agalactiae Vic (-) Oss (+) ",          
                              
                              "Streptococcus agalactiae Vic (+) Nd (-)",            
                              "Streptococcus pneumoniae Vic (-)",                    
                              "Streptococcus pneumoniae Vic (+) Mig (-)",           
                              "Streptococcus pneumoniae Vic (+) Tm (+) ",            
                              "Streptococcus pyogenes Vic (-)",                     
                              "Vibrio cholerae Vic (+)",                             
                              "Yersinia pseudotuberculosis Vic (-) Nd (+)",         
                              "Yersinia pseudotuberculosis Vic (+) (I)",             
                              "Yersinia pseudotuberculosis Vic (+) (II)",           
                              "Yersinia pseudotuberculosis Vic (+) Mig (-)"        
                              
                              
                              
                              
                              
                              
                              
                              
                            )),
                          
                          
                          
                          column(5,
                                 
                                 br(), 
                                 
                                 br(),
                                 
                                 submitButton("Submit"),
                                 
                                 br(), 
                                 
                                 br()
                                 
                          )
                          
                          
                          
                          
                          
                          
                          
                        ),
                        
                        
                        
                        
                        
                        mainPanel(
                          
                          tabsetPanel(
                            
                            type = "tabs",
                            
                            # epicurve goes here
                            
                            tabPanel(
                              
                              "Stimulon expression under different conditions",  plotOutput("Vic_boxplot"),
                              br(),
                              hr(),
                              p(" The stress abbreviations used in the plot are listed as follows:"),
                              tags$ul(
                                tags$li(tags$b("As"), " - Acidic stress"),
                                tags$li(tags$b("Bs"), " - Bile stress"),
                                tags$li(tags$b("Ctrl"), " - Control"),
                                tags$li(tags$b("Li"), " - Low iron"),
                                tags$li(tags$b("Mig"), " - Microaerophilic growth"),
                                tags$li(tags$b("Nd"), " - Nutritional downshift"),
                                tags$li(tags$b("Ns"), " - Nitrosative stress"),
                                tags$li(tags$b("Oss"), " - Osmotic stress"),
                                tags$li(tags$b("Oxs"), " - Oxidative stress"),
                                tags$li(tags$b("Sp"), " - Stationary phase"),
                                tags$li(tags$b("Tm"), " - Temperature"),
                                tags$li(tags$b("Vic"), " - Virulence inducing condition"),
                              )
                              
                            ),
                            
                            tabPanel(
                              
                              "Stimulon gene annotation data",
                              
                              dataTableOutput("mytableVic")
                              
                            )
                            
                            
                          ),
                          
                        ),
                        
                      ))
             
             
  ),
  navbarMenu("Tutorials",
           tabPanel("Tutorial 1: Step by step network construction",
           #          #h2("Tutorial 1", align="center", color="darkblue"),
                   h4("Case of study 1: Construct a gene co-expression network, how to do it?  "),
           # tags$video(id="video3", type = "video/mp4",src = "Tutorial_1.mp4", autoplay=FALSE,
           #            muted=TRUE,
           #           playsinline=TRUE,
           #           loop=FALSE,
           #           width = 1300,
           #           top="2px",
           #           controls = "controls")
           embed_url("https://youtu.be/QaN2Wfc_-do")
           #HTML('<iframe width="560" height="315" src="https://youtu.be/QaN2Wfc_-do" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
           #tags$iframe(width="560", height="315", src="https://youtu.be/QaN2Wfc_-do", frameborder="0", allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture", allowfullscreen=NA)
           ),
           #h2("Tutorial 2", align="center", color="darkblue"),
           tabPanel("Tutorial 2: Finding information of a gene of interest",
           h4(" Case of study 2: Finding genes co-expressed with streptococcal pyrogenic exotoxin B (SpeB) in Streptococcus pyogenes"),
           embed_url("https://youtu.be/hnrDZSNlYkU")
           # tags$video(id="video4", type = "video/mp4",src = "Tutorial_2.mp4", autoplay=FALSE,
           #            muted=TRUE,
           #            playsinline=TRUE,
           #            loop=FALSE,
           #            width = 1300,
           #            top="2px",
           #            controls = "controls")
           ),
           
           #h2("Tutorial 3", align="center", color="darkblue"),
           tabPanel("Tutorial 3: Core transcriptome network construction",
           h4(" Case of study 3: Co-expression of common genes in the species: E. coli, H. pylori, S. aureus "),
           
           # tags$video(id="video5", type = "video/mp4",src = "Tutorial_3.mp4", autoplay=FALSE,
           #            muted=TRUE,
           #            playsinline=TRUE,
           #            loop=FALSE,
           #            width = 1300,
           #            top="2px",
           #            controls = "controls")
           embed_url("https://youtu.be/oaOCilpGiTk")
           ),
           tabPanel("Tutorial 4: How to wotk with the stimulons?",
           #h2("Tutorial 4", align="center", color="darkblue"),
           h4(" Case of study 4: Stress specific stimulon "),
           # tags$video(id="video6", type = "video/mp4",src = "Tutorial_4.mp4", autoplay=FALSE,
           #            muted=TRUE,
           #            playsinline=TRUE,
           #            loop=FALSE,
           #            width = 1300,
           #            top="2px",
           #            controls = "controls")
           embed_url("https://youtu.be/viSzO7O7YZs")
           ) 
           
         #  img(height = 800, width = 900,src='Co-pathogenex homepage.png', align = "center"),
  ),
  
  tabPanel("General info")
    )
    
)

            

        
server <-function(input, output,session) {
    
    # Return the requested dataset ----
    datasetInput <- reactive({
        switch(input$dataset,
               "Escherichia coli EPEC" = (Escherichia_coli_EPEC_0127_H6_E2348_69_expr),
               "Escherichia coli ETEC" = (Escherichia_coli_ETEC_H10407_expr),
               "Escherichia coli UPEC" = (Escherichia_coli_UPEC_536_expr),
               "Helicobacter pylori G27" = (Helicobacter_pylori_G27_expr),
               "Salmonella enterica"=(Salmonella_expr),
               "Helicobacter pylori J99"=(Helicobacter_pylori_J99_expr),
               "Acinetobacter baumannii"=(Acinetobacter_baumannii_expr),
               "Aggregatibacter actinomycetemcomitans"=(Aggregatibacter_actinomycetemcomitans_expr),
               "Achromobacter xylosoxidans"=(Achromobacter_xylosoxidans_expr),
               "Burkholderia pseudomallei"= (Burkholderia_pseudomallei_expr),
                "Borrelia burgdorferi" = (Borrelia_burgdorferi_expr),
                "Campylobacter jejuni" = (Campylobacter_jejuni_expr),
                "Francisella tularensis" = (Francisella_tularensis_expr),
               "Haemophilus influenzae" = (Haemophilus_influenzae_expr),
               "Klebsiella pneumoniae" = (Klebsiella_pneumoniae_expr),
               "Legionella pneumophila"= (Legionella_pneumophila_expr),
               "Listeria monocytogenes"= (Listeria_monocytogenes_expr),
               "Mycobacterium tuberculosis" = (Mycobacterium_tuberculosis_expr),
               "Neisseria gonorrhoeae" = (Neisseria_gonorrhoeae_expr),
               "Neisseria meningitidis" = (Neisseria_meningitidis_expr),
               "Pseudomonas aeruginosa" = (Pseudomonas_aeruginosa_expr),
               "Staphylococcus aureus MRSA252" =(Staphylococcus_aureus_MRSA252_expr),
               "Staphylococcus aureus MSSA476" = (Staphylococcus_aureus_MSSA476_expr),
               "Staphylococcus epidermidis" = (Staphylococcus_epidermidis_expr),
               "Streptococcus agalactiae" = (Streptococcus_agalactiae_expr),
               "Streptococcus pneumoniae" = (Streptococcus_pneumoniae_expr),
               "Streptococcus pyogenes"  = (Streptococcus_pyogenes_expr),
               "Streptococcus suis"  = (Streptococcus_suis_expr),
               "Yersinia pseudotuberculosis" = (Yersinia_pseudot_TPM_expr),
               "Enterococcus faecalis" = (Enterococcus_faecalis_expr),
               "Vibrio cholerae" = (Vibrio_cholerae_expr),
               "Shigella flexneri" = (Shigella_flexneri_expr),
               
               )
    })
    
    distPlot <- reactive({
        dataset <- datasetInput()
        mat_correlation= 1-cor(t(dataset))
        sampleTree = hclust(dist(mat_correlation), method="ward.D2")
        #par(mar = c(0,4,2,0))
        plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 1)
        abline(h = 40000, col = "red")  
        
        
    })
    
    
    distPlot1 <- reactive({
        dataset <- datasetInput()
        gsg = goodSamplesGenes(dataset,verbose = 3)
        if (!gsg$allOK)
        {
            dataset=dataset[gsg$goodSamples,gsg$goodGenes]
        }
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        sft = pickSoftThreshold(dataset, powerVector = powers, verbose = 5)
        par(mfrow = c(1,2));
        cex1 = 0.6;
        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
        text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
        abline(h=0.80,col="red")
        plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
        
        
    })
    
    distPlot2 <- reactive({
        dataset <- datasetInput()
        gsg = goodSamplesGenes(dataset,verbose = 3)
        if (!gsg$allOK)
        {
            dataset=dataset[gsg$goodSamples,gsg$goodGenes]
        }
        powerin <- input$powerin
        mingene <- input$mingene
        net = blockwiseModules(dataset, power = powerin, TOMType = "unsigned", minModuleSize = mingene, reassignThreshold = 0, 
                               mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, 
                               saveTOMFileBase = "TOM", verbose = 3)
        #table(net$colors)
        mergedColors = labels2colors(net$colors)
        plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                            dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
        
    })
    
    # Return the requested analysis ----
    #graphInput <- reactive({
        #switch(input$graph,
               #"Sample overview" = distPlot(),distPlot1(),
              # "Sample overview" = distPlot1(),
              # "Network construction" = distPlot2()
                # )
  #  })
    
    # Return the requested analysis ----
   # summaryInput <- reactive({
  #      switch(input$graph,
  #             "Sample cluster" = summary1(),
  #             "Power analysis" = summary2(),
   #            "Network construction" = summary3()
  #      )
   # })
    summary1 <- reactive({
    dataset <- datasetInput()
    k <- dim(dataset)[2]
    f <- dim(dataset)[1]
    summary(t(dataset))
    paste0(input$dataset,":","Number of genes="," ",as.character(k),"-------------","Number of samples="," ",as.character(f))
    
        })
    
    summary2 <- reactive({
        dataset <- datasetInput()
        # gsg = goodSamplesGenes(dataset,verbose = 3)
        # if (!gsg$allOK)
        # {
        #     dataset=dataset[gsg$goodSamples,gsg$goodGenes]
        # }
        # 
        # powers = c(c(1:10), seq(from = 12, to=20, by=2))
        input$dataset
    })
    
    summary3 <- reactive({
        dataset <- datasetInput()
        gsg = goodSamplesGenes(dataset,verbose = 3)
        if (!gsg$allOK)
        {
            dataset=dataset[gsg$goodSamples,gsg$goodGenes]
        }
        powerin <- input$powerin
        mingene <- input$mingene
        net = blockwiseModules(dataset, power = powerin, TOMType = "unsigned", minModuleSize = mingene, reassignThreshold = 0, 
                               mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, 
                               saveTOMFileBase = "TOM", verbose = 3)
        #table(net$colors)
        as.data.frame(table(labels2colors(net$colors)))
    })
    
    summary4 <- reactive({
      dataset <- datasetInput()
      powerin<- input$powerin
      mingene<- input$mingene
      paste0(input$dataset,":","Power="," ",as.character(powerin),",","Number of genes="," ",as.character(mingene))
      
    })
    
    summary5 <- reactive({
      dataset <- datasetInput()
      powerin<- input$powerin
      mingene<- input$mingene
      paste0(input$dataset,":","Power="," ",as.character(powerin),",","Number of genes="," ",as.character(mingene))
      
    })
    
    # Generate a summary of the dataset ----
        output$summary_datasetin1 <- renderPrint({
            summary1()
        
    })
    
    output$summary_datasetin2 <- renderPrint({
      showModal(modalDialog(title = "Network construction:","It can take a minute, press dismiss if you want to change the parameters", footer=modalButton("Dismiss")))
      on.exit(removeModal())
      
      
      summary3()
      
    })
    
    output$summary_datasetin3 <- renderPrint({
      summary2()
      
    })
    
    output$summary_datasetin4 <- renderPrint({
      summary4()
      
    })
    
    output$summary_datasetin5 <- renderPrint({
      summary5()
      
    })
    
    output$clusterone =  renderPlot({distPlot()
      }) 
    output$powerplot1 =  renderPlot({distPlot1()
    }) 
    
    output$selected_graph <- renderPlot({
      showModal(modalDialog(title = "Network construction:","Almost there..", footer=NULL))
      on.exit(removeModal()) 
        distPlot2()
    })
    
    
    # Show the first "n" observations ----
    output$view <- renderTable({
        head(datasetInput(), n = input$obs)
    })
    
    dat1 <- reactive({
        dataset <- datasetInput()
        gsg = goodSamplesGenes(dataset,verbose = 3)
        if (!gsg$allOK)
        {
            dataset=dataset[gsg$goodSamples,gsg$goodGenes]
        }
        powerin <- input$powerin
        mingene <- input$mingene
        net = blockwiseModules(dataset, power = powerin, TOMType = "unsigned", minModuleSize = mingene, reassignThreshold = 0, 
                               mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, 
                               saveTOMFileBase = "TOM", verbose = 3)
        #table(net$colors)
        mergedColors = labels2colors(net$colors)
        moduleLabels = net$colors
        moduleColors = labels2colors(net$colors)
        MEs = net$MEs
        datME=moduleEigengenes(dataset,moduleColors)$eigengenes
        datME_round=round(datME,2)
        table_1_inspect=as.data.frame(table(labels2colors(net$colors)))
        colnames(table_1_inspect) <- NULL
        colnames(table_1_inspect) <- c("Module color","Number of genes")
        new_da_table = cbind(table_1_inspect,t(datME_round))
        #colnames(new_da_table) <- c("Module color","Number of genes","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3",
        #                               "Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3",
        #                              "Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(new_da_table)
    })
    
    dat2 <- reactive({
        dataset <- datasetInput()
        gsg = goodSamplesGenes(dataset,verbose = 3)
        if (!gsg$allOK)
        {
            dataset=dataset[gsg$goodSamples,gsg$goodGenes]
        }
        powerin <- input$powerin
        mingene <- input$mingene
        net = blockwiseModules(dataset, power = powerin, TOMType = "unsigned", minModuleSize = mingene, reassignThreshold = 0, 
                               mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, 
                               saveTOMFileBase = "TOM", verbose = 3)
        #table(net$colors)
        Module_Colors = labels2colors(net$colors)
        Locus_Tag= names(net$colors)
        dataset_transpose=t(scale(dataset))
        
        #rownames(table_2_inspect) <- NULL
        
        #colnames(table_2_inspect) <- NULL
        #colnames(table_2_inspect) <- c("Module color","Gene ID")
        probes = Locus_Tag
        Genes=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
        PGFam=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
        table_2_inspect= as.data.frame(cbind(Module_Colors,Locus_Tag,Genes,PGFam,dataset_transpose))
        return(as.data.frame(table_2_inspect))
    })
   
    color_from_middle <- function (data,color1,color2) 
    {
        max_val=max(abs(data))
        JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                   max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    output$mytable <- DT::renderDataTable({showModal(modalDialog(title = "Module expression:",HTML("1) Grouping genes with similar expression in modules.<br>  2) Computing module eigengene to characterize a modular expression. <br> 3) Each module is named using a color "), footer=NULL))
                                           on.exit(removeModal()) 
                                           dat1()
                                           return(datatable(dat1(), options = list(paging=FALSE, scrollX = FALSE,lengthChange = T, searching = FALSE),rownames=FALSE) %>% formatStyle(names(dat1()[-1,-2]), background= color_from_middle(dat1()[3:ncol(dat1())],'blue','red'), backgroundSize = '98% 78%',
                                                                                                                                                                backgroundRepeat = 'no-repeat', backgroundPosition = 'center'))                                                                                                                         

                                          })
    js <- c(
        "function(settings){",
        "  var instance = settings.oInstance;",
        "  var table = instance.api();",
        "  var input = instance.parent().find('.dataTables_filter input');",
        "  input.off('keyup search input').on('keyup', function(){",
        "    var keyword = '\\\\b' + input.val() + '\\\\b';",
        "    table.search(keyword, true, false).draw();",
        "  });",
        "}"
    )
    
    
    
    output$mytable2 <- DT::renderDataTable({showModal(modalDialog(title = "Gene expression:","Computing z-score per gene...", footer=NULL))
                                        on.exit(removeModal()) 
                                        dat2()
                                        return(datatable(dat2(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(dat2()[-1:-4]), backgroundColor = styleInterval(brks, clrs)))
                                          #options = list(paging=FALSE, scrollX = FALSE), 
                                          #rownames=TRUE, 
                                          #filter = "top"
                                          })   

    
    
    
       
    # observe({
    #     updateTextInput(session, inputId = "myresultscolor", value = input$textforcolor) 
    #     
    # })   
    
    
    # distPlotalluvial <- reactive({
    #     Escherichia_coli_UPEC_536_GOclear <- read_csv("Escherichia_coli_UPEC_536_GOclear.txt",col_names = FALSE)
    # 
    #     Escherichia_coli_UPEC_536_GOclear_df= as.data.frame(Escherichia_coli_UPEC_536_GOclear[1:130,9:11])
    #     colnames(Escherichia_coli_UPEC_536_GOclear_df)= c("Parent", "Child1", "Child2")
    #     return(Escherichia_coli_UPEC_536_GOclear_df)
    #           # alluvial_wide( data = Escherichia_coli_UPEC_536_GOclear_df[1:1000,2:11]
    #             #                                      , fill_by = 'all_flows'
    #                 #                                  , stratum_label_size = 2.0
    #               #                                    , stratum_width = 1/20 ) 
    #               # parcats(palluvial, marginal_histograms = TRUE, data_input = Escherichia_coli_UPEC_536_GOclear_df[1:1000,2:4])
    # 
    #     
    # 
    # 
    #     
    #     
    # })
    
        
    # output$selected_graph2 <- render_parcats({
    #     nodes <- distPlotalluvial()
    #     #req(nodes)
    #     p = alluvial_wide(nodes, max_variables = ncol(nodes), fill_by = 'first_variable')
    #     parcats(p, marginal_histograms = TRUE, data_input = nodes)
    # })
    
    graph2Input <- reactive({
        switch(input$graph2,
               "Gene expression heatmap"= distPlot3(), 
               #"Gene expression box plots"= distPlot4(), 
               "Module eigengene barplot"= distPlot5()
               
              
        )
    })
    
    
    
    distPlot3 <- reactive({
      # validate(
      #   need(input$data != "", "Please enter a module color in the box, ex: blue and hit submit")
      # )
      
        dataset <- datasetInput()
        gsg = goodSamplesGenes(dataset,verbose = 3)
        if (!gsg$allOK)
        {
            dataset=dataset[gsg$goodSamples,gsg$goodGenes]
        }
        powerin <- input$powerin
        mingene <- input$mingene
        net = blockwiseModules(dataset, power = powerin, TOMType = "unsigned", minModuleSize = mingene, reassignThreshold = 0, 
                               mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, 
                               saveTOMFileBase = "TOM", verbose = 3)
        #table(net$colors)
    
        moduleColors = labels2colors(net$colors)
        MEs = net$MEs
        datME=moduleEigengenes(dataset,moduleColors)$eigengenes
        colorh1=moduleColors
        which.module=input$textforcolor
        #par(mar=c(1, 1, 1, 1))
        col1 <- colorRampPalette(brewer.pal(7, "RdYlBu"))(256)
        col <- rev(col1)
        #plotMat(t(scale(dataset[,colorh1==which.module]) ),nrgcols=30,rlabels=NULL,rcols="black",clabels=rownames(dataset),cex.lab=0.7)
        heatmap(t(scale(dataset[,colorh1==which.module])), scale = "row", col =  col, Colv=NA, margins=c(6,10), keep.dendro = TRUE, main= as.character(input$dataset))
        legend(x="topright", legend=c("-3", "-2","-1","0","+1", "+2", "+3"),fill=c("#4575b4","#91bfdb","#e0f3f8","#ffffbf","#fee090","#fc8d59","#d73027"), title="Std. Dev.")
        #heatmaply(t(scale(dataset[,colorh1==which.module]) ))
    })
    
    
    distPlot5 <- reactive({
        dataset <- datasetInput()
        
        gsg = goodSamplesGenes(dataset,verbose = 3)
        if (!gsg$allOK)
        {
            dataset=dataset[gsg$goodSamples,gsg$goodGenes]
        }
        powerin <- input$powerin
        mingene <- input$mingene
        net = blockwiseModules(dataset, power = powerin, TOMType = "unsigned", minModuleSize = mingene, reassignThreshold = 0, 
                               mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, 
                               saveTOMFileBase = "TOM", verbose = 3)
        #table(net$colors)
        
        moduleColors = labels2colors(net$colors)
        MEs = net$MEs
        datME=moduleEigengenes(dataset, moduleColors)$eigengenes
        colorh1=moduleColors
        which.module=input$textforcolor
        #mar=c(0.3, 15.5, 8, 2)
        #plotMat(t(scale(dataset[,colorh1==which.module]) ),nrgcols=30,rlabels=NULL,rcols="black",clabels=rownames(dataset),cex.lab=0.7)
        ME=datME[, paste("ME",which.module, sep="")]
        colors_bar= ifelse(ME<0, "blue", "red")
        #par(mar=c(5, 14.2, 0, 0.85))
        par(mar=c(11,4,4,4))
        barplot(ME, col=t(colors_bar), main=as.character(input$dataset),ylab="eigengene expression", names.arg= rownames(dataset), las=2,cex.axis=1.5, cex.lab=1.5   )
        
    })
    
    
    
    
    output$myresultscolor <- renderPlot({
      
      showModal(modalDialog("Loading data...", footer=NULL))
      on.exit(removeModal())

        graph2Input() }, height = 500, width = 700 )
    
      #graph2Input() })
    
    output$download_heatmap <- downloadHandler(
      
      filename = function() { paste(input$textforcolor, 'Heatmap.png', sep='') },
      
      content = function(file) {
        showModal(modalDialog("Loading heatmap....", footer=NULL))
        on.exit(removeModal())
        png(filename = file, width = 600, height = 400, units='mm', res = 300)
        dataset <- datasetInput()
        gsg = goodSamplesGenes(dataset,verbose = 3)
        if (!gsg$allOK)
        {
          dataset=dataset[gsg$goodSamples,gsg$goodGenes]
        }
        powerin <- input$powerin
        mingene <- input$mingene
        net = blockwiseModules(dataset, power = powerin, TOMType = "unsigned", minModuleSize = mingene, reassignThreshold = 0, 
                               mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, 
                               saveTOMFileBase = "TOM", verbose = 3)
        #table(net$colors)
        
        moduleColors = labels2colors(net$colors)
        MEs = net$MEs
        datME=moduleEigengenes(dataset,moduleColors)$eigengenes
        colorh1=moduleColors
        which.module=input$textforcolor
        #par(mar=c(1, 1, 1, 1))
        col1 <- colorRampPalette(brewer.pal(7, "RdYlBu"))(256)
        col <- rev(col1)
        #plotMat(t(scale(dataset[,colorh1==which.module]) ),nrgcols=30,rlabels=NULL,rcols="black",clabels=rownames(dataset),cex.lab=0.7)
        heatmap(t(scale(dataset[,colorh1==which.module])), scale = "row", col =  col, Colv=NA, margins=c(6,10), keep.dendro = TRUE)
        legend(x="topright", legend=c("-3", "-2","-1","0","+1", "+2", "+3"),fill=c("#4575b4","#91bfdb","#e0f3f8","#ffffbf","#fee090","#fc8d59","#d73027"), title="Std. dev.")
        
               dev.off()
      }
      
    ) 
 
    output$download_barplot <- downloadHandler(
      
      filename = function() { paste(input$textforcolor, 'MEbarplot.png', sep='') },
      
      content = function(file) {
        showModal(modalDialog("Loading ME barplot...", footer=NULL))
        on.exit(removeModal())
        png(filename = file, width = 350, height = 300, units='mm', res = 300)
        dataset <- datasetInput()
        
        gsg = goodSamplesGenes(dataset,verbose = 3)
        if (!gsg$allOK)
        {
          dataset=dataset[gsg$goodSamples,gsg$goodGenes]
        }
        powerin <- input$powerin
        mingene <- input$mingene
        net = blockwiseModules(dataset, power = powerin, TOMType = "unsigned", minModuleSize = mingene, reassignThreshold = 0, 
                               mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, 
                               saveTOMFileBase = "TOM", verbose = 3)
        #table(net$colors)
        
        moduleColors = labels2colors(net$colors)
        MEs = net$MEs
        datME=moduleEigengenes(dataset, moduleColors)$eigengenes
        colorh1=moduleColors
        which.module=input$textforcolor
        #mar=c(0.3, 15.5, 8, 2)
        #plotMat(t(scale(dataset[,colorh1==which.module]) ),nrgcols=30,rlabels=NULL,rcols="black",clabels=rownames(dataset),cex.lab=0.7)
        ME=datME[, paste("ME",which.module, sep="")]
        colors_bar= ifelse(ME<0, "blue", "red")
        #par(mar=c(5, 14.2, 0, 0.85))
        par(mar=c(11,4,4,4))
        barplot(ME, col=t(colors_bar), main="",ylab="eigengene expression", names.arg= rownames(dataset), las=2, cex.axis=1.5, cex.lab=1.5  )
        
      
        dev.off()
      }
      
    )   
    
    
#####Core####
    datasetInputCore <- reactive({
        switch(input$dataset_core,
               "UPEC vs EPEC" = (UPEC_EPEC_multiExpr),
               "UPEC vs ETEC" = (UPEC_ETEC_multiExpr),
               "EPEC vs ETEC" = (EPEC_ETEC_multiExpr),
               "HP_J99 vs HP_G27" = (HP_J99_HP_G27_multiExpr),
               "MSSA476 vs MRSA252" = (MSSA476_MRSA252_multiExpr),
               
               
        )
    })
    
    distPlotCore <- reactive({
        multiExpr <- datasetInputCore()
        setLabels=setLabels_analysis[[eval(input$dataset_core)]]
        sampleTrees = list()
        for (set in 1:nSets)
        {
            sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
        }
        par(mfrow=c(2,1))
        par(mar = c(0, 4, 2, 0))
        for (set in 1:nSets)
            plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
                 xlab="", sub="", cex = 0.6)
        
        
    })
    
    
    distPlot1Core <- reactive({
        multiExpr<- datasetInputCore()
        exprSize = checkSets(multiExpr)
        setLabels=setLabels_analysis[[eval(input$dataset_core)]]
        gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
        gsg$allOK
        
        if (!gsg$allOK)
        {
          # Print information about the removed genes:
          if (sum(!gsg$goodGenes) > 0)
            printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                                      collapse = ", ")))
          for (set in 1:exprSize$nSets)
          {
            if (sum(!gsg$goodSamples[[set]]))
              printFlush(paste("In set", setLabels[set], "removing samples",
                               paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
            # Remove the offending genes and samples
            multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
          }
          # Update exprSize
          exprSize = checkSets(multiExpr)
        }
        
        powers = c(seq(4,10,by=1), seq(12,20, by=2));
        # Initialize a list to hold the results of scale-free analysis
        powerTables = vector(mode = "list", length = nSets);
        # Call the network topology analysis function for each set in turn
        for (set in 1:nSets)
          powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                             verbose = 2)[[2]]);
        collectGarbage();
        # Plot the results:
        colors = c("black", "red")
        # Will plot these columns of the returned scale free analysis tables
        plotCols = c(2,5,6,7)
        colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
                     "Max connectivity");
        # Get the minima and maxima of the plotted points
        ylim = matrix(NA, nrow = 2, ncol = 4);
        for (set in 1:nSets)
        {
          for (col in 1:length(plotCols))
          {
            ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
            ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
          }
        }
        # Plot the quantities in the chosen columns vs. the soft thresholding power
        par(mfcol = c(2,2));
        par(mar = c(4.2, 4.2 , 2.2, 0.5))
        cex1 = 0.6;
        for (col in 1:length(plotCols)) for (set in 1:nSets)
        {
          if (set==1)
          {
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                 xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
                 main = colNames[col]);
            addGrid();
          }
          if (col==1)
          {
            text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                 labels=powers,cex=cex1,col=colors[set]);
          } else
            text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
                 labels=powers,cex=cex1,col=colors[set]);
          if (col==1)
           
          {
            legend("bottomright", legend =  setLabels, col = colors, pch = 15) ;
          } else
            legend("topright", legend =  setLabels, col = colors, pch = 15) ;
        }
        
    })
    
    distPlot2Core <- reactive({
      multiExpr <- datasetInputCore()
      exprSize = checkSets(multiExpr)
      gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
      gsg$allOK
      
      if (!gsg$allOK)
      {
        # Print information about the removed genes:
        if (sum(!gsg$goodGenes) > 0)
          printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                                    collapse = ", ")))
        for (set in 1:exprSize$nSets)
        {
          if (sum(!gsg$goodSamples[[set]]))
            printFlush(paste("In set", setLabels[set], "removing samples",
                             paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
          # Remove the offending genes and samples
          multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
        }
        # Update exprSize
        exprSize = checkSets(multiExpr)
      }
      
      
      powerin_core <- input$powerin_core
      mingene_core <- input$mingene_core
      
      net = blockwiseConsensusModules(
        multiExpr, power = powerin_core, minModuleSize = mingene_core, deepSplit = 2,
        pamRespectsDendro = FALSE, 
        mergeCutHeight = 0.25, numericLabels = TRUE,
        minKMEtoStay = 0,
        saveTOMs = FALSE, verbose = 5)
      consMEs = net$multiMEs;
      moduleLabels = net$colors;
      moduleColors = labels2colors(moduleLabels)
      consTree = net$dendrograms[[1]]; 
      plotDendroAndColors(consTree, moduleColors,
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = "Consensus gene dendrogram and module colors")
        
    })
    
    # # Return the requested analysis ----
    # graphInputCore <- reactive({
    #     switch(input$graph_core,
    #            "Sample cluster" = distPlotCore(),
    #            "Power analysis" = distPlot1Core(),
    #            "Network construction" = distPlot2Core()
    #     )
    # })
    # 
    # # Return the requested analysis ----
    # summaryInputCore <- reactive({
    #     switch(input$graph_core,
    #            "Sample cluster" = summary1Core(),
    #            "Power analysis" = summary2Core(),
    #            "Network construction" = summary3Core()
    #     )
    # })
    summary1Core <- reactive({
        multiExpr <- datasetInputCore()
        exprSize = checkSets(multiExpr)
        k <- exprSize$nGenes
        #f <- exprSize$nSamples
        #summary(t(dataset))
        paste0(input$dataset_core,":","Number of common  genes="," ",as.character(k))
        
    })
    
    summary2Core <- reactive({
        multiExpr <- datasetInputCore()
        # exprSize = checkSets(multiExpr)
        # gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
        # gsg$allOK
        # 
        # if (!gsg$allOK)
        # {
        #   # Print information about the removed genes:
        #   if (sum(!gsg$goodGenes) > 0)
        #     printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
        #                                               collapse = ", ")))
        #   for (set in 1:exprSize$nSets)
        #   {
        #     if (sum(!gsg$goodSamples[[set]]))
        #       printFlush(paste("In set", setLabels[set], "removing samples",
        #                        paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        #     # Remove the offending genes and samples
        #     multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
        #   }
        #   # Update exprSize
        #   exprSize = checkSets(multiExpr)
        # }
        # 
        # powers = c(seq(4,10,by=1), seq(12,20, by=2));
        # # Initialize a list to hold the results of scale-free analysis
        # powerTables = vector(mode = "list", length = nSets);
        # # Call the network topology analysis function for each set in turn
        # for (set in 1:nSets)
        #   powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
        #                                                      verbose = 2)[[2]]);
        # collectGarbage();
        # #return(paste0(input$dataset_core,":",powerTables[[1]], powerTables[[2]]))
        return(paste0(input$dataset_core))
    })
    
    summary3Core <- reactive({
        multiExpr <- datasetInputCore()
        exprSize = checkSets(multiExpr)
        gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
        gsg$allOK
        
        if (!gsg$allOK)
        {
          # Print information about the removed genes:
          if (sum(!gsg$goodGenes) > 0)
            printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                                      collapse = ", ")))
          for (set in 1:exprSize$nSets)
          {
            if (sum(!gsg$goodSamples[[set]]))
              printFlush(paste("In set", setLabels[set], "removing samples",
                               paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
            # Remove the offending genes and samples
            multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
          }
          # Update exprSize
          exprSize = checkSets(multiExpr)
        }

     
        powerin_core <- input$powerin_core
        mingene_core <- input$mingene_core
        
        net = blockwiseConsensusModules(
          multiExpr, power = powerin_core, minModuleSize = mingene_core, deepSplit = 2,
          pamRespectsDendro = FALSE, 
          mergeCutHeight = 0.25, numericLabels = TRUE,
          minKMEtoStay = 0,
          saveTOMs = FALSE, verbose = 5)
        consMEs = net$multiMEs;
        moduleLabels = net$colors;
        moduleColors = labels2colors(moduleLabels)
        #table(net$colors)
        as.data.frame(table(labels2colors(net$colors)))
    })
    
    # Generate a summary of the dataset ----
    output$summary_datasetin_core1 <- renderPrint({
      summary1Core()
        
    })
    
    
    
    
    output$summary_datasetin_core3 <- renderPrint({
      summary2Core()
      
    })
    
    output$summary_datasetin_core2 <- renderPrint({
      showModal(modalDialog(title = "Network construction:","It can take a minute, press dismiss if you want to change the parameters", footer=modalButton("Dismiss")))
      on.exit(removeModal())
      
      
      summary3Core()
      
    }) 
    
    output$clustertwo =  renderPlot({distPlotCore()
    }) 
    output$powerplot2 =  renderPlot({distPlot1Core()
    }) 
    
    
    output$selected_graph_core <- renderPlot({
      showModal(modalDialog(title = "Network construction:","Almost there..", footer=NULL))
      on.exit(removeModal()) 
      distPlot2Core()
       
    })
    
    
    summary4Core <- reactive({
      dataset <- datasetInputCore()
      powerin<- input$powerin_core
      mingene<- input$mingene_core
      paste0(input$dataset_core,":","Power="," ",as.character(powerin),",","Number of genes="," ",as.character(mingene))
      
    })
    
    summary5Core <- reactive({
      dataset <- datasetInputCore()
      powerin<- input$powerin_core
      mingene<- input$mingene_core
      paste0(input$dataset_core,":","Power="," ",as.character(powerin),",","Number of genes="," ",as.character(mingene))
      
    })
    
    
    
    output$summary_datasetin_core4 <- renderPrint({
      summary4Core()
      
    })  
    
    output$summary_datasetin_core5 <- renderPrint({
      summary5Core()
      
    })
    
    
    # Show the first "n" observations ----
    # output$view <- renderTable({
    #     head(datasetInput(), n = input$obs)
    # })
    # 
    dat1Core <- reactive({
      multiExpr <- datasetInputCore()
      exprSize = checkSets(multiExpr)
      gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
      gsg$allOK
      
      if (!gsg$allOK)
      {
        # Print information about the removed genes:
        if (sum(!gsg$goodGenes) > 0)
          printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                                    collapse = ", ")))
        for (set in 1:exprSize$nSets)
        {
          if (sum(!gsg$goodSamples[[set]]))
            printFlush(paste("In set", setLabels[set], "removing samples",
                             paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
          # Remove the offending genes and samples
          multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
        }
        # Update exprSize
        exprSize = checkSets(multiExpr)
      }
      
      
      powerin_core <- input$powerin_core
      mingene_core <- input$mingene_core
      
      net = blockwiseConsensusModules(
        multiExpr, power = powerin_core, minModuleSize = mingene_core, deepSplit = 2,
        pamRespectsDendro = FALSE, 
        mergeCutHeight = 0.25, numericLabels = TRUE,
        minKMEtoStay = 0,
        saveTOMs = FALSE, verbose = 5)
      consMEs = net$multiMEs;
      moduleLabels = net$colors;
      moduleColors = labels2colors(moduleLabels)
      consTree = net$dendrograms[[1]]; 
      datME_1=moduleEigengenes(multiExpr[[1]]$data,moduleColors)$eigengenes
      datME_2=moduleEigengenes(multiExpr[[2]]$data,moduleColors)$eigengenes
        datME_1round=round(datME_1,2)
        datME_2round=round(datME_2,2)
        table_1_inspect=as.data.frame(table(labels2colors(net$colors)))
        colnames(table_1_inspect) <- NULL
        colnames(table_1_inspect) <- c("Module color","Number of genes")
        new_da_table = cbind(table_1_inspect,t(datME_1round),t(datME_2round))
        #colnames(new_da_table) <- c("Module color","Number of genes","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3",
        #                               "Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3",
        #                              "Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(new_da_table)
    })
    
    dat2Core <- reactive({
      multiExpr <- datasetInputCore()
      exprSize = checkSets(multiExpr)
      gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
      gsg$allOK
      
      if (!gsg$allOK)
      {
        # Print information about the removed genes:
        if (sum(!gsg$goodGenes) > 0)
          printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                                    collapse = ", ")))
        for (set in 1:exprSize$nSets)
        {
          if (sum(!gsg$goodSamples[[set]]))
            printFlush(paste("In set", setLabels[set], "removing samples",
                             paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
          # Remove the offending genes and samples
          multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
        }
        # Update exprSize
        exprSize = checkSets(multiExpr)
      }
      
      
      powerin_core <- input$powerin_core
      mingene_core <- input$mingene_core
      
      net = blockwiseConsensusModules(
        multiExpr, power = powerin_core, minModuleSize = mingene_core, deepSplit = 2,
        pamRespectsDendro = FALSE, 
        mergeCutHeight = 0.25, numericLabels = TRUE,
        minKMEtoStay = 0,
        saveTOMs = FALSE, verbose = 5)
      consMEs = net$multiMEs;
      moduleLabels = net$colors;
      moduleColors = labels2colors(moduleLabels)
      consTree = net$dendrograms[[1]]; 
      
        #table(net$colors)
        #Module_Colors = labels2colors(net$colors)
        PGFamID= names(net$colors)
      multiExpr_transpose_1=t(scale(multiExpr[[1]]$data))
      multiExpr_transpose_2=t(scale(multiExpr[[2]]$data))  
        #rownames(table_2_inspect) <- NULL
        
        #colnames(table_2_inspect) <- NULL
        #colnames(table_2_inspect) <- c("Module color","PGFam ID")
        probes = PGFamID
        PGFam_description=Annotations_core$PGFam_core[match(probes,Annotations_core$Genes_core)]
      #  PGFam=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
        table_2_inspect= as.data.frame(cbind(moduleColors,PGFamID,PGFam_description, multiExpr_transpose_1,multiExpr_transpose_2))
        return(as.data.frame(table_2_inspect))
    })
    
    color_from_middle <- function (data,color1,color2) 
    {
        max_val=max(abs(data))
        JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                   max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    output$mytable_core <- DT::renderDataTable({showModal(modalDialog(title = "Module expression:","Computing module eigengene...", footer=NULL))
      on.exit(removeModal()) 
      
        dat1Core() 
        
        return(datatable(dat1Core(),options = list(paging=FALSE, scrollX = FALSE,lengthChange = T, searching = TRUE),filter = "top",rownames=FALSE) %>% formatStyle(names(dat1Core()[-1,-2]), background= color_from_middle(dat1Core()[3:ncol(dat1())],'blue','red'), backgroundSize = '98% 78%',
                                                                                                                                                                backgroundRepeat = 'no-repeat', backgroundPosition = 'center'))
        
        
    })
    js <- c(
        "function(settings){",
        "  var instance = settings.oInstance;",
        "  var table = instance.api();",
        "  var input = instance.parent().find('.dataTables_filter input');",
        "  input.off('keyup search input').on('keyup', function(){",
        "    var keyword = '\\\\b' + input.val() + '\\\\b';",
        "    table.search(keyword, true, false).draw();",
        "  });",
        "}"
    )
    
    
    
    output$mytable2_core <- DT::renderDataTable({showModal(modalDialog(title = "Gene expression:","Computing z-score per gene...", footer=NULL))
      on.exit(removeModal()) 
        dat2Core()
        return(datatable(dat2Core(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(dat2Core()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
        #options = list(paging=FALSE, scrollX = FALSE), 
        #rownames=TRUE, 
        #filter = "top"
    })   
    
    
    
    
    
    # observe({
    #     updateTextInput(session, inputId = "myresultscolor", value = input$textforcolor) 
    #     
    # })   
    
    
    # distPlotalluvial <- reactive({
    #     Escherichia_coli_UPEC_536_GOclear <- read_csv("Escherichia_coli_UPEC_536_GOclear.txt",col_names = FALSE)
    # 
    #     Escherichia_coli_UPEC_536_GOclear_df= as.data.frame(Escherichia_coli_UPEC_536_GOclear[1:130,9:11])
    #     colnames(Escherichia_coli_UPEC_536_GOclear_df)= c("Parent", "Child1", "Child2")
    #     return(Escherichia_coli_UPEC_536_GOclear_df)
    #           # alluvial_wide( data = Escherichia_coli_UPEC_536_GOclear_df[1:1000,2:11]
    #             #                                      , fill_by = 'all_flows'
    #                 #                                  , stratum_label_size = 2.0
    #               #                                    , stratum_width = 1/20 ) 
    #               # parcats(palluvial, marginal_histograms = TRUE, data_input = Escherichia_coli_UPEC_536_GOclear_df[1:1000,2:4])
    # 
    #     
    # 
    # 
    #     
    #     
    # })
    
    
    # output$selected_graph2 <- render_parcats({
    #     nodes <- distPlotalluvial()
    #     #req(nodes)
    #     p = alluvial_wide(nodes, max_variables = ncol(nodes), fill_by = 'first_variable')
    #     parcats(p, marginal_histograms = TRUE, data_input = nodes)
    # })
    
    graph2InputCore <- reactive({
        switch(input$graph2_core,
               "Gene expression heatmap"= distPlot3Core(), 
               #"Gene expression box plots"= distPlot4(), 
               "Module eigengene barplot"= distPlot5Core()
               
               
        )
    })
    
    
    
    distPlot3Core <- reactive({
      multiExpr <- datasetInputCore()
      exprSize = checkSets(multiExpr)
      gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
      gsg$allOK
      
      if (!gsg$allOK)
      {
        # Print information about the removed genes:
        if (sum(!gsg$goodGenes) > 0)
          printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                                    collapse = ", ")))
        for (set in 1:exprSize$nSets)
        {
          if (sum(!gsg$goodSamples[[set]]))
            printFlush(paste("In set", setLabels[set], "removing samples",
                             paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
          # Remove the offending genes and samples
          multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
        }
        # Update exprSize
        exprSize = checkSets(multiExpr)
      }
      
      
      powerin_core <- input$powerin_core
      mingene_core <- input$mingene_core
      
      net = blockwiseConsensusModules(
        multiExpr, power = powerin_core, minModuleSize = mingene_core, deepSplit = 2,
        pamRespectsDendro = FALSE, 
        mergeCutHeight = 0.25, numericLabels = TRUE,
        minKMEtoStay = 0,
        saveTOMs = FALSE, verbose = 5)
      consMEs = net$multiMEs;
      moduleLabels = net$colors;
      moduleColors = labels2colors(moduleLabels)
      consTree = net$dendrograms[[1]]; 
      
        #MEs = net$MEs
        #datME=moduleEigengenes(dataset,moduleColors)$eigengenes
        colorh1=moduleColors
        which.module=input$textforcolor_core
       # par(mar=c(0.3, 15.5, 8, 2))
       # plotMat(t(scale(dataset[,colorh1==which.module]) ),nrgcols=30,rlabels=NULL,rcols="black",clabels=rownames(dataset),cex.lab=0.7)
        datME_1=moduleEigengenes(multiExpr[[1]]$data,moduleColors)$eigengenes
        ME=datME_1[, paste("ME",which.module, sep="")]
        col1 <- colorRampPalette(brewer.pal(7, "RdYlBu"))(256)
        col <- rev(col1)
        par(mar=c(1, 1, 1, 1))
        dataset1=scale(multiExpr[[1]]$data)
        dataset2=scale(multiExpr[[2]]$data)
        #dataset3= rbind(dataset1,dataset2)
        heatmap(t(rbind(dataset1[,colorh1==which.module], dataset2[,colorh1==which.module])), scale = "row", col =  col, Colv=NA, ColSideColors= rep(c("black", "grey") , each = nrow(multiExpr[[2]]$data)),margins=c(6,10),keep.dendro = TRUE, main= as.character(input$dataset_core) )
        legend(x="topright", legend=c("-3SD", "-2SD","-1SD","0","+1SD", "+2SD", "+3SD"),fill=c("#4575b4","#91bfdb","#e0f3f8","#ffffbf","#fee090","#fc8d59","#d73027"))
        #plotMat(t(scale(multiExpr[[1]]$data[,colorh1==which.module ]) ),nrgcols=30,rlabels=colnames(multiExpr[[1]]$data[,colorh1==which.module ]),rcols=which.module,clabels=rownames(multiExpr[[1]]$data), cex.main=0.7)
        #datME_2=moduleEigengenes(multiExpr[[2]]$data,moduleColors)$eigengenes
         #ME=datME_2[, paste("ME",which.module, sep="")]in
         #plotMat(t(scale(multiExpr[[2]]$data[,colorh1==which.module ]) ),nrgcols=30,rlabels=colnames(multiExpr[[2]]$data[,colorh1==which.module ]),rcols=which.module,clabels=rownames(multiExpr[[2]]$data), cex.main=0.7)
        legend(x="right", legend=c("Strain1", "Strain2"),fill=c("black","grey"))
    })
    
    
    distPlot5Core <- reactive({
      multiExpr <- datasetInputCore()
      exprSize = checkSets(multiExpr)
      gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
      gsg$allOK
      
      if (!gsg$allOK)
      {
        # Print information about the removed genes:
        if (sum(!gsg$goodGenes) > 0)
          printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                                    collapse = ", ")))
        for (set in 1:exprSize$nSets)
        {
          if (sum(!gsg$goodSamples[[set]]))
            printFlush(paste("In set", setLabels[set], "removing samples",
                             paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
          # Remove the offending genes and samples
          multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
        }
        # Update exprSize
        exprSize = checkSets(multiExpr)
      }
      
      
      powerin_core <- input$powerin_core
      mingene_core <- input$mingene_core
      
      net = blockwiseConsensusModules(
        multiExpr, power = powerin_core, minModuleSize = mingene_core, deepSplit = 2,
        pamRespectsDendro = FALSE, 
        mergeCutHeight = 0.25, numericLabels = TRUE,
        minKMEtoStay = 0,
        saveTOMs = FALSE, verbose = 5)
      consMEs = net$multiMEs;
      moduleLabels = net$colors;
      moduleColors = labels2colors(moduleLabels)
      consTree = net$dendrograms[[1]]; 
      
      #MEs = net$MEs
      #datME=moduleEigengenes(dataset,moduleColors)$eigengenes
      colorh1=moduleColors
      which.module=input$textforcolor_core
      # par(mar=c(0.3, 15.5, 8, 2))
      # plotMat(t(scale(dataset[,colorh1==which.module]) ),nrgcols=30,rlabels=NULL,rcols="black",clabels=rownames(dataset),cex.lab=0.7)
      datME_1=moduleEigengenes(multiExpr[[1]]$data,moduleColors)$eigengenes
      ME=datME_1[, paste("ME",which.module, sep="")]
      par(mfrow=c(1,2), mar=c(11,4,2,0.1))
      colors_bar= ifelse(ME<0, "blue", "red")  
    
      barplot(ME, col=t(colors_bar), main="", cex.lab=0.7,ylab="eigengene expression",xlab="",names.arg= rownames(multiExpr[[1]]$data), las=2, cex.names=0.7)   
      
    #  plotMat(t(scale(multiExpr[[1]]$data[,colorh1==which.module ]) ),nrgcols=30,rlabels=colnames(multiExpr[[1]]$data[,colorh1==which.module ]),rcols=which.module,clabels=rownames(multiExpr[[1]]$data), cex.main=0.7)
      datME_2=moduleEigengenes(multiExpr[[2]]$data,moduleColors)$eigengenes
      ME=datME_2[, paste("ME",which.module, sep="")]
  #    plotMat(t(scale(multiExpr[[2]]$data[,colorh1==which.module ]) ),nrgcols=30,rlabels=colnames(multiExpr[[2]]$data[,colorh1==which.module ]),rcols=which.module,clabels=rownames(multiExpr[[2]]$data), cex.main=0.7)
      
      colors_bar= ifelse(ME<0, "blue", "red")  
      #par(mar=c(5, 4.2, 0, 0.85))
      barplot(ME, col=t(colors_bar), main="", cex.lab=0.7,ylab="eigengene expression",xlab="",names.arg= rownames(multiExpr[[2]]$data), las=2, cex.names=0.7)   
    
        
    })
    
    
    
    
    output$myresultscolor_core <- renderPlot({
      
      showModal(modalDialog("Loading data...", footer=NULL))
      on.exit(removeModal())
        
        graph2InputCore()}, height = 500, width = 700 )
    
    output$consensus_graph_core <- renderPlot({
      multiExpr <- datasetInputCore()
      exprSize = checkSets(multiExpr)
      gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
      gsg$allOK
      
      if (!gsg$allOK)
      {
        # Print information about the removed genes:
        if (sum(!gsg$goodGenes) > 0)
          printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                                    collapse = ", ")))
        for (set in 1:exprSize$nSets)
        {
          if (sum(!gsg$goodSamples[[set]]))
            printFlush(paste("In set", setLabels[set], "removing samples",
                             paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
          # Remove the offending genes and samples
          multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
        }
        # Update exprSize
        exprSize = checkSets(multiExpr)
      }
      
      
      powerin_core <- input$powerin_core
      mingene_core <- input$mingene_core
      
      net = blockwiseConsensusModules(
        multiExpr, power = powerin_core, minModuleSize = mingene_core, deepSplit = 2,
        pamRespectsDendro = FALSE, 
        mergeCutHeight = 0.25, numericLabels = TRUE,
        minKMEtoStay = 0,
        saveTOMs = FALSE, verbose = 5)
      consMEs = net$multiMEs;
      moduleLabels = net$colors;
      moduleColors = labels2colors(moduleLabels)
      consTree = net$dendrograms[[1]]; 
      consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors);
      MET = consensusOrderMEs(consMEsC)
      par(cex = 0.6)
      #dataser_for_core= input$dataset_core
      setLabels=setLabels_analysis[[eval(input$dataset_core)]]
      plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                            zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
    })
    
    ########STIMULONS###########
    ############Acid_Stress_stimulon#############
    datasetInputStimulonAcid <- reactive({
      switch(input$select_acid_stress,
             "Burkholderia pseudomallei As (+)" = (Burkholderia_pseudomallei_expr_As),
             "Escherichia coli UPEC As (+) Sp (-)" = (Escherichia_coli_UPEC_536_expr_As),
             "Escherichia coli EPEC As (+) (I)"= (Escherichia_coli_EPEC_0127_H6_E2348_69_expr_As),
             "Escherichia coli EPEC As (+) (II)"=(Escherichia_coli_EPEC_0127_H6_E2348_69_expr_As2),
             "Enterococcus faecalis As (+)" =(Enterococcus_faecalis_expr_As),
             "Listeria monocytogenes As (+)"= (Listeria_monocytogenes_expr_As),
             "Pseudomonas aeruginosa As (+)" = (Pseudomonas_aeruginosa_expr_As),
             "Pseudomonas aeruginosa As (+) Tm (+)" = (Pseudomonas_aeruginosa_expr_As_Tm),
             "Escherichia coli ETEC As (+)"= (Escherichia_coli_ETEC_H10407_expr_As),
             "Salmonella enterica As (+)"= (Salmonella_expr_As),
             "Staphylococcus aureus MSSA476 As (+)"= (Staphylococcus_aureus_MSSA476_expr_As),
             "Staphylococcus aureus MSSA476 As (+) Ns (-)"= (Staphylococcus_aureus_MSSA476_expr_As_Ns),
             "Staphylococcus aureus MRSA252 As (+)"= (Staphylococcus_aureus_MRSA252_expr_As),
             "Klebsiella pneumoniae As (+) Sp (-)" = (Klebsiella_pneumoniae_expr_As_Sp),
             "Staphylococcus epidermidis As (+)" = (Staphylococcus_epidermidis_expr_As),
             "Shigella flexneri As (+)" = (Shigella_flexneri_expr_As),
             
             "Staphylococcus epidermidis As (+) Mig (+) "=(Staphylococcus_epidermidis_expr_Mig_up_As_up),
             "Listeria monocytogenes As (+) Oxs (+)" =(Listeria_monocytogenes_expr_Oxs_up_As_up),
             "Listeria monocytogenes As (+) Nd (-)"=(Listeria_monocytogenes_expr_As_up_Nd_down),
             "Staphylococcus epidermidis As (+) Tm (-)"=(Staphylococcus_epidermidis_expr_As_up_Tm_down)
             
      )
    })
    
    box_plot_AS <- reactive({
      AcidStress <- datasetInputStimulonAcid()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
        
        
        
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white") +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    
    output$Acidstress_boxplot <- renderPlot(
      
      box_plot_AS()
    )
    
    datAs <- reactive({
      dataset_As <- datasetInputStimulonAcid()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableAs <- DT::renderDataTable({datAs()
      
      return(datatable(datAs(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datAs()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })   
    ############SP_stimulon#############
    datasetInputStimulonSp <- reactive({
      switch(input$select_sp,
             "Neisseria gonorrhoeae Sp (-)" = (Neisseria_gonorrhoeae_expr_Sp),
             "Escherichia coli EPEC Sp (-)" = (Escherichia_coli_EPEC_0127_H6_E2348_69_expr_Sp),
             "Acinetobacter baumannii Sp (-) (I)" = (Acinetobacter_baumannii_expr_Sp1),
             "Acinetobacter baumannii Sp (-) (II)" = (Acinetobacter_baumannii_expr_Sp2),
             "Acinetobacter baumannii Sp (-) (III)" = (Acinetobacter_baumannii_expr_Sp3),
             "Acinetobacter baumannii Sp (-) (IV)"=(Acinetobacter_baumannii_expr_Sp4),
             "Streptococcus pyogenes Sp (-) Mig (-) (I)" = (Streptococcus_pyogenes_expr_Sp_Mg),
             "Burkholderia pseudomallei Sp (-)" = (Burkholderia_pseudomallei_expr_Sp_down),
             "Burkholderia pseudomallei Sp (+)"=(Burkholderia_pseudomallei_expr_Sp_up),
             "Enterococcus faecalis Sp (+)"=(Enterococcus_faecalis_expr_Sp_up),
             "Achromobacter xylosoxidans Sp (+)"=(Achromobacter_xylosoxidans_expr_Sp),
             "Borrelia burgdorferi Sp (+) (I)"=(Borrelia_burgdorferi_expr_Sp),
             "Borrelia burgdorferi Sp (+) (II)"=(Borrelia_burgdorferi_expr_Sp2),
             "Burkholderia pseudomallei Sp (+) Tm (+)"=(Burkholderia_pseudomallei_expr_Sp_Tm),
             "Campylobacter jejuni Sp (+)"=(Campylobacter_jejuni_expr_Sp),
             "Escherichia coli EPEC Sp (+)"=(Escherichia_coli_EPEC_0127_H6_E2348_69_expr_Sp),
             "Escherichia coli ETEC Sp (+)"=(Escherichia_coli_ETEC_expr_Sp),
             "Francisella tularensis Sp (+)"=(Francisella_tularensis_expr_Sp),
             "Haemophilus influenzae Sp (+) (I)"=(Haemophilus_influenzae_expr_Sp),
             "Haemophilus influenzae Sp (+) (II)"=(Haemophilus_influenzae_expr_Sp2),
             "Helicobacter pylori G27 Sp (+)"=(Helicobacter_pylori_G27_expr_Sp),
             "Klebsiella pneumoniae Sp (+)"= (Klebsiella_pneumoniae_expr_Sp),
             "Listeria monocytogenes Sp (+) Vic (+)" = (Listeria_monocytogenes_expr_Sp_Vic_up),
             "Staphylococcus aureus MRSA252 Sp (+)"=(Staphylococcus_aureus_MRSA252_expr_Sp_up),
             "Staphylococcus aureus MSSA476 Sp (+) Nd (+)"=(Staphylococcus_aureus_MSSA476_expr_Sp_up_Nd_up),
             "Staphylococcus aureus MSSA476 Sp (+)"=(Staphylococcus_aureus_MSSA476_expr_Sp_up),
             "Neisseria gonorrhoeae Sp (+) (I)"=(Neisseria_gonorrhoeae_expr_Sp1_up),
             "Neisseria gonorrhoeae Sp (+) (II)"=(Neisseria_gonorrhoeae_expr_Sp2_up),
             "Neisseria meningitidis Sp (+) Li (+)"= (Neisseria_meningitidis_expr_Sp_up_Li_up),
             "Neisseria meningitidis Sp (+)"=(Neisseria_meningitidis_Sp_up),
             "Pseudomonas aeruginosa Sp (+) Nd (+)" = (Pseudomonas_aeruginosa_expr_Sp_up_Nd_up),
             "Pseudomonas aeruginosa Sp (+)"=(Pseudomonas_aeruginosa_expr_Sp_up),
             "Salmonella enterica Sp (+) (I)"=(Salmonella_enterica_expr_Sp1_up),
             "Salmonella enterica Sp (+) (II)"=(Salmonella_enterica_expr_Sp2_up),
             "Staphylococcus epidermidis Sp (+) Tm (+)"=(Staphylococcus_epidermidis_expr_Sp_up_Tm_up),
             "Staphylococcus epidermidis Sp (+)"=(Staphylococcus_epidermidis_expr_Sp_up),
             "Shigella flexneri Sp (+)"=(Shigella_flexneri_expr_Sp_up),
             "Streptococcus pneumoniae Sp (+) (I)"=(Streptococcus_pneumoniae_expr_Sp1_up),
             "Streptococcus pneumoniae Sp (+) (II)"=(Streptococcus_pneumoniae_expr_Sp2_up),
             "Escherichia coli UPEC Sp (+)"=(Escherichia_coli_UPEC_expr_Sp_up),
             "Vibrio cholerae Sp (+) Li (+)"=(Vibrio_cholerae_expr_Sp_up_Li_up),
             "Vibrio cholerae Sp (+)" =(Vibrio_cholerae_expr_Sp_up),
             "Yersinia pseudotuberculosis Sp (+) (I)"=(Yersinia_pseudotuberculosis_expr_Sp1_up),
             "Yersinia pseudotuberculosis Sp (+) (II)"=(Yersinia_pseudotuberculosis_expr_Sp2_up),
             "Yersinia pseudotuberculosis Sp (+) Mig (-)"=(Yersinia_pseudotuberculosis_expr_Sp_up_Mig_down),
             "Aggregatibacter actinomycetemcomitans Sp (+)"=(Aggregatibacter_actinomycetemcomitans_expr_Sp_up),
             "Burkholderia pseudomallei Sp (+) Mig (+)"=(Burkholderia_pseudomallei_expr_Sp_up_Mig_up),
             "Enterococcus faecalis Sp (+) Mig (+)"=(Enterococcus_faecalis_expr_Sp_up_Mig_up),
             "Streptococcus pyogenes Sp (+) Mig (+)"=(Streptococcus_pyogenes_expr_Sp_up_Mig_up),
             "Streptococcus agalactiae Sp (+)"=(Streptococcus_agalactiae_expr_Sp_up),
             "Streptococcus agalactiae Sp (+) Mig (+)"=(Streptococcus_agalactiae_expr_Sp_up_Mig_up),
             
             
             
             "Escherichia coli UPEC Sp (-) As (+)"=(Escherichia_coli_UPEC_536_expr_As),
             "Klebsiella pneumoniae Sp (-) As (+)"=(Klebsiella_pneumoniae_expr_As_Sp),
             "Klebsiella pneumoniae Sp (-) Bs (+)"=(Klebsiella_pneumoniae_expr_Bs_Sp),
             "Pseudomonas aeruginosa Sp (-) Ns (+) (I)"=(Pseudomonas_aeruginosa_expr_Ns1_up_Sp_down),
             "Pseudomonas aeruginosa Sp (-) Ns (+) (II)"=(Pseudomonas_aeruginosa_expr_Ns2_up_Sp_down),
             "Klebsiella pneumoniae Sp (-)"=(Klebsiella_pneumoniae_expr_Sp_down),
             "Pseudomonas aeruginosa Sp (-) (I)"=(Pseudomonas_aeruginosa_expr_Sp_down),
             "Pseudomonas aeruginosa Sp (-) (II)"=(Pseudomonas_aeruginosa_expr_Sp2_down),
             "Streptococcus pyogenes Sp (-) Mig (-) (II)"=(Streptococcus_pyogenes_expr_Mig2_down_Sp_down),
             "Streptococcus pyogenes Sp (-) Mig (-) Li (+)"=(Streptococcus_pyogenes_expr_Mig_Sp_Li),
             "Escherichia coli UPEC Sp (-) Mig (-)"=(Escherichia_coli_UPEC_expr_Mig_Sp_down),
             "Shigella flexneri Sp (-) Mig (-)"=(Shigella_flexneri_expr_Mig_Sp_down)
             
             
      )
    })
    
    box_plot_SP <- reactive({
      AcidStress <- datasetInputStimulonSp()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white") +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white") +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    
    
    
    
    output$Sp_boxplot <- renderPlot(
      
      box_plot_SP()
    )
    
    datSp <- reactive({
      dataset_As <- datasetInputStimulonSp()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableSp <- DT::renderDataTable({datSp()
      
      return(datatable(datSp(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datSp()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })       
    
    
    ############Bs_stimulon#############
    datasetInputStimulonBs <- reactive({
      switch(input$select_bs,
             "Escherichia coli EPEC Bs (+)"=(Escherichia_coli_EPEC_0127_H6_E2348_69_expr_Bs),
             "Escherichia coli ETEC Bs (+)"=(Escherichia_coli_ETEC_H10407_expr_Bs),
             "Klebsiella pneumoniae Bs (+) Sp (-)"=(Klebsiella_pneumoniae_expr_Bs_Sp),
             "Klebsiella pneumoniae Bs (+)"=(Klebsiella_pneumoniae_expr_Bs),
             "Salmonella enterica Bs (+)"=(Salmonella_enterica_expr_Bs),
             "Escherichia coli UPEC Bs (+)"=(Escherichia_coli_UPEC_expr_Bs),
             "Francisella tularensis Bs (+)"=(Francisella_tularensis_expr_Bs),
             "Enterococcus faecalis Bs (+)"= (Enterococcus_faecalis_expr_Bs),
             "Listeria monocytogenes Bs (+) (I) "=(Listeria_expr_Bs1),
             "Listeria monocytogenes Bs (+) (II)"=(Listeria_expr_Bs2),
             "Staphylococcus aureus MRSA252 Bs (+)"=(Staphylococcus_aureus_MRSA252_expr_Bs),
             "Neisseria gonorrhoeae Bs (+) (I)"=(Neisseria_gonorrhoeae_expr_Bs1),
             "Neisseria gonorrhoeae Bs (+) (II)"=(Neisseria_gonorrhoeae_expr_Bs2),
             "Pseudomonas aeruginosa Bs (+)"=(Pseudomonas_aeruginosa_expr_Bs),
             "Staphylococcus epidermidis Bs (+)"=(Staphylococcus_epidermidis_expr_Bs),
             "Achromobacter xylosoxidans Bs (+)"=(Achromobacter_xylosoxidans_expr_Bs),
             "Aggregatibacter actinomycetemcomitans Bs (+) (I)"=(Aggregatibacter_actinomycetemcomitans_expr_Bs1),
             "Aggregatibacter actinomycetemcomitans Bs (+) (II)"=(Aggregatibacter_actinomycetemcomitans_expr_Bs2),
             "Aggregatibacter actinomycetemcomitans Bs (+) (III)"=(Aggregatibacter_actinomycetemcomitans_expr_Bs3),
             "Aggregatibacter actinomycetemcomitans Bs (+) (IV)"=(Aggregatibacter_actinomycetemcomitans_expr_Bs4),
             "Campylobacter jejuni Bs (+)"=(Campylobacter_jejuni_expr_Bs),
             "Streptococcus pneumoniae Bs (+) Mig (-)"=(Streptococcus_pneumoniae_expr_Bs_Mig),
             "Vibrio cholerae Bs (+)"=(Vibrio_cholerae_expr_Bs),
             "Neisseria meningitidis Bs (+)"=(Neisseria_meningitidis_expr_Bs)
             
             
      )
    })      
    
    box_plot_Bs <- reactive({
      AcidStress <- datasetInputStimulonBs()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white") +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    js <- c(
      "function(settings){",
      "  var instance = settings.oInstance;",
      "  var table = instance.api();",
      "  var input = instance.parent().find('.dataTables_filter input');",
      "  input.off('keyup search input').on('keyup', function(){",
      "    var keyword = '\\\\b' + input.val() + '\\\\b';",
      "    table.search(keyword, true, false).draw();",
      "  });",
      "}"
    )
    
    color_from_middle <- function (data,color1,color2) 
    {
      max_val=max(abs(data))
      JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                 max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    
    output$Bs_boxplot <- renderPlot(
      
      box_plot_Bs()
    )
    
    datBs <- reactive({
      dataset_As <- datasetInputStimulonBs()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableBs <- DT::renderDataTable({datBs()
      
      return(datatable(datBs(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datBs()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })       
    
    
    
    ############Nd_stimulon#############
    datasetInputStimulonNd <- reactive({
      switch(input$select_nd,
             "Enterococcus faecalis Nd (-)" =(Enterococcus_faecalis_expr_Nd),
             "Francisella tularensis Nd (-) (I)" = (Francisella_tularensis_expr_Nd),
             "Staphylococcus aureus MSSA476 Nd (-) Mig (-)" = (Staphylococcus_aureus_MSSA476_expr_Nd_Mig),
             "Staphylococcus aureus MRSA252 Nd (-) Mig (-)" = (Staphylococcus_aureus_MRSA252_expr_Nd_Mig),
             "Streptococcus agalactiae Nd (-) (I)" = (Streptococcus_agalactiae_expr_Nd),
             "Helicobacter pylori G27 Nd (-) (I)" = (Helicobacter_pylori_G27_expr_Nd),
             "Helicobacter pylori G27 Nd (-) (II)"=(Helicobacter_pylori_G27_expr_Nd2_down),
             "Helicobacter pylori J99 Nd (-)" = (Helicobacter_pylori_J99_expr_Nd),
             "Achromobacter xylosoxidans Nd (+)"=(Achromobacter_xylosoxidans_expr_Nd_up),
             "Acinetobacter baumannii Nd (+)"=(Acinetobacter_baumannii_expr_Nd_up),
             "Aggregatibacter actinomycetemcomitans Nd (+) (I)"=(Aggregatibacter_actinomycetemcomitans_expr_Nd_up),
             "Aggregatibacter actinomycetemcomitans Nd (+) (II)"=(Aggregatibacter_actinomycetemcomitans_expr_Nd2_up),
             "Borrelia burgdorferi Nd (+)"=(Borrelia_burgdorferi_expr_Nd_up),
             "Burkholderia pseudomallei Nd (+) (I)"=(Burkholderia_pseudomallei_expr_Nd1_up),
             "Burkholderia pseudomallei Nd (+) (II)"=(Burkholderia_pseudomallei_expr_Nd2_up),
             "Burkholderia pseudomallei Nd (+) (III)"=(Burkholderia_pseudomallei_expr_Nd3_up),
             "Campylobacter jejuni Nd (+)"=(Campylobacter_jejuni_expr_Nd_up),
             "Enterococcus faecalis Nd (+) (I)"=(Enterococcus_faecalis_expr_Nd1_up),
             "Enterococcus faecalis Nd (+) (II)"=(Enterococcus_faecalis_expr_Nd2_up),
             "Escherichia coli EPEC Nd (+) Oss (+)"=(Escherichia_coli_EPEC_expr_Nd_up_Oss_up),
             "Escherichia coli EPEC Nd (+)"=(Escherichia_coli_EPEC_expr_Nd_up),
             "Escherichia coli ETEC Nd (+) (I)"=(Escherichia_coli_ETEC_expr_Nd1_up),
             "Escherichia coli ETEC Nd (+) (II)"=(Escherichia_coli_ETEC_expr_Nd2_up),
             "Escherichia coli ETEC Nd (+) Mig (+)"=(Escherichia_coli_ETEC_expr_Nd_up_Mig_up),
             "Francisella tularensis Nd (+) (I)"=(Francisella_tularensis_expr_Nd1_up),
             "Francisella tularensis Nd (+) (II)"=(Francisella_tularensis_expr_Nd2_up),
             "Haemophilus influenzae Nd (+)"=(Haemophilus_influenzae_expr_Nd_up),
             "Helicobacter pylori G27 Nd (+)"=(Helicobacter_pylori_G27_expr_Nd_up),
             "Helicobacter pylori J99 Nd (+) (I)"=(Helicobacter_pylori_J99_expr_Nd_up),
             "Klebsiella pneumoniae Nd (+) (I)"=(Klebsiella_pneumoniae_expr_Nd_up),
             "Klebsiella pneumoniae Nd (+) (II)"=(Klebsiella_pneumoniae_expr_Nd2_up),
             "Listeria monocytogenes Nd (+) (I)"=(Listeria_monocytogenes_expr_Nd1_up),
             "Listeria monocytogenes Nd (+) (II)"=(Listeria_monocytogenes_expr_Nd2_up),
             "Listeria monocytogenes Nd (+) (III)"=(Listeria_monocytogenes_expr_Nd3_up),
             "Staphylococcus aureus MRSA252 Nd (+) Ns (-)"=(Staphylococcus_aureus_MRSA252_expr_Nd_up_Ns_down),
             "Staphylococcus aureus MRSA252 Nd (+)"=(Staphylococcus_aureus_MRSA252_expr_Nd_up),
             "Staphylococcus aureus MSSA476 Nd (+)"=(Staphylococcus_aureus_MSSA476_expr_Nd_up),
             "Neisseria meningitidis Nd (+)"=(Neisseria_meningitidis_expr_Nd_up),
             "Pseudomonas aeruginosa Nd (+) Ns (+)"=(Pseudomonas_aeruginosa_expr_Nd_up_Ns_up),
             "Pseudomonas aeruginosa Nd (+) Oss (+)"=(Pseudomonas_aeruginosa_expr_Nd_up_Oss_up),
             "Pseudomonas aeruginosa Nd (+)"=(Pseudomonas_aeruginosa_expr_Nd_up),
             "Salmonella enterica Nd (+)"=(Salmonella_enterica_expr_Nd_up),
             "Staphylococcus epidermidis Nd (+) (I)"=(Staphylococcus_epidermidis_expr_Nd1_up),
             "Staphylococcus epidermidis Nd (+) (II)"=(Staphylococcus_epidermidis_expr_Nd2_up),
             "Staphylococcus epidermidis Nd (+) Tm (-)"=(Staphylococcus_epidermidis_expr_Nd_up_Tm_down),
             "Shigella flexneri Nd (+)"=(Shigella_flexneri_expr_Nd_up),
             "Streptococcus pyogenes Nd (+)"=(Streptococcus_pyogenes_expr_Nd_up),
             "Streptococcus agalactiae Nd (+)"=(Streptococcus_agalactiae_expr_Nd_up),
             "Escherichia coli UPEC Nd (+) (I)"=(Escherichia_coli_UPEC_expr_Nd1_up),
             "Escherichia coli UPEC Nd (+) (II)"=(Escherichia_coli_UPEC_expr_Nd2_up),
             "Vibrio cholerae Nd (+)"=(Vibrio_cholerae_expr_Nd_up),
             "Yersinia pseudotuberculosis Nd (+)"=(Yersinia_pseudotuberculosis_expr_Nd_up),
             "Yersinia pseudotuberculosis Nd (+) Vic (-)"=(Yersinia_pseudotuberculosis_expr_Nd_up_Vic_down),
             "Enterococcus faecalis Nd (+) Oxs (+) Vic(-)"=(Enterococcus_faecalis_expr_Nd_up_Oxs_up_Vic_down),
             "Klebsiella pneumoniae Nd (+) Oxs (+)"=(Klebsiella_pneumoniae_expr_Nd_up_Oxs_up),
             "Salmonella enterica Nd (+) Oxs (+)"=(Salmonella_enterica_expr_Nd_up_Oxs_up),
             "Shigella flexneri Nd (+) Oxs (+)"=(Shigella_flexneri_expr_Nd_up_Oxs_up),
             "Staphylococcus aureus MSSA476 Nd (+) Sp (+)"=(Staphylococcus_aureus_MSSA476_expr_Sp_up_Nd_up),
             "Pseudomonas aeruginosa Nd (+) Sp (+)"=(Pseudomonas_aeruginosa_expr_Sp_up_Nd_up),
             "Listeria monocytogenes  Nd (-) Li (+)"=(Listeria_monocytogenes_expr_Li_up_Nd_down),
             "Vibrio cholerae Nd (+) Mig (+) "=(Vibrio_cholerae_expr_Mig_up_Nd_up),
             "Aggregatibacter actinomycetemcomitans  Nd (-) Mig (+)"=(Aggregatibacter_actinomycetemcomitans_expr_Mig_up_Nd_down),
             "Staphylococcus aureus MRSA252 Nd (+) Mig (+) "=(Staphylococcus_aureus_MRSA252_expr_Mig_up_Nd_up),
             "Streptococcus agalactiae Nd (-) Mig (+) "=(Streptococcus_agalactiae_expr_Mig_up_Nd_down),
             "Listeria monocytogenes Nd (-) Oss (+) "=(Listeria_monocytogenes_expr_Oss_up_Nd_down),
             "Streptococcus agalactiae Nd (-) Oss (+) "=(Streptococcus_agalactiae_expr_Oss_up_Nd_down),
             "Streptococcus pyogenes Nd (+) Tm (+) "=(Streptococcus_pyogenes_expr_Tm_up_Nd_up),
             "Streptococcus agalactiae Nd (-) Tm (+) "=(Streptococcus_agalactiae_expr_Tm_up_Nd_down),
             "Streptococcus agalactiae Nd (-) Vic (+) "=(Streptococcus_agalactiae_expr_Vic_up_Nd_down),
             "Streptococcus agalactiae Nd (-) (II)"=(Streptococcus_agalactiae_expr_Nd2_down),
             "Listeria monocytogenes Nd (-) (I)"=(Listeria_monocytogenes_expr_Nd1_down),
             
             "Listeria monocytogenes Nd (-) (II)"=(Listeria_monocytogenes_expr_Nd2_down),
             "Francisella tularensis Nd (-) (II)"=(Francisella_tularensis_expr_Nd2_down),
             "Listeria monocytogenes Nd (-) As (+)"=(Listeria_monocytogenes_expr_As_up_Nd_down),
             "Listeria monocytogenes Nd (-) (III)"=(Listeria_monocytogenes_expr_Nd3_down),
             "Escherichia coli UPEC Nd (+) (III)"=(Escherichia_coli_UPEC_expr_Nd3_up),
             "Helicobacter pylori J99 Nd (+) (II)"=(Helicobacter_pylori_J99_expr_Nd2_up)
             
             
             
      )
    })      
    
    box_plot_Nd <- reactive({
      AcidStress <- datasetInputStimulonNd()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1,color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    js <- c(
      "function(settings){",
      "  var instance = settings.oInstance;",
      "  var table = instance.api();",
      "  var input = instance.parent().find('.dataTables_filter input');",
      "  input.off('keyup search input').on('keyup', function(){",
      "    var keyword = '\\\\b' + input.val() + '\\\\b';",
      "    table.search(keyword, true, false).draw();",
      "  });",
      "}"
    )
    
    color_from_middle <- function (data,color1,color2) 
    {
      max_val=max(abs(data))
      JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                 max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    
    output$Nd_boxplot <- renderPlot(
      
      box_plot_Nd()
    )
    
    datNd <- reactive({
      dataset_As <- datasetInputStimulonNd()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableNd <- DT::renderDataTable({datNd()
      
      return(datatable(datNd(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datNd()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })       
    
    ############Li_stimulon#############
    datasetInputStimulonLi <- reactive({
      switch(input$select_li,
             "Achromobacter xylosoxidans Li (+) (I)"=(Achromobacter_xylosoxidans_expr_Li1_up),
             "Achromobacter xylosoxidans Li (+) (II)"=(Achromobacter_xylosoxidans_expr_Li2_up),
             "Acinetobacter baumannii Li (+) Oxs (+)"=(Acinetobacter_baumannii_expr_Li_up_Oxs_up),
             "Aggregatibacter actinomycetemcomitans Li (+) (I)"=(Aggregatibacter_actinomycetemcomitans_expr_Li1_up),
             "Aggregatibacter actinomycetemcomitans Li (+) (II)"=(Aggregatibacter_actinomycetemcomitans_expr_Li2_up),
             "Burkholderia pseudomallei Li (+)"=(Burkholderia_pseudomallei_expr_Li_up),
             "Campylobacter jejuni Li (+)"=(Campylobacter_jejuni_expr_Li_up),
             "Enterococcus faecalis Li (+)"=(Enterococcus_faecalis_expr_Li_up),
             "Escherichia coli EPEC Li (+)"=(Escherichia_coli_EPEC_expr_Li_up),
             "Escherichia coli ETEC Li (+)"=(Escherichia_coli_ETEC_expr_Li_up),
             "Francisella tularensis Li (+)"=(Francisella_tularensis_expr_Li_up),
             "Helicobacter pylori G27 Li (+)"=(Helicobacter_pylori_G27_expr_Li_up),
             "Helicobacter pylori J99 Li (+)"=(Helicobacter_pylori_J99_expr_Li_up),
             "Klebsiella pneumoniae Li (+)"=(Klebsiella_pneumoniae_expr_Li_up),
             "Listeria monocytogenes Li (+) Nd (-)"=(Listeria_monocytogenes_expr_Li_up_Nd_down),
             "Staphylococcus aureus MRSA252 Li (+)"=(Staphylococcus_aureus_MRSA252_expr_Li_up),
             "Neisseria gonorrhoeae Li (+)"=(Neisseria_gonorrhoeae_expr_Li_up),
             "Pseudomonas aeruginosa Li (+) (I)"=(Pseudomonas_aeruginosa_expr_Li1_up),
             "Pseudomonas aeruginosa Li (+) (II)"=(Pseudomonas_aeruginosa_expr_Li2_up),
             "Salmonella enterica Li (+) Oxs (+)"=(Salmonella_enterica_expr_Li_up_Oxs_up),
             "Salmonella enterica Li (+)"=(Salmonella_enterica_expr_Li_up),
             "Staphylococcus epidermidis Li (+)"=(Staphylococcus_epidermidis_expr_Li_up),
             "Shigella flexneri Li (+)"=(Shigella_flexneri_expr_Li_up),
             "Shigella flexneri Li (+) Oxs (+)"=(Shigella_flexneri_expr_Li_up_Oxs_up),
             "Streptococcus agalactiae Li (+) (I)"=(Streptococcus_agalactiae_expr_Li1_up),
             "Streptococcus agalactiae Li (+) (II)"=(Streptococcus_agalactiae_expr_Li2_up),
             "Streptococcus pneumoniae Li (+)"=(Streptococcus_pneumoniae_expr_Li_up),
             "Escherichia coli UPEC Li (+)"=(Escherichia_coli_UPEC_expr_Li_up),
             "Vibrio cholerae Li (+) (I)"=(Vibrio_cholerae_expr_Li1_up),
             "Vibrio cholerae Li (+) (II)"=(Vibrio_cholerae_expr_Li2_up),
             "Yersinia pseudotuberculosis Li (+)"=(Yersinia_pseudotuberculosis_expr_Li_up),
             "Neisseria meningitidis Li (+) Sp (+)"=(Neisseria_meningitidis_expr_Sp_up_Li_up),
             "Vibrio cholerae Li (+) Sp (+)"=(Vibrio_cholerae_expr_Sp_up_Li_up),
             "Escherichia coli ETEC Li (+) Tm (+) "=(Escherichia_coli_ETEC_expr_Tm_up_Li_up),
             "Streptococcus pneumoniae Li (+) Mig (-)"=(Streptococcus_pneumoniae_Li_up_Mig_down),
             "Streptococcus pyogenes Li (+) Mig (-) Sp (-) "=(Streptococcus_pyogenes_expr_Mig_Sp_Li)
             
             
             
             
             
      )
    })      
    
    box_plot_Li <- reactive({
      AcidStress <- datasetInputStimulonLi()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    js <- c(
      "function(settings){",
      "  var instance = settings.oInstance;",
      "  var table = instance.api();",
      "  var input = instance.parent().find('.dataTables_filter input');",
      "  input.off('keyup search input').on('keyup', function(){",
      "    var keyword = '\\\\b' + input.val() + '\\\\b';",
      "    table.search(keyword, true, false).draw();",
      "  });",
      "}"
    )
    
    color_from_middle <- function (data,color1,color2) 
    {
      max_val=max(abs(data))
      JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                 max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    
    output$Li_boxplot <- renderPlot(
      
      box_plot_Li()
    )
    
    datLi <- reactive({
      dataset_As <- datasetInputStimulonLi()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableLi <- DT::renderDataTable({datLi()
      
      return(datatable(datLi(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datLi()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })       
    ############Mig_stimulon#############
    datasetInputStimulonMig <- reactive({
      switch(input$select_mig,
             "Burkholderia pseudomallei Mig (+) Sp (+)"=(Burkholderia_pseudomallei_expr_Mig_up_Sp_up),
             "Enterococcus faecalis Mig (+) Vic (+)"=(Enterococcus_faecalis_expr_Mig_up_Vic_up),
             "Streptococcus pyogenes Mig (+) Sp (+)"=(Streptococcus_pyogenes_expr_Mig_up_Sp_up),
             "Streptococcus agalactiae Mig (+) Sp (+)"=(Streptococcus_agalactiae_expr_Sp_up_Mig_up),
             
             
             "Vibrio cholerae Mig (+) Nd (+)"=(Vibrio_cholerae_expr_Mig_up_Nd_up),
             "Achromobacter xylosoxidans Mig (+) (I)"=(Achromobacter_xylosoxidans_expr_Mig1_up),
             "Achromobacter xylosoxidans Mig (+) (II)"=(Achromobacter_xylosoxidans_expr_Mig2_up),
             "Aggregatibacter actinomycetemcomitans Mig (+) Nd (-)"=(Aggregatibacter_actinomycetemcomitans_expr_Mig_up_Nd_down),
             "Burkholderia pseudomallei Mig (+) (I)"=(Burkholderia_pseudomallei_expr_Mig1_up),
             "Burkholderia pseudomallei Mig (+) (II)"=(Burkholderia_pseudomallei_expr_Mig2_up),
             "Enterococcus faecalis Mig (+)"=(Enterococcus_faecalis_expr_Mig_up),
             "Escherichia coli EPEC Mig (+)"=(Escherichia_coli_EPEC_expr_Mig_up),
             "Escherichia coli ETEC Mig (+) Oss (+)"=(Escherichia_coli_ETEC_expr_Mig_up_Oss_up),
             "Escherichia coli ETEC Mig (+)"=(Escherichia_coli_ETEC_expr_Mig_up),
             "Haemophilus influenzae Mig (+)"=(Haemophilus_influenzae_expr_Mig_up),
             "Haemophilus influenzae Mig (+) Vic (+)"=(Haemophilus_influenzae_expr_Mig_up_Vic_up),
             "Klebsiella pneumoniae Mig (+) "=(Klebsiella_pneumoniae_expr_Mig_up),
             "Listeria monocytogenes Mig (+)"=(Listeria_monocytogenes_expr_Mig_up),
             "Staphylococcus aureus MRSA252 Mig (+) Nd (+)"=(Staphylococcus_aureus_MRSA252_expr_Mig_up_Nd_up),
             "Staphylococcus aureus MRSA252 Mig (+) (I)"=(Staphylococcus_aureus_MRSA252_expr_Mig1_up),
             "Staphylococcus aureus MRSA252 Mig (+) (II)"=(Staphylococcus_aureus_MRSA252_expr_Mig2_up),
             "Staphylococcus aureus MSSA476 Mig (+)"=(Staphylococcus_aureus_MSSA476_expr_Mig_up),
             "Neisseria gonorrhoeae Mig (+) (I)"=(Neisseria_gonorrhoeae_expr_Mig1_up),
             "Neisseria gonorrhoeae Mig (+) (II)"=(Neisseria_gonorrhoeae_expr_Mig2_up),
             "Neisseria meningitidis Mig (+)"=(Neisseria_meningitidis_expr_Mig_up),
             "Pseudomonas aeruginosa Mig (+) Oxs (+)"=(Pseudomonas_aeruginosa_expr_Mig_up_Oxs_up),
             "Pseudomonas aeruginosa Mig (+)"=(Pseudomonas_aeruginosa_expr_Mig_up),
             "Pseudomonas aeruginosa Mig (+) Ns (-)"=(Pseudomonas_aeruginosa_expr_Mig_up_Ns_down),
             "Salmonella enterica Mig (+)"=(Salmonella_enterica_expr_Mig_up),
             "Staphylococcus epidermidis Mig (+) (I)"=(Staphylococcus_epidermidis_expr_Mig1_up),
             "Staphylococcus epidermidis Mig (+) (II)"=(Staphylococcus_epidermidis_expr_Mig2_up),
             "Staphylococcus epidermidis Mig (+) Tm (-)"=(Staphylococcus_epidermidis_expr_Mig_up_Tm_down),
             "Staphylococcus epidermidis Mig (+) As (+)"=(Staphylococcus_epidermidis_expr_Mig_up_As_up),
             "Shigella flexneri Mig (+)"=(Shigella_flexneri_expr_Mig_up),
             "Streptococcus agalactiae Mig (+) Nd (-)"=(Streptococcus_agalactiae_expr_Mig_up_Nd_down),
             "Streptococcus pneumoniae Mig (+)"=(Streptococcus_pneumoniae_expr_Mig_up),
             "Escherichia coli UPEC Mig (+)"=(Escherichia_coli_UPEC_expr_Mig_up),
             "Escherichia coli UPEC Mig (+) Vic (+)"=(Escherichia_coli_UPEC_expr_Mig_up_Vic_up),
             "Vibrio cholerae Mig (+) (I)"=(Vibrio_cholerae_expr_Mig1_up),
             "Vibrio cholerae Mig (+) (II)"=(Vibrio_cholerae_expr_Mig2_up),
             "Yersinia pseudotuberculosis Mig (+)"=(Yersinia_pseudotuberculosis_expr_Mig_up),
             "Streptococcus pyogenes Mig (-) Sp (-) (I)"=(Streptococcus_pyogenes_expr_Sp_Mg),
             "Yersinia pseudotuberculosis Mig (-) Sp (+)"=(Yersinia_pseudotuberculosis_expr_Sp_up_Mig_down),
             
             "Enterococcus faecalis Mig (+) Sp (+)"=(Enterococcus_faecalis_expr_Sp_up_Mig_up),
             
             
             "Staphylococcus aureus MSSA476 Mig (-) Nd (-)"=(Staphylococcus_aureus_MSSA476_expr_Nd_Mig),
             "Staphylococcus aureus MRSA252 Mig (-) Nd (-)"=(Staphylococcus_aureus_MRSA252_expr_Nd_Mig),
             "Escherichia coli ETEC Mig (+) Nd (+)"=(Escherichia_coli_ETEC_expr_Nd_up_Mig_up),
             "Streptococcus pneumoniae Mig (-) Bs (+)"=(Streptococcus_pneumoniae_expr_Bs_Mig),
             "Haemophilus influenzae Mig (-) Vic (+)"=(Haemophilus_influenzae_expr_Vic_up_Mig_down),
             "Streptococcus pneumoniae Mig (-) Vic (+) "=(Streptococcus_pneumoniae_expr_Vic_up_Mig_down),
             "Yersinia pseudotuberculosis Mig (-) Vic (+) "=(Yersinia_pseudotuberculosis_expr_Vic_up_Mig_down),
             "Streptococcus pneumoniae Mig (-) Li (+)"=(Streptococcus_pneumoniae_expr_Li_up_Mig_down),
             "Streptococcus pyogenes Mig (-) Sp (-) (II)"=(Streptococcus_pyogenes_expr_Mig2_down_Sp_down),
             "Haemophilus influenzae Mig (-) (I)"=(Haemophilus_influenzae_expr_Mig_down),
             "Haemophilus influenzae Mig (-) (II)"=(Haemophilus_influenzae_expr_Mig2_down),
             "Streptococcus pyogenes Mig (-) Sp (-) Li (+)"=(Streptococcus_pyogenes_expr_Mig_Sp_Li),
             "Escherichia coli UPEC Mig (-) Sp (-)"=(Escherichia_coli_UPEC_expr_Mig_Sp_down),
             "Shigella flexneri Mig (-) Sp (-)"=(Shigella_flexneri_expr_Mig_Sp_down),
             "Listeria monocytogenes Mig (-)"=(Listeria_monocytogenes_expr_Mig_down)
             
             
             
             
             
             
      )
    })      
    
    box_plot_Mig <- reactive({
      AcidStress <- datasetInputStimulonMig()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    js <- c(
      "function(settings){",
      "  var instance = settings.oInstance;",
      "  var table = instance.api();",
      "  var input = instance.parent().find('.dataTables_filter input');",
      "  input.off('keyup search input').on('keyup', function(){",
      "    var keyword = '\\\\b' + input.val() + '\\\\b';",
      "    table.search(keyword, true, false).draw();",
      "  });",
      "}"
    )
    
    color_from_middle <- function (data,color1,color2) 
    {
      max_val=max(abs(data))
      JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                 max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    
    output$Mig_boxplot <- renderPlot(
      
      box_plot_Mig()
    )
    
    datMig <- reactive({
      dataset_As <- datasetInputStimulonMig()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableMig <- DT::renderDataTable({datMig()
      
      return(datatable(datMig(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datMig()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })
    ############Ns_stimulon#############
    datasetInputStimulonNs <- reactive({
      switch(input$select_ns,
             
             "Achromobacter xylosoxidans Ns (+)"=(Achromobacter_xylosoxidans_expr_Ns_up),
             "Aggregatibacter actinomycetemcomitans Ns (+)"=(Aggregatibacter_actinomycetemcomitans_expr_Ns_up),
             "Enterococcus faecalis Ns (+)"=(Enterococcus_faecalis_expr_Ns_up),
             "Escherichia coli EPEC Ns (+)"=(Escherichia_coli_EPEC_expr_Ns_up),
             "Escherichia coli ETEC Ns (+)"=(Escherichia_coli_ETEC_expr_Ns_up),
             "Haemophilus influenzae Ns (+)"=(Haemophilus_influenzae_expr_Ns_up),
             "Klebsiella pneumoniae Ns (+)"=(Klebsiella_pneumoniae_expr_Ns_up),
             "Listeria monocytogenes Ns (+) (I)"=(Listeria_monocytogenes_expr_Ns_up),
             "Pseudomonas aeruginosa Ns (+)"=(Pseudomonas_aeruginosa_expr_Ns_up),
             "Pseudomonas aeruginosa Ns (+) Sp (-) (I)"=(Pseudomonas_aeruginosa_expr_Ns1_up_Sp_down),
             "Pseudomonas aeruginosa Ns (+) Sp (-) (II)"=(Pseudomonas_aeruginosa_expr_Ns2_up_Sp_down),
             "Salmonella enterica Ns (+)"=(Salmonella_enterica_expr_Ns_up),
             "Staphylococcus epidermidis Ns (+)"=(Staphylococcus_epidermidis_expr_Ns_up),
             "Staphylococcus epidermidis Ns (+) Tm (-)"=(Staphylococcus_epidermidis_expr_Ns_up_Tm_down),
             "Shigella flexneri Ns (+)"=(Shigella_flexneri_expr_Ns_up),
             "Streptococcus agalactiae Ns (+)"=(Streptococcus_pneumoniae_expr_Ns_up),
             "Streptococcus pneumoniae Ns (+)"=(Streptococcus_pneumoniae_expr_Ns_up),
             "Pseudomonas aeruginosa Ns (+) Nd (+)"=(Pseudomonas_aeruginosa_expr_Nd_up_Ns_up),
             "Staphylococcus aureus MSSA476 Ns (-) As (+)"=(Staphylococcus_aureus_MSSA476_expr_As_Ns),
             "Staphylococcus aureus MRSA252 Ns (-) Nd (+)"=(Staphylococcus_aureus_MRSA252_expr_Nd_up_Ns_down),
             "Pseudomonas aeruginosa Ns (-) Mig (+) "=(Pseudomonas_aeruginosa_expr_Mig_up_Ns_down),
             "Streptococcus pyogenes Ns (+) Oss (+) "=(Streptococcus_pyogenes_expr_Oss_up_Ns_up),
             "Listeria monocytogenes Ns (+) (II)"=(Listeria_monocytogenes_expr_Ns2_up),
             "Staphylococcus aureus MSSA476 Ns (-)"=(Staphylococcus_aureus_MSSA476_expr_Ns_down),
             "Pseudomonas aeruginosa Ns (-)"=(Pseudomonas_aeruginosa_expr_Ns_down)
             
             
             
             
             
             
             
      )
    })      
    
    box_plot_Ns <- reactive({
      AcidStress <- datasetInputStimulonNs()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white") +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white") +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    js <- c(
      "function(settings){",
      "  var instance = settings.oInstance;",
      "  var table = instance.api();",
      "  var input = instance.parent().find('.dataTables_filter input');",
      "  input.off('keyup search input').on('keyup', function(){",
      "    var keyword = '\\\\b' + input.val() + '\\\\b';",
      "    table.search(keyword, true, false).draw();",
      "  });",
      "}"
    )
    
    color_from_middle <- function (data,color1,color2) 
    {
      max_val=max(abs(data))
      JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                 max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    
    output$Ns_boxplot <- renderPlot(
      
      box_plot_Ns()
    )
    
    datNs <- reactive({
      dataset_As <- datasetInputStimulonNs()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableNs <- DT::renderDataTable({datNs()
      
      return(datatable(datNs(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datNs()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })
    ############Oss_stimulon#############
    datasetInputStimulonOss <- reactive({
      switch(input$select_oss,
             
             "Achromobacter xylosoxidans Oss (+)"=(Achromobacter_xylosoxidans_expr_Oss_up),
             "Acinetobacter baumannii Oss (+)"=(Acinetobacter_baumannii_expr_Oss_up),
             "Aggregatibacter actinomycetemcomitans Oss (+)"=(Aggregatibacter_actinomycetemcomitans_expr_Oss_up),
             "Borrelia burgdorferi Oss (+)"=(Borrelia_burgdorferi_expr_Oss_up),
             "Burkholderia pseudomallei Oss (+)"=(Burkholderia_pseudomallei_expr_Oss_up),
             "Enterococcus faecalis Oss (+)"=(Enterococcus_faecalis_expr_Oss_up),
             "Escherichia coli EPEC Oss (+)"=(Escherichia_coli_EPEC_expr_Oss_up),
             "Escherichia coli ETEC Oss (+)"=(Escherichia_coli_ETEC_expr_Oss_up),
             "Francisella tularensis Oss (+)"=(Francisella_tularensis_expr_Oss_up),
             "Haemophilus influenzae Oss (+)"=(Haemophilus_influenzae_expr_Oss_up),
             "Helicobacter pylori G27 Oss (+)"=(Helicobacter_pylori_G27_expr_Oss_up),
             "Helicobacter pylori J99 Oss (+)"=(Helicobacter_pylori_J99_expr_Oss_up),
             "Klebsiella pneumoniae Oss (+) (I)"=(Klebsiella_pneumoniae_expr_Oss1_up),
             "Klebsiella pneumoniae Oss (+) (II)"=(Klebsiella_pneumoniae_expr_Oss2_up),
             "Klebsiella pneumoniae Oss (+) Tm (+)"=(Klebsiella_pneumoniae_expr_Oss_up_Tm_up),
             "Listeria monocytogenes Oss (+) Nd (-)"=(Listeria_monocytogenes_expr_Oss_up_Nd_down),
             "Staphylococcus aureus MRSA252 Oss (+)"=(Staphylococcus_aureus_MRSA252_expr_Oss_up),
             "Staphylococcus aureus MSSA476 Oss (+) (I)"=(Staphylococcus_aureus_MSSA476_expr_Oss1_up),
             "Staphylococcus aureus MSSA476 Oss (+) (II)"=(Staphylococcus_aureus_MSSA476_expr_Oss2_up),
             "Neisseria meningitidis Oss (+)"=(Neisseria_meningitidis_expr_Oss_up),
             "Pseudomonas aeruginosa Oss (+)"=(Pseudomonas_aeruginosa_expr_Oss_up),
             "Salmonella enterica Oss (+)"=(Salmonella_enterica_expr_Oss_up),
             "Shigella flexneri Oss (+) (I)"=(Shigella_flexneri_expr_Oss_up),
             "Shigella flexneri Oss (+) (II)"=(Shigella_flexneri_expr_Oss2_up),
             "Shigella flexneri Oss (+) Tm (+)"=(Shigella_flexneri_expr_Oss_up_Tm_up),
             "Streptococcus pyogenes Oss (+)"=(Streptococcus_pyogenes_expr_Oss_up),
             "Streptococcus pyogenes Oss (+) Ns (+)"=(Streptococcus_pyogenes_expr_Oss_up_Ns_up),
             "Streptococcus agalactiae Oss (+)"=(Streptococcus_agalactiae_expr_Oss_up),
             "Streptococcus agalactiae Oss (+) Nd (-)"=(Streptococcus_agalactiae_expr_Oss_up_Nd_down),
             "Streptococcus agalactiae Oss (+) Vic (-)"=(Streptococcus_agalactiae_expr_Oss_up_Vic_down),
             "Escherichia coli UPEC Oss (+)"=(Escherichia_coli_UPEC_expr_Oss_up),
             "Vibrio cholerae Oss (+)"=(Vibrio_cholerae_expr_Oss_up),
             "Yersinia pseudotuberculosis Oss (+)"=(Yersinia_pseudotuberculosis_expr_Oss_up),
             "Escherichia coli EPEC Oss (+) Nd (+)"=(Escherichia_coli_EPEC_expr_Nd_up_Oss_up),
             "Pseudomonas aeruginosa Oss (+) Nd (+)"=(Pseudomonas_aeruginosa_expr_Nd_up_Oss_up),
             "Escherichia coli ETEC Oss (+) Mig (+) "=(Escherichia_coli_ETEC_expr_Mig_up_Oss_up),
             "Burkholderia pseudomallei Oss (-)"=(Burkholderia_pseudomallei_expr_Oss_down),
             "Pseudomonas aeruginosa Oss (-)"=(Pseudomonas_aeruginosa_expr_Oss_down)
             
             
             
             
             
      )
    })      
    
    box_plot_Oss <- reactive({
      AcidStress <- datasetInputStimulonOss()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1,color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    js <- c(
      "function(settings){",
      "  var instance = settings.oInstance;",
      "  var table = instance.api();",
      "  var input = instance.parent().find('.dataTables_filter input');",
      "  input.off('keyup search input').on('keyup', function(){",
      "    var keyword = '\\\\b' + input.val() + '\\\\b';",
      "    table.search(keyword, true, false).draw();",
      "  });",
      "}"
    )
    
    color_from_middle <- function (data,color1,color2) 
    {
      max_val=max(abs(data))
      JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                 max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    
    output$Oss_boxplot <- renderPlot(
      
      box_plot_Oss()
    )
    
    datOss <- reactive({
      dataset_As <- datasetInputStimulonOss()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableOss <- DT::renderDataTable({datOss()
      
      return(datatable(datOss(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datOss()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })  
    ############Oxs_stimulon#############
    datasetInputStimulonOxs <- reactive({
      switch(input$select_oxs,
             "Acinetobacter baumannii Oxs (+)"=(Acinetobacter_baumannii_expr_Oxs_up),
             "Aggregatibacter actinomycetemcomitans Oxs (+)"=(Aggregatibacter_actinomycetemcomitans_expr_Oxs_up),
             "Burkholderia pseudomallei Oxs (+)"=(Burkholderia_pseudomallei_expr_Oxs_up),
             "Campylobacter jejuni Oxs (+) (I)"=(Campylobacter_jejuni_expr_Oxs1_up),
             "Campylobacter jejuni Oxs (+) (II)"=(Campylobacter_jejuni_expr_Oxs2_up),
             "Enterococcus faecalis Oxs (+) Vic (-)"=(Enterococcus_faecalis_expr_Oxs_up_Vic_down),
             "Escherichia coli EPEC Oxs (+) (I)"=(Escherichia_coli_EPEC_expr_Oxs1_up),
             "Escherichia coli EPEC Oxs (+) (II)"=(Escherichia_coli_EPEC_expr_Oxs2_up),
             "Escherichia coli ETEC Oxs (+) (I)"=(Escherichia_coli_ETEC_expr_Oxs_up),
             "Helicobacter pylori G27 Oxs (+)"=(Helicobacter_pylori_G27_expr_Oxs_up),
             "Helicobacter pylori J99 Oxs (+)"=(Helicobacter_pylori_J99_expr_Oxs_up),
             "Klebsiella pneumoniae Oxs (+)"=(Klebsiella_pneumoniae_expr_Oxs_up),
             "Listeria monocytogenes Oxs (+) As (+)"=(Listeria_monocytogenes_expr_Oxs_up_As_up),
             "Staphylococcus aureus MRSA252 Oxs (+)"=(Staphylococcus_aureus_MRSA252_expr_Oxs_up),
             "Staphylococcus aureus MRSA252 Oxs (+) Vic (+)"=(Staphylococcus_aureus_MRSA252_expr_Oxs_up_Vic_up),
             "Staphylococcus aureus MSSA476 Oxs (+)"=(Staphylococcus_aureus_MSSA476_expr_Oxs_up),
             "Pseudomonas aeruginosa Oxs (+)"=(Pseudomonas_aeruginosa_expr_Oxs_up),
             "Salmonella enterica Oxs (+)"=(Salmonella_enterica_expr_Oxs_up),
             "Shigella flexneri Oxs (+)"=(Shigella_flexneri_expr_Oxs_up),
             "Streptococcus pyogenes Oxs (+) (I)"=(Streptococcus_pyogenes_expr_Oxs_up),
             "Streptococcus pyogenes Oxs (+) (II)"=(Streptococcus_pyogenes_expr_Oxs2_up),
             "Escherichia coli UPEC Oxs (+)"=(Escherichia_coli_UPEC_expr_Oxs_up),
             "Yersinia pseudotuberculosis Oxs (+)"=(Yersinia_pseudotuberculosis_expr_Oxs_up),
             "Enterococcus faecalis  Oxs (+) Nd (+) Vic(-)"=(Enterococcus_faecalis_Nd_up_Oxs_up_Vic_down),
             "Klebsiella pneumoniae  Oxs (+) Nd (+)"=(Klebsiella_pneumoniae_expr_Nd_up_Oxs_up),
             "Salmonella enterica  Oxs (+) Nd (+)"=(Salmonella_enterica_expr_Nd_up_Oxs_up),
             "Shigella flexneri Oxs (+) Nd (+) "=(Shigella_flexneri_expr_Nd_up_Oxs_up),
             "Acinetobacter baumannii  Oxs (+) Li (+)"=(Acinetobacter_baumannii_expr_Li_up_Oxs_up),
             "Salmonella enterica  Oxs (+) Li (+)"=(Salmonella_enterica_expr_Li_up_Oxs_up),
             "Shigella flexneri  Oxs (+) Li (+)"=(Shigella_flexneri_expr_Li_up_Oxs_up),
             "Neisseria meningitidis Oxs (-) Tm (+) "=(Neisseria_meningitidis_expr_Tm_up_Oxs_down),
             "Escherichia coli ETEC Oxs (+) (II)"=(Escherichia_coli_ETEC_expr_Oxs2_up),
             "Achromobacter xylosoxidans Oxs (-)"=(Achromobacter_xylosoxidans_expr_Oxs_down),
             "Streptococcus agalactiae Oxs (+) (I)"=(Streptococcus_agalactiae_expr_Oxs_up),
             "Streptococcus agalactiae Oxs (+) (II)"=(Streptococcus_agalactiae_expr_Oxs2_up)
             
             
             
             
             
             
      )
    })      
    
    box_plot_Oxs <- reactive({
      AcidStress <- datasetInputStimulonOxs()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white") +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    js <- c(
      "function(settings){",
      "  var instance = settings.oInstance;",
      "  var table = instance.api();",
      "  var input = instance.parent().find('.dataTables_filter input');",
      "  input.off('keyup search input').on('keyup', function(){",
      "    var keyword = '\\\\b' + input.val() + '\\\\b';",
      "    table.search(keyword, true, false).draw();",
      "  });",
      "}"
    )
    
    color_from_middle <- function (data,color1,color2) 
    {
      max_val=max(abs(data))
      JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                 max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    
    output$Oxs_boxplot <- renderPlot(
      
      box_plot_Oxs()
    )
    
    datOxs <- reactive({
      dataset_As <- datasetInputStimulonOxs()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableOxs <- DT::renderDataTable({datOxs()
      
      return(datatable(datOxs(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datOxs()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })
    ############Tm_stimulon#############
    datasetInputStimulonTm <- reactive({
      switch(input$select_tm,
             "Achromobacter xylosoxidans Tm (+) (I)"=(Achromobacter_xylosoxidans_expr_Tm1_up),
             "Achromobacter xylosoxidans Tm (+) (II)"=(Achromobacter_xylosoxidans_expr_Tm2_up),
             "Achromobacter xylosoxidans Tm (+) (III)"=(Achromobacter_xylosoxidans_expr_Tm3_up),
             "Aggregatibacter actinomycetemcomitans Tm (+)"=(Aggregatibacter_actinomycetemcomitans_expr_Tm_up),
             "Borrelia burgdorferi Tm (+)"=(Borrelia_burgdorferi_expr_Tm_up),
             "Burkholderia pseudomallei Tm (+)"=(Burkholderia_pseudomallei_expr_Tm_up),
             "Enterococcus faecalis Tm (+)"=(Enterococcus_faecalis_expr_Tm_up),
             "Escherichia coli EPEC Tm (+)"=(Escherichia_coli_EPEC_expr_Tm_up),
             "Escherichia coli ETEC Tm (+)"=(Escherichia_coli_ETEC_expr_Tm_up),
             "Escherichia coli ETEC Tm (+) Li (+)"=(Escherichia_coli_ETEC_expr_Tm_up_Li_up),
             "Francisella tularensis Tm (+)"=(Francisella_tularensis_expr_Tm_up),
             "Haemophilus influenzae Tm (+) (I)"=(Haemophilus_influenzae_expr_Tm1_up),
             "Haemophilus influenzae Tm (+) (II)"=(Haemophilus_influenzae_expr_Tm2_up),
             "Klebsiella pneumoniae Tm (+)"=(Klebsiella_pneumoniae_expr_Tm_up),
             "Listeria monocytogenes Tm (+) (I)"=(Listeria_monocytogenes_expr_Tm1_up),
             "Listeria monocytogenes Tm (+) (II)"=(Listeria_monocytogenes_expr_Tm2_up),
             "Staphylococcus aureus MRSA252 Tm (+)"=(Staphylococcus_aureus_MRSA252_expr_Tm_up),
             "Staphylococcus aureus MSSA476 Tm (+)"=(Staphylococcus_aureus_MSSA476_expr_Tm_up),
             "Neisseria meningitidis Tm (+)"=(Neisseria_meningitidis_expr_Tm_up),
             "Neisseria meningitidis Tm (+) Oxs (-)"=(Neisseria_meningitidis_expr_Tm_up_Oxs_down),
             "Salmonella enterica Tm (+)"=(Salmonella_enterica_expr_Tm_up),
             "Staphylococcus epidermidis Tm (+)"=(Staphylococcus_epidermidis_expr_Tm_up),
             "Staphylococcus epidermidis Tm (-)"=(Staphylococcus_epidermidis_expr_Tm_down),
             "Shigella flexneri Tm (+)"=(Shigella_flexneri_expr_Tm_up),
             "Streptococcus pyogenes Tm (+) Nd (+)"=(Streptococcus_agalactiae_expr_Tm_up_Nd_down),
             "Streptococcus agalactiae Tm (+)"=(Streptococcus_agalactiae_expr_Tm_up),
             "Streptococcus agalactiae Tm (+) Nd (-)"=(Streptococcus_agalactiae_expr_Tm_up_Nd_down),
             "Streptococcus pneumoniae Tm (+) Vic (+)"=(Streptococcus_pneumoniae_expr_Tm_up_Vic_up),
             "Escherichia coli UPEC Tm (+)"=(Escherichia_coli_UPEC_expr_Tm_up),
             "Vibrio cholerae Tm (+)"=(Vibrio_cholerae_expr_Tm_up),
             "Yersinia pseudotuberculosis Tm (+)"=(Yersinia_pseudotuberculosis_expr_Tm_up),
             "Pseudomonas aeruginosa Tm (+) As (+)"=(Pseudomonas_aeruginosa_expr_As_Tm),
             "Burkholderia pseudomallei Tm (+) Sp (+)"=(Burkholderia_pseudomallei_expr_Sp_Tm),
             "Staphylococcus epidermidis Tm (+) Sp (+)"=(Staphylococcus_epidermidis_expr_Sp_up_Tm_up),
             "Staphylococcus epidermidis Tm (-) Nd (+)"=(Staphylococcus_epidermidis_expr_Nd_up_Tm_down),
             "Staphylococcus epidermidis Tm (-) Mig (+) "=(Staphylococcus_epidermidis_expr_Mig_up_Tm_down),
             "Staphylococcus epidermidis Tm (-) Ns (+) "=(Staphylococcus_epidermidis_expr_Ns_up_Tm_down),
             "Klebsiella pneumoniae  Tm (+) Oss (+)"=(Klebsiella_pneumoniae_expr_Oss_up_Tm_up),
             "Shigella flexneri Tm (+) Oss (+) " =(Shigella_flexneri_expr_Oss_up_Tm_up),
             "Shigella flexneri Tm (+) Vic (+) "=(Shigella_flexneri_expr_Vic_up_Tm_up),
             "Staphylococcus epidermidis Tm (-) As (+)"=(Staphylococcus_epidermidis_expr_As_up_Tm_down)
             
      )
    })      
    
    box_plot_Tm <- reactive({
      AcidStress <- datasetInputStimulonTm()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white") +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white" ) +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    js <- c(
      "function(settings){",
      "  var instance = settings.oInstance;",
      "  var table = instance.api();",
      "  var input = instance.parent().find('.dataTables_filter input');",
      "  input.off('keyup search input').on('keyup', function(){",
      "    var keyword = '\\\\b' + input.val() + '\\\\b';",
      "    table.search(keyword, true, false).draw();",
      "  });",
      "}"
    )
    
    color_from_middle <- function (data,color1,color2) 
    {
      max_val=max(abs(data))
      JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                 max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    
    output$Tm_boxplot <- renderPlot(
      
      box_plot_Tm()
    )
    
    datTm <- reactive({
      dataset_As <- datasetInputStimulonTm()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableTm <- DT::renderDataTable({datTm()
      
      return(datatable(datTm(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datTm()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })  
    ############Vic_stimulon#############
    datasetInputStimulonVic <- reactive({
      switch(input$select_vic,
             "Aggregatibacter actinomycetemcomitans Vic (+) (I)"=(Aggregatibacter_actinomycetemcomitans_expr_Vic1_up),
             "Aggregatibacter actinomycetemcomitans Vic (+) (II)"=(Aggregatibacter_actinomycetemcomitans_expr_Vic2_up),
             "Aggregatibacter actinomycetemcomitans Vic (+) (III)"=(Aggregatibacter_actinomycetemcomitans_expr_Vic3_up),
             "Burkholderia pseudomallei Vic (+)"=(Burkholderia_pseudomallei_expr_Vic_up),
             "Escherichia coli EPEC Vic (+)"=(Escherichia_coli_EPEC_expr_Vic_up),
             "Escherichia coli ETEC Vic (+) (I)"=(Escherichia_coli_ETEC_expr_Vic1_up),
             "Escherichia coli ETEC Vic (+) (II)"=(Escherichia_coli_ETEC_expr_Vic2_up),
             "Haemophilus influenzae Vic (+)"=(Haemophilus_influenzae_expr_Vic_up),
             "Haemophilus influenzae Vic (+) Mig (-)"=(Haemophilus_influenzae_expr_Vic_up_Mig_down),
             "Helicobacter pylori J99 Vic (+)"=(Helicobacter_pylori_J99_expr_Vic_up),
             "Klebsiella pneumoniae Vic (+)"=(Klebsiella_pneumoniae_expr_Vic_up),
             "Staphylococcus aureus MRSA252 Vic (+)"=(Staphylococcus_aureus_MRSA252_expr_Vic_up),
             "Staphylococcus aureus MSSA476 Vic (+)"=(Staphylococcus_aureus_MSSA476_expr_Vic_up),
             "Pseudomonas aeruginosa Vic (+)"=(Pseudomonas_aeruginosa_expr_Vic_up),
             "Salmonella enterica Vic (+)"=(Salmonella_enterica_expr_Vic_up),
             "Shigella flexneri Vic (+) (I)"=(Shigella_flexneri_expr_Vic1_up),
             "Shigella flexneri Vic (+) (II)"=(Shigella_flexneri_expr_Vic2_up),
             "Shigella flexneri Vic (+) Tm (+)"=(Shigella_flexneri_expr_Vic_up_Tm_up),
             "Streptococcus agalactiae Vic (+) Nd (-)"=(Streptococcus_agalactiae_expr_Vic_up_Nd_down),
             "Streptococcus pneumoniae Vic (+) Mig (-)"=(Streptococcus_pneumoniae_expr_Vic_up_Mig_down),
             "Yersinia pseudotuberculosis Vic (+) (I)"=(Yersinia_pseudotuberculosis_expr_Vic1_up),
             "Yersinia pseudotuberculosis Vic (+) (II)"=(Yersinia_pseudotuberculosis_expr_Vic2_up),
             "Yersinia pseudotuberculosis Vic (+) Mig (-)"=(Yersinia_pseudotuberculosis_expr_Vic_up_Mig_down),
             "Listeria monocytogenes Vic (+) Sp (+)"=(Listeria_monocytogenes_expr_Sp_Vic_up),
             "Yersinia pseudotuberculosis Vic (-) Nd (+)"=(Yersinia_pseudotuberculosis_expr_Nd_up_Vic_down),
             "Enterococcus faecalis Vic(-) Nd (+) Oxs (+)"=(Enterococcus_faecalis_expr_Nd_up_Oxs_up_Vic_down),
             "Enterococcus faecalis  Vic (+) Mig (+)"=(Enterococcus_faecalis_expr_Mig_up_Vic_up),
             "Haemophilus influenzae  Vic (+) Mig (+)"=(Haemophilus_influenzae_expr_Mig_up_Vic_up),
             "Escherichia coli UPEC  Vic (+) Mig (+)"=(Escherichia_coli_UPEC_expr_Mig_up_Vic_up),
             "Streptococcus agalactiae Vic (-) Oss (+) "=(Streptococcus_agalactiae_expr_Oss_up_Vic_down),
             "Enterococcus faecalis Vic (-) Oxs (+) "=(Enterococcus_faecalis_expr_Oxs_up_Vic_down),
             "Staphylococcus aureus MRSA252 Vic (+) Oxs (+) "=(Staphylococcus_aureus_MRSA252_expr_Oxs_up_Vic_up),
             "Streptococcus pneumoniae Vic (+) Tm (+) "=(Streptococcus_pneumoniae_expr_Tm_up_Vic_up),
             "Streptococcus pyogenes Vic (-)"=(Streptococcus_pyogenes_expr_Vic_down),
             "Vibrio cholerae Vic (+)"=(Vibrio_cholerae_expr_Vic_up),
             
             "Streptococcus pneumoniae Vic (-)"=(Streptococcus_pyogenes_expr_Vic_down)
             
             
             
             
             
      )
    })      
    
    box_plot_Vic <- reactive({
      AcidStress <- datasetInputStimulonVic()
      if (nrow(AcidStress) == 36 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*36)
        Module_names <- rep(colnames(AcidStress),36)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*36+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white") +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      else if (nrow(AcidStress) == 33 ) {
        
        generic_stress_sample_names=c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        #Vibrio_expr_1 = scale(Vibrio_expr[,colorh1==which.module ])
        ME_lengthOfVector=rep(NA, ncol(AcidStress)*33)
        Module_names <- rep(colnames(AcidStress),33)
        lengthOfVector_col=ncol(AcidStress)
        position_ME=1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <-      colnames(AcidStress)[i]
            position_ME = position_ME+1
          }
          position_ME = i*33+1
        }
        stress=rep(generic_stress_sample_names,ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress=stress_3
        
        ME_Kruskall_dataset = as.data.frame(cbind(ME_lengthOfVector,stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=ME_Kruskall_dataset[,2], y=as.numeric(ME_Kruskall_dataset[,1]), fill = ME_Kruskall_dataset[,2])) +  labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin( width=1, color="white") +geom_boxplot(width = .1, color= "black", alpha = .5) +scale_fill_brewer(palette = "Paired", aesthetics = "fill",guide = "none") +theme_classic(base_size = 16)
      } 
      
      
      # boxplot(t(AcidStress))
    })
    
    
    js <- c(
      "function(settings){",
      "  var instance = settings.oInstance;",
      "  var table = instance.api();",
      "  var input = instance.parent().find('.dataTables_filter input');",
      "  input.off('keyup search input').on('keyup', function(){",
      "    var keyword = '\\\\b' + input.val() + '\\\\b';",
      "    table.search(keyword, true, false).draw();",
      "  });",
      "}"
    )
    
    color_from_middle <- function (data,color1,color2) 
    {
      max_val=max(abs(data))
      JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                 max_val,color1,max_val,color1,color2,color2,max_val,max_val))
    } 
    
    
    brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
    clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
    
    
    
    
    output$Vic_boxplot <- renderPlot(
      
      box_plot_Vic()
    )
    
    datVic <- reactive({
      dataset_As <- datasetInputStimulonVic()
      
      #Module_Colors = labels2colors(net$colors)
      Locus_Tag_As= colnames(dataset_As)
      dataset_transpose_As=t(dataset_As)
      
      rownames(dataset_transpose_As) <- NULL
      
      colnames(dataset_transpose_As) <- NULL
      #colnames(table_2_inspect) <- c("Module color","Gene ID")
      #generic_stress_sample_names=c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      
      probes = Locus_Tag_As
      Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
      
      table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
      
      if (nrow(dataset_As) == 36 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
      
      else if (nrow(dataset_As) == 33 ) {
        colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        return(as.data.frame(table_2_inspect_As))
      } 
    })
    
    
    
    output$mytableVic <- DT::renderDataTable({datVic()
      
      return(datatable(datVic(), extensions = 'Buttons',options = list(paging=TRUE, lengthMenu=c(50,100,150), scrollX = FALSE,lengthChange = T, regex= FALSE, searching = TRUE,initComplete = JS(js),dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')),rownames=FALSE, class = "display") %>% formatStyle(names(datVic()[-1:-3]), backgroundColor = styleInterval(brks, clrs)))
      #options = list(paging=FALSE, scrollX = FALSE), 
      #rownames=TRUE, 
      #filter = "top"
    })   
    
    
    
        
}
# Run the application 
shinyApp(ui = ui, server = server)
