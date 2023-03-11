#####MSSA476andMRSA252#######

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
#############################



#Genes=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
PGFam_core=Annotations_all_PGFAM_new$ref_genome_function[match(probes_MRSA252,Annotations_all_PGFAM_new$comp_genome_1_patric_id)]
Genes_core=Annotations_all_PGFAM_new$ref_genome_pgfam_id[match(probes_MRSA252,Annotations_all_PGFAM_new$comp_genome_1_patric_id)]
Annotations_core_MRSA252 = as.data.frame(cbind(Genes_core,PGFam_core))

Annotations_core_all=as.data.frame(rbind(Annotations_core_ETEC,Annotations_core_EPEC,Annotations_core_UPEC,Annotations_core_HP_J99, Annotations_core_HP_G27,Annotations_core_MRSA252,Annotations_core_MSSA476 ))
Annotations_core=distinct(Annotations_core_all)

write.csv(Annotations_core,"Annotations_core.csv", row.names = FALSE)
###################MRSA252_MSSA476################################
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


