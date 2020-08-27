#hallmarks and miR maturation analysis
#andrew dhawan
#11 june 2019

#code for the paper, "A Dicer-to-Argonaute genomic switch regulates miRNA biogenesis in cancer" for analyses that I performed.

#load necessary libraries
library(ppcor)
library(reshape2)
library(penalized)
library(RankProd)
library(ggnewscale)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)

#set directory
setwd("~/Google Drive/DPhil-AndrewDhawan/Code + Data/Reprocessed GDAC Data")

cancer_types <- c('BRCA','UCEC','KIRC','LUAD','HNSC','OV','COAD','BLCA') #cancer types 

#first load in the datasets of miRs and immature miRs
#calculate all the mature:immature ratios
all_ratio_matrix <- list()
all_miRNA_data_mature <- list()
all_miRNA_data_immature <- list()
all_miRNA <- c()
for (cancer_type in cancer_types){
	if (cancer_type!='BRCA'){
		miRNA_fName <- 	 paste0(cancer_type,'/miRNA/tumour/cleaned_miRNA_mature.txt')
		mature_miRNA_data <- read.table(miRNA_fName, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
		colnames(mature_miRNA_data) <- gsub('[.]','-',colnames(mature_miRNA_data))

		miRNA_fName <- 	 paste0(cancer_type,'/miRNA/tumour/cleaned_miRNA.txt')
		immature_miRNA_data <- read.table(miRNA_fName, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
		colnames(immature_miRNA_data) <- gsub('[.]','-',colnames(immature_miRNA_data))
	}else{
		miRNA_fName <- 	 paste0(cancer_type,'/miRNA/tumour/cleaned_miRNA_mature_ductal.txt')
		mature_miRNA_data <- read.table(miRNA_fName, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
		colnames(mature_miRNA_data) <- gsub('[.]','-',colnames(mature_miRNA_data))

		miRNA_fName <- 	 paste0(cancer_type,'/miRNA/tumour/cleaned_miRNA_ductal.txt')
		immature_miRNA_data <- read.table(miRNA_fName, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
		colnames(immature_miRNA_data) <- gsub('[.]','-',colnames(immature_miRNA_data))
	}

	all_miRNA_data_mature[[cancer_type]] <- mature_miRNA_data
	all_miRNA_data_immature[[cancer_type]] <- immature_miRNA_data

	common_samples <- intersect(colnames(mature_miRNA_data),colnames(immature_miRNA_data))
	mature_miRNA_data <- mature_miRNA_data[,common_samples]
	immature_miRNA_data <- immature_miRNA_data[,common_samples]

	#first step is that we need to create the map between mature and immature  miRNA sequences
	# to do this we load first a database
	mir_database <- read.table('miRNA.txt', sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	#next we prepare a list with the index as the miRNA and the matrix inside as all immature miRNA sequences that give rise to it
 
	mirna_map <- lapply(rownames(mature_miRNA_data),function(x) intersect(unique(c(mir_database$ID[which(mir_database$Mature1_ID==x)],
	 					mir_database$ID[which(mir_database$Mature2_ID==x)])),rownames(immature_miRNA_data)))
	names(mirna_map) <- rownames(mature_miRNA_data)
	count <- 0
	mature_conv_names <- c()
	for(k in 1:length(mirna_map)){
		if (length(mirna_map[[k]]>0)){
			count <- count + 1
			mature_conv_names <- c(mature_conv_names,names(mirna_map)[k])
		}
	}

	ratio_matrix <- matrix(0, nrow=count,ncol=length(colnames(mature_miRNA_data)))
	row.names(ratio_matrix) <- mature_conv_names
	colnames(ratio_matrix) <- colnames(mature_miRNA_data)
	for ( i in 1:length(mirna_map)){
		if (length(mirna_map[[i]]) > 0){
			ratio_matrix[names(mirna_map)[i],] <- as.matrix(mature_miRNA_data[i,]/colSums(immature_miRNA_data[mirna_map[[i]],]))
		}
	}

	ratio_matrix[!(is.finite(ratio_matrix))] <- NA
	ratio_matrix <- ratio_matrix[rowSums(ratio_matrix, na.rm=T)!=0,]
	threshold <- 0.5
	ratio_matrix <- ratio_matrix[which(rowSums(is.na(ratio_matrix)) < threshold * length(colnames(ratio_matrix))),]
	all_ratio_matrix[[cancer_type]] <- ratio_matrix
	all_miRNA <- unique(c(all_miRNA,rownames(all_ratio_matrix[[cancer_type]])))
}

#print just the number of cases
for(cancer_type in cancer_types){
	print(cancer_type)
	common_cols = 	intersect(colnames(all_miRNA_data_immature[[cancer_type]]),colnames(all_miRNA_data_mature[[cancer_type]]))
	print(dim(all_miRNA_data_immature[[cancer_type]]))
	print(dim(all_miRNA_data_mature[[cancer_type]]))
	print(length(common_cols))
}

#load in the mRNA
all_mRNA_data <- list();
for (cancer_type in cancer_types){
	#print(cancer_type)
	if(cancer_type!='BRCA'){
		fname_mrna <- paste0('../Reprocessed GDAC data/',cancer_type,'/mRNA/tumour/cleaned_mRNA.txt')
	}else{
		fname_mrna <- paste0('../Reprocessed GDAC data/',cancer_type,'/mRNA/tumour/cleaned_mRNA_ductal.txt')
	}
	all_mRNA_data[[cancer_type]] <- read.table(fname_mrna, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(all_mRNA_data[[cancer_type]]) <- gsub('[.]','-',colnames(all_mRNA_data[[cancer_type]]))
	# want log2 data
	all_mRNA_data[[cancer_type]] <- log2(all_mRNA_data[[cancer_type]]+1)
	all_mRNA_data[[cancer_type]][!is.finite(as.matrix(all_mRNA_data[[cancer_type]]))] <- NA
}

#load in gene signatures
sig_fnames_list <- list();
sig_names_list <- list();
categories_of_sigs <- c('invasion','energetics','immortality','growth_suppressors','genome_instability','angiogenesis','apoptosis','proliferation','inflammation')

sig_fnames_list[['invasion']] <- c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt','invasiveness_gene_sig_entrez_marsan2014.txt')
sig_names_list[['invasion']] <- c('Hallmark: Epithelial Mesenchymal Transition','Invasiveness, Marsan 2014')

sig_fnames_list[['energetics']] <- c('HALLMARK_OXIDATIVE_PHOSPHORYLATION.txt','HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY.txt')
sig_names_list[['energetics']] <- c('Hallmark: Oxidative Phosphorylation','Hallmark: Reactive Oxygen Species Pathway')

sig_fnames_list[['immortality']] <- c('HALLMARK_G2M_CHECKPOINT.txt')
sig_names_list[['immortality']] <- c('Hallmark: G2M Checkpoint')

sig_fnames_list[['growth_suppressors']] <- c('HALLMARK_PI3K_AKT_MTOR_SIGNALING.txt','HALLMARK_XENOBIOTIC_METABOLISM.txt')
sig_names_list[['growth_suppressors']] <- c('Hallmark: PI3K AKT MTOR Signaling','Hallmark: Xenobiotic Metabolism')

sig_fnames_list[['genome_instability']] <- c('HALLMARK_DNA_REPAIR.txt','HALLMARK_P53_PATHWAY.txt')
sig_names_list[['genome_instability']] <- c('Hallmark: DNA Repair','Hallmark: p53 Pathway')

sig_fnames_list[['angiogenesis']] <- c('hypoxia_gene_sig_entrez_probes.txt','HALLMARK_ANGIOGENESIS.txt','HALLMARK_HYPOXIA.txt','angiogenesis_gene_sig_entrez_desmedt2008_pos.txt','Masiero2013angiogenesisENTREZ.txt')
sig_names_list[['angiogenesis']] <- c('Hypoxia, Buffa 2010','Hallmark: Angiogenesis','Hallmark: Hypoxia','Angiogenesis, Desmedt 2008','Angiogenesis, Masiero 2013')

sig_fnames_list[['apoptosis']] <- c('HALLMARK_APOPTOSIS.txt','apoptosis_gene_sig_entrez_desmedt2008_pos.txt')
sig_names_list[['apoptosis']] <- c('Hallmark: Apoptosis','Apoptosis, Desmedt 2008')

sig_fnames_list[['proliferation']] <-c('proliferation_gene_sig_entrez_desmedt2008_pos.txt','HALLMARK_KRAS_SIGNALING_UP.txt')
sig_names_list[['proliferation']] <-c('Proliferation, Desmedt 2008','Hallmark: KRAS Signaling Up')

sig_fnames_list[['inflammation']] <- c('HALLMARK_INFLAMMATORY_RESPONSE.txt','HALLMARK_IL2_STAT5_SIGNALING.txt','HALLMARK_IL6_JAK_STAT3_SIGNALING.txt','HALLMARK_TGF_BETA_SIGNALING.txt','HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt','immune_gene_sig_entrez_desmedt2008_pos.txt')
sig_names_list[['inflammation']] <- c('Hallmark: Inflammatory Response','Hallmark: IL2 STAT5 Signaling','Hallmark: IL6 JAK STAT3 Signaling','Hallmark: TGF Beta Signaling','Hallmark: TNFa Signaling via NFKB','Immune, Desmedt 2008')

sigs_list_by_cat <- list();
for(sig_category in categories_of_sigs){
	sigs_list_by_cat[[sig_category]] <- list();
	for(i in 1:length(sig_fnames_list[[sig_category]])){
		fname <- sig_fnames_list[[sig_category]][i]
		genes = read.csv(paste0('../miRNA_hallmarks/gene_signatures/',fname), header=F, stringsAsFactors=F, colClasses = "character")
		# print(genes)
		sigs_list_by_cat[[sig_category]][[sig_names_list[[sig_category]][i]]]<- genes
	} 
}
all_signatures <- melt(sig_names_list)$value

#------LIMIT THE MIRNA WE ARE CONSIDERING ONLY TO THOSE WHICH ARE IN THE RATIO MATRIX--------
tmp_mat <- melt(mirna_map[all_miRNA])
immature_miRNA_considered <-unique(tmp_mat[,1])
mature_miRNA_considered <- all_miRNA#unique(tmp_mat[,2])
for(cancer_type in cancer_types){
	all_miRNA_data_mature[[cancer_type]] <- all_miRNA_data_mature[[cancer_type]][intersect(rownames(all_miRNA_data_mature[[cancer_type]]),mature_miRNA_considered),]
	all_miRNA_data_immature[[cancer_type]] <- all_miRNA_data_immature[[cancer_type]][intersect(rownames(all_miRNA_data_immature[[cancer_type]]),immature_miRNA_considered),]
}
#---------------------------------------------------------------------------------------------

#calculate all the gene singature scores
# do a penalised regression
# save all the coefficients
for(model_terms in c('mature','ratio','immature')){
	all_coeffs <- list();
	all_rank_product_matrices <- list();
	for (category in categories_of_sigs){
		count <- 1
		for (gene_sig in sigs_list_by_cat[[category]]){

			sig_name <- sig_names_list[[category]][count]
			print(sig_name)
			if(model_terms != 'immature'){
				all_miR_vals <- all_miRNA
			}else{
				all_miR_vals <- immature_miRNA_considered
			}

			all_coeffs[[sig_name]] <- matrix(0,nrow=length(all_miR_vals),ncol=length(cancer_types))
			row.names(all_coeffs[[sig_name]]) <- all_miR_vals
			colnames(all_coeffs[[sig_name]]) <- cancer_types
	 		for (cancer_type in cancer_types){
	 			print(cancer_type)
				genes_present <- intersect(rownames(all_mRNA_data[[cancer_type]]),gene_sig$V1)
				#compute and score the scores
				scores <-  apply(all_mRNA_data[[cancer_type]][genes_present,], 2, function(x) median(x,na.rm=T))
				print(head(scores))
				#cross-validated linear model
				if(model_terms=='ratio'){
					miRNA_df <- all_ratio_matrix[[cancer_type]]
				}else if (model_terms == 'mature'){
					miRNA_df <- all_miRNA_data_mature[[cancer_type]]
				}else{
					miRNA_df <- all_miRNA_data_immature[[cancer_type]]
				}
				coeffs <- get_coefficients_pre_filter(cancer_type,scores,miRNA_df)
				#store the miRNA results
				all_coeffs[[sig_name]][names(coeffs),cancer_type] <- coeffs
			}
			# compute the rank-product matrix
			all_rank_product_matrices[[sig_name]] <- make_rank_prod_matrix(all_coeffs[[sig_name]])
			count <- count + 1
		}
	}

	# rank product
	rank_prod_tables <- list();
	RP_out_values <- list();
	#here we need to do the rankprod
	for (sig_name in all_signatures){
		all_coeffs[[sig_name]] <- all_coeffs[[sig_name]][which(rowSums(all_coeffs[[sig_name]]==0) < length(colnames(all_coeffs[[sig_name]]))),]
		print(dim(all_coeffs[[sig_name]]))
		RP.out <- RP(all_coeffs[[sig_name]],rep(1,length(colnames(all_coeffs[[sig_name]]))))
		RP_out_values[[sig_name]] <- RP.out
		rank_prod_tables[[sig_name]] <- topGene(RP.out,cutoff = 0.05,method="pfp",gene.names=rownames(all_coeffs[[sig_name]]))
	}

	#then save the outputs
	save(file=paste0('../../Laura AGO2 Paper/',model_terms,'_coeffs.rda'),all_coeffs)
	# save(file=paste0('../../Laura AGO2 Paper/',model_terms,'_rank_prod_output_pre_filtered.rda'),RP_out_values)
	# save(file=paste0('../../Laura AGO2 Paper/',model_terms,'_rank_prod_tables_out_pre_filtered.rda'),rank_prod_tables)
}

#now let's find which miRs are common among the three groups, first let's define the sets

up_miR_lists <- list()
down_miR_lists <- list()

for(model_terms in c('mature','immature')){ #ratio
	load(paste0('../../Laura AGO2 Paper/',model_terms,'_rank_prod_tables_out_pre_filtered.rda'))

	all_up_miRNA <- c()
	all_down_miRNA <- c()

	for(sig_name in all_signatures){
		all_down_miRNA <- c(all_down_miRNA, rownames(rank_prod_tables[[sig_name]]$Table1))
		all_up_miRNA <- c(all_up_miRNA, rownames(rank_prod_tables[[sig_name]]$Table2))
	}

	up_miR_lists[[model_terms]] <- all_up_miRNA
	down_miR_lists[[model_terms]] <- all_down_miRNA
	print(length(intersect(all_up_miRNA,all_down_miRNA)))
}

h <- list()
for(model_terms in c('mature','immature')){
	p = load(paste0('../../Laura AGO2 Paper/',model_terms,'_coeffs.rda'))

	miRs_list <- c()
	for(sig_name in names(all_coeffs)){
		miRs_list <- unique(c(miRs_list,rownames(all_coeffs[[sig_name]])))
	}
	coeffs_averaged_mat <- matrix(0,nrow=length(miRs_list),ncol=length(names(all_coeffs)))
	row.names(coeffs_averaged_mat) = miRs_list
	colnames(coeffs_averaged_mat) = names(all_coeffs)
	for(sig_name in names(all_coeffs)){
		coeffs_averaged_mat[rownames(all_coeffs[[sig_name]]),sig_name] = as.matrix(rowMedians(all_coeffs[[sig_name]]))
	}
	coeffs_averaged_mat = coeffs_averaged_mat[unique(c(up_miR_lists[[model_terms]],down_miR_lists[[model_terms]]) ),]
	col2 = rainbow(length(names(all_coeffs)))
	names(col2) =names(all_coeffs)

	# coeffs_averaged_mat <- t(apply(coeffs_averaged_mat ,1,function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)}))
	# dev.new()

	coeffs_rank_mat <- apply(coeffs_averaged_mat,2,function(x){rank(x)})
	coeffs_rank_prod = apply(coeffs_rank_mat,1,function(x){prod(x)})

	num_to_take = 12
	if(model_terms=='immature'){
		num_to_take=8
	}
	plot_mat = coeffs_averaged_mat[order(coeffs_rank_prod,decreasing=T)[1:num_to_take],]
	maxVal = max(abs(plot_mat))
	cols= colorRamp2(breaks=c(-0.1,0,0.1),colors=c('blue','white','red'))
	top_annotation = HeatmapAnnotation(Signature = names(all_coeffs),col=list(Signature=col2))
	show_leg = T
	if(model_terms=='immature'){
		top_annotation=NULL
		show_leg=F
	}
	h[[model_terms]]= Heatmap(plot_mat, col=cols,top_annotation=top_annotation,show_row_names=T,show_column_names=F, show_row_dend=F,row_names_gp = gpar(fontsize=9),cluster_columns= F,show_column_dend=F,cluster_rows=F,show_heatmap_legend=show_leg,heatmap_legend_param = list(title = "Coefficient",title_position = "topcenter"))#left_annotation=ha)
	# draw(h)
	# dev.copy(pdf,file=paste0('../../Laura AGO2 Paper/heatmap_miRNA_model_coefficients_',model_terms,'_median.pdf'),height=6,width=10)
	# dev.off()
}
dev.new()
h[['mature']] %v% h[['immature']]
dev.copy(pdf,file=paste0('../../Laura AGO2 Paper/heatmap_miRNA_model_coefficients_combined_hmap_median.pdf'),height=6,width=10)
dev.off()

#need to do the conversion in terms of immature species probably to determine the commonalities and differences as this is injective 
#first do the conversion
up_miR_lists_converted_immature <- list()
down_miR_lists_converted_immature <- list()

for(model_terms in c('mature','immature')){ 
	if(model_terms!='immature'){
		tmp_mat <- melt(mirna_map[up_miR_lists[[model_terms]]])
		up_miR_lists_converted_immature[[model_terms]] <- as.character(unique(tmp_mat[,1]))

		tmp_mat <- melt(mirna_map[down_miR_lists[[model_terms]]])
		down_miR_lists_converted_immature[[model_terms]] <- as.character(unique(tmp_mat[,1]))
	}else{
		up_miR_lists_converted_immature[[model_terms]] <- unique(up_miR_lists[[model_terms]])
		down_miR_lists_converted_immature[[model_terms]] <- unique(down_miR_lists[[model_terms]])
	}
}

print(length(intersect(up_miR_lists_converted_immature$mature,up_miR_lists_converted_immature$immature)))
print(length(setdiff(up_miR_lists_converted_immature$mature,up_miR_lists_converted_immature$immature)))
print(length(setdiff(up_miR_lists_converted_immature$immature,up_miR_lists_converted_immature$mature)))

print(length(intersect(down_miR_lists_converted_immature$mature,down_miR_lists_converted_immature$immature)))
print(length(setdiff(down_miR_lists_converted_immature$mature,down_miR_lists_converted_immature$immature)))
print(length(setdiff(down_miR_lists_converted_immature$immature,down_miR_lists_converted_immature$mature)))


#now we want concordance by signature type...
proportion_in_common_up <- c()
proportion_in_common_down <- c()
down_miRNA <- list()
up_miRNA <- list()
lengths_of_immature_up <- c()
lengths_of_immature_down <- c()
proportion_up_mature_down_immature <- c()
proportion_down_mature_up_immature <- c()
for(sig_name in all_signatures){
	down_miRNA[[sig_name]] <- list()
	up_miRNA[[sig_name]] <- list()
	for(model_terms in c('mature','immature')){ 
		load(paste0('../../Laura AGO2 Paper/',model_terms,'_rank_prod_tables_out_pre_filtered.rda'))
		down_miRNA[[sig_name]][[model_terms]] <- rownames(rank_prod_tables[[sig_name]]$Table1)
		up_miRNA[[sig_name]][[model_terms]] <-  rownames(rank_prod_tables[[sig_name]]$Table2)
	}
	tmp_mat <- melt(mirna_map[down_miRNA[[sig_name]][['mature']]])
	down_miRNA[[sig_name]][['mature_converted_immature']] <- as.character(unique(tmp_mat[,1]))
	tmp = length(intersect(down_miRNA[[sig_name]][['mature_converted_immature']],down_miRNA[[sig_name]][['immature']]))/length(down_miRNA[[sig_name]][['immature']])
	proportion_in_common_down <- c(proportion_in_common_down,tmp)
	lengths_of_immature_down <- c(lengths_of_immature_down,length(down_miRNA[[sig_name]][['immature']]))

 	tmp_mat <- melt(mirna_map[up_miRNA[[sig_name]][['mature']]])
	up_miRNA[[sig_name]][['mature_converted_immature']]  <- as.character(unique(tmp_mat[,1]))
	tmp = length(intersect(up_miRNA[[sig_name]][['mature_converted_immature']],up_miRNA[[sig_name]][['immature']]))/length(up_miRNA[[sig_name]][['immature']])
	proportion_in_common_up <- c(proportion_in_common_up,tmp)
	lengths_of_immature_up <- c(lengths_of_immature_up, length(up_miRNA[[sig_name]][['immature']]))

	tmp = length(intersect(up_miRNA[[sig_name]][['mature_converted_immature']],down_miRNA[[sig_name]][['immature']]))/length(down_miRNA[[sig_name]][['immature']])
	proportion_up_mature_down_immature <- c(proportion_up_mature_down_immature,tmp)
	if(!is.na(tmp)){
		if(tmp!=0){
			print(sig_name)
			print('up mature, down immature')
			print(intersect(up_miRNA[[sig_name]][['mature_converted_immature']],down_miRNA[[sig_name]][['immature']]))
		}
	}

	tmp = length(intersect(down_miRNA[[sig_name]][['mature_converted_immature']],up_miRNA[[sig_name]][['immature']]))/length(up_miRNA[[sig_name]][['immature']])
	proportion_down_mature_up_immature <- c(proportion_down_mature_up_immature,tmp)
	if(!is.na(tmp)){
		if(tmp!=0){
			print(sig_name)
			print('down mature, up immature')
			print(intersect(down_miRNA[[sig_name]][['mature_converted_immature']],up_miRNA[[sig_name]][['immature']]))
		}
	}
}

#-------------------------------Plotting overlap proportions-----------------------------------

col2 = rainbow(length(all_signatures))
names(col2) =all_signatures

proportions_common <- cbind(up=proportion_in_common_up,down=proportion_in_common_down)
row.names(proportions_common) <- all_signatures

melted_proportions <- melt(proportions_common)
melted_proportions <- cbind.data.frame(melted_proportions,length = lengths_of_immature_up)

dev.new()
p <- ggplot(melted_proportions[which(melted_proportions$Var2=='up'),],aes(x = reorder(Var1, value,  FUN=identity))) + geom_col(aes(y = value,fill=length) ) + coord_flip() +
 theme_minimal()+ ggtitle('Positively associated miRNA, concordance') + labs(y='Overlap proportion',x='')+
   scale_fill_gradient(low = "lightgray", high = "black",name='Number of sig.\nimmature miRNA')+#, space='Lab',name="-log p value",midpoint=median(melted_p_vals_down$prop))
	# ggplot(tmp) + 
	new_scale_fill()+
	geom_col(aes( y = -0.07,fill=Var1),position='stack' ) + 
	scale_fill_manual(values=col2,guide=F)
print(p)
dev.copy(pdf, paste0('concordance_analysis_mature_up_miRNA.pdf'),width=8,height=6)
dev.off()

proportions_common <- cbind(up=proportion_in_common_up,down=proportion_in_common_down)
row.names(proportions_common) <- all_signatures

melted_proportions <- melt(proportions_common)
melted_proportions <- cbind.data.frame(melted_proportions,length = lengths_of_immature_down)
melted_proportions <- melted_proportions[which(!is.na(melted_proportions$value)),]

dev.new()
p <- ggplot(melted_proportions[which(melted_proportions$Var2=='down'),],aes(x = reorder(Var1, value,  FUN=identity))) + geom_col(aes(y = value,fill=length) ) + coord_flip() +
 theme_minimal()+ ggtitle('Negatively associated miRNA, concordance') + labs(y='Overlap proportion',x='')+
  scale_fill_gradient(low = "lightgray", high = "black",name='Number of sig.\nimmature miRNA')+#, space='Lab',name="-log p value",midpoint=median(melted_p_vals_down$prop))
	new_scale_fill()+
	geom_col(aes( y = -0.07,fill=Var1),position='stack' ) + 
	scale_fill_manual(values=col2,guide=F)
print(p)
dev.copy(pdf, paste0('concordance_analysis_mature_down_miRNA.pdf'),width=8,height=6)
dev.off()

output_chart <- c()
for(sig_name in all_signatures){
	output_chart$immature_miRNA_down <- c(output_chart$immature_miRNA_down,paste0(down_miRNA[[sig_name]][['immature']],collapse=', '))
	output_chart$immature_miRNA_up <- c(output_chart$immature_miRNA_up,paste0(up_miRNA[[sig_name]][['immature']],collapse=', '))
	output_chart$mature_miRNA_down <- c(output_chart$mature_miRNA_down,paste0(down_miRNA[[sig_name]][['mature']],collapse=', '))
	output_chart$mature_miRNA_up <- c(output_chart$mature_miRNA_up,paste0(up_miRNA[[sig_name]][['mature']],collapse=', '))
	output_chart$mature_converted_immature_down <- c(output_chart$mature_converted_immature_down,paste0(down_miRNA[[sig_name]][['mature_converted_immature']],collapse=', '))
	output_chart$mature_converted_immature_up <- c(output_chart$mature_converted_immature_up,paste0(up_miRNA[[sig_name]][['mature_converted_immature']],collapse=', '))
}
output_chart <- as.data.frame(output_chart)
row.names(output_chart) <- all_signatures

write.table(output_chart,file='output_chart.txt',sep='\t',quote=F)

#-----------------------Which groups are each miR in---------------------------------------------------

miRNA_up_and_down <- unique(c(as.character(melt(up_miR_lists_converted_immature)[,1]),as.character(melt(down_miR_lists_converted_immature)[,1])))

membership_mat <- matrix(0,nrow=length(miRNA_up_and_down),ncol=4)
row.names(membership_mat) <- miRNA_up_and_down
count <- 0
for(model_terms in c('mature','immature')){ 
	count <- count + 1
	membership_mat[up_miR_lists_converted_immature[[model_terms]],count] <- 1
	membership_mat[down_miR_lists_converted_immature[[model_terms]],count+2] <- 1
}

colmap <- c('white','darkgrey')
names(colmap) <-as.character(0:1)
colnames(membership_mat) <- c('Up, mature','Up, immature','Down, mature','Down immature') 
dev.new()

Heatmap(membership_mat, row_names_gp = gpar(fontsize=6),row_title_gp=gpar(fontsize=9),show_column_names=T,column_names_gp=gpar(fontsize=6), show_row_dend=F,cluster_columns= F,cluster_rows=T, col=structure(c('white','darkgrey'), names=as.character(0:1)),rect_gp = gpar(col = "white", lwd = 2),gap = unit(6, "mm"),heatmap_legend_param = list(labels = c('absent','present'), at=as.character(0:1),title = "Present or absent", direction='vertical',title_position = "topcenter", legend_gp = gpar(fill = colmap)),column_title='Hallmarks associated miRNA')

dev.copy(pdf,file='membership_mat_2.pdf',width=5,height=16)
dev.off()
length(which(membership_mat[,1]==1 & membership_mat[,2]==1)) #up in mature AND immature
length(which(membership_mat[,1]==1 & membership_mat[,2]==0)) #up in mature NOT immature
length(which(membership_mat[,1]==0 & membership_mat[,2]==1)) #up in immature NOT mature
#length(which(membership_mat[,1]==1 & membership_mat[,2]==1))
length(which(membership_mat[,3]==1 & membership_mat[,4]==1)) #down in mature AND immature
length(which(membership_mat[,3]==1 & membership_mat[,4]==0)) #down in mature NOT immature
length(which(membership_mat[,3]==0 & membership_mat[,4]==1)) #down in immature NOT mature
#now let's make soem heatmaps
m1 = names(which(membership_mat[,1]==1 & membership_mat[,2]==1)) #up in mature AND immature
m2 = names(which(membership_mat[,1]==1 & membership_mat[,2]==0)) #up in mature NOT immature
m3 = names(which(membership_mat[,1]==0 & membership_mat[,2]==1)) #up in immature NOT mature
#length(which(membership_mat[,1]==1 & membership_mat[,2]==1))
m4 = names(which(membership_mat[,3]==1 & membership_mat[,4]==1)) #down in mature AND immature
m5 = names(which(membership_mat[,3]==1 & membership_mat[,4]==0)) #down in mature NOT immature
m6 = names(which(membership_mat[,3]==0 & membership_mat[,4]==1)) #down in immature NOT mature

m1_conv  = lapply(m1,function(x){intersect(mir_database[which(mir_database$ID ==x),c('Mature1_ID','Mature2_ID')],up_miR_lists[['mature']])})
m1_conv = as.character(unique(melt(m1_conv)[,1]))

m2_conv  = lapply(m2,function(x){intersect(mir_database[which(mir_database$ID ==x),c('Mature1_ID','Mature2_ID')],up_miR_lists[['mature']])})
m2_conv = as.character(unique(melt(m2_conv)[,1]))

m3_conv  = lapply(m3,function(x){mir_database[which(mir_database$ID ==x),c('Mature1_ID','Mature2_ID')]})
m3_conv = as.character(unique(c(melt(m3_conv)[,1],melt(m3_conv)[,2])))

m4_conv  = lapply(m4,function(x){intersect(mir_database[which(mir_database$ID ==x),c('Mature1_ID','Mature2_ID')],down_miR_lists[['mature']])})
m4_conv = as.character(unique(melt(m4_conv)[,1]))

m5_conv  = lapply(m5,function(x){intersect(mir_database[which(mir_database$ID ==x),c('Mature1_ID','Mature2_ID')],down_miR_lists[['mature']])})
m5_conv = as.character(unique(melt(m5_conv)[,1]))

m6_conv  = lapply(m6,function(x){mir_database[which(mir_database$ID ==x),c('Mature1_ID','Mature2_ID')]})
m6_conv = as.character(unique(c(melt(m6_conv)[,1],melt(m6_conv)[,2])))

all_samples_c <- c()
for(cancer_type in cancer_types){
	cNames = intersect(colnames(all_miRNA_data_mature[[cancer_type]]),colnames(all_miRNA_data_immature[[cancer_type]]))
	all_samples_c <- c(all_samples_c,cNames)
	print(cancer_type)
	print(length(cNames))
}

plotting_mat_immature <- matrix(NA,nrow=(length(m1)  +length(m2)+ length(m3) + length(m4) + length(m5) + length(m6)), ncol=length(all_samples_c))
row.names(plotting_mat_immature) <- c(m1,m2,m3,m4,m5,m6)
colnames(plotting_mat_immature) <- all_samples_c
split_cancer_types <- c()
for(cancer_type in cancer_types){
	cNames = intersect(colnames(all_miRNA_data_mature[[cancer_type]]),colnames(all_miRNA_data_immature[[cancer_type]]))
	cRows = intersect(c(m1,m2,m3,m4,m5,m6),rownames(all_miRNA_data_immature[[cancer_type]]))
	split_cancer_types <- c(split_cancer_types,rep(cancer_type,times=length(cNames)))
	plotting_mat_immature[cRows,cNames] <- as.matrix(all_miRNA_data_immature[[cancer_type]][cRows,cNames])
}

plotting_mat_mature <- matrix(NA,nrow=(length(m1_conv)  + length(m2_conv) + length(m3_conv) + length(m4_conv) + length(m5_conv)+ length(m6_conv)), ncol=length(all_samples_c))
row.names(plotting_mat_mature) <- c(m1_conv,m2_conv,m3_conv,m4_conv,m5_conv,m6_conv)
colnames(plotting_mat_mature) <- all_samples_c
for(cancer_type in cancer_types){
	cNames = intersect(colnames(all_miRNA_data_mature[[cancer_type]]),colnames(all_miRNA_data_immature[[cancer_type]]))
	cRows = intersect(c(m1_conv,m2_conv,m3_conv,m4_conv,m5_conv,m6_conv),rownames(all_miRNA_data_mature[[cancer_type]]))
	plotting_mat_mature[cRows,cNames] <- as.matrix(all_miRNA_data_mature[[cancer_type]][cRows,cNames])
}

plotting_mat_immature_adjusted = t(apply(plotting_mat_immature ,1,function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)}))
plotting_mat_immature_adjusted <- plotting_mat_immature_adjusted[which(rowSums(is.na(plotting_mat_immature_adjusted)) < 0.2 * length(colnames(plotting_mat_immature_adjusted))),]
cols= colorRamp2(breaks=c(-2,0,2),colors=c('blue','white','red')) 


plotting_mat_mature_adjusted = t(apply(plotting_mat_mature ,1,function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)}))
plotting_mat_mature_adjusted <- plotting_mat_mature_adjusted[which(rowSums(is.na(plotting_mat_mature_adjusted)) < 0.2 * length(colnames(plotting_mat_mature_adjusted))),]

plotting_mat <- rbind(plotting_mat_immature_adjusted,plotting_mat_mature_adjusted)

up_or_down_in_mature = as.character(((rownames(plotting_mat) %in% c(m1,m1_conv)) * 3) + ((rownames(plotting_mat) %in% c(m2,m2_conv)) * 2) + ((rownames(plotting_mat) %in% c(m3,m3_conv)) * 1))#,m2,m3,m1_conv,m2_conv,m3_conv)) * 1)
up_or_down_in_immature = as.character(((rownames(plotting_mat) %in% c(m4,m4_conv)) * 3) + ((rownames(plotting_mat) %in% c(m5,m5_conv)) * 2) + ((rownames(plotting_mat) %in% c(m6,m6_conv)) * 1))#as.character((rownames(plotting_mat) %in% c(m4,m5,m6,m4_conv,m5_conv,m6_conv)) * 1)
up_or_down_anno = cbind(Up=up_or_down_in_mature,Down=up_or_down_in_immature)
up_or_down_anno[which(up_or_down_anno=='5')] <- '3'

hMaps <- list()

for(cancer_type in cancer_types){
	cNames = intersect(colnames(all_miRNA_data_mature[[cancer_type]]),colnames(all_miRNA_data_immature[[cancer_type]]))
	if(cancer_type == 'BRCA'){

		hMaps[[cancer_type]] <- Heatmap(plotting_mat[,cNames], col=cols,show_row_names=F,split=split_rows,show_column_names=F, show_row_dend=F,cluster_columns= T,show_column_dend=F,cluster_rows=T,column_title=cancer_type,show_heatmap_legend=T,heatmap_legend_param = list(title = "Expr. Z-score",title_position = "topcenter"))#left_annotation=ha)

		}else{
			hMaps[[cancer_type]] <- Heatmap(plotting_mat[,cNames], col=cols,show_row_names=F,split=split_rows,show_column_names=F, show_row_dend=F,cluster_columns= T,show_column_dend=F,cluster_rows=T,column_title=cancer_type,show_heatmap_legend=F)

		}

}

hMaps[['membership']] <- Heatmap(up_or_down_anno, col=c('3'='snow4','2'='snow3','1'='snow2','0'='white'),show_row_names=F,split=split_rows,show_column_names=T, show_row_dend=F,cluster_columns= T,show_column_dend=F,cluster_rows=T,column_title=' ',show_heatmap_legend=T,width = unit(1, "cm"),heatmap_legend_param = list(labels=c(' ','Immature','Mature','Both'),at=c('0','1','2','3'),title = "Group",title_position = "topcenter"))#left_annotation=ha)
dev.new()
draw(hMaps[['membership']] + hMaps[["BRCA"]]+hMaps[["UCEC"]]+hMaps[["KIRC"]] + hMaps[["LUAD"]] + hMaps[["HNSC"]] + hMaps[["OV"]] + hMaps[["COAD"]] + hMaps[["BLCA"]])#column_split=split_cancer_types, 
dev.copy(pdf,'heatmaps_miRNA_mature_immature_hallmarks.pdf',width=8,height=6)
dev.off()

# load the mutation data
all_mut_data <- list()
all_mut_genes <- c()
for(cancer_type in cancer_types){
	fName_mut <- paste0('../Reprocessed GDAC data/',cancer_type,'/mutation/mutations.txt')
	all_mut_data[[cancer_type]] <- read.table(fName_mut, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(all_mut_data[[cancer_type]]) <- gsub('[.]','-',colnames(all_mut_data[[cancer_type]]))
	all_mut_data[[cancer_type]] <- ((all_mut_data[[cancer_type]] > 0) & (all_mut_data[[cancer_type]] < 16)) * 1
	all_mut_genes <- unique(c(all_mut_genes,rownames(all_mut_data[[cancer_type]])))
}

#load CNV data
all_cnv_vals <- list();
all_cnv_genes <- c()
for (cancer_type in cancer_types){
	all_cnv_vals[[cancer_type]] <- read.table(paste0(cancer_type,'/cnv/tumour/cleaned_cnv.txt'), sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(all_cnv_vals[[cancer_type]]) <- gsub('[.]','-',colnames(all_cnv_vals[[cancer_type]]))
	all_cnv_genes <- unique(c(all_cnv_genes,rownames(all_cnv_vals[[cancer_type]])))

}

p_vals_mat_mut_ratio_all <- list()
p_vals_mat_mut_all <- list()
p_vals_mat_mut_immature_all <- list()

p_vals_mat_cnv_ratio_all <- list()
p_vals_mat_cnv_all <- list()
p_vals_mat_cnv_immature_all <- list()

for(cancer_type in cancer_types){
	common_samples <- intersect(colnames(all_mut_data[[cancer_type]]),colnames(all_ratio_matrix[[cancer_type]]))
	if(length(common_samples > 9)){
		miRNA_of_interest = unique(c(up_miR_lists[['ratio']],down_miR_lists[['ratio']]))
		x = t(all_mut_data[[cancer_type]][,common_samples])
		y = t(all_ratio_matrix[[cancer_type]][intersect(rownames(all_ratio_matrix[[cancer_type]]),miRNA_of_interest),common_samples])
		n <- t(!is.na(x)) %*% (!is.na(y)) # same as count.pairwise(x,y) from psych package
		r <- cor(x, y, use = "pairwise.complete.obs",method='spearman') # MUCH MUCH faster than corr.test()
		# get a list with matrices of correlation, pvalues, standard error, etc.
		result = cor2pvalue(r,n)
		p_vals_mat_mut_ratio_all[[cancer_type]] <- sign(result$r) * result$p

	}

	common_samples <- intersect(colnames(all_mut_data[[cancer_type]]),colnames(all_miRNA_data_mature[[cancer_type]]))
	if(length(common_samples > 9)){
		miRNA_of_interest = unique(c(up_miR_lists[['mature']],down_miR_lists[['mature']]))

		x = t(all_mut_data[[cancer_type]][,common_samples])
		y = t(all_miRNA_data_mature[[cancer_type]][intersect(rownames(all_miRNA_data_mature[[cancer_type]]),miRNA_of_interest),common_samples])
		n <- t(!is.na(x)) %*% (!is.na(y)) # same as count.pairwise(x,y) from psych package
		r <- cor(x, y, use = "pairwise.complete.obs",method='spearman') # MUCH MUCH faster than corr.test()
		# get a list with matrices of correlation, pvalues, standard error, etc.
		result = cor2pvalue(r,n)
		p_vals_mat_mut_all[[cancer_type]] <- sign(result$r) * result$p
	}

	common_samples <- intersect(colnames(all_mut_data[[cancer_type]]),colnames(all_miRNA_data_immature[[cancer_type]]))
	if(length(common_samples > 9)){
		miRNA_of_interest = unique(c(up_miR_lists[['immature']],down_miR_lists[['immature']]))

		x = t(all_mut_data[[cancer_type]][,common_samples])
		y = t(all_miRNA_data_immature[[cancer_type]][intersect(rownames(all_miRNA_data_immature[[cancer_type]]),miRNA_of_interest),common_samples])
		n <- t(!is.na(x)) %*% (!is.na(y)) # same as count.pairwise(x,y) from psych package
		r <- cor(x, y, use = "pairwise.complete.obs",method='spearman') # MUCH MUCH faster than corr.test()
		# get a list with matrices of correlation, pvalues, standard error, etc.
		result = cor2pvalue(r,n)
		p_vals_mat_mut_immature_all[[cancer_type]] <- sign(result$r) * result$p
	}

	common_samples <- intersect(colnames(all_cnv_vals[[cancer_type]]),colnames(all_ratio_matrix[[cancer_type]]))
	if(length(common_samples > 9)){
		miRNA_of_interest = unique(c(up_miR_lists[['ratio']],down_miR_lists[['ratio']]))

		x = t(all_cnv_vals[[cancer_type]][,common_samples])
		y = t(all_ratio_matrix[[cancer_type]][intersect(rownames(all_ratio_matrix[[cancer_type]]),miRNA_of_interest),common_samples])
		n <- t(!is.na(x)) %*% (!is.na(y)) # same as count.pairwise(x,y) from psych package
		r <- cor(x, y, use = "pairwise.complete.obs",method='spearman') # MUCH MUCH faster than corr.test()
		# get a list with matrices of correlation, pvalues, standard error, etc.
		result = cor2pvalue(r,n)
		p_vals_mat_cnv_ratio_all[[cancer_type]] <- sign(result$r) * result$p
	}

	common_samples <- intersect(colnames(all_cnv_vals[[cancer_type]]),colnames(all_miRNA_data_mature[[cancer_type]]))
	if(length(common_samples > 9)){
		miRNA_of_interest = unique(c(up_miR_lists[['mature']],down_miR_lists[['mature']]))

		x = t(all_cnv_vals[[cancer_type]][,common_samples])
		y = t(all_miRNA_data_mature[[cancer_type]][intersect(rownames(all_miRNA_data_mature[[cancer_type]]),miRNA_of_interest),common_samples])
		n <- t(!is.na(x)) %*% (!is.na(y)) # same as count.pairwise(x,y) from psych package
		r <- cor(x, y, use = "pairwise.complete.obs",method='spearman') # MUCH MUCH faster than corr.test()
		# get a list with matrices of correlation, pvalues, standard error, etc.
		result = cor2pvalue(r,n)
		p_vals_mat_cnv_all[[cancer_type]] <- sign(result$r) * result$p
	}

	common_samples <- intersect(colnames(all_cnv_vals[[cancer_type]]),colnames(all_miRNA_data_immature[[cancer_type]]))
	if(length(common_samples > 9)){
		miRNA_of_interest = unique(c(up_miR_lists[['immature']],down_miR_lists[['immature']]))

		x = t(all_cnv_vals[[cancer_type]][,common_samples])
		y = t(all_miRNA_data_immature[[cancer_type]][intersect(rownames(all_miRNA_data_immature[[cancer_type]]),miRNA_of_interest),common_samples])
		n <- t(!is.na(x)) %*% (!is.na(y)) # same as count.pairwise(x,y) from psych package
		r <- cor(x, y, use = "pairwise.complete.obs",method='spearman') # MUCH MUCH faster than corr.test()
		# get a list with matrices of correlation, pvalues, standard error, etc.
		result = cor2pvalue(r,n)
		p_vals_mat_cnv_immature_all[[cancer_type]] <- sign(result$r) * result$p
	}
}

save(file='pValueMatrices.rda',p_vals_mat_mut_ratio_all,p_vals_mat_mut_all,p_vals_mat_mut_immature_all,p_vals_mat_cnv_ratio_all,p_vals_mat_cnv_all,p_vals_mat_cnv_immature_all)

pos_assoc_genes <- list()
neg_assoc_genes <- list()
for(j in c('ratio','all','immature')){
	if(j=='all'){
		miRNA_of_interest = unique(c(up_miR_lists[['mature']],down_miR_lists[['mature']]))
	}else if(j =='ratio'){
		miRNA_of_interest = unique(c(up_miR_lists[['ratio']],down_miR_lists[['ratio']]))

	}else{
		miRNA_of_interest = unique(c(up_miR_lists[['immature']],down_miR_lists[['immature']]))

	}
for(miR_name in miRNA_of_interest){
	pos_assoc_genes[[miR_name]] <- list()
	neg_assoc_genes[[miR_name]] <- list()

	for(k in c('cnv','mut')){
		pos_assoc_genes[[miR_name]][[k]] <- list()
		neg_assoc_genes[[miR_name]][[k]] <- list()

			if(k == 'cnv' & j == 'ratio'){
				data.set = p_vals_mat_cnv_ratio_all
			}else if(k =='cnv' & j=='all'){
				data.set = p_vals_mat_cnv_all
			}else if(k == 'mut' & j == 'ratio'){
				data.set = p_vals_mat_mut_ratio_all
			}else if(k =='mut' & j=='immature'){
				data.set = p_vals_mat_mut_immature_all
			}else if(k =='cnv' & j=='immature'){
				data.set = p_vals_mat_cnv_immature_all
			}else{
				data.set = p_vals_mat_mut_all
			}

			tmp_mat <-c()
			for(cancer_type in cancer_types){
				if(miR_name %in% colnames(data.set[[cancer_type]])){
					prev_colnames <- colnames(tmp_mat)
					tmp_mat <- cbind(tmp_mat,data.set[[cancer_type]][,miR_name])
					colnames(tmp_mat) <- c(prev_colnames,cancer_type)
				}
			}

			tmp_mat <- tmp_mat[which(rowSums(!is.na(tmp_mat))>=4),]
			RP.out <- RP(1/tmp_mat,rep(1,length(colnames(tmp_mat))))
			rank_prod_tables <- topGene(RP.out,cutoff = 0.05,method="pfp",gene.names=rownames(tmp_mat))
			neg_assoc_genes[[miR_name]][[k]][[j]] <- rownames(rank_prod_tables$Table1)
			pos_assoc_genes[[miR_name]][[k]][[j]] <- rownames(rank_prod_tables$Table2)
		}
	}
}

final_mat <- matrix('',nrow=100,ncol=24)
colnames(final_mat) <- c('CNV,miR_up,mature,pos_cor',
	'CNV,miR_up,mature,neg_cor',
	'CNV,miR_up,ratio,pos_cor',
	'CNV,miR_up,ratio,neg_cor',
	'CNV,miR_up,immature,pos_cor',
	'CNV,miR_up,immature,neg_cor',
	'CNV,miR_down,mature,pos_cor',
	'CNV,miR_down,mature,neg_cor',
	'CNV,miR_down,ratio,pos_cor',
	'CNV,miR_down,ratio,neg_cor',
	'CNV,miR_down,immature,pos_cor',
	'CNV,miR_down,immature,neg_cor',
	'Mut,miR_up,mature,pos_cor',
	'Mut,miR_up,mature,neg_cor',
	'Mut,miR_up,ratio,pos_cor',
	'Mut,miR_up,ratio,neg_cor',
	'Mut,miR_up,immature,pos_cor',
	'Mut,miR_up,immature,neg_cor',
	'Mut,miR_down,mature,pos_cor',
	'Mut,miR_down,mature,neg_cor',
	'Mut,miR_down,ratio,pos_cor',
	'Mut,miR_down,ratio,neg_cor',
	'Mut,miR_down,immature,pos_cor',
	'Mut,miR_down,immature,neg_cor')


mart <- useEnsembl('ensembl', mirror = 'useast')
mart <- useDataset("hsapiens_gene_ensembl", mart=mart)

outlists <- list();
#pos, cnv, all --> up
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['mature']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['cnv']][['all']])
}

gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[1]] <- outList
outList <- outList[1:100]
final_mat[,1] <- outList

#pos, cnv, all --> down
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['mature']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['cnv']][['all']])
}

gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[2]] <- outList

outList <- outList[1:100]
final_mat[,2] <- outList

#pos, cnv, ratio --> up
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['ratio']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['cnv']][['ratio']])
}

gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[3]] <- outList

outList <- outList[1:100]
final_mat[,3] <- outList

#pos, cnv, ratio --> down
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['ratio']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['cnv']][['ratio']])
}
gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[4]] <- outList

outList <- outList[1:100]
final_mat[,4] <- outList

#pos, cnv, ratio --> up
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['immature']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['cnv']][['immature']])
}

gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[5]] <- outList

outList <- outList[1:100]
final_mat[,5] <- outList

#pos, cnv, ratio --> down
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['immature']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['cnv']][['immature']])
}
gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[6]] <- outList

outList <- outList[1:100]
final_mat[,6] <- outList

#neg, cnv, all --> up
genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['mature']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['cnv']][['all']])
}
gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[7]] <- outList

outList <- outList[1:100]
final_mat[,7] <- outList

#neg, cnv, all --> down
genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['mature']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['cnv']][['all']])
}
gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[8]] <- outList

outList <- outList[1:100]
final_mat[,8] <- outList

#neg, cnv, ratio
genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['ratio']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['cnv']][['ratio']])
}
gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[9]] <- outList

outList <- outList[1:100]
final_mat[,9] <- outList


genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['ratio']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['cnv']][['ratio']])
}
gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[10]] <- outList

outList <- outList[1:100]
final_mat[,10] <- outList

#neg, cnv, ratio
genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['immature']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['cnv']][['immature']])
}
gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[11]] <- outList

outList <- outList[1:100]
final_mat[,11] <- outList

genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['immature']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['cnv']][['immature']])
}
gList <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes_of_interest,mart= mart)
tmp <- gList[match(names(-sort(-table(genes_of_interest))),gList[,1]),2]
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[12]] <- outList

outList <- outList[1:100]
final_mat[,12] <- outList

#-----MUT---
#pos, mut, all --> up
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['mature']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['mut']][['all']])
}

tmp <- names(-sort(-table(genes_of_interest))) 
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[13]] <- outList

outList <- outList[1:100]
final_mat[,13] <- outList

#pos, mut, all --> down
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['mature']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['mut']][['all']])
}

tmp <- names(-sort(-table(genes_of_interest))) 
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[14]] <- outList

outList <- outList[1:100]
final_mat[,14] <- outList


#pos, mut, ratio --> up
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['ratio']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['mut']][['ratio']])
}

tmp <- names(-sort(-table(genes_of_interest))) 
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[15]] <- outList

outList <- outList[1:100]
final_mat[,15] <- outList

#pos, mut, ratio --> down
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['ratio']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['mut']][['ratio']])
}
tmp <- names(-sort(-table(genes_of_interest)))
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[16]] <- outList

outList <- outList[1:100]
final_mat[,16] <- outList

#pos, mut, ratio --> up
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['immature']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['mut']][['immature']])
}

tmp <- names(-sort(-table(genes_of_interest)))
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[17]] <- outList

outList <- outList[1:100]
final_mat[,17] <- outList

#pos, mut, ratio --> down
genes_of_interest <- c()
for(miR_name in names(-sort(-table(up_miR_lists[['immature']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['mut']][['immature']])
}
tmp <- names(-sort(-table(genes_of_interest))) 
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[18]] <- outList

outList <- outList[1:100]
final_mat[,18] <- outList

#neg, mut, all --> up
genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['mature']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['mut']][['all']])
}
tmp <- names(-sort(-table(genes_of_interest))) 
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[19]] <- outList

outList <- outList[1:100]
final_mat[,19] <- outList

#neg, mut, all --> down
genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['mature']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['mut']][['all']])
}
tmp <- names(-sort(-table(genes_of_interest))) 
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[20]] <- outList

outList <- outList[1:100]
final_mat[,20] <- outList

#neg, mut, ratio
genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['ratio']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['mut']][['ratio']])
}
tmp <- names(-sort(-table(genes_of_interest))) 
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[21]] <- outList

outList <- outList[1:100]
final_mat[,21] <- outList

genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['ratio']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['mut']][['ratio']])
}
tmp <- names(-sort(-table(genes_of_interest))) 
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[22]] <- outList

outList <- outList[1:100]
final_mat[,22] <- outList

#neg, mut, ratio
genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['immature']])))){
	genes_of_interest <- c(genes_of_interest, pos_assoc_genes[[miR_name]][['mut']][['immature']])
}
tmp <- names(-sort(-table(genes_of_interest))) 
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[23]] <- outList

outList <- outList[1:100]
final_mat[,23] <- outList

genes_of_interest <- c()
for(miR_name in names(-sort(-table(down_miR_lists[['immature']])))){
	genes_of_interest <- c(genes_of_interest, neg_assoc_genes[[miR_name]][['mut']][['immature']])
}
tmp <- names(-sort(-table(genes_of_interest))) 
outList <- paste0(tmp[which(!is.na(tmp))], ' (',-sort(-table(genes_of_interest))[which(!is.na(tmp))],')')
outlists[[24]] <- outList

outList <- outList[1:100]
final_mat[,24] <- outList

write.table(final_mat,file='final_matrix_associations.txt', sep='\t',quote=F, col.names=T,row.names=F)

for(k in 1:length(outlists)){
	posn <- which(grepl('AGO2',outlists[[k]]) | grepl('EIF2C2',outlists[[k]]) | grepl('CASC7',outlists[[k]]))	
	print(paste0(colnames(final_mat)[k],': ',outlists[[k]][posn],' : ',posn,'/',length(outlists[[k]])))
}

for(k in 1:length(outlists)){
	posn <- which(grepl('DICER1',outlists[[k]]) | grepl('HERNA',outlists[[k]]) | grepl('MNG1',outlists[[k]]))	
	print(paste0(colnames(final_mat)[k],': ',outlists[[k]][posn],' : ',posn,'/',length(outlists[[k]])))
}

for(k in 1:length(outlists)){
	posn <-  which(grepl('PABPC1',outlists[[k]]) | grepl('PABPC2',outlists[[k]]) | grepl('PAB1',outlists[[k]]))
	print(paste0(colnames(final_mat)[k],': ',outlists[[k]][posn],' : ',posn,'/',length(outlists[[k]])))
}

#what miRNA in each cancer type are associated to AGO2 copy number in the strongest way? we can look at 
#correlation between AGO2 copy number and the miR expression and then compare the relative values

cor_vals_ago2_cnv_mature <- list()
cor_vals_ago2_cnv_immature <- list()
cor_vals_ago2_cnv_ratio <- list()

for(cancer_type in cancer_types){
	if('27161' %in% rownames(all_cnv_vals[[cancer_type]])){	
	cor_vals_ago2_cnv_mature[[cancer_type]] <- rep(NA,times=length(rownames(all_miRNA_data_mature[[cancer_type]])))
	names(cor_vals_ago2_cnv_mature[[cancer_type]]) <- rownames(all_miRNA_data_mature[[cancer_type]])

	com_samps <- intersect(colnames(all_miRNA_data_mature[[cancer_type]]),colnames(all_cnv_vals[[cancer_type]]))
	cnv_vals <- as.numeric(all_cnv_vals[[cancer_type]]['27161',com_samps])
	for(miR_name in rownames(all_miRNA_data_mature[[cancer_type]])){
		if(sum(!is.na(as.numeric(all_miRNA_data_mature[[cancer_type]][miR_name,com_samps])))  >= 5){

			cor_vals_ago2_cnv_mature[[cancer_type]][miR_name] <- cor(as.numeric(all_miRNA_data_mature[[cancer_type]][miR_name,com_samps]),cnv_vals,use='pairwise.complete.obs',method='spearman')
		}
	}

	cor_vals_ago2_cnv_immature[[cancer_type]] <- rep(NA,times=length(rownames(all_miRNA_data_immature[[cancer_type]])))
	names(cor_vals_ago2_cnv_immature[[cancer_type]]) <- rownames(all_miRNA_data_immature[[cancer_type]])

	com_samps <- intersect(colnames(all_miRNA_data_immature[[cancer_type]]),colnames(all_cnv_vals[[cancer_type]]))
	cnv_vals <- as.numeric(all_cnv_vals[[cancer_type]]['27161',com_samps])
	for(miR_name in rownames(all_miRNA_data_immature[[cancer_type]])){
		if(sum(!is.na(as.numeric(all_miRNA_data_immature[[cancer_type]][miR_name,com_samps])))  >= 5){

			cor_vals_ago2_cnv_immature[[cancer_type]][miR_name] <- cor(as.numeric(all_miRNA_data_immature[[cancer_type]][miR_name,com_samps]),cnv_vals,use='pairwise.complete.obs',method='spearman')
		}
	}

	cor_vals_ago2_cnv_ratio[[cancer_type]] <- rep(NA,times=length(rownames(all_ratio_matrix[[cancer_type]])))
	names(cor_vals_ago2_cnv_ratio[[cancer_type]]) <- rownames(all_ratio_matrix[[cancer_type]])

	com_samps <- intersect(colnames(all_ratio_matrix[[cancer_type]]),colnames(all_cnv_vals[[cancer_type]]))
	cnv_vals <- as.numeric(all_cnv_vals[[cancer_type]]['27161',com_samps])
	for(miR_name in rownames(all_ratio_matrix[[cancer_type]])){
		if(sum(!is.na(as.numeric(all_ratio_matrix[[cancer_type]][miR_name,com_samps])))  >= 5){

			cor_vals_ago2_cnv_ratio[[cancer_type]][miR_name] <- cor(as.numeric(all_ratio_matrix[[cancer_type]][miR_name,com_samps]),cnv_vals,use='pairwise.complete.obs',method='spearman')
		}
	}

	print(cancer_type)

	}
}

for (cancer_type in cancer_types){
	print(cancer_type)
	print('///-mature-///')
	print(head(sort(cor_vals_ago2_cnv_mature[[cancer_type]], decreasing=T)))
	print('///-immature-///')
	print(head(sort(cor_vals_ago2_cnv_immature[[cancer_type]], decreasing=T)))
	print('///-ratio-///')
	print(head(sort(cor_vals_ago2_cnv_ratio[[cancer_type]], decreasing=T)))
	print('///-next-///')
}

plot(sort(cor_vals_ago2_cnv_immature[[cancer_type]]),col='blue')
points(sort(cor_vals_ago2_cnv_mature[[cancer_type]]),col='black')
points(sort(cor_vals_ago2_cnv_ratio[[cancer_type]]),col='green')

for(cancer_type in cancer_types){

	cor_vals_ago2_cnv_mature[[cancer_type]] <- cor_vals_ago2_cnv_mature[[cancer_type]][order(cor_vals_ago2_cnv_mature[[cancer_type]])]
	cor_vals_ago2_cnv_mature[[cancer_type]] <- cor_vals_ago2_cnv_mature[[cancer_type]][which(!is.na(cor_vals_ago2_cnv_mature[[cancer_type]]))]

	plot_df <- cbind.data.frame(miRNA=names(cor_vals_ago2_cnv_mature[[cancer_type]]),vals =as.numeric(cor_vals_ago2_cnv_mature[[cancer_type]]),order=1:length(cor_vals_ago2_cnv_mature[[cancer_type]]))

	plot_df$colour <- rep('black',times=length(names(cor_vals_ago2_cnv_mature[[cancer_type]])))
	plot_df$colour[match(setdiff(up_miR_lists$mature,down_miR_lists$mature),plot_df$miRNA)] <- 'red'
	plot_df$colour[match(setdiff(down_miR_lists$mature,up_miR_lists$mature),plot_df$miRNA)] <- 'blue'

	mature_miRNA_BRCA <- ggplot(plot_df, aes(x= order, y = vals)) + 
	  geom_point(color = plot_df$colour, size = 2,alpha=0.2) + ggtitle('Mature') +xlab(' ')+ ylab('Correlation coefficient')

	p1 <- mature_miRNA_BRCA + 
	  geom_label_repel(data=rbind.data.frame(subset(plot_df,order > (length(cor_vals_ago2_cnv_mature[[cancer_type]]) -5)),plot_df[which(grepl('106b',plot_df$miRNA) | grepl('-93-',plot_df$miRNA)),]),aes(label = miRNA),
	                  box.padding   = 0.25, 
	                  point.padding = 0.25,
	                  hjust=0,direction    = "y",nudge_x=-725,
	                  segment.color = 'grey50',size=3) +
	  theme_classic()

	cor_vals_ago2_cnv_immature[[cancer_type]] <- cor_vals_ago2_cnv_immature[[cancer_type]][order(cor_vals_ago2_cnv_immature[[cancer_type]])]
	cor_vals_ago2_cnv_immature[[cancer_type]] <- cor_vals_ago2_cnv_immature[[cancer_type]][which(!is.na(cor_vals_ago2_cnv_immature[[cancer_type]]))]

	plot_df <- cbind.data.frame(miRNA=names(cor_vals_ago2_cnv_immature[[cancer_type]]),vals =as.numeric(cor_vals_ago2_cnv_immature[[cancer_type]]),order=1:length(cor_vals_ago2_cnv_immature[[cancer_type]]))
	plot_df$colour <- rep('black',times=length(names(cor_vals_ago2_cnv_immature[[cancer_type]])))
	plot_df$colour[match(setdiff(up_miR_lists$immature,down_miR_lists$immature),plot_df$miRNA)] <- 'red'
	plot_df$colour[match(setdiff(down_miR_lists$immature,up_miR_lists$immature),plot_df$miRNA)] <- 'blue'
	immature_miRNA_BRCA <- ggplot(plot_df, aes(x= order, y = vals)) + 
	  geom_point(color = plot_df$colour, size = 2,alpha=0.2)+ ggtitle('Immature') + xlab('Order')+ ylab(' ') 

	p2 <- immature_miRNA_BRCA + 
	   geom_label_repel(data=rbind.data.frame(subset(plot_df,order > (length(cor_vals_ago2_cnv_immature[[cancer_type]]) -5)),plot_df[which(grepl('106b',plot_df$miRNA) | grepl('-93$',plot_df$miRNA)),]),aes(label = miRNA),
	                  box.padding   = 0.5, 
	                  point.padding = 0.5,
	                  hjust=0,direction    = "y",nudge_x=-700,
	                  segment.color = 'grey50',size=3) +
	  theme_classic()

	cor_vals_ago2_cnv_ratio[[cancer_type]] <- cor_vals_ago2_cnv_ratio[[cancer_type]][order(cor_vals_ago2_cnv_ratio[[cancer_type]])]
	cor_vals_ago2_cnv_ratio[[cancer_type]] <- cor_vals_ago2_cnv_ratio[[cancer_type]][which(!is.na(cor_vals_ago2_cnv_ratio[[cancer_type]]))]

	plot_df <- cbind.data.frame(miRNA=names(cor_vals_ago2_cnv_ratio[[cancer_type]]),vals =as.numeric(cor_vals_ago2_cnv_ratio[[cancer_type]]),order=1:length(cor_vals_ago2_cnv_ratio[[cancer_type]]))

	plot_df$colour <- rep('black',times=length(names(cor_vals_ago2_cnv_ratio[[cancer_type]])))
	plot_df$colour[match(setdiff(up_miR_lists$ratio,down_miR_lists$ratio),plot_df$miRNA)] <- 'red'
	plot_df$colour[match(setdiff(down_miR_lists$ratio,up_miR_lists$ratio),plot_df$miRNA)] <- 'blue'

	ratio_miRNA_BRCA <- ggplot(plot_df, aes(x= order, y = vals)) + 
	  geom_point(aes(color = plot_df$colour), size = 2,alpha=0.2)+ scale_color_manual('Hallmarks association',values=c("black", "blue", "red"),labels=c("Not associated", "Negatively associated", "Positively associated"))+ ggtitle('Ratio') + xlab(' ') + ylab(' ') 

	p3 <- ratio_miRNA_BRCA + 
	  geom_label_repel(data=subset(plot_df,order > (length(cor_vals_ago2_cnv_ratio[[cancer_type]]) -5)),aes(label = miRNA),
	                  box.padding   = 0.5, 
	                  point.padding = 0.5,
	                  hjust=0,direction    = "y",nudge_x=-600,
	                  segment.color = 'grey50',size=3) +
	  geom_label_repel(data=plot_df[which(grepl('106b',plot_df$miRNA) | grepl('-93-',plot_df$miRNA)),], aes(label=miRNA), box.padding   = 0.5, 
	                  point.padding = 0.5,
	                  hjust=1,direction    = "y",nudge_x=300,nudge_y=-0.15,
	                  segment.color = 'grey50',size=3) +
	    theme_classic()

	dev.new()
	 f <-  ggarrange(p1,p2,p3,nrow=1,common.legend=T,legend='bottom')
	 f
	 annotate_figure(f,
	               top = text_grob(paste0("miRNA associated with AGO2 copy number in ",cancer_type), color = "black", face = "bold", size = 14))
	  dev.copy(pdf,paste0('miRNA_ago2_ratio_plot_',cancer_type,'.pdf'),width=10,height=6)
	  dev.off()
}

#-------we are going to repeat this analysis for AGO2 expression

cor_vals_ago2_expr_mature <- list()
cor_vals_ago2_expr_immature <- list()
cor_vals_ago2_expr_ratio <- list()

for(cancer_type in cancer_types){
	if('27161' %in% rownames(all_cnv_vals[[cancer_type]])){	
	cor_vals_ago2_expr_mature[[cancer_type]] <- rep(NA,times=length(rownames(all_miRNA_data_mature[[cancer_type]])))
	names(cor_vals_ago2_expr_mature[[cancer_type]]) <- rownames(all_miRNA_data_mature[[cancer_type]])

	com_samps <- intersect(colnames(all_miRNA_data_mature[[cancer_type]]),colnames(all_mRNA_data[[cancer_type]]))
	expr_vals <- as.numeric(all_mRNA_data[[cancer_type]]['27161',com_samps])
	for(miR_name in rownames(all_miRNA_data_mature[[cancer_type]])){
		if(sum(!is.na(as.numeric(all_miRNA_data_mature[[cancer_type]][miR_name,com_samps])))  >= 5){
			cor_vals_ago2_expr_mature[[cancer_type]][miR_name] <- cor(as.numeric(all_miRNA_data_mature[[cancer_type]][miR_name,com_samps]),expr_vals,use='pairwise.complete.obs',method='spearman')
		}
	}

	cor_vals_ago2_expr_immature[[cancer_type]] <- rep(NA,times=length(rownames(all_miRNA_data_immature[[cancer_type]])))
	names(cor_vals_ago2_expr_immature[[cancer_type]]) <- rownames(all_miRNA_data_immature[[cancer_type]])

	com_samps <- intersect(colnames(all_miRNA_data_immature[[cancer_type]]),colnames(all_mRNA_data[[cancer_type]]))
	expr_vals <- as.numeric(all_mRNA_data[[cancer_type]]['27161',com_samps])
	for(miR_name in rownames(all_miRNA_data_immature[[cancer_type]])){
		if(sum(!is.na(as.numeric(all_miRNA_data_immature[[cancer_type]][miR_name,com_samps])))  >= 5){
			cor_vals_ago2_expr_immature[[cancer_type]][miR_name] <- cor(as.numeric(all_miRNA_data_immature[[cancer_type]][miR_name,com_samps]),expr_vals,use='pairwise.complete.obs',method='spearman')
		}
	}


	cor_vals_ago2_expr_ratio[[cancer_type]] <- rep(NA,times=length(rownames(all_ratio_matrix[[cancer_type]])))
	names(cor_vals_ago2_expr_ratio[[cancer_type]]) <- rownames(all_ratio_matrix[[cancer_type]])
	com_samps <- intersect(colnames(all_ratio_matrix[[cancer_type]]),colnames(all_mRNA_data[[cancer_type]]))
	expr_vals <- as.numeric(all_mRNA_data[[cancer_type]]['27161',com_samps])
	for(miR_name in rownames(all_ratio_matrix[[cancer_type]])){
		if(sum(!is.na(as.numeric(all_ratio_matrix[[cancer_type]][miR_name,com_samps])))  >= 5){
			cor_vals_ago2_expr_ratio[[cancer_type]][miR_name] <- cor(as.numeric(all_ratio_matrix[[cancer_type]][miR_name,com_samps]),expr_vals,use='pairwise.complete.obs',method='spearman')
		}
	}
	print(cancer_type)
	}
}

for (cancer_type in cancer_types){
	print(cancer_type)
	print('///-mature-///')
	print(head(sort(cor_vals_ago2_expr_mature[[cancer_type]], decreasing=T)))
	print('///-immature-///')
	print(head(sort(cor_vals_ago2_expr_immature[[cancer_type]], decreasing=T)))
	print('///-ratio-///')
	print(head(sort(cor_vals_ago2_expr_ratio[[cancer_type]], decreasing=T)))
	print('///-next-///')
}

plot(sort(cor_vals_ago2_expr_immature[[cancer_type]]),col='blue')
points(sort(cor_vals_ago2_expr_mature[[cancer_type]]),col='black')
points(sort(cor_vals_ago2_expr_ratio[[cancer_type]]),col='green')

for(cancer_type in cancer_types){

	cor_vals_ago2_expr_mature[[cancer_type]] <- cor_vals_ago2_expr_mature[[cancer_type]][order(cor_vals_ago2_expr_mature[[cancer_type]])]
	cor_vals_ago2_expr_mature[[cancer_type]] <- cor_vals_ago2_expr_mature[[cancer_type]][which(!is.na(cor_vals_ago2_expr_mature[[cancer_type]]))]

	plot_df <- cbind.data.frame(miRNA=names(cor_vals_ago2_expr_mature[[cancer_type]]),vals =as.numeric(cor_vals_ago2_expr_mature[[cancer_type]]),order=1:length(cor_vals_ago2_expr_mature[[cancer_type]]))

	plot_df$colour <- rep('black',times=length(names(cor_vals_ago2_expr_mature[[cancer_type]])))
	plot_df$colour[match(setdiff(up_miR_lists$mature,down_miR_lists$mature),plot_df$miRNA)] <- 'red'
	plot_df$colour[match(setdiff(down_miR_lists$mature,up_miR_lists$mature),plot_df$miRNA)] <- 'blue'

	mature_miRNA_BRCA <- ggplot(plot_df, aes(x= order, y = vals)) + 
	  geom_point(color = plot_df$colour, size = 2,alpha=0.2) + ggtitle('Mature') +xlab(' ')+ ylab('Correlation coefficient')

	p1 <- mature_miRNA_BRCA + 
	  geom_label_repel(data=rbind.data.frame(subset(plot_df,order > (length(cor_vals_ago2_expr_mature[[cancer_type]]) -5)),plot_df[which(grepl('106b',plot_df$miRNA) | grepl('-93-',plot_df$miRNA)),]),aes(label = miRNA),
	                  box.padding   = 0.25, 
	                  point.padding = 0.25,
	                  hjust=0,direction    = "y",nudge_x=-725,
	                  segment.color = 'grey50',size=3) +
	  theme_classic()

	cor_vals_ago2_expr_immature[[cancer_type]] <- cor_vals_ago2_expr_immature[[cancer_type]][order(cor_vals_ago2_expr_immature[[cancer_type]])]
	cor_vals_ago2_expr_immature[[cancer_type]] <- cor_vals_ago2_expr_immature[[cancer_type]][which(!is.na(cor_vals_ago2_expr_immature[[cancer_type]]))]

	plot_df <- cbind.data.frame(miRNA=names(cor_vals_ago2_expr_immature[[cancer_type]]),vals =as.numeric(cor_vals_ago2_expr_immature[[cancer_type]]),order=1:length(cor_vals_ago2_expr_immature[[cancer_type]]))
	plot_df$colour <- rep('black',times=length(names(cor_vals_ago2_expr_immature[[cancer_type]])))
	plot_df$colour[match(setdiff(up_miR_lists$immature,down_miR_lists$immature),plot_df$miRNA)] <- 'red'
	plot_df$colour[match(setdiff(down_miR_lists$immature,up_miR_lists$immature),plot_df$miRNA)] <- 'blue'
	immature_miRNA_BRCA <- ggplot(plot_df, aes(x= order, y = vals)) + 
	  geom_point(color = plot_df$colour, size = 2,alpha=0.2)+ ggtitle('Immature') + xlab('Order')+ ylab(' ') 

	p2 <- immature_miRNA_BRCA + 
	   geom_label_repel(data=rbind.data.frame(subset(plot_df,order > (length(cor_vals_ago2_expr_immature[[cancer_type]]) -5)),plot_df[which(grepl('106b',plot_df$miRNA) | grepl('-93$',plot_df$miRNA)),]),aes(label = miRNA),
	                  box.padding   = 0.5, 
	                  point.padding = 0.5,
	                  hjust=0,direction    = "y",nudge_x=-700,
	                  segment.color = 'grey50',size=3) +
	  theme_classic()

	cor_vals_ago2_expr_ratio[[cancer_type]] <- cor_vals_ago2_expr_ratio[[cancer_type]][order(cor_vals_ago2_expr_ratio[[cancer_type]])]
	cor_vals_ago2_expr_ratio[[cancer_type]] <- cor_vals_ago2_expr_ratio[[cancer_type]][which(!is.na(cor_vals_ago2_expr_ratio[[cancer_type]]))]

	plot_df <- cbind.data.frame(miRNA=names(cor_vals_ago2_expr_ratio[[cancer_type]]),vals =as.numeric(cor_vals_ago2_expr_ratio[[cancer_type]]),order=1:length(cor_vals_ago2_expr_ratio[[cancer_type]]))

	plot_df$colour <- rep('black',times=length(names(cor_vals_ago2_expr_ratio[[cancer_type]])))
	plot_df$colour[match(setdiff(up_miR_lists$ratio,down_miR_lists$ratio),plot_df$miRNA)] <- 'red'
	plot_df$colour[match(setdiff(down_miR_lists$ratio,up_miR_lists$ratio),plot_df$miRNA)] <- 'blue'

	ratio_miRNA_BRCA <- ggplot(plot_df, aes(x= order, y = vals)) + 
	  geom_point(aes(color = plot_df$colour), size = 2,alpha=0.2)+ scale_color_manual('Hallmarks association',values=c("black", "blue", "red"),labels=c("Not associated", "Negatively associated", "Positively associated"))+ ggtitle('Ratio') + xlab(' ') + ylab(' ') 

	p3 <- ratio_miRNA_BRCA + 
	  geom_label_repel(data=rbind.data.frame(subset(plot_df,order > (length(cor_vals_ago2_expr_ratio[[cancer_type]]) -5)),plot_df[which(plot_df$miRNA == 'hsa-miR-106b-3p'),]),aes(label = miRNA),
	                  box.padding   = 0.5, 
	                  point.padding = 0.5,
	                  hjust=0,direction    = "y",nudge_x=-600,
	                  segment.color = 'grey50',size=3) +
	  geom_label_repel(data=plot_df[which(grepl('106b-5p',plot_df$miRNA) | grepl('-93-',plot_df$miRNA)),], aes(label=miRNA), box.padding   = 0.5, 
	                  point.padding = 0.5,
	                  hjust=1,direction    = "y",nudge_x=300,nudge_y=-0.1,
	                  segment.color = 'grey50',size=3) +
	  theme_classic()
	
	dev.new()
	 f <-  ggarrange(p1,p2,p3,nrow=1,common.legend=T,legend='bottom')
	 f
	 annotate_figure(f,
	               top = text_grob(paste0("miRNA associated with AGO2 expression in ",cancer_type), color = "black", face = "bold", size = 14))
	  dev.copy(pdf,paste0('miRNA_ago2_expr_plot_',cancer_type,'.pdf'),width=10,height=6)
	  dev.off()

}

#------helper functions------

get_coefficients_pre_filter <- function(cancer_type,scores, miRNA_data){
	
	#take only common subset of miRNA and scores
	common_colNames <- intersect(colnames(miRNA_data),names(scores))

	#take just the common pieces
	miRNA_data <- miRNA_data[,common_colNames]
	scores <- scores[common_colNames]

	#z-transform the scores
	scores <- as.numeric(scores) - mean(as.numeric(scores))/sd(as.numeric(scores))
	print(sum(is.na(scores)))

	#expression filter for miRNA
	expression_threshold <- 0.80 # means that at least 10% of samples must have a nonzero value of the mRNA
	miRNA_data <-miRNA_data[which((rowSums(miRNA_data==0)) < ((1-expression_threshold) * length(colnames(miRNA_data)))),]

	#remove NA values from miRNA data
	miRNA_data <- as.matrix(log2(miRNA_data))
	miRNA_data[!(is.finite(miRNA_data))] <- NA
	miRNA_data <-miRNA_data[which((rowSums(miRNA_data==0)) < ((1-expression_threshold) * length(colnames(miRNA_data)))),]

	#z-transform the miRNA data
	for (j in 1:length(rownames(miRNA_data))){
		miRNA_data[j,] <- (as.numeric(miRNA_data[j,]) - mean(as.numeric(miRNA_data[j,])))/sd(as.numeric(miRNA_data[j,]))
	}

	print(paste0("mirna " , sum(is.na(miRNA_data))))

	#first we need to subset the data into folds
	new_df <- na.omit(t(rbind(scores,miRNA_data)))
	colnames(new_df) <- c('scores',rownames(miRNA_data))

	folds <- 10
	nrows_combined_df <- 1:dim(new_df)[1]
	best_overall_error <- 99999999
	for (i in 0:(folds-1)){
		new_df_subset <- as.data.frame(new_df[!(nrows_combined_df%%folds==i),]) #takes out the 1/nth row of the data set
		#train the univaraite model
		#put these as inputs to the penalized model

		linear_models_miRNA <- matrix(,nrow=length(rownames(miRNA_data)),ncol=1)
		row.names(linear_models_miRNA) <- rownames(miRNA_data)
		for (j in 1:length(rownames(miRNA_data))){
			univariate_data <- as.data.frame(cbind(new_df_subset[,1],new_df_subset[,(j+1)]))
			colnames(univariate_data) <- c('sig_score','miRNA')
			univariate_model <- lm(formula = sig_score ~ miRNA,data = univariate_data)
	 		linear_models_miRNA[j] <- (summary(univariate_model)$coefficients)[2,4]
		}

		#significant miRNAs are those w p < 0.2:
		significant_miRNAs <- rownames(linear_models_miRNA)[which(linear_models_miRNA < 0.2 & !is.nan(linear_models_miRNA))]

		#penalised linear regression
		lambda_2_values <- c(0, 0.01, 0.1,1,10,100)
		max_likelihood <- -9999999999
		for (lambda2_val in lambda_2_values){

			cross_val_model <- optL1(response = new_df_subset[,1],penalized = new_df_subset[,significant_miRNAs], lambda2 = lambda2_val,data=as.data.frame(new_df_subset),model="linear",fold=10,trace=F)#,trace=F,maxiter=1000,tol=.Machine$double.eps^0.23)
			if ((cross_val_model$fullfit)@loglik > max_likelihood){
				best_model <<- cross_val_model
				best_lambda <- lambda2_val
			}
		}

		#now that we know the best model, let's test it on the other 1/n of the data, and record the error
		unused_df <- as.data.frame(new_df[(nrows_combined_df%%folds==i),])
		current_predictions <- predict(best_model$fullfit, penalized=unused_df[,significant_miRNAs],data=unused_df)

		cur_error <- norm((as.numeric(unused_df[,1]) - as.numeric(current_predictions)),type="2")
		if (cur_error < best_overall_error){
			best_overall_error <- cur_error
			best_overall_model <- best_model
			best_overall_lambda <- best_lambda
		}
	}

	miRNA_names_reported <- intersect(names(coef(best_overall_model$fullfit)), rownames(miRNA_data))

	#return the coefficients
	coef(best_overall_model$fullfit)[miRNA_names_reported]
}

make_rank_prod_matrix <- function(coef_matrix){
	final_answer_list <- list();
	ranks_matrix_pos <- matrix(NA,nrow=dim(coef_matrix)[1],ncol=dim(coef_matrix)[2])
	ranks_matrix_neg <- matrix(NA,nrow=dim(coef_matrix)[1],ncol=dim(coef_matrix)[2])
	row.names(ranks_matrix_pos) <- rownames(coef_matrix)
	row.names(ranks_matrix_neg) <- rownames(coef_matrix)

	for(i in 1:dim(coef_matrix)[2]){
		ranks_matrix_pos[,i] <- rank(coef_matrix[,i])/length(coef_matrix[,i]) 
		ranks_matrix_neg[,i] <- rank(-coef_matrix[,i])/length(coef_matrix[,i])
	}

	rank_prod_pos <- apply(ranks_matrix_pos,1,function(x) prod(x,na.rm=T) ^ (1/sum(!(is.na(x)))))
	rank_prod_neg <- apply(ranks_matrix_neg,1,function(x) prod(x,na.rm=T) ^ (1/sum(!(is.na(x)))))

	full_rank_matrix_pos <- cbind(ranks_matrix_pos,rank_prod_pos)
	row.names(full_rank_matrix_pos) <- rownames(ranks_matrix_pos)

	full_rank_matrix_pos <- full_rank_matrix_pos[order(full_rank_matrix_pos[,'rank_prod_pos']),]
	print(full_rank_matrix_pos[1:20,])

	full_rank_matrix_neg <- cbind(ranks_matrix_neg,rank_prod_neg)
	row.names(full_rank_matrix_neg) <- rownames(ranks_matrix_neg)
	full_rank_matrix_neg <- full_rank_matrix_neg[order(full_rank_matrix_neg[,'rank_prod_neg']),]

	print(full_rank_matrix_neg[1:20,])

 	#bootstrapping for pvalues
	N_repeats_bootstrap <- 10
	rank_prod_bootstrapped <- matrix(NA, nrow=length(row.names(full_rank_matrix_pos)),ncol=N_repeats_bootstrap)
	row.names(rank_prod_bootstrapped) <- rownames(full_rank_matrix_pos)

	final_mat_boot <- full_rank_matrix_pos[,1:length(cancer_types)]
	row.names(final_mat_boot) <- rownames(full_rank_matrix_pos)

	for (i in 1:N_repeats_bootstrap){
		for (col in 1:length(cancer_types)){
			resampled <- sample.int(length(full_rank_matrix_pos[,col]),length(full_rank_matrix_pos[,col]),replace=F)
			final_mat_boot[,col] <- as.numeric(full_rank_matrix_pos[resampled,col])
		}
		rank_prod <- apply(final_mat_boot,1,function(x) (prod(na.omit(x))^(1/length(na.omit(x)))))
		rank_prod_bootstrapped[rownames(final_mat_boot),i] <- rank_prod

	}
	p_values <- sapply(1:dim(rank_prod_bootstrapped)[1],function(x) sum(rank_prod_bootstrapped[x,] < full_rank_matrix_pos[x,(length(cancer_types)+1)]) / length(rank_prod_bootstrapped[x,]) )
	final_mat_pos <- cbind(full_rank_matrix_pos,p_values)
	final_mat_pos <- final_mat_pos[order(final_mat_pos[,'p_values']),]
	final_answer_list[['pos']] <- final_mat_pos

	N_repeats_bootstrap <- 10
	rank_prod_bootstrapped <- matrix(NA, nrow=length(row.names(full_rank_matrix_neg)),ncol=N_repeats_bootstrap)
	row.names(rank_prod_bootstrapped) <- rownames(full_rank_matrix_neg)

	final_mat_boot <- full_rank_matrix_neg[,1:length(cancer_types)]
	row.names(final_mat_boot) <- rownames(full_rank_matrix_neg)

	for (i in 1:N_repeats_bootstrap){
		for (col in 1:length(cancer_types)){
			resampled <- sample.int(length(full_rank_matrix_neg[,col]),length(full_rank_matrix_neg[,col]),replace=F)
			final_mat_boot[,col] <- as.numeric(full_rank_matrix_neg[resampled,col])
		}
		rank_prod <- apply(final_mat_boot,1,function(x) (prod(na.omit(x))^(1/length(na.omit(x)))))
		rank_prod_bootstrapped[rownames(final_mat_boot),i] <- rank_prod
	}
	p_values <- sapply(1:dim(rank_prod_bootstrapped)[1],function(x) sum(rank_prod_bootstrapped[x,] < full_rank_matrix_neg[x,(length(cancer_types)+1)]) / length(rank_prod_bootstrapped[x,]) )
	final_mat_neg <- cbind(full_rank_matrix_neg,p_values)
	final_mat_neg <- final_mat_neg[order(final_mat_neg[,'p_values']),]
	final_answer_list[['neg']] <- final_mat_neg

	final_answer_list
}


cor2pvalue = function(r, n) {
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  se <- sqrt((1-r*r)/(n-2))
  out <- list(r, n, t, p, se)
  names(out) <- c("r", "n", "t", "p", "se")
  return(out)
}

