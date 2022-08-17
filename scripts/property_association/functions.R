#---------------------------------------------------------------------------------------------------------------------------

##### FUNCTIONS FOR MICROBIOME_RECONSTRUCTION - DIVERSITY #####

#---------------------------------------------------------------------------------------------------------------------------

##### LIBRARIES #####

#---------------------------------------------------------------------------------------------------------------------------

# library(vegan)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
library(ggsci)
# library(RColorBrewer)

#---------------------------------------------------------------------------------------------------------------------------

# OTHER FUNCTIONS

#---------------------------------------------------------------------------------------------------------------------------

source("../general_functions.R")

#---------------------------------------------------------------------------------------------------------------------------

# Find outliers from distance matrix

#---------------------------------------------------------------------------------------------------------------------------

find_outliers = function(dist_mat, thr){
	
	max_dist_mat=sapply(1:ncol(dist_mat), function(x){
		pos=which(dist_mat[x,]==max(dist_mat[x,]))
		return(colnames(dist_mat)[pos])
		}
	)
	max_dist_mat=unlist(max_dist_mat)
	if(length(unique(max_dist_mat))==1){
		out=c(unique(max_dist_mat))
		return(out)
	}
	if(any(table(max_dist_mat)>dim(dist_mat)[1]*thr)){
		out=c(names(table(max_dist_mat)[table(max_dist_mat)>dim(dist_mat)[1]*thr]))
		return(out)
	} else {
		return(NA)
	}

}

#---------------------------------------------------------------------------------------------------------------------------

# PCA

#---------------------------------------------------------------------------------------------------------------------------

PCA=function(taxa_values, num_features=1000, scaling=TRUE, out_removal=NA){
	
	### Checks

		# taxa_values must be a data frame
		if(class(taxa_values)!="data.frame"){
			stop("taxa_values must be a data frame")
		}
				
	### PCA function

		# Remove not values columns
		pos=c()
		for(i in 1:ncol(taxa_values)){
			types=class(taxa_values[,i])
			if(types=="integer" | types=="numeric" | types=="double"){
				pos[[length(pos)+1]]=TRUE
			} else {
				pos[[length(pos)+1]]=FALSE
			}
		}
		tab=taxa_values[,pos]
	
		# Convert NA to zeros
		tab[is.na(tab)]=0

		# Remove samples with all zeros (can happen if I select the bacteria to be used)
		p=apply(tab, 1, function(x){
				all(x==0)
			}
		)
		tab=tab[!p,]

		# Feature selection by standard deviation 
		std_dev=apply(tab, 2, sd)
		if(num_features=="all" | num_features>ncol(tab)){
			num_features=ncol(tab)
		}
		most_variable=sort(unlist(std_dev), decreasing=TRUE)[1:num_features]
		most_variable_features=unique(names(most_variable))
		sel_col=sapply(colnames(tab), function(x){
				if(any(x==most_variable_features)){
					return(TRUE)
				} else {
					return(FALSE)
				}
			}
		)
		mat_pca=tab[,sel_col]

		# to use scale.=TRUE, remove columns with all zeros and costant
		if(scaling==TRUE){
			pos_col=sapply(1:ncol(mat_pca), function(y){
					if(all(mat_pca[,y]==0)){
						return(FALSE)
					} else if(var(mat_pca[,y])==0){
						return(FALSE)
					} else {
						return(TRUE)
					}
				}
			)
			mat_pca=mat_pca[,pos_col]
		}

		# PCA table
		if(scaling==TRUE){
			pca_data=prcomp(mat_pca, scale.=TRUE)
		} else {
			pca_data=prcomp(mat_pca, scale.=FALSE)	
		}
	
	### Removal of outliers

		if(!is.na(out_removal)){
			# Detect outliers
			dist_tab=as.matrix(dist(pca_data$x, method="euclidean"))
			
			all_outliers=c()
			tmp_dist=dist_tab
			repeat{
				out=find_outliers(dist_mat=tmp_dist, thr=out_removal)
				if(is.na(out)){
					break
				} else {
					all_outliers=c(all_outliers, out)
					p=sapply(rownames(tmp_dist), function(x){
							if(any(x==out)){
								return(FALSE)
							} else {
								return(TRUE)
							}
						}
					)
					tmp_dist=tmp_dist[p, p]
				}
			}

			# Remove outliers and rerun the function
			if(length(all_outliers)!=0){
				pos=sapply(rownames(taxa_values), function(x){
						if(any(x==all_outliers)){
							return(FALSE)
						} else {
							return(TRUE)
						}
					}
				)
				new_taxa_values=taxa_values[pos,]
				new_res=PCA(taxa_values=new_taxa_values, 
						num_features=num_features, 
						out_removal=out_removal,
						scaling=scaling
						)

				return(new_res)
			}
		}

	### Return 
	
		return(pca_data)
	
}

#---------------------------------------------------------------------------------------------------------------------------

# PCA plots

#---------------------------------------------------------------------------------------------------------------------------

PCA_plots=function(pca_data, metadata, cat_properties=NA, cont_properties=NA, PCs, val_removed, pal="default"){
	# Features
		features=c(cat_properties, cont_properties)
		features=features[!is.na(features)]
		
	# Palette
		if(all(names(pal)=="default")){
			pal=rep("default", length(features))
		}
		names(pal)=features

	# Match samples metadata and pca_data
		metadata=match_metadata(metadata=metadata, taxa_values=pca_data$x, colname_metadata="file_id")

	# PCA values
		plotting_table=data.frame(pca_data$x, metadata[,features])
		# Control colname (if length(features)==1)
		if(any(colnames(plotting_table)=="metadata...features.")){
			colnames(plotting_table)[colnames(plotting_table)=="metadata...features."]=features
		}
		# Force to numeric or categorical values
		if(!all(is.na(cat_properties))){
			if(length(cat_properties)>1){
				plotting_table[,cat_properties]=lapply(plotting_table[,cat_properties], as.factor)
			} else {
				plotting_table[,cat_properties]=as.factor(plotting_table[,cat_properties])
			}
		}
		if(!all(is.na(cont_properties))){
			for(i in cont_properties){
				plotting_table[,i]=as.numeric(plotting_table[,i])
			}
		}
							
	# Scree plot
		# Percentage of variability explained
			eigs=pca_data$sdev^2
			perc=eigs/sum(eigs)*100

		# plot
			df=data.frame(PCs=paste("PC", 1:50, sep=''), perc=perc[1:50])
			fact_order=paste("PC", 1:50, sep='')
			df$PCs=factor(df$PCs, levels=fact_order)
			df$cum_perc=cumsum(df$perc)
					
			scree_plot1=ggplot(df, aes(x=PCs, y=perc))+
							geom_bar(fill = "#0073C2FF", stat="identity")+
							theme_light()+
							theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10), 
							axis.title.x=element_blank())
			scree_plot2=ggplot(df, aes(x=PCs, y=cum_perc))+
							geom_bar(fill = "#0073C2FF", stat="identity")+
							theme_light()+
							theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10), 
							axis.title.x=element_blank())

	# Combinations of PCs
		comb=combn(PCs, 2)
		n_PCs=max(comb)

	# Scatter plot 2D
		if(any(!is.na(cat_properties))){
			cat_plots=lapply(cat_properties, function(x){
					# Legend
						# Count observations for each group (to be added to legend)
							tmp=plotting_table[,c("PC1","PC2",x)]
							colnames(tmp)=c("x_values", "y_values", "color_group")
							num_observations=table(tmp[,ncol(tmp)])
										
						# Remove labels for which we have zero samples
							pos=which(num_observations!=0 & num_observations!="logical")
							num_observations=num_observations[pos]
							new_labels=paste(names(num_observations), "(", num_observations, ")", sep='')
										
						# Extract the legend and convert it to a ggplot
							PCA_tmp=ggplot(tmp)+
									geom_point(size=2,aes(x=color_group, y=x_values,
													col=factor(color_group, labels=new_labels)))+
									theme_light()+
									labs(color=x)+
									theme(legend.title=element_text(size=20), legend.text=element_text(size=20))
							if(all(pal[[x]]!="default")){ 
								PCA_tmp=PCA_tmp+scale_color_manual(values=pal[[x]])
							}			
							leg=get_legend(PCA_tmp)
							leg_plot=as_ggplot(leg)
					# Scatter
						PCA_scatter_plots=apply(comb, 2, function(y){
								PCA_plot=ggplot(plotting_table)+
											geom_point(size=3,alpha=0.5, aes_string(x=paste("PC", y[[1]], sep=''), y=paste("PC", y[[2]], sep=''), col=x))+
											xlab(paste("PC", y[[1]], " (", signif(perc[[y[[1]]]], 3), "%)", sep=''))+
											ylab(paste("PC", y[[2]], " (", signif(perc[[y[[2]]]], 3), "%)", sep=''))+
											theme_light()+
											theme(legend.position="none")+
											labs(color=x)
								if(all(pal[[x]]!="default")){ 
									PCA_plot=PCA_plot+scale_color_manual(values=pal[[x]])
								}
								return(PCA_plot)
							}
						)
					# Boxplots
						pos=sapply(plotting_table[,x], function(y){
								!any(y==val_removed[[x]])
							}
						)
						plotting_table2=plotting_table[pos,]
						if(length(unique(plotting_table2[,x]))>2){
							stat_test="kruskal.test"
						} else {
							stat_test="wilcox.test"
						}
						PCA_boxplots=lapply(PCs, function(y){
								PCA_boxplot=ggplot(plotting_table2, aes_string(x=x, y=paste("PC", y, sep=''), color=x))+
													geom_boxplot(notch=TRUE, outlier.size=-1)+
													geom_jitter(position=position_jitter(0.2), size=2)+
													ylab(paste("PC", y, sep=''))+
													xlab(new_labels)+
													theme_light()+
													theme(legend.position="none", 
															axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank())+
												stat_compare_means(method=stat_test, #size=7,
																	label.y.npc="bottom", label.x.npc="center")
								if(all(pal[[x]]!="default")){ 
									PCA_boxplot=PCA_boxplot+scale_color_manual(values=pal[[x]])
								}											
								return(PCA_boxplot)
							}
						)
					# Return
					return(list(legend=leg_plot, scatters=PCA_scatter_plots, by_PCs=PCA_boxplots))
				}
			)
			names(cat_plots)=cat_properties
		} else {
			cat_plots=NA
		}
		if(!is.na(cont_properties)){
			cont_plots=lapply(cont_properties, function(x){
					# Legend
						# Count observations for each group (to be added to legend)
							tmp=plotting_table[,c("PC1","PC2",x)]
							colnames(tmp)=c("x_values", "y_values", "color_group")
										
						# Extract the legend and convert it to a ggplot
							PCA_tmp=ggplot(tmp)+
									geom_point(size=2,alpha=0.5, aes(x=x_values, y=y_values, col=color_group))+
										labs(color=x)+
										theme_light()
							if(pal[[x]]=="default"){
								PCA_tmp=PCA_tmp+scale_color_gradient(low="blue", high="red")
							} else {
								PCA_tmp=PCA_tmp+scale_color_gradient(low=pal[[x]][1], high=pal[[x]][2])
							}
							leg=get_legend(PCA_tmp)
							leg_plot=as_ggplot(leg)
					# Scatter
						PCA_scatter_plots=apply(comb, 2, function(y){
								PCA_plot=ggplot(plotting_table)+
											geom_point(size=3,alpha=0.5, aes_string(x=paste("PC", y[[1]], sep=''), y=paste("PC", y[[2]], sep=''), col=x))+
											xlab(paste("PC", y[[1]], sep=''))+
											ylab(paste("PC", y[[2]], sep=''))+
											theme_light()+
											theme(legend.position="none")+
											labs(color=x)
								if(pal[[x]]=="default"){
									PCA_plot=PCA_plot+scale_color_gradient(low="blue", high="red")
								} else {
									PCA_plot=PCA_plot+scale_color_gradient(low=pal[[x]][1], high=pal[[x]][2])
								}
							}
						)
					# Correlations
						PCA_corplots=lapply(PCs, function(y){
								my_data=cbind(as.data.frame(plotting_table[,x]), as.data.frame(plotting_table[,paste("PC", y, sep='')]))
								colnames(my_data)=c(x, paste("PC", y, sep=''))
								cor_plot=ggscatter(my_data, x=x, y=paste("PC", y, sep=''), 
															add="reg.line", conf.int=FALSE,
															cor.coef=TRUE, 
															cor.coeff.args=list(method="pearson", aes(color="blue"), 
																						show.legend=FALSE, size=7))+
											stat_cor(method="spearman", cor.coef.name="rho", 
														label.x.npc="left", label.y.npc="bottom", 
														aes(color="red"), show.legend=FALSE, size=7)+
											font("xy.text", size=25)+
											font("xy.title", size=30) 
								return(cor_plot)
							}
						)
					# Return
					return(list(legend=leg_plot, scatters=PCA_scatter_plots, by_PCs=PCA_corplots))
				}
			)
			names(cont_plots)=cont_properties
		} else {
			cont_plots=NA
		}
	# Return
	return(list(scree=scree_plot1, cum_scree=scree_plot2, cat_plot=cat_plots, cont_plot=cont_plots))
	
}

#---------------------------------------------------------------------------------------------------------------------------

# PCA plots 2D, combination of 2 features

#---------------------------------------------------------------------------------------------------------------------------

PCA_plots_2properties=function(pca_data, metadata, PCs, prop1=NA, prop2=NA, pal="default"){

	# Match samples metadata and pca_data
		metadata=match_metadata(metadata=metadata, taxa_values=pca_data$x, 
									colname_metadata="file_id")

	# PCA values
		features=c(prop1, prop2)
		names(features)=c(names(prop1), names(prop2))
		plotting_table=data.frame(pca_data$x, metadata[,unique(features)])
		# Force to numeric or categorical values
		for(i in 1:length(prop1)){
			if(names(prop1)[[i]]=="cat"){
				plotting_table[,prop1[[i]]]=as.factor(plotting_table[,prop1[[i]]])
			} else if(names(prop1)[[i]]=="cont"){
				plotting_table[,prop1[[i]]]=as.numeric(plotting_table[,prop1[[i]]])
			}
		}
		if(length(unique(prop2))==1){
			plotting_table[,unique(prop2)]=as.factor(plotting_table[,unique(prop2)])
		} else {
			plotting_table[,unique(prop2)]=lapply(plotting_table[,unique(prop2)], as.factor)
		}
							
	# Percentage of variability explained
		eigs=pca_data$sdev^2
		perc=eigs/sum(eigs)*100

	# Combinations of PCs
		comb=combn(PCs, 2)
		n_PCs=max(comb)

	# Scatter plot 2D
		scatter_plots=lapply(1:length(prop1), function(x){
				if(names(prop1)[[x]]=="cat"){
					# Legend
						# Count observations for each group (to be added to legend)
							tmp=plotting_table[,c("PC1", "PC2", prop1[[x]], prop2[[x]])]
							colnames(tmp)=c("x_values", "y_values", "color_group1", "color_group2")
							num_observations1=table(tmp[,"color_group1"])
							num_observations2=table(tmp[,"color_group2"])
										
						# Remove labels for which we have zero samples
							pos=which(num_observations1!=0 & num_observations1!="logical")
							num_observations1=num_observations1[pos]
							new_labels1=paste(names(num_observations1), "(", num_observations1, ")", sep='')
							pos=which(num_observations2!=0 & num_observations2!="logical")
							num_observations2=num_observations2[pos]
							new_labels2=paste(names(num_observations2), "(", num_observations2, ")", sep='')
										
						# Extract the legend and convert it to a ggplot
							PCA_tmp=ggplot(tmp)+
									geom_point(size=2,aes(x=color_group1, y=x_values, shape=color_group2,
													col=factor(color_group1, labels=new_labels1)))+
									scale_shape_manual(values=c(19, 17))+
									labs(color=prop1[[x]], shape=prop2[[x]])+
									theme_light()
							if(pal[[x]]=="brewer"){ 
								PCA_tmp=PCA_tmp+scale_colour_brewer(palette="brew")
							} else if(pal[[x]]=="default"){ 
								PCA_tmp=PCA_tmp
							} else if(pal[[x]]!="default"){ 
								PCA_tmp=PCA_tmp+scale_color_manual(values=pal[[x]])
							}
							leg=get_legend(PCA_tmp)
							leg_plot=as_ggplot(leg)
					# Scatter
						PCA_scatter_plots=apply(comb, 2, function(y){
								PCA_plot=ggplot(plotting_table)+
											geom_point(size=2.5, alpha=0.7, 
											aes_string(shape=prop2[[x]],
												x=paste("PC", y[[1]], sep=''), 
												y=paste("PC", y[[2]], sep=''), col=prop1[[x]]))+
											scale_shape_manual(values=c(19, 17))+
											xlab(paste("PC", y[[1]], " (", signif(perc[[y[[1]]]], 3), "%)", sep=''))+
											ylab(paste("PC", y[[2]], " (", signif(perc[[y[[2]]]], 3), "%)", sep=''))+
											theme_light()+
											theme(legend.position="none")+
											labs(color=prop1[[x]])
								if(pal[[x]]=="brewer"){ 
									PCA_plot=PCA_plot+scale_colour_brewer(palette="brew")
								} else if(pal[[x]]=="default"){ 
									PCA_plot=PCA_plot
								} else if(pal[[x]]!="default"){ 
									PCA_plot=PCA_plot+scale_color_manual(values=pal[[x]])
								}
								return(PCA_plot)
							}
						)
				} else if(names(prop1)[[x]]=="cont"){
					# Legend
						# Count observations for each group (to be added to legend)
							tmp=plotting_table[,c("PC1", "PC2", prop1[[x]], prop2[[x]])]
							colnames(tmp)=c("x_values", "y_values", "color_group1", "color_group2")
						
						# Remove labels for which we have zero samples
							pos=which(num_observations2!=0 & num_observations2!="logical")
							num_observations2=num_observations2[pos]
							new_labels2=paste(names(num_observations2), "(", num_observations2, ")", sep='')
										
						# Extract the legend and convert it to a ggplot
							PCA_tmp=ggplot(tmp)+
									geom_point(size=2.5, alpha=0.5, aes(x=x_values, y=y_values, 
												col=color_group1, shape=color_group2))+
										scale_shape_manual(values=c(19, 17))+
										labs(color=prop1[[x]], shape=prop2[[x]])+
										theme_light()
							if(pal[[x]]=="default"){
								PCA_tmp=PCA_tmp+scale_color_gradient(low="blue", high="red")
							} else {
								PCA_tmp=PCA_tmp+scale_color_gradient(low=pal[[x]][1], high=pal[[x]][2])
							}
							leg=get_legend(PCA_tmp)
							leg_plot=as_ggplot(leg)
					# Scatter
						PCA_scatter_plots=apply(comb, 2, function(y){
								PCA_plot=ggplot(plotting_table)+
											geom_point(size=2.5, alpha=0.5, 
												aes_string(shape=prop2[[x]],
													x=paste("PC", y[[1]], sep=''), 
													y=paste("PC", y[[2]], sep=''), col=prop1[[x]]))+
											scale_shape_manual(values=c(19, 17))+
											xlab(paste("PC", y[[1]], sep=''))+
											ylab(paste("PC", y[[2]], sep=''))+
											theme(legend.position="none")+
											theme_light()+
											labs(color=prop1[[x]])
								if(pal[[x]]=="default"){
									PCA_plot=PCA_plot+scale_color_gradient(low="blue", high="red")
								} else {
									PCA_plot=PCA_plot+scale_color_gradient(low=pal[[x]][1], high=pal[[x]][2])
								}
								return(PCA_plot)
							}
						)
				}

				# Return
				return(list(legend=leg_plot, scatter=PCA_scatter_plots))
			}
		)
		
	# Return
	return(scatter_plots)
	
}

#---------------------------------------------------------------------------------------------------------------------------

# Test associations and correlations of PCs to features

#---------------------------------------------------------------------------------------------------------------------------

test_dimreduction_to_features=function(dimred_data, metadata, cat_features=NA, cont_features=NA, val_removed){

	# Match samples metadata and dimred_data
		metadata=match_metadata(metadata=metadata, taxa_values=dimred_data, 
									colname_metadata="file_id")

	# Table of PCA and features
		features=c(cat_features, cont_features)
		features=features[!is.na(features)]
		plotting_table=data.frame(dimred_data, metadata[,features])
		# Control colname (if length(features)==1)
		if(any(colnames(plotting_table)=="metadata...features.")){
			colnames(plotting_table)[colnames(plotting_table)=="metadata...features."]=features
		}
	
	# Force to numeric or categorical values
		if(!any(is.na(cat_features))){
			for(i in cat_features){
				plotting_table[,i]=as.character(plotting_table[,i])
				if(any(is.na(plotting_table[,i]))){
					plotting_table[is.na(plotting_table[,i]),i]="unknown"
				}
			}
		}
		if(!any(is.na(cont_features))){
			for(i in cont_features){
				plotting_table[,i]=as.numeric(plotting_table[,i])
			}
		}
	
	# Wilcoxon/Kruskal test
		if(any(!is.na(cat_features))){
			cat_test=lapply(cat_features, function(feat){
					pos=sapply(plotting_table[,feat], function(x){
							!any(x==val_removed[[feat]])
						}
					)
					test_table=plotting_table[pos,]
					interesting_feat=unique(plotting_table[,feat])
					pos=sapply(interesting_feat, function(x){
							!any(x==val_removed[[feat]])
						}
					)
					interesting_feat=interesting_feat[pos]
					# Apply test (wilcox or kruskal)
					if(length(unique(test_table[,feat]))==2){
						p_val_vect=sapply(colnames(dimred_data), function(z){
								wilcox.test(formula(paste(z, " ~ ", feat, sep='')), data=test_table)$p.value
							}
						)
						test_name=c("wilc_p", "wilc_q")
					} else {
						p_val_vect=sapply(colnames(dimred_data), function(z){
								kruskal.test(formula(paste(z, " ~ ", feat, sep='')), data=test_table)$p.value
							}
						)
						test_name=c("krusk_p", "krusk_q")
					}
					# Correct with fdr
					q_val_vect=p.adjust(p_val_vect, method="fdr")
					# make table
					test_tab=matrix(c(p_val_vect, q_val_vect), nrow=ncol(dimred_data), byrow=FALSE)
					rownames(test_tab)=colnames(dimred_data)
					colnames(test_tab)=test_name
					return(test_tab)
				}
			)
			names(cat_test)=cat_features
		} else {
			cat_test=NA
		}
				
	# Correlation tests
		if(!is.na(cont_features)){
			cont_test=lapply(cont_features, function(feat){
					test_table=plotting_table
					# Pearson test
					pearson_test=lapply(colnames(dimred_data), function(z){
							r=cor.test(test_table[,z], test_table[,feat], method="spearman")$estimate
							p=cor.test(test_table[,z], test_table[,feat], method="spearman")$p.value
							return(list(r=r, p=p))
						}
					)
					pearson_r=sapply(pearson_test, function(x){x$r})
					pearson_p=sapply(pearson_test, function(x){x$p})
					pearson_q=p.adjust(pearson_p, method="fdr")
				
					# Spearman test
					spearman_test=lapply(colnames(dimred_data), function(z){
							r=cor.test(test_table[,z], test_table[,feat], method="spearman")$estimate
							p=cor.test(test_table[,z], test_table[,feat], method="spearman")$p.value
							return(list(r=r, p=p))
						}
					)
					spearman_r=sapply(spearman_test, function(x){x$r})
					spearman_p=sapply(spearman_test, function(x){x$p})
					spearman_q=p.adjust(spearman_p, method="fdr")
				
					# Tab
					test_tab=data.frame(pearson_r, pearson_p, pearson_q,
											spearman_r, spearman_p, spearman_q)
					colnames(test_tab)=c("pearson_r", "pearson_p", "pearson_q",
											"spearman_r", "spearman_p", "spearman_q")

					# Return
					return(test_tab)
				}
			)
			names(cont_test)=cont_features
		} else {
			cont_test=NA
		}

	# Return
		return(list(cat_test=cat_test, cont_test=cont_test))

}

#---------------------------------------------------------------------------------------------------------------------------

# Heatmap of features association/correlation to PCs

#---------------------------------------------------------------------------------------------------------------------------

feat_ass_heatmap=function(PCA_ass_cor, PCs=1:6, stat_val="q"){

	# Table of q values
		feat_test=unlist(PCA_ass_cor, recursive=FALSE)
		new_labs=c(names(PCA_ass_cor[["cat_test"]]), names(PCA_ass_cor[["cont_test"]]))
		feat_test=feat_test[!is.na(feat_test)]
		plotting_list=lapply(feat_test, function(x){
				pos=grep("_p", colnames(x))
				if(length(pos)==2){
					pos=grep("spearman_p", colnames(x))
				}
				p_val=x[PCs, pos]
				return(p_val)
			}
		)

	# Heatmap
	if(stat_val=="q"){
		q_val=p.adjust(unlist(plotting_list), method="fdr")
		plotting_table=data.frame(matrix(q_val, nrow=length(plotting_list), byrow=TRUE))
		rownames(plotting_table)=new_labs
		colnames(plotting_table)=paste("PC", PCs, sep='')

		plotting_tab=melt(data.frame(plotting_table, rownames(plotting_table)))
		colnames(plotting_tab)=c("features", "PCs", "values")
		# labels order
		is_num=can_be_numeric(levels(plotting_tab[,1]))
		if(is_num){
			plotting_tab[,1]=factor(plotting_tab[,1], levels=sort(unique(as.numeric(as.character(plotting_tab[,1]))), decreasing=TRUE))
		} else {
			plotting_tab[,"features"]=factor(plotting_tab[,"features"], levels=rev(new_labs))
		}
		
		heat_map=ggplot(plotting_tab, aes(PCs, features, width=1, height=1)) +
					geom_tile(aes(fill=values)) + 
					geom_text(aes(label=round(values, 2)), size=4) +
					theme_light()+
					coord_fixed() +
					scale_fill_gradient2(name=paste(stat_val, " values", sep=""),
											low=pal_nejm("default")(8)[1],
											mid=pal_nejm("default")(8)[6],
											high=pal_nejm("default")(8)[2],
											midpoint = 0.4) +
					theme(axis.title.x = element_blank(), axis.title.y = element_blank())#,
	} else if(stat_val=="p"){
		p_val=unlist(plotting_list)
		plotting_table=data.frame(matrix(p_val, nrow=length(plotting_list), byrow=TRUE))
		rownames(plotting_table)=new_labs
		colnames(plotting_table)=paste("PC", PCs, sep='')

		# remove not significant p (otherwise I cannot plot different color for significant/not significant)
		plotting_tab=melt(data.frame(plotting_table, rownames(plotting_table)))
		colnames(plotting_tab)=c("features", "PCs", "values")
		tmp=sapply(plotting_tab$values, function(x){
				if(is.na(x) | x>=0.05){
					return(NA)
				} else {
					return(x)
				}
			}
		)
		plotting_tab=data.frame(plotting_tab, tmp)
		# labels order
		is_num=can_be_numeric(levels(plotting_tab[,1]))
		if(is_num){
			plotting_tab[,1]=factor(plotting_tab[,1], levels=sort(unique(as.numeric(as.character(plotting_tab[,1]))), decreasing=TRUE))
		} else {
			plotting_tab[,"features"]=factor(plotting_tab[,"features"], levels=rev(new_labs))
		}
		
		if(all(is.na(plotting_tab$tmp))){
			heat_map=ggplot(plotting_tab, aes(PCs, features, width=1, height=1)) +
							geom_tile() + 
							geom_text(aes(label=signif(values, 2)), size=3) +
							theme_light()+
							coord_fixed() +
							scale_fill_gradient2(name=paste(stat_val, " values", sep="")) +
							theme(axis.title.x = element_blank(), axis.title.y = element_blank())
		} else {
			heat_map=ggplot(plotting_tab, aes(PCs, features, width=1, height=1)) +
							geom_tile(aes(fill=tmp)) + 
							geom_text(aes(label=signif(values, 2)), size=3) +
							theme_light() + 
							coord_fixed() +
							scale_fill_gradient2(name=paste(stat_val, " values", sep="")) +
							theme(axis.title.x = element_blank(), axis.title.y = element_blank())		
		}
	}
	# Return
		return(heat_map)
	
}

#---------------------------------------------------------------------------------------------------------------------------

# To add grids to scatterplot3d

#---------------------------------------------------------------------------------------------------------------------------

#' Add grids to a scatterplot3d
#' 
#' @description The goal of this function is to add grids on an existing
#'  plot created using the package scatterplot3d
#' @param x,y,z numeric vectors specifying the x, y, z coordinates of points.
#'  x can be a matrix or a data frame containing 3 columns corresponding to
#'  the x, y and z coordinates. In this case the arguments y and z are optional
#' @param grid specifies the facet(s) of the plot on which grids should be drawn.
#'  Possible values are the combination of "xy", "xz" or "yz".
#'  Example: grid = c("xy", "yz"). The default value is TRUE to add grids only on xy facet.
#' @param col.grid,lty.grid color and line type to be used for grids
#' @param lab a numerical vector of the form c(x, y, len).
#'  The values of x and y give the (approximate) number of tickmarks on the x and y axes.
#' @param lab.z the same as lab, but for z axis
#' @param scale.y of y axis related to x- and z axis
#' @param angle angle between x and y axis
#' @param "xlim, ylim, zlim" the x, y and z limits (min, max) of the plot.
#' 
#' @note
#' Users who want to extend an existing scatterplot3d graphic with the
#'  function addgrids3d, should consider to set the arguments scale.y, angle, ...,
#'  to the value used in scatterplot3d.
#' 
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' 
#' @example
#' library(scatterplot3d)
#' data(iris)
#' scatterplot3d(iris[, 1:3], pch = 16, grid=T, box=F)
#' addgrids3d(iris[, 1:3], grid = c("xy", "xz", "yz"))
addgrids3d <- function(x, y=NULL, z=NULL, grid = TRUE,
                    col.grid = "grey", lty.grid = par("lty"),
                    lab = par("lab"), lab.z = mean(lab[1:2]),
                    scale.y = 1, angle = 40,
                    xlim=NULL, ylim=NULL, zlim=NULL){
  
  
  if(inherits(x, c("matrix", "data.frame"))){
    x <- as.data.frame(x)
    y <- unlist(x[,2])
    z <- unlist(x[,3])
    x <- unlist(x[,1])
  }
  
  p.lab <- par("lab")
  
  angle <- (angle%%360)/90
  yz.f <- scale.y * abs(if (angle < 1) angle else if (angle >3) angle - 4 else 2 - angle)
  yx.f <- scale.y * (if (angle < 2) 1 - angle else angle - 3)
  
  
  # x axis range
  x.range <- range(x[is.finite(x)], xlim)
  x.prty <- pretty(x.range, n = lab[1], min.n = max(1, min(0.5 *lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  x <- x/x.scal
  x.range <- range(x.prty)/x.scal
  x.max <- ceiling(x.range[2])
  x.min <- floor(x.range[1])
  if (!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2]/x.scal))
    x.min <- min(x.min, floor(xlim[1]/x.scal))
  }
  x.range <- range(x.min, x.max)
  
  # y axis range
  y.range <- range(y[is.finite(y)], ylim)
  y.prty <- pretty(y.range, n = lab[2], min.n = max(1, min(0.5 *lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  y <- (y - y.add)/y.scal
  y.max <- (max(y.prty) - y.add)/y.scal
  if (!is.null(ylim))
    y.max <- max(y.max, ceiling((ylim[2] - y.add)/y.scal))
  
  # Z axis range
  z.range <- range(z[is.finite(z)], zlim)
  z.prty <- pretty(z.range, n = lab.z, min.n = max(1, min(0.5 *lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  z <- z/z.scal
  z.range <- range(z.prty)/z.scal
  z.max <- ceiling(z.range[2])
  z.min <- floor(z.range[1])
  if (!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2]/z.scal))
    z.min <- min(z.min, floor(zlim[1]/z.scal))
  }
  z.range <- range(z.min, z.max)
  
  # Add grid
  if ("xy" %in% grid || grid == TRUE) {
    i <- x.min:x.max
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + 
               z.min, col = col.grid, lty = lty.grid)
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min, x.max + 
               (i * yx.f), i * yz.f + z.min, col = col.grid, lty = lty.grid)
  }
   
  if ("xz" %in% grid) {
    i <- x.min:x.max
    segments(i + (yx.f * y.max), yz.f * y.max + z.min, 
             i + (yx.f * y.max), yz.f * y.max + z.max, 
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.max + temp,temp1 + i , col = col.grid, lty = lty.grid)
    
  }
  
  if ("yz" %in% grid) {
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min,  
             x.min + (i * yx.f) ,i * yz.f + z.max,  
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.min, i , col = col.grid, lty = lty.grid)
    }
  
}