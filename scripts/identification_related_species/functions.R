#---------------------------------------------------------------------------------------------------------------------------

##### FUNCTIONS FOR MICROBIOME_RECONSTRUCTION - DIVERSITY #####

#---------------------------------------------------------------------------------------------------------------------------

##### LIBRARIES #####

#---------------------------------------------------------------------------------------------------------------------------

library(coin)

#---------------------------------------------------------------------------------------------------------------------------

# OTHER FUNCTIONS

#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------

# Indep test from coin

#---------------------------------------------------------------------------------------------------------------------------

indep_test=function(taxa_values, metadata, property, to_be_removed=NA, correction, alternative="two.sided"){
	
	# Since the test give problems with continuous values, I bin them
		is_cont=can_be_numeric(metadata[,property])
		if(is_cont){
			decision=property_binning_decision(metadata=metadata, columns_to_be_tested=property, perc_thr=30, perc_NA=30)
			if(is.na(decision)){
				decision=summary(metadata[,property])[c(2,3,5)]
				binned_tab=sapply(metadata[,property], function(y){
						if(is.na(y)){
							return("unknown")
						} else if(y<=decision[1]){
							return("low")
						} else if(y<=decision[2] & y>decision[1]){
							return("medium_low")
						} else if(y<=decision[3] & y>decision[2]){
							return("medium_high")
						} else if(y>decision[3]){
							return("high")
						} 
					}
				)
				metadata[,property]=binned_tab
			} else {
				binned_tab=bin_properties(metadata=metadata, binning_decision=decision)
				metadata[,property]=binned_tab
			}
		}

	# Sample removal
		if(!is.na(to_be_removed)){
			pos=sapply(metadata[,property], function(y){
					!any(y==to_be_removed)
				}
			)
			taxa_values=taxa_values[pos,]
			metadata=metadata[pos,]
		}		
	
	# Remove columns with all zeros
		pos=apply(taxa_values, 2, function(x){
				!all(x==0)
			}
		)
		taxa_values=taxa_values[,pos]

	# Test
		name_test=c()
		p_val=sapply(1:ncol(taxa_values), function(x){
				col_lab=as.character(colnames(taxa_values)[x])
				test_table=data.frame(taxa_values[,col_lab], metadata[,property], as.character(metadata[,correction]))
				colnames(test_table)=c("val", property, correction)
				# remove subgroups with 1 sample only
				is_cont_corr=can_be_numeric(test_table[,correction])
				if(!is_cont_corr){
					joined_subgroups=apply(test_table, 1, function(y){
							paste(y[2:3], collapse="_")
						}
					)
					tmp=table(joined_subgroups)
					if(any(tmp==1)){
						one_subgroup=names(tmp)[tmp==1]
						p=sapply(joined_subgroups, function(y){
								if(any(y==one_subgroup)){
									return(FALSE)
								} else {
									return(TRUE)
								}
							}
						)
						test_table=test_table[p,]
					}
					if(length(unique(test_table[,2]))==1){
						return(NA)
					}
				}
				# tryCatch: if error or warnings it returns NA
				tmp=tryCatch({independence_test(formula(paste("val ~ ", property, " | ", correction, sep='')), 
									alternative=alternative, data=test_table)}, 
								error=function(x){return(NA)}, warning=function(x){return(NA)})
				if(is.na(tmp)){
					return(NA)
				}
				if(alternative!="two.sided"){
					full_descr=coin::show(tmp)$data.name
					name_test=strsplit(full_descr, split="\n")[[1]][2]
					name_test=strsplit(name_test, split="\\(")[[1]][2]
					name_test=gsub(") ", "", name_test)
					name_test=gsub(", ", "_", name_test)
					name_test<<-name_test
				}
				if(is.na(statistic(tmp))){
					p_val=NA
				} else {
					p_val=coin::pvalue(tmp)
				}
				return(p_val)
			}
		)
		q_val=p.adjust(p_val, method="fdr")

	# Create table
		final_tab=cbind(p_val, q_val)
		if(alternative!="two.sided"){
			colnames(final_tab)=paste("greater_", name_test, "_", c("indep_p", "indep_q"), sep="")
		} else {
			colnames(final_tab)=c("indep_p", "indep_q")
		}
		rownames(final_tab)=colnames(taxa_values)
	
	# Return
	return(final_tab)				

}

#---------------------------------------------------------------------------------------------------------------------------

# Wilcoxon or Kruskal-Wallis test

#---------------------------------------------------------------------------------------------------------------------------

wilc_krusk_test=function(taxa_values, metadata, property, to_be_removed=NA){
	
	# Sample removal
		if(!is.na(to_be_removed)){
			pos=sapply(metadata[,property], function(y){
					!any(y==to_be_removed)
				}
			)
			taxa_values=taxa_values[pos,]
			metadata=metadata[pos,]
		}		
		
	# Choose statistic method
		samples_used=table(metadata[,property])
		if(length(samples_used)>2){
			stat_method="krusk"
		} else {
			stat_method="wilcox"
		}
	
	# Remove columns with all zeros
		pos=apply(taxa_values, 2, function(x){
				!all(x==0)
			}
		)
		taxa_values=taxa_values[,pos]

	# Test
	if(stat_method=="wilcox"){
		# Wilcoxon test for all the variables
		p_val=sapply(1:ncol(taxa_values), function(x){
				col_lab=as.character(colnames(taxa_values)[x])
				test_table=data.frame(taxa_values[,col_lab], metadata[,property])
				colnames(test_table)=c("val", property)
				p_val=wilcox.test(formula(paste("val ~ ", property, sep='')), data=test_table)$p.value
				return(p_val)
			}
		)
		q_val=p.adjust(p_val, method="fdr")
	} else if(stat_method=="krusk"){
		# Kruskal test for all the variables
		p_val=sapply(1:ncol(taxa_values), function(x){
				col_lab=as.character(colnames(taxa_values)[x])
				test_table=data.frame(taxa_values[,col_lab], metadata[,property])
				colnames(test_table)=c("val", property)
				p_val=kruskal.test(formula(paste("val ~ ", property, sep='')), data=test_table)$p.value
				return(p_val)
			}
		)
		q_val=p.adjust(p_val, method="fdr")
	}

	# Create table
		final_tab=cbind(p_val, q_val)
		colnames(final_tab)=c(paste(stat_method, "_p", sep=''), paste(stat_method, "_q", sep=''))
		rownames(final_tab)=colnames(taxa_values)
	
	# Return
	return(final_tab)	

}

#---------------------------------------------------------------------------------------------------------------------------

# correlation test

#---------------------------------------------------------------------------------------------------------------------------

corr_test=function(taxa_values, metadata, property, to_be_removed=NULL, met="spearman"){
	
	# Sample removal
		if(!is.null(to_be_removed)){
			pos=sapply(metadata[,property], function(y){
					!any(y==to_be_removed)
				}
			)
			taxa_values=taxa_values[pos,]
			metadata=metadata[pos,]
		}		
	
	# Remove columns with all zeros
		pos=apply(taxa_values, 2, function(x){
				!all(x==0)
			}
		)
		taxa_values=taxa_values[,pos]

	# Test
		r_p_val=sapply(1:length(taxa_values), function(x){
				col_lab=as.character(colnames(taxa_values)[x])
				pos=which(is.na(metadata[,property]))
				if(length(pos)==0){
					r=cor.test(taxa_values[,col_lab], metadata[,property], method=met)$estimate
					p=cor.test(taxa_values[,col_lab], metadata[,property], method=met)$p.value
				} else {
					r=cor.test(taxa_values[-pos,col_lab], metadata[-pos,property], method=met)$estimate
					p=cor.test(taxa_values[-pos,col_lab], metadata[-pos,property], method=met)$p.value
				}
				return(c(r=r, p=p))
			}
		)
		r_p_val=t(r_p_val)
		q_val=p.adjust(r_p_val[,"p"], method="fdr")

	# Create table
		final_tab=cbind(r_p_val, q_val)
		colnames(final_tab)=c(paste(met, "_r", sep=''), paste(met, "_p", sep=''), paste(met, "_q", sep=''))
		rownames(final_tab)=colnames(taxa_values)
		
	# Return
	return(final_tab)				

}

#---------------------------------------------------------------------------------------------------------------------------

# Mean and median

#---------------------------------------------------------------------------------------------------------------------------

measures_of_central_tendency=function(taxa_values, metadata, property=NULL, to_be_removed=NA, met="median"){
	
	# Sample removal
		if(!is.na(to_be_removed)){
			pos=sapply(metadata[,property], function(y){
					!any(y==to_be_removed)
				}
			)
			taxa_values=taxa_values[pos,]
			metadata=metadata[pos,]
		}		
	
	# Remove columns with all zeros
		pos=apply(taxa_values, 2, function(x){
				!all(x==0)
			}
		)
		taxa_values=taxa_values[,pos]

	# info
		if(!is.null(property)){
			groups=sort(as.character(unique(metadata[,property])))
		}
				
	# Mean
		if(met=="mean" & !is.null(property)){
			tot=apply(taxa_values, 2, mean)
			mean_tabs=lapply(groups, function(x){
					pos=which(metadata[,property]==x)
					val=apply(taxa_values[pos,], 2, mean)
				}
			)
			names(mean_tabs)=groups

			tot_nonzeros=apply(taxa_values, 2, function(x){
					mean(x[x!=0])
				}
			)
			mean_nonzeros_tabs=lapply(groups, function(x){
					pos=which(metadata[,property]==x)
					val=apply(taxa_values[pos,], 2, function(x){
							mean(x[x!=0])
						}
					)
				}
			)
			names(mean_nonzeros_tabs)=paste("nonzeros_", groups,sep='')
			
			df=data.frame(tot, mean_tabs, tot_nonzeros, mean_nonzeros_tabs, check.names=FALSE) #check.names to avoid the adding of X to numeric labels
			colnames(df)=paste("mean_", colnames(df), sep='')
		} else if(met=="mean" & is.null(property)){
			tot=apply(taxa_values, 2, mean)
			
			tot_nonzeros=apply(taxa_values, 2, function(x){
					mean(x[x!=0])
				}
			)
			
			df=data.frame(tot, tot_nonzeros, check.names=FALSE) #check.names to avoid the adding of X to numeric labels
			colnames(df)=paste("mean_", colnames(df), sep='')
		}

	# Median
		if(met=="median" & !is.null(property)){
			tot=apply(taxa_values, 2, median)
			median_tabs=lapply(groups, function(x){
					pos=which(metadata[,property]==x)
					val=apply(taxa_values[pos,], 2, median)
				}
			)
			names(median_tabs)=groups

			tot_nonzeros=apply(taxa_values, 2, function(x){
					median(x[x!=0])
				}
			)
			median_nonzeros_tabs=lapply(groups, function(x){
					pos=which(metadata[,property]==x)
					val=apply(taxa_values[pos,], 2, function(x){
							median(x[x!=0])
						}
					)
				}
			)
			names(median_nonzeros_tabs)=paste("nonzeros_", groups,sep='')
			
			df=data.frame(tot, median_tabs, tot_nonzeros, median_nonzeros_tabs, check.names=FALSE)
			colnames(df)=paste("median_", colnames(df), sep='')
		} else if(met=="median" & is.null(property)){
			tot=apply(taxa_values, 2, median)
			
			tot_nonzeros=apply(taxa_values, 2, function(x){
					median(x[x!=0])
				}
			)
				
			df=data.frame(tot, tot_nonzeros, check.names=FALSE)
			colnames(df)=paste("median_", colnames(df), sep='')
		}

	# Presence
		if(met=="presence" & !is.null(property)){
			tot=apply(taxa_values, 2, function(x){
					(length(x[x!=0])/length(x))*100
				}
			)
			presence_tabs=lapply(groups, function(x){
					pos=which(metadata[,property]==x)
					val=taxa_values[pos,]
					pres=apply(val, 2, function(y){
							(length(y[y!=0])/length(y))*100
						}
					)
					return(pres)
				}
			)
			names(presence_tabs)=groups

			df=data.frame(tot, presence_tabs, check.names=FALSE)
			colnames(df)=paste("presence_", colnames(df), sep='')
		} else if(met=="presence" & is.null(property)){
			tot=apply(taxa_values, 2, function(x){
					(length(x[x!=0])/length(x))*100
				}
			)
			
			df=data.frame(tot, check.names=FALSE)
			colnames(df)=paste("presence_", colnames(df), sep='')
		}

	# Generalized log fold change
	# (see Zeller, "Meta-analysis of fecal metagenomes reveals global microbial signatures that are specific for colorectal cancer", 
	# section "Univariate meta-analysis for the identification of CRC-associated gut microbial species")
		if(met=="genlogFC"){
			if(length(groups)!=2){
				return(NULL)
			}
			if(any(taxa_values>1)){
				print("Values must be between 0 and 1 (=relative abundances!)")
				return(NULL)
			}
			min_score=min(taxa_values[taxa_values!=0])
			no_zeros=taxa_values+min_score
			log_taxa_values=log(no_zeros, base=10)
			colnames(log_taxa_values)=colnames(taxa_values)

			gen_log_FC_values=apply(log_taxa_values, 2, function(x){
					log_values=lapply(groups, function(z){
							pos=metadata[,property]==z
							log_score_values=x[pos]
							perc_log_score_values=quantile(log_score_values, probs=seq(0.1, 0.9, 0.1), na.rm=TRUE)
							return(perc_log_score_values)
						}
					)
					log_FC_values=mean(log_values[[1]]-log_values[[2]])
					return(log_FC_values)
				}
			)

			df=data.frame(gen_log_FC_values, check.names=FALSE)
			colnames(df)=paste("genlogFC_", groups[[1]], "_", groups[[2]], sep='')
		}
		
	# Return
	return(df)				

}
