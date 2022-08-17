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
