#---------------------------------------------------------------------------------------------------------------------------

##### GENERAL FUNCTIONS FOR MICROBIOME_RECONSTRUCTION #####

#---------------------------------------------------------------------------------------------------------------------------

##### LIBRARIES #####

#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------

##### PALETTES #####

#---------------------------------------------------------------------------------------------------------------------------

# palette default ggplot2
gg_color_hue=function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# my palettes
my_palettes=list(
				default="default",
				COAD_IEO_rest=c("#D43F3AFF", "#E48B32", "#BFBFBF", "#BFBFBF", "#BFBFBF", "#BFBFBF", "#BFBFBF", "#BFBFBF", "#BFBFBF"),
				locuszoom2=c("#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF", "#357EBDFF", "#9632B8FF", "#B8B8B8FF", "#EE4C97FF"),
				CMS=c("#FFA62B", "#709AE1FF", "#FD7446FF", "#1A9993FF", "#A4B07D", "#91331FFF"),
				Left_right=c("#1F77B4FF", "#FF7F0EFF", "#BFBFBF"),
				MSI=c("#2CA02CFF", "#D62728FF", "#BFBFBF"),
				viridis7=c("#440154FF", "#443A83FF", "#31688EFF", "#21908CFF", "#35B779FF", "#8FD744FF", "#FDE725FF"),
				nejm=c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF"),
				jama=c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF"),
				brewer_Paired_Dark2_Set2=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", 
																		"#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
																		"#A6761D", "#666666", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", 
																		"#B3B3B3")
			)

#---------------------------------------------------------------------------------------------------------------------------

##### FUNCTIONS #####

#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------

# Check if vector can be converted to numeric

#---------------------------------------------------------------------------------------------------------------------------

can_be_numeric=function(x, n_to_be_cat=3){
	# Remove factor because they can be converted to numeric
	if(class(x)=="factor"){
		val=as.character(x)
	} else {
		val=x
	}
	numNAs=sum(is.na(val))
	numNAs_new=suppressWarnings(sum(is.na(as.numeric(val))))
	res=numNAs_new==numNAs
	# if, when I convert to numeric, I get the same number of NAs, 
	# it means that I am not entroducing new NAs so the values can be converted to NA
	if(sum(is.na(val))==length(val)){ 
		res=FALSE
	}
	# If the value is numeric but I have only 3 entries, I consider it categorical
	if(dim(table(val))<=n_to_be_cat){
		res=FALSE
	}
	return(res)
}

#---------------------------------------------------------------------------------------------------------------------------

# Match rows of table to metadata

#---------------------------------------------------------------------------------------------------------------------------

match_to_metadata=function(metadata, other_tab, colname_metadata="rownames", colname_other_tab_values="rownames"){
	# if need to match to rownames, make a new column
	if(colname_metadata=="rownames"){
		metadata=data.frame(rownames(metadata), metadata)
		colnames(metadata)[1]=colname_metadata
		colname_metadata=colname_metadata
	}
	# if need to match to rownames, make a new column
	if(colname_other_tab_values=="rownames"){
		other_tab=data.frame(rownames(other_tab), other_tab)
		colnames(other_tab)[1]=colname_metadata
		colname_other_tab_values=colname_metadata
	}
	# order tabs
	other_tab_new_list=lapply(as.character(metadata[,colname_metadata]), function(x){
			p=which(other_tab[,colname_other_tab_values]==x)
			if(length(p)==0){
				tab=matrix(rep(NA, ncol(other_tab)), ncol=ncol(other_tab))
				colnames(tab)=colnames(other_tab)
				return(tab)
			} else {
				return(other_tab[p,])
			}
		}
	)
	other_tab_new_tmp=do.call(rbind, other_tab_new_list)
	other_tab_new=other_tab_new_tmp[,colnames(other_tab_new_tmp)!=colname_other_tab_values]
	# if there is only one column, it is returning a vector, so I convert it to matrix
	if(class(other_tab_new)!="matrix" & class(other_tab_new)!="data.frame"){
		other_tab_new=as.matrix(other_tab_new)
		colnames(other_tab_new)=colnames(other_tab_new_tmp)[colnames(other_tab_new_tmp)!=colname_other_tab_values]
	}
	
	return(other_tab_new)
}

#---------------------------------------------------------------------------------------------------------------------------

# Match rows and columns of two tables, set columns to be selected

#---------------------------------------------------------------------------------------------------------------------------

match_tables=function(tab1, tab2){
	# Check if dimensions are the same
		if(any(dim(tab1)!=dim(tab2))){
			print("Tables with different dimensions")
			col_names=intersect(colnames(tab1), colnames(tab2))
			tab1=tab1[,col_names]
			tab2=tab2[,col_names]
		}

	# Check same columns
		if(length(intersect(colnames(tab1), colnames(tab2)))!=ncol(tab1)){
			stop("Not same column labels")
		}
		if(length(intersect(colnames(tab2), colnames(tab1)))!=ncol(tab2)){
			stop("Not same column labels")
		}
		if(ncol(tab1)==1){ # with 1 column it returns a vector
			row_names=rownames(tab2)
			tab2=matrix(tab2[,colnames(tab1)])
			colnames(tab2)=colnames(tab1)
			rownames(tab2)=row_names
		} else if(nrow(tab1)==1){ # with 1 row it returns a vector
			row_names=rownames(tab2)
			tab2=t(matrix(tab2[,colnames(tab1)]))
			colnames(tab2)=colnames(tab1)
			rownames(tab2)=row_names
		} else {
			tab2=tab2[,colnames(tab1)]
		}	

	# Chech same rows
		if(length(intersect(rownames(tab1), rownames(tab2)))!=nrow(tab1)){
			stop("Not same row labels")
		}
		if(length(intersect(rownames(tab2), rownames(tab1)))!=nrow(tab2)){
			stop("Not same row labels")
		}
		if(ncol(tab1)==1){ # with 1 column it returns a vector
			tab2=matrix(tab2[rownames(tab1),])
			rownames(tab2)=rownames(tab1)
			colnames(tab2)=colnames(tab1)
		} else if(ncol(tab1)==1){ # with 1 row it returns a vector
			tab2=t(matrix(tab2[rownames(tab1),]))
			rownames(tab2)=rownames(tab1)
			colnames(tab2)=colnames(tab1)
		} else {
			tab2=tab2[rownames(tab1),]
		}

	# Return
	return(list(tab1, tab2))
}

#---------------------------------------------------------------------------------------------------------------------------

# Match metadata to bacteria values tables

#---------------------------------------------------------------------------------------------------------------------------

match_metadata=function(metadata, taxa_values, colname_metadata, colname_taxa_values="rownames"){
	# if need to match to rownames, make a new column
	if(colname_taxa_values=="rownames"){
		taxa_values=data.frame(rownames(taxa_values), taxa_values)
		colnames(taxa_values)[1]=colname_metadata
		colname_taxa_values=colname_metadata
	}
	# Match the tables
	pos=sapply(taxa_values[,colname_taxa_values], function(x){
			p=which(as.character(metadata[,colname_metadata])==as.character(x))
			if(length(p)==0){
				print(paste("Sample ", x, " missing", sep=''))
				return(NA)
			}
			return(p)
		}
	)
	metadata_new=metadata[pos,]
	metadata_new[,colname_metadata]=taxa_values[,colname_taxa_values]
	return(metadata_new)
}

#---------------------------------------------------------------------------------------------------------------------------

# Join metadata

#---------------------------------------------------------------------------------------------------------------------------

join_metadata=function(metadata_paths, join_by, matching1, matching2, metadata_is_list=FALSE, separator="\t", keep_column=""){
	# If the metadata to be joined are a list of matrices or data.frames, or a list of paths
	if(!metadata_is_list){
		metadatas=lapply(metadata_paths, function(x){
				read.csv(x, sep=separator, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
			}
		)
	} else {
		metadatas=metadata_paths
	}
	# join tables by columns or rows
	if(join_by=="columns"){
		new_full_metadata=metadatas[[1]]
		for(i in 2:length(metadatas)){
			metadata_to_be_added=match_to_metadata(metadata=new_full_metadata, other_tab=metadatas[[i]], 
													colname_metadata=matching1, colname_other_tab_values=matching2)
			new_full_metadata=data.frame(new_full_metadata, metadata_to_be_added)
		}
		full_metadata=new_full_metadata
	} else if(join_by=="rows"){
		# Common columns
		col_names_tot=unique(unlist(lapply(metadatas, colnames)))
		col_pres=sapply(col_names_tot, function(x){
				tmp=sapply(metadatas, function(y){
						if(any(colnames(y)==x)){
							return(1)
						} else {
							return(0)
						}
					}
				)
				return(sum(tmp))
			}
		)
		sel_colnames=names(col_pres)[col_pres==length(metadatas)]
		new_metadatas=lapply(metadatas, function(x){
				x[,sel_colnames]
			}
		)
		full_metadata=do.call(rbind, new_metadatas)
		# to force a column to be kept: keep_column is a list with 2 vectors
		# the vector named val lists the column labels that must be kept
		# the vector named fill lists the values that must be added to the tables that miss the column
		if(all(keep_column!="")){
		  for(i in 1:length(keep_column[[1]])){
		    tmp=lapply(metadatas, function(x){
		        if(any(colnames(x)==keep_column$val[[i]])){
		          return(x[,keep_column$val[[i]]])
		        } else {
		          return(rep(keep_column$fill[[i]], nrow(x)))
		        }
		      }
		    )
		    new_val=unlist(tmp)
		    full_metadata=data.frame(full_metadata, new_val)
		    colnames(full_metadata)[ncol(full_metadata)]=keep_column[[i]]
		  }
		}
	}

	# Convert to numeric
	cont_column=NULL
	for(i in 1:ncol(full_metadata)){
		cont_column[[i]]=can_be_numeric(full_metadata[,i])
	}
	for(i in which(cont_column)){
		full_metadata[,i]=as.numeric(as.character(full_metadata[,i]))
	}
	# Add unknown to categorical
	cat_column=!cont_column
	# Remove
	for(i in which(cat_column)){
		for(y in 1:nrow(full_metadata)){
			if(is.na(full_metadata[y,i])){
				full_metadata[y,i]="unknown"
			}
		}
	}

	return(full_metadata)
}

#---------------------------------------------------------------------------------------------------------------------------

# Join taxa tables

#---------------------------------------------------------------------------------------------------------------------------

join_taxa_table=function(taxa_paths){
	# upload the taxa tables
	taxa=lapply(taxa_paths, function(x){
			read.csv(x, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
		}
	)
	names(taxa)=NULL
	# Common columns
	col_names_tot=unique(unlist(sapply(taxa, colnames)))
	new_taxa_val=lapply(taxa, function(x){
			missing_taxa=setdiff(col_names_tot, colnames(x))
			if(length(missing_taxa)!=0){
				tab_missing=matrix(rep(0, length(missing_taxa)*nrow(x)), ncol=length(missing_taxa))
				rownames(tab_missing)=rownames(x)
				colnames(tab_missing)=missing_taxa
				res=data.frame(x, tab_missing)
				colnames(res)=gsub("X", "", colnames(res))
				res=res[,col_names_tot]
			} else {
				res=x
			}	
			return(res)
		}
	)
	taxa_values=do.call(rbind, new_taxa_val)
	return(taxa_values)
}

#---------------------------------------------------------------------------------------------------------------------------

# Create a new property from another one 

#---------------------------------------------------------------------------------------------------------------------------

create_new_property=function(old_feat, met, metadata){
	old_val=metadata[,old_feat]
	old_val=as.data.frame(old_val)
	colnames(old_val)=old_feat
	# choose the method to create the new property
	if(met=="corr_plate_id"){
		tab=table(old_val)
		thr=nrow(metadata)*0.1
		# number of samples per plate_id must be more than 10% of all samples but
		# consider 10% of the samples is too much when I consider a lot of samples
		# so I set that if the plate id has 5 samples, it has not to be joined
		if(thr>5){
			thr=5
		}
		repeat{
			to_be_joined=names(tab[tab<thr])
			new_val=sapply(old_val[,1], function(x){
					if(any(as.character(x)==to_be_joined)){
						return("other")
					} else {
						return(as.character(x))
					}
				}
			)
			new_tab=table(new_val)
			if(all(new_tab>=thr)){
				break
			} else {
				thr=thr+1
			}
		}
	} else if(met=="sum_feat"){
		for(i in 1:ncol(old_val)){
			old_val[,i]=as.numeric(as.character(old_val[,i]))
		}
		new_val=apply(old_val, 1, sum, na.rm=TRUE)
	} else if(met=="mult_feat"){
		for(i in 1:ncol(old_val)){
			old_val[,i]=as.numeric(as.character(old_val[,i]))
		}
		new_val=apply(old_val, 1, function(x){
				val=1
				for(i in 1:length(x)){
					val=val*as.numeric(as.character(x[[i]]))
				}
				return(val)
			}
		)
	} else if(met=="perc_feat"){
		for(i in 1:ncol(old_val)){
			old_val[,i]=as.numeric(as.character(old_val[,i]))
		}
		new_val=apply(old_val, 1, function(x){
				if(any(is.na(x))){
					return(NA)
				} else {
					(as.numeric(as.character(x[1]))/as.numeric(as.character(x[2])))*100
				}
			}
		)
	} else if(met=="range_0to1"){
		for(i in 1:ncol(old_val)){
			old_val[,i]=as.numeric(as.character(old_val[,i]))
		}
		new_val=apply(old_val, 2, function(x){
				(x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
			}
		)
	} else if(met=="zero_to_NA"){
		new_val=old_val
		for(i in 1:ncol(old_val)){
			new_val[,i]=as.numeric(as.character(new_val[,i]))
			p=which(!is.na(new_val[,i]) & new_val[,i]==0)
			new_val[p,i]=NA
		}
	} 
	
	return(new_val)
}

#---------------------------------------------------------------------------------------------------------------------------

# Decide to bin or keep continuous the properties

#---------------------------------------------------------------------------------------------------------------------------

property_binning_decision=function(metadata, columns_to_be_tested, perc_thr=30, perc_NA=30){
	
	# no NA
		perc_nas=sapply(columns_to_be_tested, function(x){
				n_nas=length(which(is.na(metadata[,x])))
				perc=(n_nas/nrow(metadata))*100
				return(perc)
			}
		)
		nas_cont_feat=names(perc_nas[which(perc_nas>=perc_NA)])
		nonas_cont_feat=names(perc_nas[which(perc_nas<perc_NA)])

	# no zero
		perc_zeros=sapply(nonas_cont_feat, function(x){
				n_zeros=length(which(metadata[!is.na(metadata[,x]),x]==0))
				perc=(n_zeros/nrow(metadata[!is.na(metadata[,x]),]))*100
				return(perc)
			}
		)
		zeros_cont_feat=names(perc_zeros[which(perc_zeros>=perc_thr)])
		nozeros_cont_feat=names(perc_zeros[which(perc_zeros<perc_thr)])

	# normality check
		is_norm_p=sapply(nozeros_cont_feat, function(x){
				if(var(metadata[!is.na(metadata[,x]),x])!=0){
					p_norm=shapiro.test(metadata[!is.na(metadata[,x]),x])$p.value
					return(p_norm)
				} else {
					return(NA)
				}
			}
		)
		is_norm_q=p.adjust(is_norm_p[!is.na(is_norm_p)], method="fdr")
		norm_cont_feat=nozeros_cont_feat[is_norm_q>=0.05]
		not_norm_cont_feat=nozeros_cont_feat[is_norm_q<0.05]

	# normality
		is_normal=sapply(norm_cont_feat, function(x){
				mean(metadata[,x], na.rm=TRUE)
			}
		)

	# bimodality
		is_bimodal=sapply(not_norm_cont_feat, function(x){
				is_bi=is.bimodal(metadata[,x])
				if(is_bi){
					modes_values=sort(Modes(metadata[,x])[[1]], decreasing=FALSE)
					den=density(metadata[,x], na.rm=TRUE)
					pos1=which(den$x==modes_values[[1]])
					pos2=which(den$x==modes_values[[2]])
					if(length(pos1)==0){
						pos1=which(abs(den$x-modes_values[[1]])==min(abs(den$x-modes_values[[1]])))
					}
					if(length(pos2)==0){
						pos2=which(abs(den$x-modes_values[[2]])==min(abs(den$x-modes_values[[2]])))
					}
					min_val=min(den$y[pos1:pos2])
					pos3=which(den$y==min_val)
					if(length(pos3)>1){
						pos3=pos3[round(length(pos3)/2, 1)]
					}
					return(den$x[pos3])
				} else {
					return(NA)
				}
			}
		)
		bimodal_not_norm_cont_feat=not_norm_cont_feat[!is.na(is_bimodal)]
		not_bimodal_not_norm_cont_feat=not_norm_cont_feat[is.na(is_bimodal)]

	# trimodality
		is_trimodal=lapply(not_bimodal_not_norm_cont_feat, function(x){
				is_tri=is.trimodal(metadata[,x])
				if(is_tri){
					modes_values=sort(Modes(metadata[,x])[[1]], decreasing=FALSE)
					den=density(metadata[,x], na.rm=TRUE)
					pos1=which(den$x==modes_values[[1]])
					pos2=which(den$x==modes_values[[2]])
					pos3=which(den$x==modes_values[[3]])
					if(length(pos1)==0){
						pos1=which(abs(den$x-modes_values[[1]])==min(abs(den$x-modes_values[[1]])))
					}
					if(length(pos2)==0){
						pos2=which(abs(den$x-modes_values[[2]])==min(abs(den$x-modes_values[[2]])))
					}
					if(length(pos3)==0){
						pos2=which(abs(den$x-modes_values[[3]])==min(abs(den$x-modes_values[[3]])))
					}
					min_val1=min(den$y[pos1:pos2])
					pos4=which(den$y==min_val1)
					min_val2=min(den$y[pos2:pos3])
					pos5=which(den$y==min_val2)
					return(c(den$x[pos4], den$x[pos5]))
				} else {
					return(NA)
				}
			}
		)
		names(is_trimodal)=not_bimodal_not_norm_cont_feat
		trimodal_not_bimodal_not_norm_cont_feat=not_bimodal_not_norm_cont_feat[!is.na(is_trimodal)]
		not_trimodal_not_bimodal_not_norm_cont_feat=not_bimodal_not_norm_cont_feat[is.na(is_trimodal)]

	# rest
	 	rest=lapply(not_trimodal_not_bimodal_not_norm_cont_feat, function(x){
				dec=summary(metadata[,x])[c(2,3,5)]
				dec=dec[!duplicated(dec)]
				return(dec)
			}
		)
		names(rest)=not_trimodal_not_bimodal_not_norm_cont_feat

	# final thresholds
		res=lapply(columns_to_be_tested, function(x){
				if(any(x==nas_cont_feat)){
					return(NA)
				} else if(any(x==zeros_cont_feat)){
					return(0)
				} else if(any(x==norm_cont_feat)){
					return(is_normal[[x]])
				} else if(any(x==bimodal_not_norm_cont_feat)){
					return(is_bimodal[[x]])
				} else if(any(x==trimodal_not_bimodal_not_norm_cont_feat)){
					return(is_trimodal[[x]])
				} else {
					return(rest[[x]])
				} 
			}
		)
		names(res)=columns_to_be_tested

	# return
		return(res)

}

#---------------------------------------------------------------------------------------------------------------------------

# Add binned columns

#---------------------------------------------------------------------------------------------------------------------------

bin_properties=function(metadata, binning_decision){
	
	to_be_binned=names(binning_decision[!is.na(binning_decision)])

	new_feat_list=lapply(to_be_binned, function(x){
			if(length(binning_decision[[x]])==1){
				if(binning_decision[x]==0){
					new_col=sapply(metadata[,x], function(y){
							if(is.na(y)){
								return("unknown")
							} else if(y==0){
								return("absent")
							} else {
								return("present")
							}
						}
					)
				} else {
					new_col=sapply(metadata[,x], function(y){
							if(is.na(y)){
								return("unknown")
							} else if(y<=binning_decision[x]){
								return("low")
							} else {
								return("high")
							}
						}
					)
				}
			} else if(length(binning_decision[[x]])==2){
				new_col=sapply(metadata[,x], function(y){
						if(is.na(y)){
							return("unknown")
						} else if(y<=binning_decision[[x]][1]){
							return("low")
						} else if(y>binning_decision[[x]][1] & y<=binning_decision[[x]][2]){
							return("medium")
						} else if(y>binning_decision[[x]][2]){
							return("high")
						}
					}
				)
			} else if(length(binning_decision[[x]])==3){
				new_col=sapply(metadata[,x], function(y){
						if(is.na(y)){
							return("unknown")
						} else if(y<=binning_decision[[x]][1]){
							return("low")
						} else if(y>binning_decision[[x]][1] & y<=binning_decision[[x]][2]){
							return("medium_low")
						} else if(y>binning_decision[[x]][2] & y<=binning_decision[[x]][3]){
							return("medium_high")
						} else if(y>binning_decision[[x]][3]){
							return("high")
						}
					}
				)
			} 
			
			return(new_col)
		}
	)
	new_feat_tab=do.call(cbind, new_feat_list)
	colnames(new_feat_tab)=paste("bin_", to_be_binned, sep='')

	# return
		return(new_feat_tab)

}
