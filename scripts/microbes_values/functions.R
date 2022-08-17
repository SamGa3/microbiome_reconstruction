#---------------------------------------------------------------------------------------------------------------------------

##### FUNCTIONS FOR MICROBIOME_RECONSTRUCTION - MICROBES_VALUES #####

#---------------------------------------------------------------------------------------------------------------------------

##### LIBRARIES #####

#---------------------------------------------------------------------------------------------------------------------------

library(data.table)
library(plyr)

#---------------------------------------------------------------------------------------------------------------------------

##### FUNCTIONS #####

#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------

# Extract tables from pathseq outputs

#---------------------------------------------------------------------------------------------------------------------------
		
create_reads_table=function(pathseq_output_folder, samples_id, taxon, kingdom, reads_type){
	# Create Pathseq outputs paths
		# list file folders
		folders=dir(pathseq_output_folder) 
		# create paths
		paths=paste(pathseq_output_folder, "/", folders, "/score_out.txt", sep='')
		# extract Pathseq outputs
		pathseq_outputs=lapply(paths, function(x){
				if(file.exists(x)){
					table=fread(file=x, sep="\t", dec=".", quote="\"",
								nrows=Inf, header=TRUE, data.table=options(datatable.fread.datatable=FALSE),
								na.strings=getOption("datatable.na.strings","NA"), 
								select = c("tax_id","type","name","kingdom",reads_type),
								stringsAsFactors=FALSE)
					return(table)
				} else {
					return(NA)
				}
			}
		)
		# remove NA
		pos_na=!is.na(pathseq_outputs)
		pathseq_outputs=pathseq_outputs[pos_na]
		folders=folders[pos_na]
		# file_id
		samples_id_tmp=gsub(".bam", "_out", samples_id[,"File Name"])
		file_id=sapply(folders, function(x){
				pos=which(samples_id_tmp==x)
				samples_id[pos,"File ID"]
			}
		)
		names(pathseq_outputs)=file_id
	# Select by kingdom and taxon
		selected_list=lapply(pathseq_outputs, function(x){
				# select by kingdom and taxon
				pos=which(x[,"kingdom"]==kingdom & x[,"type"]==taxon)
				selection=as.data.frame(matrix(x[pos,reads_type], nrow=1))
				colnames(selection)=as.character(x[pos,"tax_id"])
				return(selection)
			}
		)
		names(selected_list)=names(pathseq_outputs)
		# Make the list a data.frame
		selected_table=do.call(rbind.fill, selected_list)
		rownames(selected_table)=names(selected_list)
		# remove NAs
		selected_table[is.na(selected_table)]=0
	# Return
		return(selected_table)
}

#---------------------------------------------------------------------------------------------------------------------------

# Merge taxa in table

#---------------------------------------------------------------------------------------------------------------------------

merge_taxa_in_tab=function(tab, merging_tax){
	# Check if there are taxa that must be merged
		to_be_merged=intersect(colnames(tab), as.character(merging_tax[,1]))
		if(length(to_be_merged)==0){
			return(tab)
		}

	# Create a table without the taxa that must be merged
		common_merging=setdiff(colnames(tab), as.character(merging_tax[,1]))
		new_tab=tab[,common_merging]

	# Sum values of taxa that must be merged
		for(i in 1:length(to_be_merged)){
			merging_tax_id=to_be_merged[[i]]
			pos_merging_tax=which(merging_tax[,1]==merging_tax_id)
			real_tax=as.character(merging_tax[pos_merging_tax,2])
			# If pathseq reports only the value of the bacteria that must be merged
			if(!any(colnames(new_tab)==real_tax)){
				new_tab2=cbind.data.frame(new_tab, tab[,merging_tax_id])
				colnames(new_tab2)[ncol(new_tab2)]=real_tax
				new_tab=new_tab2
			} else {
				new_tab[,real_tax]=new_tab[,real_tax]+tab[,merging_tax_id]
			}
		}

	# Return
		return(new_tab)
}
