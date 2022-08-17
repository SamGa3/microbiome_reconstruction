#---------------------------------------------------------------------------------------------------------------------------

##### FUNCTIONS FOR MICROBIOME_RECONSTRUCTION - BATCH_CORRECTION #####

#---------------------------------------------------------------------------------------------------------------------------

##### LIBRARIES #####

#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------

# OTHER FUNCTIONS

#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------

# Transformation of microbial data to batch correct them

#---------------------------------------------------------------------------------------------------------------------------

transform_data_to_batchcorrect=function(taxa_values){
	# Remove microbes with all zeros
		tmp_matrix=t(as.matrix(taxa_values))
		pos=apply(tmp_matrix, 1, function(x){
				if(sum(x)==0){
					return(FALSE)
				} else {
					return(TRUE)
				}
			}
		)
		taxa_values_matrix=tmp_matrix[pos,]

	# Scale values
		scale_taxa_values_matrix=scale(scale=TRUE, x=taxa_values_matrix, center=TRUE)

	# To log transform data, I must remove zeros. To do it, I add to all the values the absolute value of the minimum value 
	# present in the table
		min_score=min(scale_taxa_values_matrix[scale_taxa_values_matrix!=0])
		min_score_up=abs(min_score)+abs(min_score)/2
		no_zeros=scale_taxa_values_matrix+min_score_up
		scale_log_taxa_values_matrix=log(no_zeros, base=2)	

	# Return
	return(scale_log_taxa_values_matrix)
}
