args = commandArgs(T)

samples = args[1]
lane_dir = args[2]
output_dir_root = args[3]

sample_list = read.delim(samples)

data = NULL
for(i in 1:nrow(sample_list))		
{
	lane = as.character(sample_list$Lane[i])	
	cfrna_sample = as.character(sample_list$SequencingID[i])
	short_sample_id = gsub("Sample_", "", gsub("_cfrna", "", cfrna_sample))
	cfrna_ID = as.character(sample_list$ID[i])
	
	lane_path = file.path(lane_dir, lane)
	
	freq_vars = file.path(lane_path, "exp_analysis/deduped/output/barcode-deduped",
		cfrna_sample, paste0(cfrna_sample, ".sorted.dualindex-deduped.sorted.freq.paired.Q30.txt"))
	
	freq_indels = file.path(lane_path, "exp_analysis/deduped/output/barcode-deduped",
		cfrna_sample, paste0(cfrna_sample, ".sorted.dualindex-deduped.sorted.indels.paired.Q30.txt"))
	
	star_bam = file.path(lane_path, "exp_analysis/STAR/second_pass/",
		paste0(cfrna_sample, ".Aligned.sortedByCoord.out.bam"))

	star_fusion = file.path(lane_path, "genotyping/STAR_Fusion", 
		cfrna_sample, "star-fusion.preliminary/star-fusion.fusion_candidates.preliminary.wSpliceInfo.wAnnot")			

	brk = F
	for(path in c(freq_vars, freq_indels, star_fusion))
	{
		if(!file.exists(path))
		{
			#stop(paste("File", path, "doesn't exist!"))
			cat(paste("File", path, "doesn't exist!\n"))
			brk = T			
		}
	}
	
	if(brk)
	{				
		stop("Some input files are missing. Make sure you are providinng the correct input!")
	}

	if(length(list.files(file.path(lane_path), pattern = "*ample2barcode*", full.names = T)) == 0)
	{
		cat(paste0("No sample2barcode file found for lane: ", lane, ". Files in the lane folder: ", paste(list.files(lane_path), collapse = "\n\t"), "\n"))
		next
	}
	
	s2b = read.delim(rev(sort(list.files(file.path(lane_path), pattern = "*ample2barcode*", full.names = T)))[1], header = F)
	
	new_file = gsub(".sorted.dualindex-deduped.sorted.freq.paired.Q30.txt", "", basename(freq_vars))
	output_dir = file.path(output_dir_root, paste0(lane,"_",cfrna_ID))
	dir.create(output_dir, recursive = T, showWarning = F)
	
	system(paste("ln -sf ", gsub(".freq.paired.Q30.txt", ".bam", freq_vars), file.path(output_dir, basename(gsub(".freq.paired.Q30.txt", ".bam", freq_vars)))))
	system(paste("ln -sf ", gsub(".freq.paired.Q30.txt", ".bam.bai", freq_vars), file.path(output_dir, basename(gsub(".freq.paired.Q30.txt", ".bam.bai", freq_vars)))))	
	system(paste("ln -sf ", freq_vars, file.path(output_dir, basename(freq_vars))))
	system(paste("ln -sf ", freq_indels, file.path(output_dir, basename(freq_indels))))	
	system(paste("ln -sf ", star_bam, file.path(output_dir, basename(star_bam))))	
	system(paste("cp -f ", star_fusion, file.path(output_dir, basename(star_fusion))))
	
	sample_type = "Control"		
	if(!grepl("^CTR", short_sample_id) && !grepl("^SU2", short_sample_id))
	{
		sample_type = strsplit(short_sample_id, "\\.")[[1]][1]
		sample_type = gsub("\\d", "", sample_type)		
		sample_type = strsplit(sample_type, "\\-")[[1]][1]
	}else{
		if(grepl("CS", short_sample_id))
		{
			sample_type = "Spike"
		}
		if(grepl("^SU2", short_sample_id))
		{
			sample_type = "LDCT"
		}
	}		

	if(nrow(s2b[s2b$V1 == short_sample_id, ]) == 0)
	{
		stop(paste("No sample2barcode record for lane/file:", lane, "/", short_sample_id))
	}

	df = data.frame(Lane = lane, ID = cfrna_ID, Selector = as.character(s2b[s2b$V1 == short_sample_id, "V5"]),
		SampleType = sample_type, Freq = file.path(output_dir, basename(freq_vars)), Indel = file.path(output_dir, basename(freq_indels)),
		Fusion = file.path(output_dir, basename(star_fusion)))
	data = rbind(data, df)
		
}
data$Donor = unlist(lapply(as.character(data$ID), function(x) strsplit(x, ".", fixed = T)[[1]][1]))
data$Donor = unlist(lapply(as.character(data$Donor), function(x) strsplit(x, "-", fixed = T)[[1]][1]))

write.table(data, file.path(output_dir_root, "sample_info.txt"), sep = "\t", row.names = F)
