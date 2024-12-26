library(reshape2)

args = commandArgs(T)
input_dir = args[1]
print(input_dir)
output_dir_root = args[2]
output_dir = args[2]
dir.create(output_dir, recursive = T, showWarning = F)

samples = read.delim(file.path(input_dir, "sample_info.txt"))
training_samples = read.delim(file.path(input_dir, "background_files",  "training_samples.txt"))

all_data = NULL
not_called = NULL
for(i in 1:nrow(samples))
{
	sample = samples[i,]
	print(paste(sample$Lane, sample$ID))
	res_dir = file.path(output_dir_root, "call_data", paste0(sample$Lane, "_", sample$ID))
	
	if(!file.exists(file.path(res_dir, "calls_SNV_idels.txt")))
	{
		not_called = rbind(not_called, sample[,1:3])
		next
	}
	df = read.delim(file.path(res_dir, "calls_SNV_idels.txt"), stringsAsFactors = F) 
	if(nrow(df) == 0)
	{
		not_called = rbind(not_called, sample[,1:3])
		next
	}

	df$Lane = sample$Lane
	df$ID = sample$ID
	df$SampleType = sample$SampleType	
	
	df$Selector = gsub("_withERCC.bed", "", sample$Selector)	
	df$Selector = gsub("_with_ERCC.bed", "", df$Selector)	
	df$Selector = gsub("_TWIST.bed", "", df$Selector)
	df$Selector = gsub("_", " ", df$Selector)

	if(sample$ID %in% training_samples$ID)
	{
		df$SampleType = "Training"
	}
	
	df$REF = ifelse(df$REF == TRUE, "T", as.character(df$REF))
	df$VAR = ifelse(df$VAR == TRUE, "T", as.character(df$VAR))
	
	df$AF = df$NR / df$DEPTH
	df$CallType = ifelse(grepl("-", df$VAR, fixed = T) | grepl("+", df$VAR, fixed = T), "Indel", "SNV")
	all_data = rbind(all_data, df)
}

for(i in 1:nrow(samples))
{
	sample = samples[i,]
	print(paste(sample$Lane, sample$ID))
	res_dir = file.path(output_dir_root, "call_data", paste0(sample$Lane, "_", sample$ID))
	
	if(!file.exists(file.path(res_dir, "calls_fusions.txt")))
	{
		not_called = rbind(not_called, sample[,1:3])
		next
	}
	df1 = read.delim(file.path(res_dir, "calls_fusions.txt"), stringsAsFactors = F) 
	
	if(nrow(df1) == 0)
	{
		not_called = rbind(not_called, sample[,1:3])
		next
	}
	if(nrow(df1) > 1)
	{
		df1 = df1[which.max(df1$JunctionReadCount),]
	}
	
	template = all_data[1:nrow(df1),1:which(colnames(all_data) == "BKG_MAX_AF")] 
	template[,] = NA
	template$GENE = paste(gsub("--", "/ ", df1[,1]))
	
	template$NR = df1$JunctionReadCount
	template$DEPTH = df1$JunctionReadCount + df1$SpanningFragCount
		
	df = template
	df$Lane = sample$Lane
	df$ID = sample$ID
	df$SampleType = sample$SampleType
	df$Selector = gsub("_withERCC.bed", "", sample$Selector)	
	df$Selector = gsub("_with_ERCC.bed", "", df$Selector)	
	df$Selector = gsub("_", " ", df$Selector)
	if(sample$ID %in% training_samples$ID)
	{
		df$SampleType = "Training"
	}
		
	df$AF = df1$FFPM
	df$CallType = "Fusion"
	all_data = rbind(all_data, df)
}

for(i in 1:nrow(samples))
{
	sample = samples[i,]
	print(paste(sample$Lane, sample$ID))
	res_dir = dirname(as.character(sample$Freq))
	
	if(!file.exists(file.path(res_dir, "called_snv.cfRNA_METex14_skipping.sam")))
	{
		not_called = rbind(not_called, sample[,1:3])
		next
	}

	df1 = read.delim(file.path(res_dir, "called_snv.cfRNA_METex14_skipping.sam"), stringsAsFactors = F, header = F) 
	
	if(nrow(df1) == 0)
	{
		not_called = rbind(not_called, sample[,1:3])
		next
	}
	if(df1[2,1] == 0)
	{
		next
	}
	template = all_data[1,1:which(colnames(all_data) == "BKG_MAX_AF")] 
	template[,] = NA
	template$GENE = "MET"
	template$NR = df1[2,1]
	template$DEPTH = df1[1,1]
	template$AF = template$NR * 1.0 / template$DEPTH
	
	df = template
	df$Lane = sample$Lane
	df$ID = sample$ID
	df$SampleType = sample$SampleType
	df$Selector = gsub("_withERCC.bed", "", sample$Selector)	
	df$Selector = gsub("_with_ERCC.bed", "", df$Selector)	
	df$Selector = gsub("_", " ", df$Selector)
	if(sample$ID %in% training_samples$ID)
	{
		df$SampleType = "Training"
	}
	
	df$CallType = "METex14_skipping"
	all_data = rbind(all_data, df)
}

write.table(all_data, file.path(output_dir, "variant_call_table.preliminary.txt"), sep = "\t", row.names = F)

subset_data = all_data[all_data$NR > 0 & all_data$CallType == "SNV" & all_data$BKG_MAX_AF > 0,]
subset_data$Q = p.adjust(subset_data$P_NR, method = "BH")
subset_data2 = all_data[!(all_data$NR > 0 & all_data$CallType == "SNV"),]
subset_data2$Q = 1
subset_data_3 = all_data[all_data$NR > 0 & all_data$CallType == "SNV" & all_data$BKG_MAX_AF == 0,]
subset_data_3$Q = 0

all_data = rbind(subset_data2, subset_data_3, subset_data)

tmp = paste(all_data$CHR, all_data$POS, all_data$REF, all_data$VAR)
splits = split(all_data, tmp)


all_data = do.call(rbind, lapply(splits, function(spl){
	if(spl$CallType[1] == "SNV")
	{		
		q_thres = 0.1
		spl$Called = spl$Q < q_thres
	}else{
		spl$Called = spl$CallType == "Fusion" | 
				  spl$CallType == "METex14_skipping" | 				  
				  grepl("-", spl$VAR, fixed = T) & nchar(spl$VAR) > 6 | 
				  grepl("+", spl$VAR, fixed = T) & nchar(spl$VAR) > 6
	}
	spl
	}))

write.table(all_data, file.path(output_dir, "variant_call_table.txt"), sep = "\t", row.names = F)
