args = commandArgs(T)

seed = as.integer(args[1])
unique_traning_samples = as.integer(args[2])
input_dir = args[3]
output_dir = args[4]
dir.create(output_dir, recursive =  T, showWarning = F)

set.seed(seed)

data = read.delim(file.path(input_dir, "sample_info.txt"))
 
samples = data[order(data$ID),]
#candidate_samples = samples[samples$SampleType == "Control" & samples$Selector == "CFRNA_V1_withERCC.bed",]
candidate_samples = samples[samples$SampleType == "Control",]

control_set = unique(sort(as.character(candidate_samples$Donor)))
training_controls = sample(control_set, unique_traning_samples)

candidate_samples = candidate_samples[candidate_samples$Donor %in% training_controls,]
write.table(candidate_samples, file.path(output_dir, "training_samples.txt"), row.names = F, sep = "\t")
 