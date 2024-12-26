args = commandArgs(T)

seed = as.integer(args[1])
unique_traning_samples = args[2]
input_dir = args[3]
output_dir = args[4]
dir.create(output_dir, recursive =  T, showWarning = F)

set.seed(seed)

data = read.delim(file.path(input_dir, "sample_info.txt"))

samples = read.delim(unique_traning_samples)
samples$TMP = paste(samples$Lane, samples$ID)
data$TMP = paste(data$Lane, data$ID)

candidate_samples = data[order(data$ID),]
candidate_samples = candidate_samples[candidate_samples$TMP %in% samples$TMP,]
write.table(candidate_samples[,colnames(candidate_samples) != "TMP"], file.path(output_dir, "training_samples.txt"), row.names = F, sep = "\t")
 