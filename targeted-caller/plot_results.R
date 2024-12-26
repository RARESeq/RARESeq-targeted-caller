library(ggplot2)
library(cowplot) 
library(reshape2)
library(viridis)
library(ComplexHeatmap)
source("targeted-caller/misc.R")

args = commandArgs(T)
input_dir = args[1]

output_dir_root = args[2]
output_dir = args[2]
dir.create(output_dir, recursive = T, showWarning = F)

samples = read.delim(file.path(input_dir, "sample_info.txt"))
training_samples = read.delim(file.path(input_dir, "background_files", "training_samples.txt"))
all_data = read.delim(file.path(output_dir, "variant_call_table.txt"))

lvls = unique(sort(as.character(all_data$SampleType)))
lvls = c("Training", "Control", lvls[!lvls %in% c("Training", "Control")])
all_data$SampleType = factor(as.character(all_data$SampleType), levels = lvls)

gene_calls = aggregate(Called ~ Lane+ID+SampleType+CallType+GENE+Selector, all_data, FUN = sum)
gene_calls$gene_variant = paste0(gene_calls$CallType, "--", gene_calls$GENE)

casted = dcast(gene_calls, Lane + ID + SampleType + Selector  ~ gene_variant, value.var = "Called")

colorbar_df = casted[,1:4]
mat = casted[,-c(1:4)]
mat[is.na(mat)] = 0
rownames(mat) = paste(colorbar_df$Lane, colorbar_df$ID, sep = "_")
rownames(colorbar_df) = rownames(mat)
mat = t(mat)
if(nrow(mat) == 0 || ncol(mat) == 0)
{
	stop("No mutation calls! Please make sure you provided enough samples to the pipeline.")	
}
type = unlist(lapply(rownames(mat), function(x) {strsplit(x, "--")[[1]][1]}))
genes = unique(unlist(lapply(rownames(mat), function(x) {strsplit(x, "--")[[1]][2]})))
mat_list= split(as.data.frame(mat), type)
mat_list = lapply(mat_list, function(x){
	rownames(x) = unlist(lapply(rownames(x), function(y) {strsplit(y, "--")[[1]][2]}))
	x =x[match(genes, rownames(x)),]
	rownames(x) = genes
	x[is.na(x)] = 0
	as.matrix(x)
	})

heatmap_color_annotation <- function(data, columns, palettes = NULL, dataset = "Carcinoma", name = "ann", which = "column", annotation_legend_param = list(), show_annotation_name = T)
{
	#print(columns)
	df = data[,match(columns, colnames(data)),drop = F]
	rownames(df) = rownames(data)
	colnames(df) = columns
	
	top_colors = list()
	for(col in columns)
	{
		if(!is.factor(df[,col]))
		{
			df[,col] = as.factor(as.character(df[,col]))
		}
		if(all(is.na(df[,col])))
		{
			df[,col] = as.factor(as.character(paste('random_', rep(sample(1:2, 1), nrow(df)))))
		}
		if(is.null(palettes))
		{
			type_cols = get_colors(length(levels(df[,col])), column = col, dataset = dataset)
		}else{
			type_cols = get_colors(length(levels(df[,col])), palette = palettes[[col]], dataset = dataset)
		}
		names(type_cols) = levels(df[,col])
		top_colors[[col]] = type_cols
	}
	if(which == "row")
	{
		annotation_name_side = "top"
	}else{
		annotation_name_side = "left"
	}
	HeatmapAnnotation(df = df, name = name, col = top_colors, annotation_name_side = annotation_name_side, annotation_legend_param = annotation_legend_param,
		show_annotation_name = show_annotation_name, which = which, border = T)
}

palettes = list(Selector = 1)

top_annotation = heatmap_color_annotation(colorbar_df, c("Selector"),
			palettes = palettes,
 			name = "top_ann", dataset = "", 
			show_annotation_name = T, #annotation_name_side = "left",
			annotation_legend_param = list(legend_direction = "horizontal", nrow = 1, by_row = F,
				title_position = "topcenter", title_gp = gpar(fontsize = 11)))

col = c(SNV = "#D6372E", Indel = "#5189BB", Fusion = "#70B460", METex14_skipping = "#985EA8")

tbt = as.data.frame(table(as.character(colorbar_df$SampleType)))
colorbar_df$SampleTypePretty = paste0(colorbar_df$SampleType, " (n=", tbt[match(colorbar_df$SampleType, tbt[,1]), 2], ")")
lvls = paste0(levels(colorbar_df$SampleType), " (n=", tbt[match(levels(colorbar_df$SampleType), tbt[,1]), 2], ")")
colorbar_df$SampleTypePretty = factor(colorbar_df$SampleTypePretty, levels = lvls)
	
h <- oncoPrint(mat_list,	
	column_split = colorbar_df$SampleTypePretty, 
	top_annotation = top_annotation, border = T,
	row_names_side = "left", pct_side = "right",
	row_names_gp = gpar(fontface = "italic"),
	width = unit(4, "in"), height = unit(1/7 * length(genes), "in"),
	column_title_rot = 90,
    alter_fun = list(
        SNV = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["SNV"], col = NA)),
        Fusion = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Fusion"], col = NA)),
        METex14_skipping = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["METex14_skipping"], col = NA)),
        Indel = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Indel"], col = NA)),
        background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#FEFEFE", col = NA))
    ), col = col)


pdf(file.path(output_dir, paste0("oncoprint.pdf")), height = 8, width = 15)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side  = "bottom")
dev.off()

png(file.path(output_dir, paste0("oncoprint.png")), height = 8, width = 15, units = "in", res = 200)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side  = "bottom")
dev.off()

#Debug outputs

stat = "AF"
ylab = "AF (%)"
all_data$DEPTH_DISCRETE = factor(as.character(all_data$DEPTH %/% 250 * 250), levels = as.character(unique(sort(all_data$DEPTH %/% 250 * 250))))
all_data$shape = ifelse(all_data$SampleType %in% c("Control", "Training"), as.character(all_data$SampleType), "Cancer") 
summary = aggregate(as.formula(paste(stat, "~DEPTH_DISCRETE+SampleType")), data = all_data, FUN = mean)
summarysd = aggregate(as.formula(paste(stat, "~DEPTH_DISCRETE+SampleType")), data = all_data, FUN = function(x) sqrt(var(x) / length(x)))
summary$Low = summary$AF - summarysd$AF
summary$High = summary$AF + summarysd$AF
pdf(file.path(output_dir, paste0("depth_vs_AF_aggregated.pdf")), height = 8, width = 8)	
g <- ggplot(summary,aes_string(x="DEPTH_DISCRETE", y=stat)) +  
	geom_bar(stat = "identity") + 
	geom_errorbar(aes(ymin=Low, ymax=High), width=.2) + 	
	theme_bw() + 	
	theme(aspect.ratio= 1/1)+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
	scale_color_manual(values = rev(get_colors(2, col = "State"))) +	
	labs(x = "Depth", y = ylab, color = "Called?")  + 
	facet_wrap(~SampleType, scales = "free", ncol = 6)
plot(g)
g <- ggplot(all_data,aes_string(x="DEPTH_DISCRETE", y=stat)) +  		
	geom_point(aes(color = Called, shape = shape)) + 
	geom_boxplot() +  	
	theme_bw() + 	
	theme(aspect.ratio= 1/1)+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
	scale_color_manual(values = rev(get_colors(2, col = "State"))) +	
	labs(x = "Depth", y = ylab, color = "Called?")  + 
	facet_wrap(~SampleType, scales = "free", ncol = 6)
plot(g)
dev.off()


subset_data = all_data[all_data$CallType == "SNV",]	
subset_data$variant = paste(subset_data$GENE, subset_data$POS, subset_data$VAR)
n = ceiling(sqrt(length(table(subset_data$variant))))
pdf(file.path(output_dir, paste0("depth_vs_AF.pdf")), height = n * 2, width = n * 2.5 * 2)

g <- ggplot(subset_data,aes_string(x="DEPTH", y=stat)) +  
	#geom_violin() +  
	geom_point(aes(color = Called, shape = shape), size = 3) + 
	
	theme_bw() + 
	#scale_y_continuous(limits = c(0, max(snvs_called$Calls, na.rm = T))) + 
	theme(aspect.ratio= 1/1.2)+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
	scale_color_manual(values = rev(get_colors(2, col = "State"))) +
	#stat_compare_means(hide.ns = F, label = "p.signif") + 
		#comparisons = list(c("0","0.5"), c("0","1"), c("0","2.5"), c("0","5"), c("0","10"), c("0","25"), c("0","50"), c("0","100")))+ 
	labs(x = "Depth", y = ylab, color = "Called?")  + 
	facet_wrap(~variant, scales = "free")		

g2 <- ggplot(subset_data,aes_string(x="SampleType", y=stat)) +  
	geom_violin() +  
	geom_jitter(aes(color = Called), width = 0.1, height =0, size = 3) + 
	theme_bw() + 
	#scale_y_continuous(limits = c(0, max(snvs_called$Calls, na.rm = T))) + 
	theme(aspect.ratio= 1/1.5)+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
	scale_color_manual(values = rev(get_colors(2, col = "State"))) +
	#stat_compare_means(hide.ns = F, label = "p.signif") + 
		#comparisons = list(c("0","0.5"), c("0","1"), c("0","2.5"), c("0","5"), c("0","10"), c("0","25"), c("0","50"), c("0","100")))+ 
	labs(x = "Sample type", y = ylab, color = "Called?")  + 
	facet_wrap(~variant, scales = "free")
plot(plot_grid(g2, g, ncol = 2))
dev.off()

