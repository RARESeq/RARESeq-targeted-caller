library(RColorBrewer)
get_colors <- function(count, column = NULL, palette = NULL, dataset = "Carcinoma")
{
	if(!is.null(column))
	{
		if(column %in% c("Cluster", "State"))
		{
			palette = 1
		}
		if(column %in% c("Tissue", "SourceTissue"))
		{
			palette = 2
		}
		if(column %in% c("Dataset"))
		{
			palette = 3
		}
		if(column %in% c("acronym", "Histology"))
		{
			palette = 4
		}
		if(column %in% c("CellType", "Cell Type",  "Cell type", "CallType"))
		{
			palette = 5
		}
		if(column %in% c("Filtered"))
		{
			palette = 6
		}
		if(column %in% c("CAF", "InGroundTruth"))
		{
			palette = 7
		}
		if(tolower(column) %in% c("metacluster", "ecotype"))
		{
			palette = 8 
		}
		if(tolower(column) %in% c("coo"))
		{
			palette = 9
		}
	}
	
	if(is.null(palette))
	{
		palette = 8
	}

	if(palette == 1)
	{
		cols = c(brewer.pal(9, "Set1"), brewer.pal(12, "Paired")[-c(2, 4, 6, 8, 10, 11, 12)], brewer.pal(8, "Dark2")[-c(5,7,8)] , brewer.pal(8, "Set3"))
		cols[6] = "gold"
		a1 <- col2rgb(cols)
		a2 <- rgb2hsv(a1)
		cols=hsv(a2[1,] * 0.93, a2[2,], a2[3,])
		if(dataset == "Lymphoma")
		{
			cols = c("#EB7D5B", "#FED23F", "#B5D33D","#6CA2EA", "#442288", cols) 
		}
	}else{
		if(palette == 2)
		{
			cols = rep(c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"), 50)	
		}else{
			if(palette == 3)
			{
				cols = rep(c(colorRampPalette(brewer.pal(8, "Accent"))(8)[c(3,5,6,4,7,8,1,2)], colorRampPalette(brewer.pal(8, "Set2"))(8), colorRampPalette(brewer.pal(8, "Pastel1"))(8)), 10000)
				cols[6] = "lightblue"
			}else{ 
				if(palette == 4) 
				{
					cols = rep(c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2")[-c(5,6,7)] , rev(brewer.pal(8, "Set2")), rev(brewer.pal(8, "Accent"))), 50)
				}else{
					if(palette == 5) 
					{
						cols = rep(rev(c(brewer.pal(8, "Set3"), brewer.pal(9, "Set1"), brewer.pal(12, "Paired")[-c(2, 4, 6, 8, 10, 11, 12)], brewer.pal(8, "Dark2")[-c(5,7,8)])), 50)
					}else{
						if(palette == 6) 
						{
							cols = rev(c("#016450", "#fff7fb"))
						}else{
							if(palette == 7) 
							{
								cols = rev(c("blue", "hotpink4","salmon", "darkorchid", "deeppink1", "lightpink","orange2", "aquamarine", "brown1", "aquamarine4"))
							}else{
								if(palette == 8)  
								{
									cols = rep(c(brewer.pal(12, "Paired")[-c(2, 4, 6, 8, 10, 11, 12)], brewer.pal(8, "Dark2")[-c(5,7,8)] , brewer.pal(8, "Set3")), 50)
									cols[6] = "gold"
									aux = cols[10]
									cols[10] = cols[2]
									cols[2] = aux
									aux = cols[5]
									cols[5] = cols[9]
									cols[9] = aux
									aux = cols[3]
									cols[3] = cols[8]
									cols[8] = aux
									aux = cols[3]
									cols[3] = cols[7]
									cols[7] = aux
									aux = cols[4]
									cols[4] = cols[7]
									cols[7] = aux
									cols_copy = cols 

									cols = c(rgb(214, 55, 46, maxColorValue = 255), rgb(250, 221, 75, maxColorValue = 255), rgb(112, 180, 96, maxColorValue = 255), rgb(230, 144, 193, maxColorValue = 255), rgb(152, 94, 168, maxColorValue = 255), rgb(163, 163, 163, maxColorValue = 255), rgb(183, 211, 229, maxColorValue = 255), rgb(230, 216, 194, maxColorValue = 255), rgb(240, 143, 53, maxColorValue = 255), rgb(81, 137, 187, maxColorValue = 255))
									cols = c(cols, cols_copy)
									cols = cols[-c(11:17)]
									if(dataset == "Lymphoma")
									{
										cols = c("#3E2883", "#76A2E4", "#8E549F", "#458833", "#BBD058", "#FFFD61", "#F8D25D", "#F3A83B", "#EC5428", cols)
									}
								}else{
									if(palette == 9)  
									{
										cols = rep(c("#73BBCE", "#DD7A1E", "#7F7F7F"), 100)
									}else{
										if(palette == 10)  
										{
											cols = brewer.pal(4, "YlGnBu")
										}else{
											cols = rep(c(colorRampPalette(brewer.pal(8, "Accent"))(limit), colorRampPalette(brewer.pal(8, "Set2"))(limit), colorRampPalette(brewer.pal(8, "Pastel1"))(limit)), 10000)
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	cols[1:count] 
}
