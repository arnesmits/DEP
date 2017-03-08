single_prot_plot <- function(data, protein, type) {
  if(type == "centered") {
    df <- exprs(data) - fData(data)$mean
    df %<>% data.frame(.) %>% rownames_to_column(.)
    df %<>% filter(rowname == protein) %>% gather(ID, val, 2:ncol(.)) %>% left_join(., pData(data), by = "ID") %>% mutate(sample = gsub("_", " ", sample))
    df$replicate <- as.factor(df$replicate)
    p1 <- ggplot(df, aes(sample, val, col = replicate)) + geom_hline(yintercept = 0) + theme_bw() +
      stat_summary(fun.y = "mean", colour = "black", size = 0, geom = "bar", fill = "black") +
      geom_point(shape = 17, size = 4) + labs(title = unique(df$name), x = "Baits", y = "Enrichment (log2)") +
      theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 1), axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.title = element_text(size=14,face="bold"), legend.position="top")
  }
  if(type == "contrast") {
    df <- fData(data) %>% .[,grep("_diff$", colnames(.))] %>% rownames_to_column(.)
    colnames(df)[2:ncol(df)] %<>%  gsub("_diff", "", .) %>% gsub("[_]", " ", .) %>% gsub("[-]", "- ", .)
    df %<>% filter(rowname == protein) %>% gather(sample, LFC, 2:ncol(.))
    p1 <- ggplot(df, aes(sample, LFC)) + geom_hline(yintercept = 0) + theme_bw() +
      geom_bar(stat = "unique", size = 0, col = "black", fill = "black") + labs(title = unique(df$name), x = "Baits", y = "Enrichment (log2)") +
      theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 1), axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.title = element_text(size=14,face="bold"), legend.position="top")
  }
  p1
}

plot_heatmap <- function(data, type, k = 6, col_limit = 6) {
  data <- data[fData(data)$sign == "+"]

  if(type == "centered") {
    df <- exprs(data) - fData(data)$mean

    set.seed(1)
    kmeans <- kmeans(df,k)
    order <- df %>% data.frame() %>% cbind(., cluster = kmeans$cluster) %>% mutate(row = apply(., 1, function(x) max(x))) %>% group_by(cluster) %>% summarize(index=sum(row)/n()) %>% arrange(desc(index)) %>% collect %>% .[[1]] %>% match(seq(1:k),.)
    kmeans$cluster <- order[kmeans$cluster]
  }
  if(type == "contrast") {
    df <- fData(data)
    df <- df[,grep("_diff", colnames(df))]
    colnames(df) <- gsub("_diff","" , colnames(df))

    set.seed(1)
    kmeans <- kmeans(df,k)
    order <- cbind(df, cluster = kmeans$cluster) %>% gather(sample, diff, 1:(ncol(.)-1)) %>% group_by(cluster) %>% summarize(row = mean(diff)) %>% arrange(desc(row)) %>% collect %>% .[[1]] %>% match(seq(1:k),.)
    kmeans$cluster <- order[kmeans$cluster]
  }

  ht1 = Heatmap(df, col = colorRamp2(seq(-col_limit, col_limit, (col_limit/5)), rev(brewer.pal(11, "RdBu"))), split = kmeans$cluster,
                cluster_columns = T, cluster_rows = T, na_col = "grey", clustering_distance_rows = "pearson",clustering_distance_columns = "pearson",
                row_names_side = "left", column_names_side = "top", show_row_names = T, show_column_names = T,
                heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal", legend_width = unit(5, "cm"), title_position = "lefttop"),
                name = "Enrichment (Log2)", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 16))
  draw(ht1, heatmap_legend_side = "top")
}

volcano <- function(data, contrast, labelsize = 3, add_names = TRUE) {
  feat_data <- fData(data)
  diff <- grep(paste(contrast, "_diff", sep = ""), colnames(feat_data))
  padj <- grep(paste(contrast, "_p.adj", sep = ""), colnames(feat_data))
  sign <- grep(paste(contrast, "_sign", sep = ""), colnames(feat_data))
  feat_data$x <- feat_data[,diff]
  feat_data$y <- -log10(feat_data[,padj])
  feat_data$z <- feat_data[,sign]
  extra <- feat_data %>% filter(z == "+")
  name1 <- gsub("[.].*", "", contrast)
  name2 <- gsub(".*[.]", "", contrast)
  if(add_names) {
    p1 <- ggplot(feat_data, aes(x, y)) + geom_vline(xintercept = 0) + geom_point(col = "grey") + geom_point(data = extra, col = "black") + geom_text_repel(data = extra, aes(label = name), size = labelsize, point.padding = unit(0.3, "lines")) + theme_bw() +
      theme(legend.position="none", axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + labs(x = "Log2 Fold Change", y = "-log10 adjusted P value") +
      geom_text(data = data.frame(), aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1, 0), vjust = c(-1, -1), label = c(name1, name2), size = 5, fontface = "bold"))
  } else {
    p1 <- ggplot(feat_data, aes(x, y)) + geom_vline(xintercept = 0) + geom_point(col = "grey") + geom_point(data = extra, col = "black") + theme_bw() +
      theme(legend.position="none", axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + labs(x = "Log2 Fold Change", y = "-log10 adjusted P value") +
      geom_text(data = data.frame(), aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1, 0), vjust = c(-1, -1), label = c(name1, name2), size = 5, fontface = "bold"))
  }
  p1
}

plot_norm <- function(raw, norm) {
  df1 <- exprs(raw) %>% data.frame() %>% rownames_to_column(.) %>% gather(ID, val, 2:ncol(.)) %>% left_join(., pData(raw), by = "ID") %>% mutate(var = "original")
  df2 <- exprs(norm) %>% data.frame() %>% rownames_to_column(.) %>% gather(ID, val, 2:ncol(.)) %>% left_join(., pData(norm), by = "ID") %>% mutate(var = "normalized")
  df <- rbind(df1, df2)
  ggplot(df, aes(x = ID, y = val, fill = sample)) + geom_boxplot(notch = T) + coord_flip() + facet_wrap(~var, ncol = 1) + theme_bw() + labs(x = "", y = "Log2 Intensity") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.title = element_text(size=14,face="bold"), legend.position="right")
}

plot_missval <- function(data) {
  df <- exprs(data) %>% data.frame(.)
  missval <- df[apply(df, 1, function(x) any(is.na(x))),]
  missval <- ifelse(is.na(missval), 0, 1)
  ht2 = Heatmap(missval, col = c("white","black"), row_names_side = "left", column_names_side = "top", show_row_names = F, show_column_names = T,
                name = "Missing values (0)", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 16))
  draw(ht2, heatmap_legend_side = "top")
}

plot_numbers <- function(data) {
  df <- exprs(data) %>% data.frame() %>% rownames_to_column() %>% gather(ID, bin, 2:ncol(.)) %>% mutate(bin = ifelse(is.na(bin), 0, 1))
  stat <- df %>% group_by(ID) %>% summarize(n = n(), sum = sum(bin)) %>% left_join(., pData(data), by = "ID")
  ggplot(stat, aes(x = ID, y = sum, fill = sample)) + geom_histogram(stat = "identity") + theme_bw() + labs(title = "Proteins per sample", x = "", y = "Number of ProteinGroups") + geom_hline(yintercept = unique(stat$n)) +
    theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 1, size=10), axis.title=element_text(size=12,face="bold"), legend.text=element_text(size=10), legend.title = element_text(size=12,face="bold"), legend.position="right")
}

plot_frequency <- function(data) {
  df <- exprs(data) %>% data.frame() %>% rownames_to_column() %>% gather(ID, bin, 2:ncol(.)) %>% mutate(bin = ifelse(is.na(bin), 0, 1))
  stat <- df %>% group_by(rowname) %>% summarize(sum = sum(bin))
  table <- table(stat$sum) %>% data.frame()
  ggplot(table, aes(x = Var1, y = Freq)) + geom_histogram(stat = "identity") + theme_bw() + labs(title = "Protein identifications overlap between samples", x = "Identified in number of samples", y = "Number of ProteinGroups") +
    theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 1, size=12), axis.title=element_text(size=12,face="bold"), legend.text=element_text(size=12), legend.title = element_text(size=12,face="bold"), legend.position="right")

}
