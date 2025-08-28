#' Plot synteny comparison between two genomes using GFF and BLAST
#'
#' @param query_gff Path to query GFF file
#' @param subject_gff Path to subject GFF file
#' @param blast_tsv Path to BLAST tabular file (outfmt 6)
#' @param output Path to output PNG file, or NULL to skip saving
#' @param target_genes Vector of gene labels to highlight
#' @param target_labels Vector of labels to display
#' @param target_colors Vector of colors for highlighted BLAST bands
#' @param query_color Color for query gene arrows
#' @param subject_color Color for subject gene arrows
#' @param query_label Label to display at the left of the query line
#' @param subject_label Label to display at the left of the subject line
#' @param plot_width Width of output PNG (default 10)
#' @param plot_height Height of output PNG (default 6)
#'
#' @return Saves a PNG plot (if output is specified) and returns ggplot object
#' @export
VisBlast <- function(query_gff,
                     subject_gff,
                     blast_tsv,
                     output = "synteny.png",
                     target_genes = c("geneA"),
                     target_labels = c("geneA"),
                     target_colors = c("#ec6800"),
                     query_color = "#028760",
                     subject_color = "#028760",
                     query_label = "Query",
                     subject_label = "Subject",
                     plot_width = 10,
                     plot_height = 6) {

  # 依存パッケージ自動インストール
  packages <- c("ggplot2", "dplyr", "grid")
  for(p in packages){
    if(!requireNamespace(p, quietly = TRUE)){
      install.packages(p, dependencies = TRUE)
    }
    library(p, character.only = TRUE)
  }

  # ファイル存在チェック
  if(!file.exists(query_gff)) stop("query_gff file does not exist.")
  if(!file.exists(subject_gff)) stop("subject_gff file does not exist.")
  if(!file.exists(blast_tsv)) stop("blast_tsv file does not exist.")

  # 引数チェック
  len_genes <- length(target_genes)
  if(length(target_labels) != len_genes){
    warning("target_genes and target_labels have different lengths. Labels will be matched by position where possible.")
    target_labels <- rep(target_labels[1], len_genes)
  }
  if(length(target_colors) != len_genes){
    warning("target_genes and target_colors have different lengths. Colors will be matched by position where possible.")
    target_colors <- rep(target_colors[1], len_genes)
  }

  # GFF読み込み
  query_gff <- tryCatch(
    read.table(query_gff, sep="\t", header=FALSE, stringsAsFactors=FALSE),
    error = function(e) stop("Failed to read query GFF: ", e$message)
  )
  subject_gff <- tryCatch(
    read.table(subject_gff, sep="\t", header=FALSE, stringsAsFactors=FALSE),
    error = function(e) stop("Failed to read subject GFF: ", e$message)
  )
  colnames(query_gff) <- colnames(subject_gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

  query_genes <- query_gff %>% filter(type=="gene") %>% mutate(y_center=100)
  subject_genes <- subject_gff %>% filter(type=="gene") %>% mutate(y_center=0)

  if(nrow(query_genes)==0) warning("No genes found in query GFF.")
  if(nrow(subject_genes)==0) warning("No genes found in subject GFF.")

  # BLAST読み込み
  blast <- tryCatch(
    read.table(blast_tsv, sep="\t", header=FALSE),
    error = function(e) stop("Failed to read BLAST TSV: ", e$message)
  )
  if(nrow(blast)==0){
    warning("BLAST file is empty. No synteny bands will be drawn.")
    blast <- data.frame(qseqid=character(0), sseqid=character(0), pident=numeric(0),
                        length=numeric(0), mismatch=numeric(0), gapopen=numeric(0),
                        qstart=numeric(0), qend=numeric(0), sstart=numeric(0),
                        send=numeric(0), evalue=numeric(0), bitscore=numeric(0))
  }
  colnames(blast) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                       "qstart","qend","sstart","send","evalue","bitscore")

  # GFF遺伝子領域に重なるBLASTだけ抽出
  target_blast <- tryCatch({
    do.call(rbind, lapply(1:nrow(query_genes), function(i){
      gene <- query_genes[i, ]
      hits <- blast %>% filter(
        qseqid == gene$seqid &
          !(qend < gene$start | qstart > gene$end)
      )
      if(nrow(hits) > 0) hits
    }))
  }, error=function(e){data.frame()})

  # 矢印データ作成（色固定）
  arrow_df <- function(df, y_center, color_default, target_genes, target_labels){
    if(nrow(df)==0) return(data.frame())
    n <- nrow(df)
    labels <- rep("", n)
    for(i in seq_len(n)){
      idx <- match(df$attributes[i], target_genes)
      if(!is.na(idx) && idx <= length(target_labels)){
        labels[i] <- target_labels[idx]
      }
    }
    df %>% mutate(
      xmin = start,
      xmax = end,
      y_center = y_center,
      direction = ifelse(strand == "+", 1, -1),
      color = color_default,
      label = labels,
      xstart = ifelse(direction==1, xmin, xmax),
      xend   = ifelse(direction==1, xmax, xmin)
    )
  }

  query_arrow <- arrow_df(query_genes, 100, query_color, target_genes, target_labels)
  subject_arrow <- arrow_df(subject_genes, 0, subject_color, target_genes, target_labels)

  label_df <- query_arrow %>% filter(label != "")

  # 帯データ作成（マルチカラー対応、矢印色は変更せず）
  if(nrow(target_blast) > 0){
    band_df <- do.call(rbind, lapply(1:nrow(target_blast), function(i){
      query_match <- query_genes %>% filter(start <= target_blast$qend[i] & end >= target_blast$qstart[i])
      attr_value <- if(nrow(query_match) > 0) query_match$attributes[1] else NA
      idx <- match(attr_value, target_genes)
      is_target <- !is.na(attr_value) && attr_value %in% target_genes
      color <- if(is_target) target_colors[idx] else "#7b7c7d"
      label <- if(is_target) target_labels[idx] else NA
      data.frame(
        x = c(target_blast$qstart[i], target_blast$qend[i], target_blast$send[i], target_blast$sstart[i]),
        y = c(100,100,0,0),
        group = i,
        label = label,
        color = color
      )
    }))
  } else {
    band_df <- data.frame(x=numeric(0), y=numeric(0), group=integer(0), label=character(0), color=character(0))
  }

  # 背景領域
  xmin_all <- min(c(target_blast$qstart, target_blast$qend,
                    target_blast$sstart, target_blast$send,
                    query_gff$start, query_gff$end,
                    subject_gff$start, subject_gff$end), na.rm = TRUE)
  xmax_all <- max(c(target_blast$qstart, target_blast$qend,
                    target_blast$sstart, target_blast$send,
                    query_gff$start, query_gff$end,
                    subject_gff$start, subject_gff$end), na.rm = TRUE)
  query_bg <- data.frame(y=100)
  subject_bg <- data.frame(y=0)

  # プロット作成
  gg <- ggplot() +
    geom_rect(data=query_bg, aes(xmin=xmin_all, xmax=xmax_all, ymin=y, ymax=y+2), fill="#595857") +
    geom_rect(data=subject_bg, aes(xmin=xmin_all, xmax=xmax_all, ymin=y, ymax=y+2), fill="#595857") +
    geom_text(aes(x = xmin_all - 10, y = 110, label = query_label), hjust = 1, vjust = 0, size = 5, fontface="bold") +
    geom_text(aes(x = xmin_all - 10, y = 10,   label = subject_label), hjust = 1, vjust = 0, size = 5, fontface="bold") +
    geom_polygon(data=band_df, aes(x=x, y=y, group=group, fill=color), alpha=0.8) +
    geom_segment(data=query_arrow, aes(x=xstart, xend=xend, y=y_center+1, yend=y_center+1, color=color),
                 arrow=arrow(length=unit(0.4,"cm"), type="closed"), size=2) +
    geom_segment(data=subject_arrow, aes(x=xstart, xend=xend, y=y_center+1, yend=y_center+1, color=color),
                 arrow=arrow(length=unit(0.4,"cm"), type="closed"), size=2) +
    geom_text(data=label_df, aes(x=(xmin+xmax)/2, y=y_center+8, label=label),
              color="red", size=4, fontface="bold") +
    scale_fill_identity() +
    scale_color_identity() +
    xlab("Position (bp)") + ylab("") +
    theme_minimal() +
    theme(panel.grid=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

  # ファイル保存は NULL でない場合のみ
  if(!is.null(output)){
    out_dir <- dirname(output)
    if(!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)
    ggsave(output, gg, width = plot_width, height = plot_height)
  }

  return(gg)
}
