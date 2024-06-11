# Fig S1A ------
prepareDataFromScRNA <- function (object = NULL, diffData = NULL, showAverage = TRUE, 
                                  cells = NULL, group.by = "ident", assays = "RNA", slot = "data", 
                                  scale.data = TRUE, cluster.order = NULL, keep.uniqGene = TRUE, 
                                  sep = "_") 
{
  markerGene <- unique(diffData$gene)
  if (showAverage == TRUE) {
    mean_gene_exp <- Seurat::AverageExpression(object, features = markerGene, 
                                               group.by = group.by, assays = assays, slot = slot) %>% 
      data.frame() %>% as.matrix()
    name1 <- gsub(pattern = paste0(assays, ".", sep = ""), 
                  replacement = "", colnames(mean_gene_exp))
    colnames(mean_gene_exp) <- gsub(pattern = "\\.", replacement = " ", 
                                    name1)
    colnames(mean_gene_exp) <- levels(Seurat::Idents(object))
    if (scale.data == TRUE) {
      mean_gene_exp <- t(scale(t(mean_gene_exp)))
    }
    if (!is.null(cluster.order)) {
      mean_gene_exp <- mean_gene_exp[, cluster.order]
    }
    geneMode = "average"
  }
  else {
    cell.order <- data.frame(cell.id = names(Seurat::Idents(object)), 
                             cell.ident = Seurat::Idents(object))
    if (is.null(cluster.order)) {
      cell.order$cell.ident <- factor(cell.order$cell.ident, 
                                      levels = levels(Seurat::Idents(object)))
    }
    else {
      cell.order$cell.ident <- factor(cell.order$cell.ident, 
                                      levels = cluster.order)
    }
    cell.order <- cell.order[order(cell.order$cell.ident), 
                             ]
    getassy <- Seurat::GetAssayData(object = object, slot = slot)[features = markerGene, 
                                                                  cells = NULL, drop = FALSE] %>% as.matrix()
    id.order <- match(cell.order$cell.id, colnames(getassy))
    getassy <- getassy[, id.order]
    colnames(getassy) <- paste(colnames(getassy), cell.order$cell.ident, 
                               sep = "|")
    mean_gene_exp <- getassy
    if (scale.data == TRUE) {
      mean_gene_exp <- t(scale(t(mean_gene_exp)))
    }
    geneMode = "all"
  }
  merMat <- data.frame(mean_gene_exp, check.names = FALSE) %>% 
    tibble::rownames_to_column(., var = "gene")
  cinfo.gene <- diffData[, c("cluster", "gene")]
  cn <- unique(cinfo.gene$cluster)
  wide.res <- purrr::map_df(seq_along(cn), function(x) {
    tmp <- cinfo.gene[which(cinfo.gene$cluster == cn[x]), 
                      ]
    tmp2 <- merMat[which(merMat$gene %in% tmp$gene), ] %>% 
      dplyr::mutate(cluster = as.character(x))
    return(tmp2)
  })
  if (keep.uniqGene == TRUE) {
    wide.res <- wide.res %>% dplyr::distinct(., gene, .keep_all = TRUE)
    geneType = paste("unique", sep, sep = "|")
  }
  else {
    wide.res <- wide.res %>% dplyr::mutate(., gene = make.unique(gene, 
                                                                 sep = sep))
    geneType = paste("nounique", sep, sep = "|")
  }
  df <- reshape2::melt(wide.res, id.vars = c("cluster", "gene"), 
                       variable.name = "cell_type", value.name = "norm_value")
  df$cluster_name <- paste("cluster ", df$cluster, sep = "")
  if (showAverage == FALSE) {
    df$cell_type <- sapply(strsplit(as.character(df$cell_type), 
                                    split = "\\|"), "[", 2)
  }
  cl.info <- data.frame(table(wide.res$cluster)) %>% dplyr::mutate(Var1 = as.numeric(as.character(Var1))) %>% 
    dplyr::arrange(Var1)
  id <- unique(df$cluster_name)
  df <- purrr::map_df(seq_along(id), function(x) {
    tmp <- df %>% dplyr::filter(cluster_name == id[x])
    tmp %>% dplyr::mutate(cluster_name = paste(cluster_name, 
                                               " (", cl.info$Freq[x], ")", sep = ""))
  })
  df$cluster_name <- factor(df$cluster_name, levels = paste("cluster ", 
                                                            cl.info$Var1, " (", cl.info$Freq, ")", sep = ""))
  return(list(wide.res = wide.res, long.res = df, type = "scRNAdata", 
              geneMode = geneMode, geneType = geneType))
}
ct.markers.11ct.top50 <- ct.markers.11ct.top50[c(1:450,501:550,451:500),]

sn.data <- prepareDataFromScRNA(combined.sct,
                                diffData = ct.markers.11ct,
                                showAverage = T,
                                group.by = 'large_ct',keep.uniqGene = F,
                                scale.data = T)

tmp.vis <-function (object = NULL, ht.col.list = list(col_range = c(-2, 
                                                                    0, 2), col_color = c("#08519C", "white", "#A50F15")), border = TRUE, 
                    plot.type = c("line", "heatmap", "both"), ms.col = c("#0099CC", "grey90", "#CC3333"), line.size = 0.1, line.col = "grey90",  
                    add.mline = TRUE, mline.size = 2, mline.col = "#CC3333", 
                    ncol = 4, ctAnno.col = NULL, set.md = "median", textbox.pos = c(0.5, 0.8), textbox.size = 8, panel.arg = c(2, 0.25, 4, "grey90", 
                                                                                                                          NA), ggplot.panel.arg = c(2, 0.25, 4, "grey90", NA), 
                    annoTerm.data = NULL, annoTerm.mside = "right", termAnno.arg = c("grey95", "grey50"), add.bar = FALSE, bar.width = 8, textbar.pos = c(0.8, 
                                                                                             0.8), go.col = NULL, go.size = NULL, by.go = "anno_link", 
                    annoKegg.data = NULL, annoKegg.mside = "right", keggAnno.arg = c("grey95",  "grey50"), add.kegg.bar = FALSE, kegg.col = NULL, kegg.size = NULL, 
                    by.kegg = "anno_link", word_wrap = TRUE, add_new_line = TRUE, 
                    add.box = FALSE, boxcol = NULL, box.arg = c(0.1, "grey50"), 
                    add.point = FALSE, point.arg = c(19, "orange", "orange", 
                                                     1), add.line = TRUE, line.side = "right", markGenes = NULL, 
                    markGenes.side = "right", genes.gp = c("italic", 10, NA), 
                    term.text.limit = c(10, 18), mulGroup = NULL, lgd.label = NULL, 
                    show_row_names = FALSE, subgroup.anno = NULL, annnoblock.text = TRUE, 
                    annnoblock.gp = c("white", 8), add.sampleanno = TRUE, sample.group = NULL, 
                    sample.col = NULL, sample.order = NULL, cluster.order = NULL, 
                    sample.cell.order = NULL, HeatmapAnnotation = NULL, column.split = NULL, 
                    cluster_columns = FALSE, pseudotime_col = NULL, gglist = NULL, 
                    ...) 
{
  ComplexHeatmap::ht_opt(message = FALSE)
  if (is.null(ht.col.list[["col_range"]])) {
    col_range = c(-2, 0, 2)
  }
  else {
    col_range = ht.col.list[["col_range"]]
  }
  if (is.null(ht.col.list[["col_color"]])) {
    col_color = c("#08519C", "white", "#A50F15")
  }
  else {
    col_color = ht.col.list[["col_color"]]
  }
  col_fun = circlize::colorRamp2(col_range, col_color)
  plot.type <- match.arg(plot.type)
  if (plot.type == "line") {
    data <- data.frame(object$long.res)
    if (!is.null(sample.order)) {
      data$cell_type <- factor(data$cell_type, levels = sample.order)
    }
    line <- ggplot2::ggplot(data, ggplot2::aes(x = cell_type, 
                                               y = norm_value))
    if (object$type == "mfuzz") {
      line <- line + ggplot2::geom_line(ggplot2::aes(color = membership, 
                                                     group = gene), size = line.size) + ggplot2::scale_color_gradient2(low = ms.col[1], 
                                                                                                                       mid = ms.col[2], high = ms.col[3], midpoint = 0.5)
    }
    else {
      line <- line + ggplot2::geom_line(ggplot2::aes(group = gene), 
                                        color = line.col, size = line.size)
    }
    if (add.mline == TRUE) {
      if (object$type == "wgcna") {
        linec <- unique(data$modulecol)
        names(linec) <- linec
        line <- line + ggplot2::geom_line(stat = "summary", 
                                          fun = "median", size = mline.size, ggplot2::aes(group = 1, 
                                                                                          color = modulecol)) + ggplot2::scale_color_manual(values = linec)
      }
      else {
        line <- line + ggplot2::geom_line(stat = "summary", 
                                          fun = "median", colour = mline.col, size = mline.size, 
                                          ggplot2::aes(group = 1))
      }
    }
    else {
      line <- line
    }
    line1 <- line + ggplot2::theme_classic(base_size = 14) + 
      ggplot2::ylab("Normalized expression") + ggplot2::xlab("") + 
      ggplot2::theme(axis.ticks.length = ggplot2::unit(0.1, 
                                                       "cm"), axis.text.x = ggplot2::element_text(angle = 45, 
                                                                                                  hjust = 1, color = "black"), strip.background = ggplot2::element_blank()) + 
      ggplot2::facet_wrap(~cluster_name, ncol = ncol, 
                          scales = "free")
    return(line1)
  }
  else {
    data <- data.frame(object$wide.res, check.names = FALSE) %>% 
      dplyr::arrange(as.numeric(as.character(cluster)))
    if (object$type == "mfuzz") {
      mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
        dplyr::select(-gene, -cluster, -membership)
    }
    else if (object$type == "wgcna") {
      mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
        dplyr::select(-gene, -cluster, -modulecol)
    }
    else if (object$type == "scRNAdata") {
      mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
        dplyr::select(-gene, -cluster)
    }
    else if (object$type == "monocle") {
      mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
        dplyr::select(-gene, -cluster)
    }
    else {
      mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
        dplyr::select(-gene, -cluster)
    }
    rownames(mat) <- data$gene
    if (object$geneMode == "all" | ncol(mat) > 20) {
      use_raster = TRUE
    }
    else {
      use_raster = FALSE
    }
    if (!is.null(sample.order)) {
      mat <- mat[, sample.order]
    }
    cl.info <- data.frame(table(data$cluster)) %>% dplyr::mutate(Var1 = as.numeric(as.character(Var1))) %>% 
      dplyr::arrange(Var1)
    cluster.num <- nrow(cl.info)
    subgroup <- lapply(1:nrow(cl.info), function(x) {
      nm <- rep(as.character(cl.info$Var1[x]), cl.info$Freq[x])
      paste("C", nm, sep = "")
    }) %>% unlist()
    if (!is.null(cluster.order)) {
      subgroup <- factor(subgroup, levels = paste("C", 
                                                  cluster.order, sep = ""))
      cluster_row_slices = FALSE
    }
    else {
      cluster_row_slices = TRUE
    }
    if (object$geneMode == "all" & object$type == "scRNAdata") {
      celltype <- sapply(strsplit(colnames(mat), split = "\\|"), 
                         "[", 2)
      cell.num.info <- table(celltype)[unique(celltype)]
      if (is.null(sample.cell.order)) {
        column_split = factor(rep(names(cell.num.info), 
                                  cell.num.info), levels = unique(celltype))
      }
      else {
        column_split = factor(rep(names(cell.num.info), 
                                  cell.num.info), levels = sample.cell.order)
      }
      if (is.null(sample.col)) {
        block.col = 1:length(cell.num.info)
      }
      else {
        block.col = sample.col
      }
    }
    else {
      if (is.null(sample.group)) {
        sample.info = colnames(mat)
        if (is.null(HeatmapAnnotation)) {
          if (object$geneType == "branched") {
            if (ncol(mat) == 200) {
              column_split = rep(c("branch1", "branch2"), 
                                 each = 100)
            }
            else {
              column_split = rep(levels(object$pseudotime), 
                                 rep(100, ncol(mat)/100))
            }
          }
          else {
            column_split = NULL
          }
        }
        else {
          column_split = column.split
        }
      }
      else {
        sample.info = sample.group
        column_split = sample.group
      }
      if (object$type != "monocle") {
        sample.info <- factor(sample.info, levels = unique(sample.info))
        if (is.null(sample.col)) {
          scol <- circlize::rand_color(n = length(sample.info))
          names(scol) <- sample.info
        }
        else {
          scol <- sample.col
          names(scol) <- sample.info
        }
      }
      else {
        sample.info <- factor(object$pseudotime, levels = unique(object$pseudotime))
        if (is.null(pseudotime_col)) {
          if (object$geneType == "branched") {
            if (length(unique(object$pseudotime)) == 
                3) {
              pseudotime_col <- c("red", "grey80", "blue")
            }
            else {
              pseudotime_col <- circlize::rand_color(n = length(unique(object$pseudotime)))
            }
          }
          else {
            pseudotime_col <- c("blue", "red")
          }
        }
        else {
          pseudotime_col <- pseudotime_col
        }
        if (is.null(sample.col)) {
          if (object$type != "monocle") {
            scol <- circlize::rand_color(n = length(sample.info))
            names(scol) <- sample.info
          }
          else {
            if (object$geneType == "branched") {
              if (length(unique(object$pseudotime)) == 
                  3) {
                scol <- rep(pseudotime_col, table(object$pseudotime)[unique(object$pseudotime)])
                names(scol) <- sample.info
              }
              else {
                scol <- rep(pseudotime_col, table(object$pseudotime)[unique(object$pseudotime)])
                names(scol) <- sample.info
              }
            }
            else {
              scol <- (grDevices::colorRampPalette(pseudotime_col))(100)
              names(scol) <- sample.info
            }
          }
        }
        else {
          scol <- sample.col
          names(scol) <- sample.info
        }
      }
    }
    if (add.sampleanno == TRUE) {
      if (object$geneMode == "all" & object$type == "scRNAdata") {
        topanno = ComplexHeatmap::HeatmapAnnotation(cluster = ComplexHeatmap::anno_block(gp = grid::gpar(fill = block.col), 
                                                                                         labels = NULL), show_annotation_name = FALSE)
      }
      else {
        if (is.null(HeatmapAnnotation)) {
          topanno = ComplexHeatmap::HeatmapAnnotation(sample = sample.info, 
                                                      col = list(sample = scol), gp = grid::gpar(col = ifelse(object$type == 
                                                                                                                "monocle", NA, "white")), show_legend = ifelse(object$type == 
                                                                                                                                                                 "monocle", FALSE, TRUE), show_annotation_name = FALSE)
        }
        else {
          topanno = HeatmapAnnotation
        }
      }
    }
    else {
      topanno = NULL
    }
    if (is.null(ctAnno.col)) {
      colanno <- jjAnno::useMyCol("stallion", n = cluster.num)
    }
    else {
      colanno <- ctAnno.col
    }
    names(colanno) <- 1:cluster.num
    align_to = split(1:nrow(mat), subgroup)
    anno.block <- ComplexHeatmap::anno_block(align_to = align_to, 
                                             panel_fun = function(index, nm) {
                                               npos = as.numeric(unlist(strsplit(nm, split = "C"))[2])
                                               grid::grid.rect(gp = grid::gpar(fill = colanno[npos], 
                                                                               col = NA))
                                               if (annnoblock.text == TRUE) {
                                                 grid::grid.text(label = paste("n:", length(index), 
                                                                               sep = ""), rot = 90, gp = grid::gpar(col = annnoblock.gp[1], 
                                                                                                                    fontsize = as.numeric(annnoblock.gp[2])))
                                               }
                                             }, which = "row")
    if (!is.null(markGenes)) {
      rowGene <- rownames(mat)
      annoGene <- markGenes
      gene.col <- data %>% dplyr::select(gene, cluster) %>% 
        dplyr::filter(gene %in% annoGene)
      gene.col <- purrr::map_df(1:cluster.num, function(x) {
        tmp <- gene.col %>% dplyr::filter(as.numeric(cluster) == 
                                            x) %>% dplyr::mutate(col = colanno[x])
      })
      gene.col <- gene.col[match(annoGene, gene.col$gene), 
      ]
      if (is.na(genes.gp[3])) {
        gcol = gene.col$col
      }
      else {
        gcol = genes.gp[3]
      }
      index <- match(annoGene, rowGene)
      geneMark = ComplexHeatmap::anno_mark(at = index, 
                                           labels = annoGene, which = "row", side = markGenes.side, 
                                           labels_gp = grid::gpar(fontface = genes.gp[1], 
                                                                  fontsize = as.numeric(genes.gp[2]), col = gcol))
    }
    else {
      geneMark = NULL
    }
    right_annotation = ComplexHeatmap::rowAnnotation(gene = geneMark, 
                                                     cluster = anno.block)
    if (object$type == "monocle" | object$geneMode == "all" | 
        ncol(mat) > 20) {
      show_column_names = FALSE
    }
    else {
      show_column_names = TRUE
    }
    if (object$geneType == "non-branched") {
      rg <- range(as.numeric(as.character(sample.info)))
      col_fun2 = circlize::colorRamp2(c(rg[1], rg[2]), 
                                      pseudotime_col)
      lgd = ComplexHeatmap::Legend(col_fun = col_fun2, 
                                   title = "pseudotime")
      lgd_list = list(lgd)
    }
    else if (object$geneType == "branched") {
      if (length(levels(sample.info)) == 3) {
        lgd = ComplexHeatmap::Legend(labels = levels(sample.info), 
                                     legend_gp = grid::gpar(fill = pseudotime_col), 
                                     title = "branch")
      }
      else {
        lgd = ComplexHeatmap::Legend(labels = levels(sample.info), 
                                     legend_gp = grid::gpar(fill = pseudotime_col), 
                                     title = "branch")
      }
      lgd_list = list(lgd)
    }
    else {
      lgd_list <- NULL
    }
    if (plot.type == "heatmap") {
      htf <- ComplexHeatmap::pheatmap(as.matrix(mat), name = "Relative \n expression", 
                                     cluster_columns = cluster_columns, show_row_names = show_row_names, 
                                     border = border, column_split = column_split, 
                                     row_split = subgroup, cluster_row_slices = cluster_row_slices, 
                                     column_names_side = "top", show_column_names = show_column_names, 
                                     top_annotation = topanno, right_annotation = right_annotation, cluster_rows = F,
                                     breaks = seq(-4,4,length.out = 256),
                                     col = col_fun, use_raster = use_raster, ...)
      ComplexHeatmap::draw(htf, merge_legend = TRUE, annotation_legend_list = lgd_list)
    }
    else {
      rg = range(mat)
      if (!is.null(gglist)) {
        anno_ggplot2 = ComplexHeatmap::anno_zoom(align_to = align_to, 
                                                 which = "row", panel_fun = function(index, 
                                                                                     nm) {
                                                   g <- gglist[[nm]]
                                                   g <- grid::grid.grabExpr(print(g))
                                                   grid::pushViewport(grid::viewport())
                                                   grid::grid.rect()
                                                   grid::grid.draw(g)
                                                   grid::popViewport()
                                                 }, size = grid::unit(as.numeric(ggplot.panel.arg[1]), 
                                                                      "cm"), gap = grid::unit(as.numeric(ggplot.panel.arg[2]), 
                                                                                              "cm"), width = grid::unit(as.numeric(ggplot.panel.arg[3]), 
                                                                                                                        "cm"), side = "right", link_gp = grid::gpar(fill = ggplot.panel.arg[4], 
                                                                                                                                                                    col = ggplot.panel.arg[5]))
      }
      else {
        anno_ggplot2 = NULL
      }
      panel_fun = function(index, nm) {
        if (add.box == TRUE & add.line != TRUE) {
          xscale = c(-0.1, 1.1)
        }
        else {
          xscale = c(0, 1)
        }
        grid::pushViewport(grid::viewport(xscale = xscale, 
                                          yscale = c(0, 1)))
        grid::grid.rect()
        if (object$geneMode == "all" & object$type == 
            "scRNAdata") {
          mulGroup <- cell.num.info
          grid::grid.lines(x = c(0, 1), y = rep(0.5, 
                                                2), gp = grid::gpar(col = "black", lty = "dashed"))
          cu <- cumsum(mulGroup)
          seqn <- data.frame(st = c(1, cu[1:(length(cu) - 
                                               1)] + 1), sp = c(cu[1], cu[2:length(cu)]))
        }
        else {
          if (is.null(mulGroup)) {
            mulGroup <- ncol(mat)
            seqn <- data.frame(st = 1, sp = ncol(mat))
          }
          else {
            mulGroup <- mulGroup
            grid::grid.lines(x = c(0, 1), y = rep(0.5, 
                                                  2), gp = grid::gpar(col = "black", lty = "dashed"))
            cu <- cumsum(mulGroup)
            seqn <- data.frame(st = c(1, cu[1:(length(cu) - 
                                                 1)] + 1), sp = c(cu[1], cu[2:length(cu)]))
          }
        }
        if (object$geneMode == "all" & object$type == 
            "scRNAdata") {
          cell.ave <- purrr::map_dfr(1:nrow(seqn), function(x) {
            tmp <- seqn[x, ]
            tmpmat <- mat[index, c(tmp$st:tmp$sp)]
            rg <- range(mat[index, ])
            if (set.md == "mean") {
              mdia <- mean(rowMeans(tmpmat))
            }
            else if (set.md == "median") {
              mdia <- stats::median(apply(tmpmat, 1, 
                                          stats::median))
            }
            else {
              message("supply mean/median !")
            }
            res <- data.frame(x = x, val = mdia)
            return(res)
          })
          if (add.line == TRUE) {
            grid::grid.lines(x = scales::rescale(cell.ave$x, 
                                                 to = c(0, 1)), y = scales::rescale(cell.ave$val, 
                                                                                    to = c(0.1, 0.9)), gp = grid::gpar(lwd = 3, 
                                                                                                                       col = mline.col))
          }
        }
        else {
          lapply(1:nrow(seqn), function(x) {
            tmp <- seqn[x, ]
            tmpmat <- mat[index, c(tmp$st:tmp$sp)]
            if (set.md == "mean") {
              mdia <- colMeans(tmpmat)
            }
            else if (set.md == "median") {
              mdia <- apply(tmpmat, 2, stats::median)
            }
            else {
              message("supply mean/median !")
            }
            pos = scales::rescale(1:ncol(tmpmat), to = c(0, 
                                                         1))
            if (is.null(boxcol)) {
              boxcol <- rep("grey90", ncol(tmpmat))
            }
            else {
              boxcol <- boxcol
            }
            if (add.box == TRUE) {
              lapply(1:ncol(tmpmat), function(x) {
                ComplexHeatmap::grid.boxplot(scales::rescale(tmpmat[, 
                                                                    x], to = c(0, 1), from = c(rg[1] - 
                                                                                                 0.5, rg[2] + 0.5)), pos = pos[x], 
                                             direction = "vertical", box_width = as.numeric(box.arg[1]), 
                                             outline = FALSE, gp = grid::gpar(col = box.arg[2], 
                                                                              fill = boxcol[x]))
              })
            }
            if (add.point == TRUE) {
              grid::grid.points(x = scales::rescale(1:ncol(tmpmat), 
                                                    to = c(0, 1)), y = scales::rescale(mdia, 
                                                                                       to = c(0, 1), from = c(rg[1] - 0.5, 
                                                                                                              rg[2] + 0.5)), pch = as.numeric(point.arg[1]), 
                                gp = grid::gpar(fill = point.arg[2], 
                                                col = point.arg[3]), size = grid::unit(as.numeric(point.arg[4]), 
                                                                                       "char"))
            }
            if (add.line == TRUE) {
              grid::grid.lines(x = scales::rescale(1:ncol(tmpmat), 
                                                   to = c(0, 1)), y = scales::rescale(mdia, 
                                                                                      to = c(0, 1), from = c(rg[1] - 0.5, 
                                                                                                             rg[2] + 0.5)), gp = grid::gpar(lwd = 3, 
                                                                                                                                            col = mline.col[x]))
            }
          })
        }
        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                                                "ComplexHeatmap")
        text <- paste("Gene Size:", nrow(mat[index, 
        ]), sep = " ")
        grid.textbox(text, x = textbox.pos[1], y = textbox.pos[2], 
                     gp = grid::gpar(fontsize = textbox.size, fontface = "italic", 
                                     ...))
        grid::popViewport()
      }
      if (!is.null(subgroup.anno)) {
        align_to = split(1:nrow(mat), subgroup)
        align_to = align_to[subgroup.anno]
      }
      else {
        align_to = subgroup
      }
      anno = ComplexHeatmap::anno_link(align_to = align_to, 
                                       which = "row", panel_fun = panel_fun, size = grid::unit(as.numeric(panel.arg[1]), 
                                                                                               "cm"), gap = grid::unit(as.numeric(panel.arg[2]), 
                                                                                                                       "cm"), width = grid::unit(as.numeric(panel.arg[3]), 
                                                                                                                                                 "cm"), side = line.side, link_gp = grid::gpar(fill = panel.arg[4], 
                                                                                                                                                                                               col = panel.arg[5]))
      if (!is.null(annoTerm.data)) {
        termanno <- annoTerm.data
        if (ncol(termanno) == 2) {
          colnames(termanno) <- c("id", "term")
        }
        else if (ncol(termanno) == 3) {
          colnames(termanno) <- c("id", "term", "pval")
        }
        else if (ncol(termanno) == 4) {
          colnames(termanno) <- c("id", "term", "pval", 
                                  "ratio")
        }
        else {
          message("No more than 4 columns!")
        }
        if (is.null(go.col)) {
          gocol <- circlize::rand_color(n = nrow(termanno))
        }
        else {
          gocol <- go.col
        }
        if (is.null(go.size)) {
          gosize <- rep(12, nrow(termanno))
        }
        else {
          if (go.size == "pval") {
            termanno.tmp <- purrr::map_df(unique(termanno$id), 
                                          function(x) {
                                            tmp <- termanno %>% dplyr::filter(id == 
                                                                                x) %>% dplyr::mutate(size = scales::rescale(-log10(pval), 
                                                                                                                            to = term.text.limit))
                                          })
            gosize <- termanno.tmp$size
          }
          else {
            gosize <- go.size
          }
        }
        termanno <- termanno %>% dplyr::ungroup() %>% 
          dplyr::mutate(col = gocol, fontsize = gosize)
        term.list <- lapply(1:length(unique(termanno$id)), 
                            function(x) {
                              tmp = termanno[which(termanno$id == unique(termanno$id)[x]), 
                              ]
                              df <- data.frame(text = tmp$term, col = tmp$col, 
                                               fontsize = tmp$fontsize)
                              return(df)
                            })
        names(term.list) <- unique(termanno$id)
        if (!is.null(subgroup.anno)) {
          align_to2 = split(seq_along(subgroup), subgroup)
          align_to2 = align_to2[subgroup.anno]
          term.list = term.list[subgroup.anno]
        }
        else {
          align_to2 = subgroup
          term.list = term.list
        }
        textbox = ComplexHeatmap::anno_textbox(align_to2, 
                                               term.list, word_wrap = word_wrap, add_new_line = add_new_line, 
                                               side = annoTerm.mside, background_gp = grid::gpar(fill = termAnno.arg[1], 
                                                                                                 col = termAnno.arg[2]), by = by.go)
        if (ncol(termanno) - 2 > 2) {
          anno_gobar <- function(data = NULL, bar.width = 0.1, 
                                 align_to = NULL, panel.arg = panel.arg, 
                                 ...) {
            if (ncol(data) - 2 == 3) {
              data <- data %>% dplyr::mutate(bary = -log10(pval))
            }
            else {
              data <- data %>% dplyr::mutate(bary = ratio)
            }
            ComplexHeatmap::anno_zoom(align_to = align_to, 
                                      which = "row", panel_fun = function(index, 
                                                                          nm) {
                                        grid::pushViewport(grid::viewport(xscale = c(0, 
                                                                                     1), yscale = c(0, 1)))
                                        grid::grid.rect()
                                        tmp <- data %>% dplyr::filter(id == 
                                                                        nm)
                                        grid::grid.segments(x0 = rep(0, nrow(tmp)), 
                                                            x1 = scales::rescale(rev(tmp$bary), 
                                                                                 to = c(0.1, 0.9)), y0 = scales::rescale(1:nrow(tmp), 
                                                                                                                         to = c(0.1, 0.9)), y1 = scales::rescale(1:nrow(tmp), 
                                                                                                                                                                 to = c(0.1, 0.9)), gp = grid::gpar(lwd = bar.width, 
                                                                                                                                                                                                    col = rev(tmp$col), lineend = "butt"))
                                        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                                                                                "ComplexHeatmap")
                                        text <- nm
                                        grid.textbox(text, x = textbar.pos[1], 
                                                     y = textbar.pos[2], gp = grid::gpar(fontsize = textbox.size, 
                                                                                         fontface = "italic", col = unique(tmp$col), 
                                                                                         ...))
                                        grid::popViewport()
                                      }, size = grid::unit(as.numeric(panel.arg[1]), 
                                                           "cm"), gap = grid::unit(as.numeric(panel.arg[2]), 
                                                                                   "cm"), width = grid::unit(as.numeric(panel.arg[3]), 
                                                                                                             "cm"), side = "right", link_gp = grid::gpar(fill = termAnno.arg[1], 
                                                                                                                                                         col = termAnno.arg[2]), ...)
          }
          baranno = anno_gobar(data = termanno, align_to = align_to2, 
                               panel.arg = panel.arg, bar.width = bar.width)
        }
        if (add.bar == TRUE) {
          baranno
        }
        else {
          baranno = NULL
        }
      }
      else {
        textbox = NULL
        baranno = NULL
      }
      if (!is.null(annoKegg.data)) {
        termanno <- annoKegg.data
        if (ncol(termanno) == 2) {
          colnames(termanno) <- c("id", "term")
        }
        else if (ncol(termanno) == 3) {
          colnames(termanno) <- c("id", "term", "pval")
        }
        else if (ncol(termanno) == 4) {
          colnames(termanno) <- c("id", "term", "pval", 
                                  "ratio")
        }
        else {
          message("No more than 4 columns!")
        }
        if (is.null(kegg.col)) {
          gocol <- circlize::rand_color(n = nrow(termanno))
        }
        else {
          gocol <- kegg.col
        }
        if (is.null(kegg.size)) {
          gosize <- rep(12, nrow(termanno))
        }
        else {
          if (kegg.size == "pval") {
            termanno.tmp <- purrr::map_df(unique(termanno$id), 
                                          function(x) {
                                            tmp <- termanno %>% dplyr::filter(id == 
                                                                                x) %>% dplyr::mutate(size = scales::rescale(-log10(pval), 
                                                                                                                            to = term.text.limit))
                                          })
            gosize <- termanno.tmp$size
          }
          else {
            gosize <- kegg.size
          }
        }
        termanno <- termanno %>% dplyr::ungroup() %>% 
          dplyr::mutate(col = gocol, fontsize = gosize)
        term.list <- lapply(1:length(unique(termanno$id)), 
                            function(x) {
                              tmp = termanno[which(termanno$id == unique(termanno$id)[x]), 
                              ]
                              df <- data.frame(text = tmp$term, col = tmp$col, 
                                               fontsize = tmp$fontsize)
                              return(df)
                            })
        names(term.list) <- unique(termanno$id)
        if (!is.null(subgroup.anno)) {
          align_to2 = split(seq_along(subgroup), subgroup)
          align_to2 = align_to2[subgroup.anno]
          term.list = term.list[subgroup.anno]
        }
        else {
          align_to2 = subgroup
          term.list = term.list
        }
        textbox.kegg = ComplexHeatmap::anno_textbox(align_to2, 
                                                    term.list, word_wrap = word_wrap, add_new_line = add_new_line, 
                                                    side = annoKegg.mside, background_gp = grid::gpar(fill = keggAnno.arg[1], 
                                                                                                      col = keggAnno.arg[2]), by = by.kegg)
        if (ncol(termanno) - 2 > 2) {
          anno_keggbar <- function(data = NULL, bar.width = 0.1, 
                                   align_to = NULL, panel.arg = panel.arg, 
                                   ...) {
            if (ncol(data) - 2 == 3) {
              data <- data %>% dplyr::mutate(bary = -log10(pval))
            }
            else {
              data <- data %>% dplyr::mutate(bary = ratio)
            }
            ComplexHeatmap::anno_zoom(align_to = align_to, 
                                      which = "row", panel_fun = function(index, 
                                                                          nm) {
                                        grid::pushViewport(grid::viewport(xscale = c(0, 
                                                                                     1), yscale = c(0, 1)))
                                        grid::grid.rect()
                                        tmp <- data %>% dplyr::filter(id == 
                                                                        nm)
                                        grid::grid.segments(x0 = rep(0, nrow(tmp)), 
                                                            x1 = scales::rescale(rev(tmp$bary), 
                                                                                 to = c(0.1, 0.9)), y0 = scales::rescale(1:nrow(tmp), 
                                                                                                                         to = c(0.1, 0.9)), y1 = scales::rescale(1:nrow(tmp), 
                                                                                                                                                                 to = c(0.1, 0.9)), gp = grid::gpar(lwd = bar.width, 
                                                                                                                                                                                                    col = rev(tmp$col), lineend = "butt"))
                                        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                                                                                "ComplexHeatmap")
                                        text <- nm
                                        grid.textbox(text, x = textbar.pos[1], 
                                                     y = textbar.pos[2], gp = grid::gpar(fontsize = textbox.size, 
                                                                                         fontface = "italic", col = unique(tmp$col), 
                                                                                         ...))
                                        grid::popViewport()
                                      }, size = grid::unit(as.numeric(panel.arg[1]), 
                                                           "cm"), gap = grid::unit(as.numeric(panel.arg[2]), 
                                                                                   "cm"), width = grid::unit(as.numeric(panel.arg[3]), 
                                                                                                             "cm"), side = "right", link_gp = grid::gpar(fill = keggAnno.arg[1], 
                                                                                                                                                         col = keggAnno.arg[2]), ...)
          }
          baranno.kegg = anno_keggbar(data = termanno, 
                                      align_to = align_to2, panel.arg = panel.arg, 
                                      bar.width = bar.width)
        }
        if (add.kegg.bar == TRUE) {
          baranno.kegg
        }
        else {
          baranno.kegg = NULL
        }
      }
      else {
        textbox.kegg = NULL
        baranno.kegg = NULL
      }
      if (line.side == "right") {
        if (markGenes.side == "right") {
          right_annotation2 = ComplexHeatmap::rowAnnotation(gene = geneMark, 
                                                            cluster = anno.block, line = anno, anno_ggplot2 = anno_ggplot2, 
                                                            textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                                                            baranno.kegg = baranno.kegg)
          left_annotation = NULL
        }
        else {
          right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block, 
                                                            line = anno, anno_ggplot2 = anno_ggplot2, 
                                                            textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                                                            baranno.kegg = baranno.kegg)
          left_annotation = ComplexHeatmap::rowAnnotation(gene = geneMark)
        }
      }
      else {
        if (markGenes.side == "right") {
          right_annotation2 = ComplexHeatmap::rowAnnotation(gene = geneMark, 
                                                            cluster = anno.block, anno_ggplot2 = anno_ggplot2, 
                                                            textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                                                            baranno.kegg = baranno.kegg)
          left_annotation = ComplexHeatmap::rowAnnotation(line = anno)
        }
        else {
          right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block, 
                                                            anno_ggplot2 = anno_ggplot2, textbox = textbox, 
                                                            bar = baranno, textbox.kegg = textbox.kegg, 
                                                            baranno.kegg = baranno.kegg)
          left_annotation = ComplexHeatmap::rowAnnotation(line = anno, 
                                                          gene = geneMark)
        }
      }
      if (object$type == "monocle" | object$geneMode == 
          "all" | ncol(mat) > 20) {
        show_column_names = FALSE
      }
      else {
        show_column_names = TRUE
      }
      htf <- ComplexHeatmap::Heatmap(as.matrix(mat), name = "Z-score", 
                                     cluster_columns = cluster_columns, show_row_names = show_row_names, 
                                     border = border, column_split = column_split, 
                                     top_annotation = topanno, right_annotation = right_annotation2, 
                                     left_annotation = left_annotation, column_names_side = "top", 
                                     show_column_names = show_column_names, row_split = subgroup, 
                                     cluster_row_slices = cluster_row_slices, col = col_fun, 
                                     use_raster = use_raster, ...)
      if (is.null(mulGroup)) {
        ComplexHeatmap::draw(htf, merge_legend = TRUE, 
                             annotation_legend_list = lgd_list)
      }
      else {
        if (is.null(lgd.label)) {
          lgd.label <- paste("group", 1:length(mulGroup), 
                             sep = "")
        }
        else {
          lgd.label <- lgd.label
        }
        lgd_list2 = ComplexHeatmap::Legend(labels = lgd.label, 
                                           type = "lines", legend_gp = grid::gpar(col = mline.col, 
                                                                                  lty = 1))
        if (!is.null(lgd_list)) {
          lgd_list_com <- ComplexHeatmap::packLegend(lgd_list, 
                                                     lgd_list2)
        }
        else {
          lgd_list_com <- lgd_list2
        }
        ComplexHeatmap::draw(htf, annotation_legend_list = lgd_list_com, 
                             merge_legend = TRUE)
      }
    }
  }
}

enrich.50go <- enrichCluster(object = sc.data,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        organism = "mmu",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 1234)

pdf('sc.top50.50go.pdf',height = 10,width = 14,onefile = F)

tmp.vis(object = sc.data,
        ht.col.list = list(col_range = c(-4, 0, 4)),
        plot.type = "both",
        column_names_rot = 45,
        show_row_dend = F,
        markGenes = markGenes,
        markGenes.side = "left",
        annoTerm.data = enrich.50go,
        line.side = "left",
        cluster.order = c(1:11),
        go.col = rep(cr.col[1:11],each = 5),
        sample.col = cr.col.use,
        ctAnno.col = cr.col.use,
        go.size = 8,
        add.bar = T)
dev.off()

# fig S1B -------
library(sfsmisc)
library(MASS)
celltypes <- unique(combined.sct@meta.data$large_ct)
set.seed(1234)

getEuclideanDistance <- function(celltype, lowcv = T){
  print(paste("Working on", celltype))
  library(hopach)
  tmp <- subset(combined.sct, cells = WhichCells(combined.sct, idents = celltype))
  
  counts <- GetAssayData(object = tmp, slot = "counts",assay = 'RNA')
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) > 0
  expr <- counts[keep_genes, ]
  
  ifelse(min(table(tmp$condition)) > 300, 
         expr <- expr[,c(rownames(tmp@meta.data[tmp$condition == 'MC',])[sample(1:nrow(tmp@meta.data[tmp$condition == 'MC',]),300)],
                         rownames(tmp@meta.data[tmp$condition == 'MO',])[sample(1:nrow(tmp@meta.data[tmp$condition == 'MO',]),300)])
                       ],
         expr <- expr)
  tmp <- subset(tmp,cells = colnames(expr))
  
  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  
  nsample <- min(table(tmp@meta.data$condition)[c("MC", "MO")])
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
  print(nsample)
  MO_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$condition == "MO")], nsample)
  MC_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$condition == "MC")], nsample)
  ds_expr_r <- ds_expr[, c(MC_r, MO_r)]
  
  if(lowcv){
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
  }
  else{
    genes <- rownames(ds_expr_r)
  }
  
  calcEuclDist <- function(matr, MC, MO){
    tmp <- data.matrix(sqrt(matr[genes, MC]))
    mean <- rowMeans(sqrt(matr[genes, MC]))
    d_MC <- distancevector(t(tmp), mean , d="euclid")
    names(d_MC) <- MC
    
    tmp <- data.matrix(sqrt(matr[genes, MO]))
    mean <- rowMeans(sqrt(matr[genes, MO]))
    d_MO <- distancevector(t(tmp), mean , d="euclid")
    names(d_MO) <- MO
    
    list(MC = d_MC, MO = d_MO)
  }
  ds <- calcEuclDist(matr = ds_expr_r, MO = MO_r, MC = MC_r)
  ds
}
res <- lapply(celltypes, function(x) getEuclideanDistance(x, lowcv = F))
names(res) <- celltypes
res_original <- res
options(scipen = 3)
diffs <- unlist(lapply(res_original, function(x) log2(mean(x[[2]]) / mean(x[[1]]))))
pvals <- unlist(lapply(res_original, function(x) wilcox.test(x[[1]], x[[2]])$p.value))
adj_pvals <- p.adjust(pvals, method = "BH")

ord <- c(1,7,3,4,10,5,6,11,8,2,9)
par(mar = c(15,5,5,5))
boxplot(do.call(c, res[ord]),
            gap = c(0.155,0.17),
        las = 2, outline = F,
        col = c("#0B996F", "#D6570D"),
        ylab = "Transcriptional noise", 
        xaxt = 'n', ylim = c(0.11,0.19),
        axes = F,
        boxplot.width = 0.5)

axis(1, at = seq(from = 1.5, to = 22, by = 2), names(res)[ord], las = 2)
axis(side = 2, at = c(0.11, 0.14, 0.16), labels = FALSE)

# fig S1C -----
adipo.sub <- subset(combined.sct, cells = WhichCells(combined.sct, idents = c('Brown adipocytes','White adipocytes','Prolif.adipocytes','FAPs')))
DimPlot(adipo.sub,raster = T)

adipo.sub <- RunUMAP(adipo.sub,dims = 1:30)
DimPlot(adipo.sub,raster = T,label = F, cols = c(npg.col), group.by = 'large_ct',raster.dpi = c(512,512),pt.size = 1) + ggtitle('Cell types') +
  NoAxes() + 
  theme(panel.background=element_rect(fill='transparent', color='black'), 
        title = element_text(size = 10),
        legend.text = element_text(size = 10), legend.key.height=unit(.8,"line"))


adipo.sub <- CreateSeuratObject(counts = GetAssayData(object = adipo.sub, slot = "counts"),
                                meta.data = adipo.sub@meta.data)
# re-integrated adipocytes -------
adipo.split <- SplitObject(adipo.sub,split.by = 'samples')

adipo.split <- lapply(X = adipo.split, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst',nfeatures = 1500)
})
features <- SelectIntegrationFeatures(adipo.split,nfeatures = 1500)

plan(multisession, workers=40)
adipo.split <- lapply(X = adipo.split, FUN = function(x){
  x <- ScaleData(x, vars.to.regress = c("mitoRatio",'G2M.Score','S.Score'),features = features)
  x <- RunPCA(x,features = features)
})
save.image()
anchors <- FindIntegrationAnchors(object.list = adipo.split,
                                  # normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction = 'rpca',k.anchor = 20)
adipo.inte <- IntegrateData(anchorset = anchors)
DefaultAssay(adipo.inte)

adipo.inte$condition <- NA
adipo.inte$condition[which(str_detect(adipo.inte@meta.data$cells, "^HDF_"))] <- "HFD"
adipo.inte$condition[which(str_detect(adipo.inte@meta.data$cells, "^ND_"))] <- "ND"
adipo.inte$gander <- factor(c(rep('female',21684), rep('male',24825), rep('female',28220),rep('male',31450)))

adipo.inte <- adipo.inte %>% 
  ScaleData() %>%
  RunPCA( verbose = FALSE) %>% 
  RunUMAP( reduction = "pca", dims = 1:30, verbose  = FALSE) %>% 
  FindNeighbors( reduction = "pca", dims = 1:30) %>% 
  FindClusters( resolution = c(0.4, 0.6, 0.8, 1.0,1.2))                    

adipo.inte$condition <- NA
adipo.inte$condition[which(str_detect(adipo.inte@meta.data$cells, "^HDF_"))] <- "MO"
adipo.inte$condition[which(str_detect(adipo.inte@meta.data$cells, "^ND_"))] <- "MC"

adipo.inte$gander <- NA
adipo.inte$gander[which(str_detect(adipo.inte@meta.data$samples, "HFD_F1"))] <- "female"
adipo.inte$gander[which(str_detect(adipo.inte@meta.data$samples, "HFD_F2"))] <- "female"
adipo.inte$gander[which(str_detect(adipo.inte@meta.data$samples, "ND_F1"))] <- "female"
adipo.inte$gander[which(str_detect(adipo.inte@meta.data$samples, "ND_F2"))] <- "female"
adipo.inte$gander[which(str_detect(adipo.inte@meta.data$samples, "HFD_M1"))] <- "male"
adipo.inte$gander[which(str_detect(adipo.inte@meta.data$samples, "HFD_M2"))] <- "male"
adipo.inte$gander[which(str_detect(adipo.inte@meta.data$samples, "ND_M1"))] <- "male"
adipo.inte$gander[which(str_detect(adipo.inte@meta.data$samples, "ND_M2"))] <- "male"
adipo.inte$con_gan <- paste0(adipo.inte$condition,'_',adipo.inte$gander)
adipo.inte$con_gan <- factor(adipo.inte$con_gan, levels = c('MC_male','MO_male',
                                                          'MC_female','MO_female'))

adipo.inte$large_ct[which(str_detect(adipo.inte@meta.data$large_ct,"Prolif.adipocytes"))] <- "Brown adipocytes"
adipo.inte$large_ct <- factor(adipo.inte$large_ct, levels = c('FAPs','Brown adipocytes','White adipocytes'))
adipo.inte$large_ct <- adipo.inte$large_ct 

dittoDimPlot(adipo.inte, "large_ct", size = 1)
DimPlot(adipo.inte,raster = T,label = T, cols = npg.col, group.by = 'large_ct',raster.dpi = c(4096,4096),pt.size = 6) + ggtitle('Cell types') +
  NoAxes() + 
  # scale_color_igv() + 
  theme(panel.background=element_rect(fill='transparent', color='black'), 
        title = element_text(size = 10),
        legend.text = element_text(size = 10), legend.key.height=unit(.8,"line"))

# fig S1D -----
adipo.inte.5000 <- subset(adipo.inte,cells = c(rownames(adipo.inte@meta.data[sample(c(1:106179),5000),])))
library(monocle3)
library(Seurat)
library(viridis)
library(ggsci)

mono <- GetAssayData(adipo.inte.5000,assay = 'RNA', layer = 'counts')
cell_meta <- adipo.inte.5000@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(mono))
rownames(gene_annotation) <- rownames(mono)

cds <- new_cell_data_set(mono,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 10)
plot_pc_variance_explained(cds)

# mono umap
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds, reduction_method="UMAP", color_cells_by="large_ct", group_label_size = 5,cell_size = 1) + ggtitle('cds.umap')

# seurat umap
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(adipo.inte.5000, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="large_ct") + ggtitle('bov.umap')

# cds <- cluster_cells(cds,reduction_method = 'UMAP',resolution = .001)
cds <- cluster_cells(cds,reduction_method = 'UMAP',resolution = 1)
# get trajectory
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = F,
           label_branch_points = FALSE, color_cells_by="large_ct", group_label_size = 0,cell_size = 1) + 
  scale_color_npg() +
  ggtitle("Pseeudotime") +
  NoAxes() + 
  theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 10),
        legend.text = element_text(size = 8), legend.key.height=unit(.8,"line"))


test <- order_cells(cds)
plot_cells(test, color_cells_by = "pseudotime", 
           label_cell_groups = F, 
           label_leaves = FALSE,  
           label_branch_points = F,
           group_label_size = 1,
           cell_size = .1,
           label_roots = F,
           cell_stroke = 1,
           trajectory_graph_segment_size = 1
) + NoAxes() + 
  theme(panel.background=element_rect(fill='transparent', color='black'), 
        title = element_text(size = 10),
        legend.text = element_text(size = 10), legend.key.height=unit(.8,"line"))


# fig S1E -----
library(CytoTRACE2)
library(patchwork)
cytotrace2_result <- cytotrace2(adipo.inte.5000,is_seurat = T)
annotation <- data.frame(phenotype = adipo.inte.5000@meta.data$large_ct) %>% 
  set_rownames(., colnames(adipo.inte.5000))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = T)                      
library(ggpubr)
cytotrace2_result$large_ct <- factor(cytotrace2_result$large_ct,levels = c('FAPs','Brown Adipocytes','White Adipocytes'))

p1 <- ggboxplot(cytotrace2_result@meta.data, x="large_ct", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",
                fill="large_ct",
                palette = "npg",
                xlab = F, 
                bxp.errorbar=T,
                bxp.errorbar.width=0.5, 
                size=1, 
                outlier.shape=NA,
                legend = "right",
                alpha = 0.8) + 
  ylab('Potency score')  + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),legend.position = 'none')

my_comparisons <- list(c("FAPs", "Brown Adipocytes"), c("Brown Adipocytes", "White Adipocytes"))
p1+stat_compare_means(comparisons = my_comparisons,
                      method = "wilcox.test")

# fig S1F & G & H------
w.adipo <- subset(combined.sct, cells = WhichCells(combined.sct, idents = c('White adipocytes')))
DimPlot(w.adipo,raster = F)

w.adipo <- RunUMAP(w.adipo,dims = 1:30)
w.adipo <- CreateSeuratObject(counts = GetAssayData(object = w.adipo, slot = "counts"),
                                 meta.data = w.adipo@meta.data)
# re-integrated adipocytes -------
adipo.split <- SplitObject(w.adipo,split.by = 'samples')
plan(multisession, workers=1)
adipo.split <- lapply(X = adipo.split, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst',nfeatures = 1500)
})
features <- SelectIntegrationFeatures(adipo.split,nfeatures = 1500)

plan(multisession, workers=40)
adipo.split <- lapply(X = adipo.split, FUN = function(x){
  x <- ScaleData(x, vars.to.regress = c("mitoRatio",'G2M.Score','S.Score'),features = features)
  x <- RunPCA(x,features = features)
})
save.image()
anchors <- FindIntegrationAnchors(object.list = adipo.split,
                                  # normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction = 'rpca',k.anchor = 20)
w.adipo.inte <- IntegrateData(anchorset = anchors)
rm(w.adipo,anchors,adipo.split)
DefaultAssay(w.adipo.inte) <- 'integrated'
w.adipo.inte <- w.adipo.inte %>% 
  ScaleData() %>%
  RunPCA( verbose = FALSE) %>% 
  RunUMAP( reduction = "pca", dims = 1:30, verbose = FALSE) 

w.adipo.inte$con_gan <- factor(w.adipo.inte$con_gan, levels = c('ND_male','HFD_male','ND_female','HFD_female'))
DimPlot(w.adipo.inte,raster = F,label = F, cols = simpson.col, group.by = 'con_gan',split.by = 'con_gan',raster.dpi = c(4096,4096),pt.size = .5) + ggtitle('') +
  NoAxes() + 
  # scale_color_igv() + 
  theme(panel.background=element_rect(fill='transparent', color='black'), 
        title = element_text(size = 10),
        legend.text = element_text(size = 10), legend.key.height=unit(.8,"line"))

comlist <- t(combn(unique(w.adipo.inte$condition), 2))

for (ct in unique(w.adipo.inte$gander)) {
  print(ct)
  Idents(w.adipo.inte) <- 'gander'
  sub.seu.1 <- subset(w.adipo.inte, idents = ct)
  for (gp in nrow(comlist)) {
    age1 <- comlist[gp,1]
    age2 <- comlist[gp,2]
    DEGs <- FindMarkers(sub.seu.1,
                        ident.1 = age2,
                        ident.2 = age1, 
                        logfc.threshold = 0.01,
                        group.by = 'condition')
    DEGs.name <- paste('w.adipo.all',ct,age2,age1,sep = '_')
    assign(DEGs.name,DEGs)
  }
}
w.adipo.all.deg <- list(w.adipo.all_male_MC_MO = w.adipo.all_male_ND_HFD,
                        w.adipo.all_female_MC_MO = w.adipo.all_female_ND_HFD)
rio::export(w.adipo.all.deg, file = '../01.analysis/063fig.final//Table6.white adipocytes wilcox.test.lfc0.01.xlsx',row.names = T)

keyvals.colour <- ifelse(
  w.adipo.all_male_ND_HFD$avg_log2FC < -0.25, 'royalblue',
  ifelse(w.adipo.all_male_ND_HFD$avg_log2FC > 0.25, 'red3',
         'gray80'))
keyvals.colour[is.na(keyvals.colour)] <- 'gray80'
names(keyvals.colour)[keyvals.colour == 'red3'] <- 'MC high'
names(keyvals.colour)[keyvals.colour == 'gray80'] <- 'No sig'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'MO high'

EnhancedVolcano(w.adipo.all_male_ND_HFD,
                lab = rownames(w.adipo.all_male_ND_HFD),
                selectLab = 'Ucp1',
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 0.25,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 4,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                # col = c('gray80', 'gray80', 'gray80', 'red3'),
                colCustom = keyvals.colour,
                colAlpha = 2/5,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Male',
                subtitle = '',
                xlim = c(-1.5,1.5)
) 


keyvals.colour <- ifelse(
  w.adipo.all_female_ND_HFD$avg_log2FC < -0.25, 'royalblue',
  ifelse(w.adipo.all_female_ND_HFD$avg_log2FC > 0.25, 'red3',
         'gray80'))
keyvals.colour[is.na(keyvals.colour)] <- 'gray80'
names(keyvals.colour)[keyvals.colour == 'red3'] <- 'mND high'
names(keyvals.colour)[keyvals.colour == 'gray80'] <- 'No sig'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'mHFD high'

EnhancedVolcano(w.adipo.all_female_ND_HFD,
                lab = rownames(w.adipo.all_female_ND_HFD),
                selectLab = 'Ucp1',
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 0.25,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 4,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                # col = c('gray80', 'gray80', 'gray80', 'red3'),
                colCustom = keyvals.colour,
                colAlpha = 2/5,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Female',
                subtitle = '',
                xlim = c(-1.5,1.5)
) 

# fig S1I & J& K ----
plan(multisession, workers=1)
nbrOfWorkers()
# subtype from all cells -------
fap.only <- subset(combined.sct, cells = WhichCells(combined.sct, idents = c('FAPs')))
DimPlot(fap.only,raster = F)

fap.only <- RunUMAP(fap.only,dims = 1:30)
fap.only <- CreateSeuratObject(counts = GetAssayData(object = fap.only, slot = "counts"),
                                 meta.data = fap.only@meta.data)
# re-integrated fapcytes -------
fap.split <- SplitObject(fap.only,split.by = 'samples')

fap.split <- lapply(X = fap.split, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst',nfeatures = 1500)
})
features <- SelectIntegrationFeatures(fap.split,nfeatures = 1500)

plan(multisession, workers=40)
fap.split <- lapply(X = fap.split, FUN = function(x){
  x <- ScaleData(x, vars.to.regress = c("mitoRatio",'G2M.Score','S.Score'),features = features)
  x <- RunPCA(x,features = features)
})
save.image()
anchors <- FindIntegrationAnchors(object.list = fap.split,
                                  # normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction = 'rpca',k.anchor = 20)
fap.only.inte <- IntegrateData(anchorset = anchors)
fap.only.inte <- fap.only.inte %>% 
  ScaleData() %>%
  RunPCA( verbose = FALSE) %>% 
  RunUMAP( reduction = "pca", dims = 1:30, verbose = FALSE) 

DimPlot(fap.only.inte,raster = F,label = F, cols = simpson.col, group.by = 'con_gan',split.by = 'con_gan',raster.dpi = c(4096,4096),pt.size = .5) + ggtitle('') +
  NoAxes() + 
  # scale_color_igv() + 
  theme(panel.background=element_rect(fill='transparent', color='black'), 
        title = element_text(size = 10),
        legend.text = element_text(size = 10), legend.key.height=unit(.8,"line"))

comlist <- t(combn(unique(fap.only.inte$condition), 2))

for (ct in unique(fap.only.inte$gander)) {
  print(ct)
  Idents(fap.only.inte) <- 'gander'
  sub.seu.1 <- subset(fap.only.inte, idents = ct)
  for (gp in nrow(comlist)) {
    age1 <- comlist[gp,1]
    age2 <- comlist[gp,2]
    DEGs <- FindMarkers(sub.seu.1,
                        ident.1 = age2,
                        ident.2 = age1, 
                        logfc.threshold = 0.01,
                        group.by = 'condition')
    DEGs.name <- paste('fap.all',ct,age2,age1,sep = '_')
    assign(DEGs.name,DEGs)
  }
}
fap.all.deg <- list(fap.all_male_MC_MO = fap.all_male_ND_HFD,
                      fap.all_female_MC_MO = fap.all_female_ND_HFD)
rio::export(fap.all.deg, file = '../01.analysis/063fig.final//Table5.FAPs wilcox.test.lfc0.01.xlsx',row.names = T)

keyvals.colour <- ifelse(
  fap.all_male_ND_HFD$avg_log2FC < -0.25, 'royalblue',
  ifelse(fap.all_male_ND_HFD$avg_log2FC > 0.25, 'red3',
         'gray80'))
keyvals.colour[is.na(keyvals.colour)] <- 'gray80'
names(keyvals.colour)[keyvals.colour == 'red3'] <- 'MC high'
names(keyvals.colour)[keyvals.colour == 'gray80'] <- 'No sig'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'MO high'

EnhancedVolcano(fap.all_male_ND_HFD,
                lab = rownames(fap.all_male_ND_HFD),
                selectLab = 'Ucp1',
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 0.25,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 4,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                # col = c('gray80', 'gray80', 'gray80', 'red3'),
                colCustom = keyvals.colour,
                colAlpha = 2/5,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Male',
                subtitle = '',
                xlim = c(-1.5,1.5)
) 


keyvals.colour <- ifelse(
  fap.all_female_ND_HFD$avg_log2FC < -0.25, 'royalblue',
  ifelse(fap.all_female_ND_HFD$avg_log2FC > 0.25, 'red3',
         'gray80'))
keyvals.colour[is.na(keyvals.colour)] <- 'gray80'
names(keyvals.colour)[keyvals.colour == 'red3'] <- 'mND high'
names(keyvals.colour)[keyvals.colour == 'gray80'] <- 'No sig'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'mHFD high'

EnhancedVolcano(fap.all_female_ND_HFD,
                lab = rownames(fap.all_female_ND_HFD),
                selectLab = 'Ucp1',
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 0.25,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 4,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                # col = c('gray80', 'gray80', 'gray80', 'red3'),
                colCustom = keyvals.colour,
                colAlpha = 2/5,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Female',
                subtitle = '',
                xlim = c(-1.5,1.5)
) 
                       
