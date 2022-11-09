library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(here)

here::i_am("PlotCellTypeDist.R")

### Load the cX raw_data
cX_data <- read_csv(here("data", paste0("Tonsil_Tissues.csv")))

cX_subpop <- cX_data %>%
    mutate(
        celltype = case_when(
            `subpop:Epithelial_cells` == 1 ~ "Epithelial",
            `subpop:Endothelial_cells` == 1 ~ "Endothelial",
            `subpop:Lymphatic_cells` == 1 ~ "Lymphatic",
            `subpop:Dendritic_cells` == 1 ~ "DC",
            `subpop:M1_Macrophage` == 1 ~ "M1",
            `subpop:M2_Macrophage` == 1 ~ "M2",
            `subpop:Neutrophils` == 1 ~ "Neutrophil",
            `subpop:CD8_T_cells` == 1 ~ "CD8",
            `subpop:CD4_T_cells` == 1 ~ "CD4",
            `subpop:Treg_cells` == 1 ~ "Treg",
            `subpop:B_cells` == 1 ~ "B",
            `subpop:Other_cells` == 1 ~ "Other",
            TRUE ~ "Other"
        )
    )

### Load Flowsom cX2 data
cX_flowsom <- read_csv(here("data", paste0("flowsom_cX2_celltypes.csv"))) %>%
    rename(
        `Region` = Well,
        `label:mask:cell_region` = cellID,
        `pos_x:mask:cell_region` = pos_x,
        `pos_y:mask:cell_region` = pos_y,
        `celltype` = CX2_full
    ) %>%
    select(
        Region,
        `label:mask:cell_region`,
        #`pos_x:mask:cell_region`,
        #`pos_y:mask:cell_region`,
        celltype
    )

### Join the data
cX_data <- full_join(
    cX_data, cX_flowsom,
    by = c("Region", "label:mask:cell_region")
)

cX_data <- cX_data %>% rename(
    Cytokeratin = `mean_intensity:CK:cell_region`,
    CD31  = `mean_intensity:CD31:cell_region`,
    Podoplanin = `mean_intensity:Podoplanin:cell_region`,
    CD11b = `mean_intensity:CD11b:cell_region`,
    CD11c = `mean_intensity:CD11c:cell_region`,
    CD15  = `mean_intensity:CD15:cell_region`,
    CD68  = `mean_intensity:CD68:cell_region`,
    CD163 = `mean_intensity:CD163:cell_region`,
    CD206 = `mean_intensity:CD206:cell_region`,
    CD20  = `mean_intensity:CD20:cell_region`,
    CD8   = `mean_intensity:CD8:cell_region`,
    CD4   = `mean_intensity:CD4:cell_region`,
    FoxP3 = `mean_intensity:FOXP3:cell_region`,
    a.SMA = `mean_intensity:aSMA:cell_region`,
    b.catenin = `mean_intensity:bcatenin:cell_region`,
    PD.1  = `mean_intensity:PD1:cell_region`,
    PD.L1 = `mean_intensity:PDL1:cell_region`,
    CD45  = `mean_intensity:CD45:cell_region`,
    EGFR = `mean_intensity:EGFR:cell_region`,
    Vimentin = `mean_intensity:Vimentin:cell_region`
) %>% mutate(
    software = "cellXpress 2"
)


### Read the mesmer data
mesmer_data <- read_csv(
    here("data", paste0("Tonsil_mesmer_features_raw.csv"))
) %>%
    rename(
        `Region` = PointNum,
        `label:mask:cell_region` = cellID,
        `pos_x:mask:cell_region` = pos_x,
        `pos_y:mask:cell_region` = pos_y,
        `area:mask:cell_region`  = cell_Area
    ) %>% mutate(
        Region = recode(
            Region,
            "A1" = "A01",
            "A2" = "A02",
            "A3" = "A03",
            "A4" = "A04",
            "A5" = "A05",
            "A6" = "A06"
        ),
        software = "DeepCell"
    )

### Add the ROI annotation
mesmer_is_ROI <- NULL

for (region_name in unique(cX_data$Region)) {

    message("Process ", region_name, "...")

    ### Generate a matrix for faster computations
    cX_data_pos <- cX_data %>%
        filter(Region == region_name) %>%
        select("pos_x:mask:cell_region", "pos_y:mask:cell_region") %>%
        as.matrix()

    cX_is_ROI_all <- cX_data %>%
        filter(Region == region_name) %>%
        pull(`is_member:ROI:All`)

    cX_celltypes_all <- cX_data %>%
        filter(Region == region_name) %>%
        pull(`celltype`)

    mesmer_data_pos <- mesmer_data %>%
        filter(Region == region_name) %>%
        select("pos_x:mask:cell_region", "pos_y:mask:cell_region") %>%
        as.matrix()

    ### Loop through all the mesmer cells
    mesmer_is_ROI_all <- rep(1, nrow(mesmer_data_pos))
    #cX_celltypes <- rep("Other", nrow(mesmer_data_pos))

    for (c_idx in seq_len(nrow(mesmer_data_pos))) {
        ### Find the distance of the T cells to all PDL1 hi cells
        cell_pos_cur <- mesmer_data_pos[c_idx, ]
        dist_all <- sweep(cX_data_pos, 2, cell_pos_cur)
        dist_all <- sqrt(dist_all[, 1]^2 + dist_all[, 2]^2)

        ### Make sure own will not be selected
        dist_all[dist_all == 0] <- Inf

        mesmer_is_ROI_all[c_idx] <- cX_is_ROI_all[which.min(dist_all)]
        #cX_celltypes[c_idx] <- cX_celltypes_all[which.min(dist_all)]
    }

    mesmer_is_ROI_tmp <- tibble(
        `Region` = region_name,
        `label:mask:cell_region` = mesmer_data %>%
            filter(Region == region_name) %>%
            pull(`label:mask:cell_region`),
        `is_member:ROI:All` = mesmer_is_ROI_all
        #`cX_celltype` = cX_celltypes
    )

    mesmer_is_ROI <- bind_rows(mesmer_is_ROI, mesmer_is_ROI_tmp)
}

### Join the data
mesmer_data <- full_join(
    mesmer_data, mesmer_is_ROI,
    by = c("Region", "label:mask:cell_region")
)


remove_outliers <- function(x, na.rm = TRUE, ...) {
    if (length(x) < 10) {
        rep(NA, length(x))

    } else {
        qnt <- quantile(x, probs = c(.01, .99), na.rm = na.rm, ...)
        #H <- 1.5 * IQR(x, na.rm = na.rm)
        y <- x
        y[x < qnt[1]] <- NA
        y[x > qnt[2]] <- NA
        y
    }
}

### Join the two datsets togeher
mkr_names <- c(
    "Cytokeratin" = "cytokeratin",
    "CD31" = "CD31",
    "Podoplanin" = "podoplanin",
    "CD11b" = "CD11b",
    "CD11c" = "CD11c",
    "CD15" = "CD15",
    "CD68" = "CD68",
    "CD163" = "CD163",
    "CD206" = "CD206",
    "CD20" = "CD20",
    "CD8" = "CD8",
    "CD4" = "CD4",
    "FoxP3" = "FOXP3",
    "PD.1" = "PD-1",
    "PD.L1" = "PD-L1",
    "CD45" = "CD45",
    "EGFR" = "EGFR",
    "a.SMA" = "αSMA",
    "b.catenin" = "β-catenin",
    "Vimentin" = "vimentin"
)

celltype_names <- c(
    ## Lymphoid cells
    "B"="B", "CD4"="Th", "CD8"="Tc", "Treg"="Treg",
    ## Myeloid cells
    "DC"="DC", "M1"="M1", "M2"="M2", "Neutrophil"="Neutro.",
    ## Others
    "Lymphatic"="Lymph.", "Endothelial"="Endo.", "Epithelial"="Epi."
)

raw_data <-
    bind_rows(
        cX_data %>% select(
            Region, celltype, software,
            `is_member:ROI:All`,
            all_of(names(mkr_names))
        ) %>% pivot_longer(
            cols = all_of(names(mkr_names))
        ),
        mesmer_data %>% select(
            Region, celltype, software,
            `is_member:ROI:All`,
            all_of(names(mkr_names))
        ) %>% pivot_longer(
            cols = all_of(names(mkr_names))
        )
    ) %>% mutate(
        value = log2(value + 1),
        celltype = recode(celltype, !!!celltype_names),
        celltype = factor(celltype, levels = celltype_names),
        name = recode(name, !!!mkr_names)
    ) %>% group_by(
        name
    ) %>%
    mutate(
        value = remove_outliers(value)
    ) %>% ungroup()

### Plot composition
celltype_cols <- c(
    "B" = "#96C3DC",
    "Th" = "#1B63A5",
    "Tc" = "#A4DB77",
    "Treg" = "#FFFF88",
    "DC" = "#2C9321",
    "M1" = "#FC6908",
    "M2" = "#BEA0CC",
    "Neutro." = "#562888",
    "Lymph." = "#FBB25c",
    "Endo." = "#F88587",
    "Epi." = "#D90017"
)

celldist_data <- raw_data %>%
    select(
        software, celltype
    ) %>%
    group_by(
        software, celltype
    ) %>%
    summarise(
        n = n()
    ) %>% 
    ungroup() %>%
    group_by(
        software
    ) %>%
    summarise(
        celltype = celltype,
        n = n,
        pct = n / sum(n) * 100
    ) %>%
    ungroup()

celltypes <- celldist_data %>%
    filter(software == "cellXpress 2") %>%
    arrange(-pct) %>%
    pull(celltype) %>%
    as.character()

celldist_data <- celldist_data %>%
    mutate(celltype = factor(celltype, levels = celltypes))

write_csv(celldist_data, file = here("results", "tonsil_celldist.csv"))

celldist_plot <- celldist_data %>%
    ggplot(
        aes(x = software, y = pct, fill = celltype)
    ) +
    geom_col() +
    scale_fill_manual(
        name = "Cell types",
        values = celltype_cols,
        breaks = names(celltype_cols)
    ) +
    scale_y_continuous(
        name = "% of all cells",
        limits = c(0, 101), expand = c(0, 0)
    ) +
    coord_flip() +
    theme_classic(base_size = 12) +
    theme(
        axis.text = element_text(colour = "black"),
        axis.title.y = element_blank()
    )

ggsave(
    here("figures", "CellTypeDist_all.pdf"),
    plot = celldist_plot,
    device = cairo_pdf,
    width = 4,
    height = 1,
    units = "in",
    dpi = 600
)


### Plot PD-1 and PD-L1 inside and outside GC
PD1_dist_plot <- raw_data %>%
    filter(
        name == "PD-1",
        celltype %in% c("Th", "Tc", "Treg"),
        Region %in% c("A03", "A04", "A05", "A06")
    ) %>%
    mutate(
        is_GC = factor(
            `is_member:ROI:All`,
            labels = c("Outside GC", "Inside GC")
        )
    ) %>%
    select(
        Region, software, is_GC, value
    ) %>%
    group_by(
        software, is_GC, Region
    ) %>%
    summarise(
        PD1_mean = mean(value, na.rm = TRUE)
    ) %>% ungroup(
    ) %>%
    ggplot(
        aes(x = Region, y = PD1_mean, fill = is_GC)
    ) +
    geom_bar(
        stat = "identity", width = .75, position = "dodge"
    ) +
    scale_y_continuous(
        name = "mean T-cell PD-1\nintensity (log2 f.u.)",
        limits = c(0, 15),
        breaks = c(0, 5, 10, 15),
        expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    facet_grid(cols = vars(software))  +
    theme(
        strip.background = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "darkgray")
    )

ggsave(
    here("figures", "pd1_dist_plot.pdf"),
    plot = PD1_dist_plot,
    device = cairo_pdf,
    width = 5,
    height = 1.7,
    units = "in",
    dpi = 600
)


## Plot M1 distribution
m1_dist_plot <- raw_data %>%
    filter(
        name == "CD68", ## Any marker will do
        Region %in% c("A03", "A04", "A05", "A06")
    ) %>%
    mutate(
        is_GC = factor(
            `is_member:ROI:All`,
            labels = c("Outside GC", "Inside GC")
        )
    ) %>%
    select(
        Region, software, is_GC, celltype
    ) %>%
    group_by(
        software, Region, is_GC, celltype
    ) %>%
    summarise(
        n = n()
    ) %>% ungroup() %>%
    complete(
        software, Region, is_GC, celltype,  fill = list(n = 0)
    ) %>%
    group_by(
        software, Region, is_GC
    ) %>%
    summarise(
        celltype = celltype,
        n = n,
        pct = n / sum(n) * 100
    ) %>%
    ungroup() %>%
    filter(
        celltype == "M1"
    ) %>%
    ggplot(
        aes(x = Region, y = pct, fill = is_GC)
    ) +
    geom_bar(
        stat = "identity", width = 0.75, position = "dodge") +
    scale_y_continuous(
        name = "M1 cell %",
        limits = c(0, 6),
        breaks = c(0, 2, 4, 6),
        expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    facet_grid(cols = vars(software)) +
    theme(
        strip.background = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "darkgray")
    )

ggsave(
    here("figures", "m1_dist_plot.pdf"),
    plot = m1_dist_plot,
    device = cairo_pdf,
    width = 5,
    height = 1.7,
    units = "in",
    dpi = 600
)

## Plot M2 distribution
m2_dist_plot <- raw_data %>%
    filter(
        name == "CD68", ## Any marker will do
        Region %in% c("A03", "A04", "A05", "A06")
    ) %>%
    mutate(
        is_GC = factor(
            `is_member:ROI:All`,
            labels = c("Outside GC", "Inside GC")
        )
    ) %>%
    select(
        Region, software, is_GC, celltype
    ) %>%
    group_by(
        software, Region, is_GC, celltype
    ) %>%
    summarise(
        n = n()
    ) %>% ungroup() %>%
    complete(
        software, Region, is_GC, celltype,  fill = list(n = 0)
    ) %>%
    group_by(
        software, Region, is_GC
    ) %>%
    summarise(
        celltype = celltype,
        n = n,
        pct = n / sum(n) * 100
    ) %>%
    ungroup() %>%
    filter(
        celltype == "M2"
    ) %>%
    ggplot(
        aes(x = Region, y = pct, fill = is_GC)
    ) +
    geom_bar(
        stat = "identity", width = 0.75, position = "dodge") +
    scale_y_continuous(
        name = "M2 cell %",
        limits = c(0, 8),
        breaks = c(0, 2, 4, 6, 8),
        expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    facet_grid(cols = vars(software)) +
    theme(
        strip.background = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "darkgray")
    )

ggsave(
    here("figures", "m2_dist_plot.pdf"),
    plot = m2_dist_plot,
    device = cairo_pdf,
    width = 5,
    height = 1.7,
    units = "in",
    dpi = 600
)



GeomSplitViolin <- ggproto(
    "GeomSplitViolin", GeomViolin,
    draw_group = function(
        self, data, ..., draw_quantiles = NULL
    ) {
        data <- transform(
            data,
            xminv = x - violinwidth * (x - xmin),
            xmaxv = x + violinwidth * (xmax - x)
        )

        grp <- data[1, "group"]

        newdata <- plyr::arrange(
            transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
            if (grp %% 2 == 1) y else -y
        )

        newdata <- rbind(
            newdata[1, ],
            newdata,
            newdata[nrow(newdata), ],
            newdata[1, ]
        )

        newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <-
            round(newdata[1, "x"])

        if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
            stopifnot(
                all(draw_quantiles >= 0),
                all(draw_quantiles <= 1)
            )
            quantiles <- ggplot2:::create_quantile_segment_frame(
                data, draw_quantiles
            )

            aesthetics <- data[
                rep(1, nrow(quantiles)),
                setdiff(names(data), c("x", "y")),
                drop = FALSE
            ]

            aesthetics$alpha <- rep(1, nrow(quantiles))
            both <- cbind(quantiles, aesthetics)
            quantile_grob <- GeomPath$draw_panel(both, ...)
            ggplot2:::ggname(
                "geom_split_violin",
                grid::grobTree(
                    GeomPolygon$draw_panel(newdata, ...),
                    quantile_grob
                )
            )
        } else {
                ggplot2:::ggname(
                    "geom_split_violin",
                    GeomPolygon$draw_panel(newdata, ...)
                )
        }
    }
)

geom_split_violin <- function(
    mapping = NULL, data = NULL,
    stat = "ydensity", position = "identity", ..., 
    draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
    show.legend = NA, inherit.aes = TRUE
) {
    layer(
        data = data, mapping = mapping, stat = stat,
        geom = GeomSplitViolin,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            trim = trim, scale = scale,
            draw_quantiles = draw_quantiles, na.rm = na.rm,
            ...
        )
    )
}


### Get the color
all_seg_labels <- c(
    cellprofiler = "Propagation (CellProfiler)",
    qupath = "StarDist (QuPath)",
    deepcell = "Mesmer (DeepCell)",
    cellXpress1 = "Watershed (cellXpress 2)",
    cX2_merged = "CellShape AI (cellXpress 2)",
    cellXpress2 = "CellShape AI (cellXpress 2)"
)

all_pt_colors <- rev(brewer.pal(length(all_seg_labels) - 1, "Set1"))
all_pt_colors <- c(all_pt_colors, all_pt_colors[length(all_pt_colors)])
names(all_pt_colors) <- names(all_seg_labels)

software_colors <- c(all_pt_colors[["cellXpress2"]], all_pt_colors[["deepcell"]])

PlotCellTypeDist <- function(
    mkr_name,
    y_lim = c(0,18),
    y_breaks = seq(0, 18, 2),
    region_names = c("A03", "A04", "A05", "A06"),
    is_GC_only = FALSE
) {
    ### Prepare plot data
    plot_data <- if (is_GC_only) {
        raw_data %>% filter(
            Region %in% region_names,
            `is_member:ROI:All` > 0,
            name == mkr_name
        )
    } else {
        raw_data %>% filter(
            Region %in% region_names,
            name == mkr_name
        )
    }

    ### Find the difference
    celltypes <- levels(plot_data$celltype)

    p_vals <- sapply(
        celltypes,
        function(celltype_cur) {
            a <- plot_data %>%
                filter(
                    celltype == celltype_cur,
                    software == "cellXpress 2"
                ) %>%
                pull(value)

            b <- plot_data %>%
                filter(
                    celltype == celltype_cur,
                    software == "DeepCell"
                ) %>%
                pull(value)

            cX2_median    <- median(a, na.rm = TRUE)
            mesmer_median <- median(b, na.rm = TRUE)

            if ((sum(!is.na(a)) < 5) || (sum(!is.na(b)) < 5)) {
                return(c(
                    cX = sum(!is.na(a)), m = sum(!is.na(b)), p = NA,
                    cX2_median = cX2_median, mesmer_median = mesmer_median
                ))

            } else if (abs(cX2_median - mesmer_median) < log2(1.1)) {
                return(c(
                    cX = sum(!is.na(a)), m = sum(!is.na(b)), p = NaN,
                    cX2_median = cX2_median, mesmer_median = mesmer_median
                ))
            }

            tryCatch(
                {
                    test_res <- wilcox.test(a, b, na.rm = TRUE)
                    c(
                        cX = sum(!is.na(a)), m = sum(!is.na(b)),
                        p = test_res$p.value,
                        cX2_median = cX2_median, mesmer_median = mesmer_median
                    )
                },
                error = function(cond) {
                    # Choose a return value in case of error
                    c(
                        cX = sum(!is.na(a)), m = sum(!is.na(b)), p = NA,
                        cX2_median = cX2_median, mesmer_median = mesmer_median
                    )
                }
            )
        }
    )

    p_adj <- p.adjust(p_vals["p",], "fdr")
    p_labels <- rep("n.s.", length(p_adj))
    p_labels[p_adj < 0.01] <- "*"
    p_labels[p_adj < 0.001] <- "**"
    p_labels[p_adj < 0.0001] <- "***"
    p_labels[is.na(p_adj)] <- "n.a."
    p_labels[is.nan(p_adj)] <- "s.d."

    cellnum_labels <- paste0("n=", p_vals["cX", ], ",\n", p_vals["m", ])
    p_labels <- paste0(p_labels, "\n", cellnum_labels)

    p_table <- tibble(
        celltype = factor(celltypes, levels = celltypes),
        p_labels = p_labels,
        software = "cellXpress 2"
    )

    ### Find max cell types
    message(
        "cX2 max cell type = ",
        celltypes[which.max(p_vals["cX2_median", ])]
    )
    message(
        "Mesmer max cell type = ",
        celltypes[which.max(p_vals["mesmer_median", ])]
    )

    ### Generate the plot
    p <- plot_data %>%
        ggplot(
            aes(x = celltype, y = value, fill = software)
        ) +
        geom_split_violin(
            width = 1.1,
            alpha = .4, trim = FALSE, scale = "area",
            show.legend = TRUE, na.rm = TRUE
        ) +
        geom_boxplot(
            width = .25, alpha = .6, size = 0.4,
            outlier.shape = NA, #21,
            show.legend = FALSE, na.rm = TRUE
        ) +
        geom_text(
            data = p_table, aes(y = y_lim[2], label = p_labels),
            vjust = 1, lineheight = 0.8, size = 3.2
        ) +
        scale_y_continuous(
            name = paste0(
                "mean cellular ", mkr_name, "\nintensity (log2 f.u.)"
            ),
            limits = y_lim,
            breaks = y_breaks
        ) +
        scale_fill_manual(
            name = "Tools",
            values = c(
                all_pt_colors[["cellXpress2"]],
                all_pt_colors[["deepcell"]]
            )
            #breaks = c("cellXpress 2", "DeepCell")
            #labels = seg_labels
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 1))) +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 12) +
        theme(
            axis.text = element_text(colour = "black"),
            axis.title.x = element_blank(),
            panel.grid.major.y =  element_line(colour = "darkgray")
        )

    print(p)

    ggsave(
        here(
            "figures",
            paste0(
                "CellTypeDist_", mkr_name, "_",
                paste0(region_names, collapse = ","), ".pdf"
            )
        ),
        plot = p,
        device = cairo_pdf,
        width = 8.7,
        height = 2,
        units = "in",
        dpi = 600
    )
}

PlotCellTypeDist(
    "cytokeratin", y_lim = c(2, 19), y_breaks = c(4, 8, 12),
    region_names = c("A01", "A02")
) # Epithelial cells - ns

PlotCellTypeDist(
    "CD31", y_lim = c(-3, 19), y_breaks = seq(0, 12, 4),
    region_names = c("A01", "A02")
) # Endothelial cells - cX2

PlotCellTypeDist(
    "podoplanin", y_lim = c(1, 16), y_breaks = seq(2, 10, 4),
    region_names = c("A01", "A02")
) # Lymphatic OK - ns

PlotCellTypeDist(
    "CD11c", y_lim = c(3, 19), y_breaks = seq(6, 14, 4),
) # DC - mesmer

PlotCellTypeDist(
    "CD11b", y_lim = c(3, 17), y_breaks = seq(4, 12, 4),
) # Neutrophil - ns

PlotCellTypeDist(
    "CD68", y_lim = c(2, 17), y_breaks = seq(4, 12, 4),
) # M1 - cX2. M2 - ns

PlotCellTypeDist(
    "CD163", y_lim = c(-1.5, 18), y_breaks = seq(0, 12, 4)
) # M2 - cX2

PlotCellTypeDist(
    "CD20", y_lim = c(-1, 20), y_breaks = seq(0, 15, 5),
) # B-cell - mesmer. DeepCell lymphatic very high

PlotCellTypeDist(
    "CD8", y_lim = c(1.5, 14.5), y_breaks = seq(2, 10, 4)
) # CD8 - ns

PlotCellTypeDist(
    "CD4", y_lim = c(-1, 21), y_breaks = seq(0, 15, 5)
) # CD4 - cX2, CD8 - mesmer, Treg - cX2

PlotCellTypeDist(
    "FOXP3", y_lim = c(5.5, 14.5), y_breaks = seq(6, 12, 2)
) # Treg - mesmer

#PlotCellTypeDist(
#    "CD45", y_lim = c(-1, 20), y_breaks = seq(0, 16, 4)
#)

PlotCellTypeDist(
    "PD-1",
    y_lim = c(1.5, 19), y_breaks = seq(2, 14, 4),
    is_GC_only = TRUE
)
