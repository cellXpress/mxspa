###' ===========================================================================
###' R script for ploting the benchmark results
###' Written by Lit-Hsin Loo, 2022
###' ===========================================================================
library(tidyverse)
library(RColorBrewer)
library(here)

PlotPerf <- function(
    dataset_name,
    xlim, ylim,
    filename
) {
    ### Common definitions
    all_seg_labels <- c(
        cellprofiler = "Propagation (CellProfiler)",
        qupath = "StarDist (QuPath)",
        deepcell = "Mesmer (DeepCell)",
        cellXpress1 = "Watershed (cellXpress 2)",
        cX2_merged = "CellShape AI (cellXpress 2)",
        cellXpress2 = "CellShape AI (cellXpress 2)"
    )

    all_pt_shapes <- c(
        cellprofiler = 25,
        qupath = 24,
        deepcell = 23,
        cellXpress1 = 21,
        cX2_merged = 22,
        cellXpress2 = 22
    )

    all_pt_sizes <- c(
        cellprofiler = 2,
        qupath = 2,
        deepcell = 2.5,
        cellXpress1 = 2.75,
        cX2_merged = 2.75,
        cellXpress2 = 2.75
    )

    all_pt_colors <- rev(brewer.pal(length(all_pt_sizes) - 1, "Set1"))
    all_pt_colors <- c(all_pt_colors, all_pt_colors[length(all_pt_colors)])
    names(all_pt_colors) <- names(all_pt_sizes)

    ### Process the perf_data
    plot_data <- perf_data %>%
        filter(Dataset_name == dataset_name) %>%
        group_by(Seg_method) %>%
        summarise(
            F1_score = mean(F1_score),
            MMS = mean(MMS)
        )
    print(plot_data)

    #### Get the methods (the order is very important)
    seg_methods <- names(all_seg_labels)[
        names(all_seg_labels) %in% plot_data$Seg_method
    ]

    ### Make sure the level is the same order
    plot_data <- plot_data %>%
        mutate(
            Seg_method = fct_relevel(Seg_method, seg_methods)
        )

    ### Make sure all other configs are in the same order
    seg_labels <- all_seg_labels[seg_methods]
    pt_shapes <- all_pt_shapes[seg_methods]
    pt_sizes <- all_pt_sizes[seg_methods]
    pt_colors <- all_pt_colors[seg_methods]

    ### Generate the plot
    scatter_plot <- plot_data %>%
        ggplot(
            aes(x = MMS, y = F1_score)
        ) +
        geom_point(
            aes(shape = Seg_method, fill = Seg_method, size = Seg_method)
        ) +
        scale_fill_manual(
            values = pt_colors,
            breaks = rev(seg_methods),
            labels = seg_labels
        ) +
        scale_shape_manual(
            values = pt_shapes,
            guide = "none"
        ) +
        scale_size_manual(
            values = pt_sizes,
            guide = "none"
        ) +
        guides(
            fill = guide_legend(
                override.aes = list(
                    shape = rev(pt_shapes),
                    size = rev(pt_sizes)
                )
            )
        ) +
        scale_x_continuous(
            name = "Mean matched score",
            limits = xlim
        ) +
        scale_y_continuous(
            name = "F1 score",
            limits = ylim
        ) +
        coord_fixed() +
        theme_classic() +
        theme(
            axis.text = element_text(color = "black")
        )

    ggsave(
        filename,
        plot = scatter_plot,
        width = 4, height = 2, units = "in", dpi = 400
    )
}

perf_data <- read_csv(
    here(
        "Projects", "2022_Segmentation_Benchmark", "final_results",
        "win_execution_time.csv"
    )
)
