###' ===========================================================================
###' R script for ploting the benchmark results
###' Written by Lit-Hsin Loo, 2022
###' ===========================================================================
library(tidyverse)
library(RColorBrewer)
library(here)

### Common definitions
all_seg_labels <- c(
    cellprofiler = "Propagation (CellProfiler)",
    qupath = "StarDist (QuPath)",
    deepcell = "Mesmer (DeepCell)",
    cellXpress1 = "Watershed (cellXpress 1/2)",
    cX2_merged = "CellShape AI (cellXpress 2)",
    cellXpress2 = "CellShape AI (cellXpress 2)"
)

all_pt_colors <- rev(brewer.pal(length(all_seg_labels) - 1, "Set1"))
all_pt_colors <- c(all_pt_colors, all_pt_colors[length(all_pt_colors)])
names(all_pt_colors) <- names(all_seg_labels)

### Read the data
plot_data <- read_csv(
    here(
        "Projects", "2022_Segmentation_Benchmark", "final_results",
        "win_execution_time.csv"
    )
) %>% select(
    Dataset, Method, `Speed_(pixels/min)`
)

### Get the valid seg_methods
#seg_methods <- names(all_seg_labels)[
#    names(all_seg_labels) %in% plot_data$Method
#]
seg_methods <- c("cellXpress2", "cellprofiler", "qupath",  "cellXpress1")

### Make sure all other configs are in the same order
seg_labels <- all_seg_labels[seg_methods]
pt_colors <- all_pt_colors[seg_methods]

time_plot <- plot_data %>%
    mutate(
        Method = fct_relevel(Method, seg_methods)
    ) %>%
    arrange(Dataset, `Speed_(pixels/min)`) %>%
    ggplot(
        aes(x = Dataset, y = `Speed_(pixels/min)`, fill = Method)
    ) + geom_bar(
        stat = "identity",
        position = "dodge",
        color = "black"
    ) +
    scale_fill_manual(
        name = NULL,
        values = pt_colors,
        breaks = c("cellXpress2", "cellXpress1", "qupath", "cellprofiler"),
        labels = seg_labels
    ) +
    scale_x_discrete(
        labels = function(x) stringr::str_wrap(x, width = 8)
    ) +
    xlab("Image sets") +
    ylab("Processing speed\n(million pixels/min)") +
    theme_classic() +
    theme(
        legend.key.size = unit(0.15, 'in'),
        legend.background=element_blank(),
        legend.text=element_text(size = 8),
        legend.position = c(0.4, 0.7),
        axis.text = element_text(color = "black")
    )

ggsave(
    here(
        "Projects", "2022_Segmentation_Benchmark", "figures",
        "win_execution_time.png"
    ),
    plot = time_plot,
    width = 3.1, height = 2, units = "in", dpi = 400
)

ggsave(
    here(
        "Projects", "2022_Segmentation_Benchmark", "figures",
        "win_execution_time.pdf"
    ),
    plot = time_plot,
    width = 3.1, height = 2, units = "in", dpi = 400
)