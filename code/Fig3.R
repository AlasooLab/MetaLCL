library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(gridExtra)


data <- read_tsv("USP18_transcriptional_targets.tsv")
data_transduction <- read_tsv("USP18_transduction_df.tsv")


# define colors
get_color <- function(value) {
  ifelse(is.na(value), "white",
         ifelse(value == 1, "gray",
                ifelse(value > 0, "red", "blue")))
}

data$"USP18 targets" <- data$"USP18 targets" * -1
data_transduction$"USP18 targets" <- data_transduction$"USP18 targets" * -1


#new column for colors
data_long <- data %>%
  gather(variable, value, -gene_name) %>%
  mutate(color = sapply(value, get_color))

data_long_transduction <- data_transduction %>%
  gather(variable, value, -gene_name) %>%
  mutate(color = sapply(value, get_color))

data_long$variable <- factor(data_long$variable, levels = c("PID", "ChEMBL","GWAS","DE","USP18 targets"))
data_long_transduction$variable <- factor(data_long_transduction$variable, levels = c("PID", "ChEMBL","GWAS","DE","USP18 targets"))


#divide data into subsets for separate plots
genes1 <- unique(data_long$gene_name)[1:19]
genes2 <- unique(data_long$gene_name)[20:38]

plot1_data <- data_long %>% filter(gene_name %in% genes1)
plot2_data <- data_long %>% filter(gene_name %in% genes2)

plot1 <- ggplot(plot1_data, aes(x = gene_name, y = variable, fill = color)) +
  geom_tile(color = "white") +
  scale_fill_identity() +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16,face = "bold"),
  )

plot2 <- ggplot(plot2_data, aes(x = gene_name, y = variable, fill = color)) +
  geom_tile(color = "white") +
  scale_fill_identity() +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16,face = "bold"),
  )

plot_transduction <- ggplot(data_long_transduction, aes(x = gene_name, y = variable, fill = color)) +
  geom_tile(color = "white") +
  scale_fill_identity() +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16,face = "bold"),
        )

combined_plot <- grid.arrange(plot_transduction,plot1, plot2, ncol = 1)
ggsave(plot = combined_plot, filename = "all_plots_Fig3.jpg", width = 12, height = 12, units = "in", dpi = 600)
