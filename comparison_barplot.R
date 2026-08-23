library(data.table)

setwd("/Users/ekim14/Downloads/Pleiotropy_Project/Final")

# Read the CSV files
file_t1 <- fread("Pleiotropy_K4_ssh_t1.csv")
file_t2 <- fread("Pleiotropy_K4_ssh_t2.csv")


library(ggplot2)

plot_data <- data.frame(
  Category = c("≥2 cancers", "≥3 cancers", "All cancers"),
  Count = c(
    length(unique(file_t1$rsid)),
    length(unique(file_t2$rsid)),
    0
  )
)

ggplot(plot_data, aes(x = Category, y = Count)) +
  geom_col(width = 0.6, fill = "steelblue") +
  geom_text(aes(label = Count), vjust = -0.3, size = 5) +
  labs(
    x = NULL,
    y = "Number of unique SNPs"
  ) +
  theme_classic(base_size = 14)

