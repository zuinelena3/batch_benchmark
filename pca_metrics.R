library(ggplot2)

files <- list.files(path = "performance_pca_metrics/", pattern = ".rds", full.names = T)

metrics <- lapply(files, function(x) readRDS(x))
metrics <- lapply(metrics, function(x) {
  x$wasserstein <- log(x$wasserstein)
  x$wasserstein <- 1 - ((x$wasserstein - min(x$wasserstein)) / (max(x$wasserstein) - min(x$wasserstein)))
  x
})

bool <- sapply(metrics, function(x) any(grepl("ilisi", colnames(x))))
metrics <- metrics[bool]
metrics <- lapply(metrics, function(x) x[!is.na(x$ilisi), ])
metrics <- do.call(rbind, metrics)

metrics <- metrics[, sapply(metrics, is.numeric)]

saveRDS(metrics, file = "perfomance_algo_pca_metrics.rds")

pca <- prcomp(metrics, scale. = TRUE)

loadings <- data.frame(
  PC1 = pca$rotation[,1] * 5, 
  PC2 = pca$rotation[,2] * 5,
  Variable = rownames(pca$rotation),
  stringsAsFactors = FALSE
)

coords <- pca$x

p <- ggplot(coords, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, color = "darkgrey") +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "#920000", linewidth = 0.8) +
  geom_text(data = loadings,
            aes(x = PC1, y = PC2, label = Variable),
            size = 8, vjust = -0.5, hjust = 0.1, color = "#920000", fontface = "bold") +
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)"))
ggsave("figures/figure3.pdf", p, width = 30, height = 25, units = "cm")

