library(dplyr)
library(psych)
library(ggplot2)
library(ggfortify)
library(patchwork)

load_params <- function(file) {
  env <- new.env()
  load(file = file, envir = env)
  as.list(env)
}

files <- list.files(path = "features/params", pattern = ".RData", full.names = TRUE)
ll <- lapply(files, load_params)
ll <- unlist(ll, recursive = FALSE)
names(ll) <- tools::file_path_sans_ext(basename(files))

params <- lapply(ll, function(x) {
  median <- as.data.frame(t(apply(x[, -1], 2, median, na.rm = TRUE)))
  df <- cbind(data.frame(dataset_name = unique(x$dataset_name)), median)
})
params <- do.call(rbind, params)

n_cells <- as.numeric(params$n_cells)
n_genes <- as.numeric(params$n_genes)

params$dropout_rate <- params$lib_loc <-params$lib_scale <- NULL
params$mean_rate <- params$mean_shape <- NULL
params$out_facloc[is.na(params$out_facloc)] <- params$out_facscale[is.na(params$out_facscale)] <- 0
params <- params[,!grepl("q10th", colnames(params))]
params <- params[,!grepl("q50th", colnames(params))]
params <- params[,!grepl("q90th", colnames(params))]
params <- params[,!grepl("out_f", colnames(params))]

# colnames(params)[7] <- "out_prob"

# Log transformation ####
skewness <- function(x) {
  n <- length(x)
  m3 <- sum((x - mean(x))^3) / n
  s3 <- sd(x)^3
  m3 / s3
}
skews <- apply(params[, -1], 2, skewness)
print(skews)

# The log-transformation is necessary
log_transf <- function(mat, base = exp(1)) {
  apply(mat, 2, function(x) {
    sapply(x, function(val) {
      if (val < 0) {
        -log((-val) + 1, base)
      } else if (val == 0) {
        0
      } else {
        log(val + 1, base)
      }
    })
  })
}

feat <- log_transf(params[, -1])

cor_mat <- abs(cor(feat))
dist_mat <- as.dist(1 - abs(cor_mat))

hc <- hclust(dist_mat)
plot(hc)

clusters <- cutree(hc, k = 5)

# For each cluster, select feature with highest variance
selected_features <- sapply(unique(clusters), function(cl) {
  slc_feat <- names(clusters[clusters == cl])
  variances <- apply(feat[, slc_feat, drop = FALSE], 2, var)
  slc_feat[which.max(variances)]
})

factoextra::fviz_nbclust(feat, FUNcluster = factoextra::hcut, method = "wss") +
  geom_vline(xintercept = 5, linetype = 2)

feat <- feat[, selected_features]
params <- params[, c("dataset_name", selected_features)]
save(params, file = "results/params.RData")

# PCA ####
pca <- prcomp(feat, scale. = TRUE, center = TRUE)
summary(pca)
save(pca, file = "results/pca.RData")
# save(pca, file = "/home/elena/Documents/WP1/packages/BatChef/inst/extdata/pca.rda")

# visualization
coords <- as.data.frame(pca$x)[, 1:2]
hull <- chull(coords[, 1], coords[, 2])
hull <- c(hull, hull[1]) 
polygon_coords <- coords[hull, ]

ggplot(coords, aes(x = PC1, y = PC2)) + geom_point(size = 3) + theme_classic() +
  geom_path(data = polygon_coords, aes(x = PC1, y = PC2), color = "black", linewidth = 0.8)


loadings <- data.frame(
  PC1 = pca$rotation[,1] * 5, 
  PC2 = pca$rotation[,2] * 5,
  Variable = rownames(pca$rotation),
  stringsAsFactors = FALSE
)

loadings$Variable[5] <- "overexpr_prop"

p <- ggplot(coords, aes(x = PC1, y = PC2)) +
  geom_point(size = 5, color = "darkgrey") +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "#920000", linewidth = 0.8) +
  geom_text(data = loadings, 
            aes(x = PC1, y = PC2, label = Variable),
            size = 10, vjust = -0.5, hjust = 0.1, color = "#920000", fontface = "bold") +
  theme_classic() + geom_path(data = polygon_coords, aes(x = PC1, y = PC2), color = "black", linewidth = 0.8) +
  theme(text = element_text(size = 25)) +
  labs(x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)"))
ggsave("figures/figure2.pdf", p, width = 30, height = 25, units = "cm")

# Supplementary Figure 1

info <- readxl::read_xlsx(path = "features/wp1_benchmark_datasets.xlsx", sheet = "Datasets")

coords$organism <- "Homo sapiens"
coords$organism[grepl("_m_", rownames(coords))] <- "Mus musculus"
coords$organism[grepl("_hm_", rownames(coords))] <- "Multiple-species"
coords$organism[grepl("_h_p_", rownames(coords))] <- "Multiple-species"
coords$organism[grepl("tm_", rownames(coords))] <- "Mus musculus"

coords$technology <- NA
coords$technology[grepl("tm_drop", rownames(coords))] <- info$Technology[grepl("tm_drop", info$`Dataset Name`)]
coords$technology[grepl("tm_brain", rownames(coords))] <- info$Technology[grepl("tm_fac", info$`Dataset Name`)]
coords$technology[grepl("tm_fac", rownames(coords))] <- info$Technology[grepl("tm_fac", info$`Dataset Name`)]
coords$technology[grepl("ts_", rownames(coords))] <- info$Technology[grepl("ts", info$`Dataset Name`)]
coords$technology[is.na(coords$technology)] <- info$Technology[!c(info$`Dataset Name` %in% c("tm_drop", "tm_fac", "ts"))]

coords$technology[grepl("Multiple", coords$technology)] <- "Multiple-platform"

coords$type <- NA
coords$type[grepl("tm_drop", rownames(coords))] <- info$Type[grepl("tm_drop", info$`Dataset Name`)]
coords$type[grepl("tm_brain", rownames(coords))] <- info$Type[grepl("tm_fac", info$`Dataset Name`)]
coords$type[grepl("tm_fac", rownames(coords))] <- info$Type[grepl("tm_fac", info$`Dataset Name`)]
coords$type[grepl("ts_", rownames(coords))] <- info$Type[grepl("ts", info$`Dataset Name`)]

coords$type[is.na(coords$type)] <- info$Type[!c(info$`Dataset Name` %in% c("tm_drop", "tm_fac", "ts"))]
coords$type[grepl("nuclei", coords$type)] <- "Single-cell/nucleus"


org_palette <- MetBrewer::met.brewer(name = "Java", n = 3)[c(1, 3, 2)]

p1 <- ggplot(coords, aes(x = PC1, y = PC2, color = organism)) + 
  geom_point(size = 4) + theme_classic() +
  scale_color_manual(values = org_palette)+
  geom_path(data = polygon_coords, aes(x = PC1, y = PC2), color = "black", linewidth = 0.8) +
  theme(text = element_text(size = 25)) + labs(color = "Organism") +
  labs(x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)"))

tech_palette <- c(MetBrewer::met.brewer(name = "Derain", n = 7)[5],
                  MetBrewer::met.brewer(name = "Archambault", n = 7)[2:7],
                  MetBrewer::met.brewer(name = "Kandinsky", n = 4)[1],
                  MetBrewer::met.brewer(name = "Cassatt2", n = 4)[3], 
                  MetBrewer::met.brewer(name = "OKeeffe1", n = 11)[c(8, 10)])

levels <- levels(as.factor(coords$technology))
p2 <- ggplot(coords, aes(x = PC1, y = PC2, color = technology)) + 
  geom_point(size = 4) + theme_classic() +
  scale_color_manual(values = tech_palette)+
  geom_path(data = polygon_coords, aes(x = PC1, y = PC2), color = "black", linewidth = 0.8) +
  theme(text = element_text(size = 25)) + labs(color = "Technology") +
  labs(x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)"))

type_palette <- c("#65CCFF", "#ebd56c", "#FF9932")

p3 <- ggplot(coords, aes(x = PC1, y = PC2, color = type)) + 
  geom_point(size = 4) + theme_classic() +
  scale_color_manual(values = type_palette)+
  geom_path(data = polygon_coords, aes(x = PC1, y = PC2), color = "black", linewidth = 0.8) +
  theme(text = element_text(size = 25)) + labs(color = "Protocol") +
  labs(x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)"))


p <- p1/p2/p3 +
  plot_annotation(tag_levels = 'A') +
  plot_layout(nrow = 3) 
ggsave(p, file = "figures/suppl_figure1.pdf", width = 25, height = 41, units = "cm")


