library(ggplot2)
library(dplyr)
library(sf)
library(patchwork)
library(tidyr)
library(e1071)
library(psych)
library(purrr)

# Algorithms performance #####

files <- c(list.files(path = "performance_algo/performance_real", pattern = ".rds", full.names = TRUE),
           list.files(path = "performance_algo/performance_sims", pattern = ".rds", full.names = TRUE))

ll <- lapply(files, readRDS)
names(ll) <- tools::file_path_sans_ext(basename(files))

ll <- lapply(ll, function(x) {
  x <- data.frame(method = x$method, wasserstein = x$wasserstein, ari = x$ari)
  x <- x[x$method != "PCA", ]
  x
})

methods <- c("limma", "combat", "seuratv3_cca", "seuratv3_rpca",
             "seuratv5_cca", "seuratv5_rpca", "fastmnn", "harmony", "harmony_gpu",
             "scanorama", "scvi", "scmerge2_un", "scmerge2_semi",
             "liger", "bbknn", "concord", "scgen")

# Add the missing approaches in each dataset
diff <- lapply(ll, function(x) setdiff(methods, x$method))
diff <- lapply(diff, function(x){
  x <- as.data.frame(x)
  x$method <- x$x
  x$x <- NULL
  if (dim(x)[1] == 0){
    x <- x
  }
  else {
    x$wasserstein <- 1000
    x$ari <- 0
    x
  }
})

ll <- lapply(1:length(ll), function(x) {
  ll[[x]] <- rbind(ll[[x]], diff[[x]])
})

ll <- lapply(ll, function(x) x[order(match(x$method, methods)), ])
ll <- lapply(ll, function(x){
  x$wasserstein <- log(x$wasserstein)
  x$wasserstein <- 1 - ((x$wasserstein - min(x$wasserstein)) / (max(x$wasserstein) - min(x$wasserstein)))
  x$ari[x$ari < 0] <- 0
  method <- x$method
  m <- x[, 2:3]
  geo_mean <- data.frame(apply(m, 1, geometric.mean))
  colnames(geo_mean) <- "geom_mean"
  x <- data.frame(method = method, wasserstein = x$wasserstein, ari = x$ari,
                  geo_mean = geo_mean)
})

perf <- lapply(ll, function(x) {
  method <- x$method
  geom <- t(as.matrix(x$geom_mean))
  geom <- as.data.frame(geom)
  colnames(geom) <- method
  geom
})

perf <- as.data.frame(do.call(rbind, perf))
rownames(perf) <- sub("_metrics", "", tools::file_path_sans_ext(basename(files)))
rownames(perf)[grepl("tabula_h_", rownames(perf))] <- "tabula_sapiens_eye_ovary"

# Metadata
load("results/params.RData")
load("results/params_sim.RData")

metadata <- params[params$dataset_name %in% rownames(perf), ]
metadata <- cbind(metadata, perf[metadata$dataset_name, ])

colnames(metadata)[6] <- "overexpr_prop"

metadata_sims <-  params_sim[params_sim$dataset_name %in% rownames(perf), ]
colnames(metadata_sims)[6] <- "overexpr_prop"
metadata_sims <- cbind(metadata_sims, perf[metadata_sims$dataset_name, ])
metadata_sims <- metadata_sims[, colnames(metadata)]

metadata <- rbind(metadata, metadata_sims)

# Sort batch correction methods performance
order_perf <- function(input) {
  oder <- t(apply(input, 1, function(x) {
    names(sort(x, decreasing = TRUE))
  }))
  oder <- matrix(as.factor(oder), nrow = nrow(oder),
                 ncol = ncol(oder))
  rownames(oder) <- rownames(input)
  return(oder)
}

sort_perf <- function(input) {
  perf <- t(apply(input, 1, function(x) {
    sort(x, decreasing = TRUE)
  }))
  perf <- matrix(as.factor(perf), nrow = nrow(perf),
                 ncol = ncol(perf))
  rownames(perf) <- rownames(input)
  return(perf)
}

# with scmerge2-semi #####
out1 <- sort_perf(input = perf)
order1 <- order_perf(input = perf)

opt_meth <- as.numeric(out1[, 1])
names(opt_meth) <- order1[, 1]
names(opt_meth) <- ifelse(opt_meth < 0.55, "unreliable", names(opt_meth))

# without scmerge2-semi ####

perf <- perf[, !grepl("scmerge2_semi", colnames(perf))]

out2 <- sort_perf(input = perf)
order2 <- order_perf(input = perf)

opt_meth1 <- as.numeric(out2[, 1])
names(opt_meth1) <- order2[, 1]
names(opt_meth1) <- ifelse(opt_meth1 < 0.55, "unreliable", names(opt_meth1))

opt_meth2 <- as.numeric(out2[, 2])
names(opt_meth2) <- order2[, 2]

# Visualization ####
load("results/pca.RData")
rownames(pca$rotation)[5] <- "overexpr_prop"

load("results/pca_sims.RData")
pca <- pca$x[, 1:2]
pca <- as.data.frame(pca)

pca <- rbind(pca, pca_sim)

pca <- pca[rownames(out1), ]

coords <- pca
coords$best_method <- names(opt_meth)
coords$best1_method <- names(opt_meth1)
coords$best2_method <- names(opt_meth2)
coords$worst_method <- order2[, ncol(order2)]

# Difference between the first two optimal methods
col1 <- as.numeric(out2[, 1])
col2 <- as.numeric(out2[, 2])

max <- pmax(col1, col2)

diff <- abs((col1 - col2)/max)

coords$dist <- diff
coords$sim <- grepl("sim_", rownames(pca))

pltt <- c(MetBrewer::met.brewer(name = "Hiroshige")[5],
          MetBrewer::met.brewer(name = "Hokusai1")[2],
          MetBrewer::met.brewer(name = "Hiroshige")[c(6,8)],
          MetBrewer::met.brewer(name = "Hokusai3")[2],
          MetBrewer::met.brewer(name = "Hokusai1")[4],
          MetBrewer::met.brewer(name = "Isfahan2")[1:4],
          MetBrewer::met.brewer(name = "Cassatt2")[3], "white",  "darkgreen", "#6551CC", "#CCBFFF", "blue", "#00A087")
names(pltt) <- c("bbknn", "scvi", "liger", "fastmnn", "scmerge2_un", "scanorama",
                 "seuratv5_rpca", "seuratv5_cca", "seuratv3_cca",
                 "seuratv3_rpca", "combat", "unreliable", "scmerge2_semi", "limma", "concord",
                 "harmony", "scgen")

hull <- chull(coords[, 1], coords[, 2])
hull <- c(hull, hull[1])
polygon_coords <- coords[hull, ]
poly <- polygon_coords[-11, ]

points_sf <- st_as_sf(coords, coords = c("PC1", "PC2"))

polygon_sf <- st_sfc(st_polygon(list(as.matrix(polygon_coords[, 1:2]))))
polygon_sf <- st_sf(polygon_sf)

voronoi <- st_voronoi(st_union(points_sf))
vor_sf <- st_collection_extract(voronoi)
vor_clip <- st_intersection(st_as_sf(vor_sf), polygon_sf)

vor_colored <- st_join(vor_clip, points_sf, join = st_nearest_feature)

p1 <- ggplot() +
  geom_sf(data = vor_colored, aes(fill = dist,  shape = sim), color = "black") +
  geom_sf(data = polygon_sf, fill = NA, color = "black", linewidth = 0.1) +
  geom_sf(data = points_sf, aes(shape = sim), fill = NA, color = "black", size =3) +
  scale_fill_gradient(low = "white", high = "#264DFF") +
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(fill = "Distance", shape = "Simulation") + xlab("PC1") + ylab("PC2")

p3 <-  ggplot(coords, aes(x = PC1, y = PC2, fill = best1_method, shape = sim)) +
  geom_voronoi(outline = poly, color = 1, linewidth = 0.1) +
  geom_point(size = 4) +
  scale_fill_manual(values = pltt) +
  labs(fill = "Method", shape = "Simulation") +
  theme_classic() +
  theme(text = element_text(size = 25)) + guides(fill = "none")

p5 <- ggplot(coords, aes(x = PC1, y = PC2, fill = best_method, shape = sim)) +
  geom_voronoi(outline = poly, color = 1, linewidth = 0.1) +
  geom_point(size = 4) +
  scale_fill_manual(values = pltt) +
  labs(fill = "Method", shape = "Simulation") +
  theme_classic() +
  theme(text = element_text(size = 25))

p <- p5 / p3 + plot_layout(guides = "collect")+
  plot_annotation(tag_levels = 'A')
ggsave(p, file = "figures/suppl_figure2.pdf", width = 25, height = 31, units = "cm")

## SVM #####

x <- metadata[, 2:6]
y <- as.factor(coords$best1_method)
dat <- data.frame(x, y = as.factor(y))

set.seed(369)
tune_result <- tune(svm, y ~ ., data = dat,
                    kernel = "radial", type = "C-classification",
                    ranges = list(gamma = 2^(-5:3), cost = 2^(-3:5)),
                    tunecontrol = tune.control(sampling = "cross"))

svm_best <- svm(y ~ ., data = dat, kernel = "radial",
                cost = tune_result$best.parameters$cost, type = "C-classification", gamma = 1)

table(svm_best$fitted)

opt <- coords$best1_method
names(opt) <- rownames(coords)
opt <- as.factor(opt)

names <- names(opt)[grepl("unr", opt)]

svm_best$fitted[names(svm_best$fitted) %in% names] <- "unreliable"
coords$svm_best <- svm_best$fitted
1 - sum(diag(table(coords$svm_best, coords$best1_method)))/sum(table(coords$svm_best))

res <- list(coords = coords, color = pltt, pca = pca, svm = svm_best)

save(res, file = "/home/elena/Documents/WP1/packages/BatChef/inst/extdata/svm.rda")

p <- p1 + p2 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )
ggsave("figures/figure4.pdf", p, width = 31, height =15, units = "cm")

#### with sf ##############
hull <- chull(coords[, 1], coords[, 2])
hull <- c(hull, hull[1])
polygon_coords <- coords[hull, ]
poly <- polygon_coords[-11, ]

points_sf <- st_as_sf(coords, coords = c("PC1", "PC2"))

polygon_sf <- st_sfc(st_polygon(list(as.matrix(polygon_coords[, 1:2]))))
polygon_sf <- st_sf(polygon_sf)

voronoi <- st_voronoi(st_union(points_sf))
vor_sf <- st_collection_extract(voronoi)
vor_clip <- st_intersection(st_as_sf(vor_sf), polygon_sf)

vor_colored <- st_join(vor_clip, points_sf, join = st_nearest_feature)

p2 <- ggplot() +
  geom_sf(data = vor_colored, aes(fill = svm_best), color = "black") +
  geom_sf(data = polygon_sf, fill = NA, color = "black", linewidth = 0.1) +
  geom_sf(data = points_sf, aes(shape = sim), fill = NA, color = "black", size =3) +
  scale_fill_manual(values = pltt) +
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(fill = "Method", shape = "Simulation") + xlab("PC1") + ylab("PC2")

ggsave("figures/figure_svm.jpg", p2, width = 25, height =15, units = "cm")

p1 <- ggplot() +
  geom_sf(data = vor_colored, aes(fill = dist), color = "black") +
  geom_sf(data = polygon_sf, fill = NA, color = "black", linewidth = 0.1) +
  geom_sf(data = points_sf, aes(shape = sim), fill = NA, color = "black", size =3) +
  scale_fill_gradient(low = "white", high = "#264DFF", breaks = c(0.1, 0.3, 0.5),
                      labels = c("0.1", "0.3", "0.5"))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(fill = "Delta", shape = "Simulation") + xlab("PC1") + ylab("PC2")

p <- p1 + p2 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )
ggsave("figures/figure4.pdf", p, width = 31, height =15, units = "cm")

# only real world dataset

coords_real <- coords[coords$sim == FALSE, ]

hull <- chull(coords_real[, 1], coords_real[, 2])
hull <- c(hull, hull[1])
polygon_coords <- coords_real[hull, ]
poly <- polygon_coords[-11, ]

points_sf <- st_as_sf(coords_real, coords = c("PC1", "PC2"))

polygon_sf <- st_sfc(st_polygon(list(as.matrix(polygon_coords[, 1:2]))))
polygon_sf <- st_sf(polygon_sf)

voronoi <- st_voronoi(st_union(points_sf))
vor_sf <- st_collection_extract(voronoi)
vor_clip <- st_intersection(st_as_sf(vor_sf), polygon_sf)

vor_colored <- st_join(vor_clip, points_sf, join = st_nearest_feature)

p2 <- ggplot() +
  geom_sf(data = vor_colored, aes(fill = svm_best), color = "black") +
  geom_sf(data = polygon_sf, fill = NA, color = "black", linewidth = 0.1) +
  geom_sf(data = points_sf, aes(shape = sim), fill = NA, color = "black", size =3) +
  scale_fill_manual(values = pltt) +
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(fill = "Method", shape = "Simulation") + xlab("PC1") + ylab("PC2")

ggsave("figures/figure_svm_real.jpg", p2, width = 25, height =15, units = "cm")


################################

# library(sf)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(purrr)
#
# # pivot longer
# coords_long <- pivot_longer(
#   coords,
#   cols = c(svm_best, svm_worst),
#   names_to = "name",
#   values_to = "value"
# )
#
# coords_long$name <- factor(coords_long$name,
#                            levels = c("svm_best", "svm_worst"),
#                            labels = c("Best", "Worst"))
#
# # convert to sf points
# points_sf <- st_as_sf(
#   coords_long,
#   coords = c("PC1", "PC2"),
#   crs = NA
# )
#
# polygon_sf <- st_sfc(st_polygon(list(as.matrix(polygon_coords[, 1:2]))))
# polygon_sf <- st_sf(polygon_sf)
#
# # function to build voronoi for one facet
# make_voronoi <- function(dat){
#
#   # union geometry
#   vor <- st_voronoi(st_union(dat))
#
#   # extract polygons
#   vor_sf <- st_collection_extract(vor)
#
#   # convert to sf
#   vor_sf <- st_as_sf(vor_sf)
#
#   # clip to boundary
#   vor_clip <- st_intersection(vor_sf, polygon_sf)
#
#   # assign nearest point attributes
#   vor_join <- st_join(
#     vor_clip,
#     dat,
#     join = st_nearest_feature
#   )
#
#   vor_join
# }
#
# # build per facet
# voronoi_sf <- points_sf %>%
#   split(.$name) %>%
#   map(make_voronoi) %>%
#   bind_rows()
#
# # plot
# p <- ggplot() +
#   geom_sf(
#     data = voronoi_sf,
#     aes(fill = value),
#     color = "black",
#     linewidth = 0.1
#   ) +
#   geom_sf(
#     data = points_sf,
#     aes(shape = sim),
#     size = 4
#   ) +
#   facet_wrap(~name) +
#   scale_fill_manual(values = pltt) +
#   theme_classic() +
#   theme(
#     legend.position = "bottom",
#     text = element_text(size = 30)
#   )+
#   # theme(text = element_text(size = 30)) +
#   labs(fill = "Method", shape = "Simulation")
#
# p
#
