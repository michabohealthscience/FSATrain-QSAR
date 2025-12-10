# ==============================
# Install and Load Required Packages
# ==============================

packages <- c("fingerprint", "Rtsne", "plotly", "readr", "dplyr", "data.table", "viridisLite")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) install.packages(p)
}
lapply(packages, library, character.only = TRUE)

# # ==============================
# # Load Data and Filter Valid SMILES
# # ==============================
# 
# data <- read_csv("all-molecules-modelling_and_pred.csv")
# data <- data %>% filter(!is.na(SMILES))
# 
# mols <- parse.smiles(data$SMILES)
# valid_idx <- sapply(mols, Negate(is.null))
# data <- data[valid_idx, ]
# mols <- mols[valid_idx]
# 
# # ==============================
# # Compute Fingerprints and Matrix
# # ==============================
# 
# fps <- lapply(mols, get.fingerprint, type = "extended")
# fp_mat <- fp.to.matrix(fps)
# fp_mat_unique <- unique(fp_mat)
# unique_idx <- !duplicated(fp_mat)
# data_unique <- data[unique_idx, ]
# 
# saveRDS(fp_mat_unique, "fp_mat_unique.RDS")
# saveRDS(data_unique, "data_unique.RDS")


# ==============================
# Compute t-SNE
# ==============================
fp_mat_unique <- readRDS('pre_generated/fp_mat_unique.RDS')
data_unique <- readRDS('pre_generated/data_unique.RDS')

n_samples <- nrow(fp_mat_unique)
if (n_samples < 4) stop("Too few unique molecules for t-SNE")
perplexity <- min(30, floor((n_samples - 1) / 3))
cat("Using perplexity =", perplexity, "for", n_samples, "unique samples\n")

set.seed(42)
tsne_result <- Rtsne(fp_mat_unique, dims = 2, perplexity = perplexity)

# ==============================
# Merge t-SNE with Metadata
# ==============================

tsne_df <- data.frame(
  X = tsne_result$Y[, 1],
  Y = tsne_result$Y[, 2],
  SMILES = data_unique$SMILES,
  Name = data_unique$name,
  Activity = data_unique$Activity
)

# ==============================
# Compute MW and logP for MW-logP Plot
# ==============================

# df <- fread("all-molecules-modelling_and_pred.csv")
# mols_all <- parse.smiles(df$SMILES)
# 
# 
# df$MW <- sapply(mols_all, get.exact.mass)
# df$logP <- sapply(mols_all, get.xlogp)
# df_clean <- df[!is.na(MW) & !is.na(logP)]
# 
# saveRDS(df_clean, 'df_clean.RDS')

df_clean <- readRDS('pre_generated/df_clean.RDS')

# ==============================
# Standardized Color Vector (Based on Activity)
# ==============================

# Collect unique activity values (excluding NA)
all_activities <- na.omit(c(tsne_df$Activity, df_clean$Activity))
unique_vals <- sort(unique(all_activities))
viridis_colors <- viridis(length(unique_vals))
activity_map <- setNames(viridis_colors, unique_vals)

get_color_vector <- function(activity_vec) {
  sapply(seq_along(activity_vec), function(i) {
    if (is.na(activity_vec[i])) "red" else activity_map[as.character(activity_vec[i])]
  })
}

tsne_colors <- get_color_vector(tsne_df$Activity)
mwlogp_colors <- get_color_vector(df_clean$Activity)

# ==============================
# Plot 1: t-SNE Projection (with Activity in Tooltip)
# ==============================

plot1 <- plot_ly(
  tsne_df,
  x = ~X, y = ~Y,
  type = 'scatter',
  mode = 'markers',
  text = ~paste(
    "Name:", Name,
    "<br>SMILES:", SMILES,
    "<br>Activity:", ifelse(is.na(Activity), "NA", Activity)
  ),
  marker = list(color = tsne_colors, size = 10),
  hoverinfo = 'text',
  showlegend = FALSE
) %>%
  layout(
    title = "t-SNE of Molecular Fingerprints",
    xaxis = list(title = "t-SNE 1"),
    yaxis = list(title = "t-SNE 2"),
    showlegend = FALSE
  )

# ==============================
# Plot 2: MW vs logP (with Activity in Tooltip)
# ==============================

plot2 <- plot_ly(
  df_clean,
  x = ~MW,
  y = ~logP,
  type = "scatter",
  mode = "markers",
  text = ~paste(
    "Name:", name,
    "<br>SMILES:", SMILES,
    "<br>Activity:", ifelse(is.na(Activity), "NA", Activity)
  ),
  marker = list(color = mwlogp_colors, size = 10),
  hoverinfo = 'text',
  showlegend = FALSE
) %>%
  layout(
    title = "MW vs logP",
    xaxis = list(title = "Molecular Weight"),
    yaxis = list(title = "logP"),
    showlegend = FALSE
  )

# ==============================
# Combine Plots into Interactive Dashboard
# ==============================

subplot(
  plot1, plot2,
  nrows = 1,
  margin = 0.05,
  titleX = TRUE,
  titleY = TRUE
) %>%
  layout(
    title = "Chemical Space Comparison (t-SNE vs MW-logP)",
    showlegend = FALSE
  )
