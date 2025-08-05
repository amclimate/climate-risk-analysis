# Climate Risk Index Analysis 

This repository provides an R-based pipeline for calculating climate risk indices (Heat, Drought, Flood) from NetCDF data, using both weighted average and PCA methods.

## 🔧 Dependencies

Install the required packages:

```r
install.packages(c("ncdf4", "abind", "terra", "ClimProjDiags", "s2dv", 
                   "multiApply", "RColorBrewer", "sf", "rnaturalearth",
                   "ggplot2", "dplyr", "gt", "webshot2"))
webshot2::install_phantomjs()
# === 0. Load Required Libraries ===
library(ncdf4)
library(abind)
library(terra)
library(ClimProjDiags)
library(s2dv)
library(multiApply)
library(RColorBrewer)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(dplyr)
library(gt)
library(webshot2)

# === 1. Load NetCDF Files ===
tmax_nc <- nc_open("TMAX_2010_2024.nc")
ppt_nc  <- nc_open("PR_2010_2024.nc")

lon   <- ncvar_get(tmax_nc, "longitude")
lat   <- ncvar_get(tmax_nc, "latitude")
time  <- ncvar_get(tmax_nc, "valid_time")
dates <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")
attr(dates, "calendar") <- "proleptic_gregorian"

# === 2. Load Variables ===
tmax_data <- aperm(ncvar_get(tmax_nc, "mx2t"), c(3, 1, 2))
ppt_data  <- aperm(ncvar_get(ppt_nc,  "tp"),   c(3, 1, 2))
stopifnot(length(dates) == dim(tmax_data)[1])

# === 2.5 Convert Precipitation Units ===
ppt_data <- ppt_data * 86400  # Convert kg/m²/s to mm/day

# === 3. Detrend Temperature ===
ts_matrix <- matrix(tmax_data, nrow = dim(tmax_data)[1])
anomaly_matrix <- apply(ts_matrix, 2, function(ts) {
  DailyAno(data = ts, dates = dates, calendar = "proleptic_gregorian")
})
anomaly_data <- array(anomaly_matrix, dim = dim(tmax_data))
dim(anomaly_data) <- setNames(dim(anomaly_data), c("time", "lon", "lat"))
trend <- Trend(anomaly_data, time_dim = "time")
detrended_data <- tmax_data - (anomaly_data - trend$detrended)

# === 4. Wrap into 5D Arrays ===
wrap_5D <- function(x, var_name) {
  arr <- array(
    x,
    dim = c(model = 1, var = 1, dim(x)),
    dimnames = list(
      model = "model1",
      var   = var_name,
      time  = as.character(dates),
      lon   = as.character(lon),
      lat   = as.character(lat)
    )
  )
  dim(arr) <- setNames(dim(arr), c("model", "var", "time", "lon", "lat"))
  attr(arr, "Variables")$dat1$time <- dates
  return(arr)
}

tmax_5D      <- wrap_5D(tmax_data, "mx2t")
detrended_5D <- wrap_5D(detrended_data, "mx2t")
ppt_5D       <- wrap_5D(ppt_data, "tp")

# === 5. Climate Indices ===
heat_thresh <- Threshold(detrended_5D, qtiles = 0.9, dates = dates, calendar = "proleptic_gregorian")

heat_index <- Climdex(detrended_5D, metric = "t90p", threshold = heat_thresh,
                      dates = dates, calendar = "proleptic_gregorian", timedim = "time")

drought_index <- Climdex(ppt_5D, metric = "cdd",
                         dates = dates, calendar = "proleptic_gregorian", timedim = "time")

flood_index <- Climdex(ppt_5D, metric = "rx5day",
                       dates = dates, calendar = "proleptic_gregorian", timedim = "time")

# === 6. Standardize Indices ===
standardize_index <- function(index_result) {
  mean_vals <- apply(index_result, 1, mean, na.rm = TRUE)
  sd_vals   <- apply(index_result, 1, sd,   na.rm = TRUE)
  mean_arr  <- array(mean_vals, dim = dim(index_result))
  sd_arr    <- array(sd_vals,   dim = dim(index_result))
  (index_result - mean_arr) / sd_arr
}

HeatExtremeIndex     <- standardize_index(heat_index$result)
DroughtExtremeIndex  <- standardize_index(drought_index$result)
FloodingExtremeIndex <- standardize_index(flood_index$result)

# === 7. Combine Indices via Weighted Average ===
combined_index <- CombineIndices(
  indices = list(HeatExtremeIndex, DroughtExtremeIndex, FloodingExtremeIndex),
  weights = c(0.5, 0.3, 0.2)
)

# === 8. Save to GeoTIFF ===
save_index_geotiff <- function(index, filename) {
  mean_map <- MeanDims(index, dim = "year")[1, 1, , ]
  r <- rast(t(mean_map))
  ext(r) <- c(min(lon), max(lon), min(lat), max(lat))
  crs(r) <- "EPSG:4326"
  writeRaster(r, filename, overwrite = TRUE)
}

save_index_geotiff(HeatExtremeIndex,     "ExtremeHeatIndex.tif")
save_index_geotiff(DroughtExtremeIndex,  "DroughtIndex.tif")
save_index_geotiff(FloodingExtremeIndex, "FloodingIndex.tif")
save_index_geotiff(combined_index,       "CombinedRiskIndex.tif")

# === 9. Crop to Iran and Save Vector ===
iran_shape <- ne_countries(country = "Iran", returnclass = "sf")
r_comb     <- rast("CombinedRiskIndex.tif")
iran_vect  <- project(vect(iran_shape), crs(r_comb))
r_iran     <- mask(r_comb, iran_vect)
writeRaster(r_iran, "CombinedRiskIndex_Iran.tif", overwrite = TRUE)

# === 10. Plot Combined Map ===
breaks <- seq(-3, 3, 0.5)
colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks) - 1)

PlotEquiMap(
  MeanDims(combined_index, dim = "year"),
  lon = lon, lat = lat,
  toptitle = "Combined Climate Risk Index (Iran)",
  brks = breaks,
  filled.continents = FALSE,
  col_inf = colors[1],
  col_sup = colors[length(colors)],
  col_mid = colors[2:(length(colors) - 1)],
  colNA = "gray",
  fileout = "CombinedRiskIndex_RdYlBu.png"
)

# === 11. PCA-Based Risk Index ===
reshape_for_pca <- function(index) {
  mean_index <- MeanDims(index, dim = "year")[1, 1, , ]
  as.vector(mean_index)
}

heat_flat    <- reshape_for_pca(HeatExtremeIndex)
drought_flat <- reshape_for_pca(DroughtExtremeIndex)
flood_flat   <- reshape_for_pca(FloodingExtremeIndex)

climate_matrix <- cbind(heat_flat, drought_flat, flood_flat)
valid_rows <- complete.cases(climate_matrix) & apply(climate_matrix, 1, function(x) sd(x) > 0)
climate_matrix <- climate_matrix[valid_rows, ]

pca <- prcomp(climate_matrix, scale. = TRUE)
pca_scores <- pca$x[, 1]

pca_map <- matrix(NA, nrow = length(lon), ncol = length(lat))
flat_index <- reshape_for_pca(HeatExtremeIndex)
pca_map[!is.na(flat_index)][valid_rows] <- pca_scores

r_pca <- rast(t(pca_map))
ext(r_pca) <- c(min(lon), max(lon), min(lat), max(lat))
crs(r_pca) <- "EPSG:4326"
writeRaster(r_pca, "ClimateRisk_PCA.tif", overwrite = TRUE)

PlotEquiMap(
  array(pca_map, dim = c(1, 1, length(lon), length(lat))),
  lon = lon, lat = lat,
  toptitle = "PCA-Based Climate Risk Index",
  brks = breaks,
  filled.continents = FALSE,
  col_inf = colors[1],
  col_sup = colors[length(colors)],
  col_mid = colors[2:(length(colors) - 1)],
  colNA = "gray",
  fileout = "ClimateRisk_PCA_RdYlBu.png"
)

# === 12. Time Series Preparation ===
ts_index <- function(index, name) {
  years <- 2010:2024
  mean_ts <- sapply(1:15, function(i) {
    mean(index[i, 1, 1, , ], na.rm = TRUE)
  })
  data.frame(year = years, value = mean_ts, type = name)
}

df_heat     <- ts_index(HeatExtremeIndex, "Heat")
df_drought  <- ts_index(DroughtExtremeIndex, "Drought")
df_flood    <- ts_index(FloodingExtremeIndex, "Flood")
df_combined <- ts_index(combined_index, "Combined")

# === 13. Regression Plot Helper ===
add_regression <- function(df) {
  fit <- lm(value ~ year, data = df)
  coeff <- summary(fit)$coefficients
  intercept <- format(coeff[1, 1], digits = 2)
  slope     <- format(coeff[2, 1], digits = 2)
  r2        <- format(summary(fit)$r.squared, digits = 3)
  p_val     <- format.pval(coeff[2, 4], digits = 3, eps = 0.001)
  substitute(
    italic(y) == a + b %.% italic(x) * "," ~ 
    ~ R^2 ~ "=" ~ r2 ~ "," ~ italic(p) ~ "=" ~ pval,
    list(a = intercept, b = slope, r2 = r2, pval = p_val)
  )
}

make_ts_plot <- function(df, color, title) {
  eq <- add_regression(df)
  ggplot(df, aes(x = year, y = value)) +
    geom_line(color = color, linewidth = 1.1) +
    geom_point(color = color, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = color, linetype = "dashed") +
    annotate("text", x = min(df$year), y = max(df$value, na.rm = TRUE),
             label = as.character(as.expression(eq)), parse = TRUE, hjust = 0, size = 5) +
    labs(title = title, x = "Year", y = "Standardized Index") +
    theme_minimal(base_size = 14)
}

# === 14. Export Individual Trend Plots ===
p1 <- make_ts_plot(df_heat,    "red",    "Heat Index Trend")
p2 <- make_ts_plot(df_drought, "orange", "Drought Index Trend")
p3 <- make_ts_plot(df_flood,   "blue",   "Flood Index Trend")
p4 <- make_ts_plot(df_combined,"black",  "Combined Risk Index Trend")

ggsave("Heat_Index_Trend.png",     plot = p1, width = 8, height = 5, dpi = 300)
ggsave("Drought_Index_Trend.png",  plot = p2, width = 8, height = 5, dpi = 300)
ggsave("Flood_Index_Trend.png",    plot = p3, width = 8, height = 5, dpi = 300)
ggsave("Combined_Index_Trend.png", plot = p4, width = 8, height = 5, dpi = 300)

# === 15. Combined Time Series ===
df_all <- rbind(df_heat, df_drought, df_flood, df_combined)

ggplot(df_all, aes(x = year, y = value, color = type)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("Heat" = "red", "Drought" = "orange", "Flood" = "blue", "Combined" = "black")) +
  labs(
    title = "Yearly Standardized Climate Indices with Trend Lines (2010–2024)",
    x = "Year", y = "Standardized Index", color = "Index Type"
  ) +
  theme_minimal(base_size = 14)

ggsave("Climate_Indices_TimeSeries_With_Trend.png", width = 8, height = 5, dpi = 300)

# === 16. Export PDF of All Plots ===
pdf("Climate_Index_Trends_4Plots.pdf", width = 8, height = 5)
print(p1); print(p2); print(p3); print(p4)
dev.off()

# === 17. Export CSV Summary ===
write.csv(df_all, "Climate_Indices_TimeSeries.csv", row.names = FALSE)

# === 18. Summary Table with gt ===
summarize_index <- function(df) {
  fit <- lm(value ~ year, data = df)
  data.frame(
    type = unique(df$type),
    mean = mean(df$value, na.rm = TRUE),
    slope = coef(fit)[2],
    intercept = coef(fit)[1],
    r_squared = summary(fit)$r.squared,
    p_value = coef(summary(fit))[2, 4]
  )
}

summary_df <- bind_rows(
  summarize_index(df_heat),
  summarize_index(df_drought),
  summarize_index(df_flood),
  summarize_index(df_combined)
)

summary_df %>%
  mutate(p_value = format.pval(p_value, digits = 3, eps = 0.001)) %>%
  gt() %>%
  tab_header(title = "Summary Statistics of Climate Indices (2010–2024)") %>%
  fmt_number(columns = c(mean, slope, intercept, r_squared), decimals = 3) %>%
  cols_label(
    type = "Index Type",
    mean = "Mean",
    slope = "Trend Slope",
    intercept = "Intercept",
    r_squared = "R²",
    p_value = "p-value"
  ) %>%
  gtsave("Climate_Indices_SummaryTable.png")
