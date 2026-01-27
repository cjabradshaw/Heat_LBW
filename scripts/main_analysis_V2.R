#####################################
# Stage 0: Loading libraries (v2)
#####################################

# Reproducible startup
set.seed(123)
options(stringsAsFactors = FALSE, scipen = 999)
Sys.setenv(TZ = "Australia/Sydney")

# Packages
pkgs <- c(
  "data.table","dplyr","tidyr","purrr","tibble","stringr",
  "lubridate","ggplot2","patchwork",
  "mice","splines","dlnm","glmmTMB","multcomp",
  "sf","viridis",
  "DescTools","car","grid",
  "dagitty","ggdag","ggrepel",
  "mvtnorm","readr","forcats","scales"
)

to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))


## =========================
## Stage 0.1: CONFIG (v2)
## =========================
base_dir <- "/path/to/your/project"  

dir_data_raw   <- file.path(base_dir, "data_raw")
dir_data_proc  <- file.path(base_dir, "data_processed")
dir_outputs    <- file.path(base_dir, "outputs")
dir_figures    <- file.path(base_dir, "figures")

dir.create(dir_data_raw,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_data_proc, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_outputs,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figures,   showWarnings = FALSE, recursive = TRUE)

path_lbw        <- file.path(dir_data_raw, "LBW_data.csv")
path_population <- file.path(dir_data_raw, "population_data.csv")
path_dist_births <- file.path(dir_data_raw, "dist_births_mor_mpi.csv")

stopifnot(file.exists(path_lbw), file.exists(path_population))


#####################################
# Stage 1: DAG figure (R)
#####################################


# ===== DAG (PM2.5 mediator; Wealth -> Place solid; Precip -> Tmean dotted) =====
g_primary <- dagitty('dag {
Tmean [exposure]
LBW   [outcome]
SES
Wealth
Place
Province
Trend
PM25
Humidity
Precip
Urban
U

Province -> Tmean
Province -> PM25
Province -> LBW
Province -> Urban

Trend -> Tmean
Trend -> PM25
Trend -> LBW

SES -> Urban
SES -> Place
SES -> LBW

Wealth -> Urban
Wealth -> LBW
Wealth -> Place

Urban -> PM25

U -> Tmean
U -> PM25

Tmean -> PM25
PM25  -> LBW
Tmean -> LBW
}')

# Coordinates
coordinates(g_primary) <- list(
  x = c(Tmean=0.35, LBW=0.00, SES=-1.45, Urban=-0.35, Province=1.45, Trend=-0.90,
        Humidity=-0.70, Precip=1.00, PM25=0.95, Place=-1.20, U=0.10, Wealth=-1.75),
  y = c(Tmean=0.55, LBW=-0.15, SES=-0.95, Urban=-0.70, Province=-0.05, Trend=0.95,
        Humidity=0.95,  Precip=0.80, PM25=-0.75, Place=-0.40, U=0.10, Wealth=-0.35)
)

# Build node/edge data frames
node_pos <- as.data.frame(dagitty::coordinates(g_primary)) %>%
  tibble::rownames_to_column("name")

confounders <- c("SES","Wealth","Urban","Province","Trend","Humidity","Precip")

nodes <- node_pos %>%
  filter(name != "U") %>%
  mutate(
    role = case_when(
      name == "Tmean"       ~ "Exposure",
      name == "LBW"         ~ "Outcome",
      name == "PM25"        ~ "Mediator",
      name %in% confounders ~ "Confounder",
      TRUE                  ~ "Other"
    )
  ) %>%
  left_join(tibble(
    name  = c("Tmean","LBW","SES","Wealth","Urban","Province","Trend",
              "Humidity","Precip","PM25","Place"),
    label = c("Mean Temperature (tmean)","Low Birth Weight","Maternal SES (education)",
              "Wealth Index","Urbanicity","Province","Calendar Trend (Year)",
              "Relative Humidity","Precipitation","PM2.5","Place of Delivery")
  ), by = "name")

edges_raw <- as.data.frame(dagitty::edges(g_primary)) %>%
  transmute(.from = as.character(v), .to = as.character(w))
pos_from <- node_pos %>% dplyr::select(name, x_from = x, y_from = y)
pos_to   <- node_pos %>% dplyr::select(name, x_to   = x, y_to   = y)

# Primary (solid) edges = all DAG edges except from/to U
edges_solid <- edges_raw %>%
  filter(.from != "U", .to != "U") %>%
  left_join(pos_from, by = c(".from" = "name")) %>%
  left_join(pos_to,   by = c(".to"   = "name")) %>%
  transmute(x = x_from, y = y_from, xend = x_to, yend = y_to,
            style = "Primary (solid)")

# Secondary (dotted) edges — add ONLY what you want dotted
sens_pairs <- tibble(
  .from = c("Place","Humidity","Precip","Humidity"),
  .to   = c("LBW",  "LBW",     "Tmean","Tmean")  # Precip -> Tmean requested
) %>%
  distinct(.from, .to)   # avoid duplicates
sens_pairs <- add_row(sens_pairs, .from = "Precip", .to = "LBW") %>% distinct(.from,.to)

edges_dotted <- sens_pairs %>%
  left_join(pos_from, by = c(".from" = "name")) %>%
  left_join(pos_to,   by = c(".to"   = "name")) %>%
  transmute(x = x_from, y = y_from, xend = x_to, yend = y_to,
            style = "Secondary (dotted)")

edges_all <- bind_rows(edges_solid, edges_dotted)

# ---- Plot ----
p <- ggplot() +
  geom_curve(
    data = edges_all,
    aes(x = x, y = y, xend = xend, yend = yend, linetype = style),
    curvature = 0.08, linewidth = 0.6, lineend = "round",
    arrow = arrow(type = "closed", length = unit(0.15, "in"))
  ) +
  scale_linetype_manual(
    values = c("Primary (solid)" = "solid", "Secondary (dotted)" = "dashed"),
    name = "Path type"
  ) +
  geom_point(
    data = nodes,
    aes(x = x, y = y, fill = role),
    shape = 21, size = 6, color = "black", stroke = 0.6
  ) +
  ggrepel::geom_label_repel(
    data = nodes,
    aes(x = x, y = y, label = label),
    size = 3.2, fill = "white",
    label.r = unit(0.1, "lines"),
    box.padding = 0.25, label.padding = unit(0.12, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.size = 0.2, max.overlaps = 40
  ) +
  scale_fill_manual(
    values = c(
      "Exposure"   = "#1b9e77",
      "Mediator"   = "#e78ac3",
      "Confounder" = "#7570b3",
      "Other"      = "#BDBDBD",
      "Outcome"    = "#d95f02"
    ),
    breaks = c("Exposure","Mediator","Confounder","Other","Outcome"),
    name = "Node role"
  ) +
  labs(title = "Causal DAG: Temperature and Low Birth Weight") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

# Save
ggsave("Supplementary_Figure_S1_DAG_temperature_LBW.png", p, width = 7.6, height = 4.8, dpi = 300)
ggsave("Supplementary_Figure_S1_DAG_temperature_LBW.pdf",  p, width = 7.6, height = 4.8)

#######################
#Stage 2: Data analysis
#######################

set.seed(123)

# Load data
LBW <- fread(path_lbw) %>%
  rename(province = Region)

# Ensure the date column is in Date format
LBW$date  <- as.Date(LBW$date)
LBW$Year  <- year(LBW$date)
LBW$Month <- month(LBW$date)

# Load population data and merge with the main dataset
population <- fread(path_population) %>%
  rename(province = Province)

LBW <- left_join(LBW, population, by = c("province", "Month", "Year"))

# Exclude specific provinces and filter years
LBW <- LBW %>%
  filter(!province %in% c("AJK")) %>%    # Exclude AJK
  filter(Year >= 2008 & Year <= 2017)    # Keep 2008–2017

# Convert relevant variables to factors
LBW_clean <- LBW %>%
  mutate(
    Education          = as.factor(Education),
    `Place of delivery`= as.factor(`Place of delivery`),
    `Wealth Index`     = as.factor(`Wealth Index`),
    province           = as.factor(province),
    `Mother Age`       = as.factor(`Mother Age`)
  ) %>%
  # Rename to avoid spaces/special chars before imputation
  rename(
    Place_of_delivery = `Place of delivery`,
    Wealth_Index      = `Wealth Index`,
    Mother_Age        = `Mother Age`
  )

# Quick NA counts for Urban/Rural (robust to data.table)
LBW_clean %>%
  dplyr::select(Urban_prop, Rural_prop) %>%
  is.na() %>%
  colSums()

#####################################
# Stage 3: Multiple imputation
#####################################

# Build imputation frame
impute_data <- LBW_clean %>%
  dplyr::select(Education, Place_of_delivery, Wealth_Index,
                Education_prop, Wealth_Index_prop, Urban_prop, Rural_prop)

# Optional: inspect missingness pattern
md.pattern(impute_data)

# Methods: active for Urban_prop, passive for Rural_prop
imputation_methods <- make.method(impute_data)
imputation_methods["Education"]          <- "polyreg"
imputation_methods["Place_of_delivery"]  <- "polyreg"
imputation_methods["Wealth_Index"]       <- "polyreg"
imputation_methods["Education_prop"]     <- "pmm"
imputation_methods["Wealth_Index_prop"]  <- "pmm"
imputation_methods["Urban_prop"]         <- "pmm"
# Props are 0–100 here; if you switch to 0–1 later, change to "~I(1 - Urban_prop)".
imputation_methods["Rural_prop"]         <- "~I(100 - Urban_prop)"

# Predictor matrix: passive variable shouldn't predict or be predicted
pred <- make.predictorMatrix(impute_data)
pred[, "Rural_prop"] <- 0
pred["Rural_prop", ] <- 0

imputed <- mice(
  impute_data,
  method          = imputation_methods,
  predictorMatrix = pred,
  m               = 5,
  printFlag       = FALSE
)

# Show any logged events
if (!is.null(imputed$loggedEvents) && nrow(imputed$loggedEvents) > 0) {
  print(imputed$loggedEvents)
}

# Take one completed dataset (you also pool in modeling below)
LBW_clean_imputed <- complete(imputed, 1)

# Safety net: recompute Rural deterministically & clamp
LBW_clean_imputed <- LBW_clean_imputed %>%
  mutate(
    Rural_prop = 100 - Urban_prop,
    Rural_prop = pmin(pmax(Rural_prop, 0), 100)
  )

# Replace back into LBW_clean
LBW_clean <- LBW_clean %>%
  mutate(
    Education          = LBW_clean_imputed$Education,
    Place_of_delivery  = LBW_clean_imputed$Place_of_delivery,
    Wealth_Index       = LBW_clean_imputed$Wealth_Index,
    Education_prop     = LBW_clean_imputed$Education_prop,
    Wealth_Index_prop  = LBW_clean_imputed$Wealth_Index_prop,
    Urban_prop         = LBW_clean_imputed$Urban_prop,
    Rural_prop         = LBW_clean_imputed$Rural_prop
  )

# Verify no NAs & Urban+Rural identity
LBW_clean %>%
  dplyr::select(Urban_prop, Rural_prop) %>%
  is.na() %>%
  colSums()

LBW_clean %>%
  summarise(
    min_sum = min(Urban_prop + Rural_prop, na.rm = TRUE),
    max_sum = max(Urban_prop + Rural_prop, na.rm = TRUE),
    max_dev = max(abs(Urban_prop + Rural_prop - 100), na.rm = TRUE)
  )

#############
# Stage 4: Merge precipitation and scale variables
############


# Province-level precipitation
path_precip <- file.path(dir_data_raw, "province_precipitation.csv")  
stopifnot(file.exists(path_precip))
province_precipitation <- fread(path_precip)
province_precipitation$date <- as.Date(province_precipitation$date)
province_precipitation$ADM1_EN <- as.factor(province_precipitation$ADM1_EN)

province_precipitation <- province_precipitation %>%
  rename(province = ADM1_EN) %>%
  mutate(
    province = case_when(
      province == "Azad Kashmir" ~ "AJK",
      province == "Balochistan"  ~ "Baluchistan",
      province == "Gilgit Baltistan" ~ "GB",
      province == "Khyber Pakhtunkhwa" ~ "KPK",
      TRUE ~ province
    )
  )

LBW_clean <- LBW_clean %>%
  left_join(province_precipitation, by = c("province", "date"))

# If your CSV uses 'precipitation' instead of 'total_precipitation', fall back:
if (!"total_precipitation" %in% names(LBW_clean) && "precipitation" %in% names(LBW_clean)) {
  LBW_clean$total_precipitation <- LBW_clean$precipitation
}

# Final single scaling block (keeps attributes for DLNM projections)
LBW_clean <- LBW_clean %>%
  arrange(province, date) %>%
  mutate(
    scaled_tmean         = scale(tmean_C, center = TRUE, scale = TRUE),
    scaled_humidity      = scale(humidity, center = TRUE, scale = TRUE),
    scaled_precipitation = scale(total_precipitation, center = TRUE, scale = TRUE),
    PM25_scaled          = scale(PM25, center = TRUE, scale = TRUE),
    Education_prop_scaled    = scale(Education_prop, center = TRUE, scale = TRUE),
    Wealth_Index_prop_scaled = scale(Wealth_Index_prop, center = TRUE, scale = TRUE),
    Urban_prop_scaled        = scale(Urban_prop, center = TRUE, scale = TRUE),
    Rural_prop_scaled        = scale(Rural_prop, center = TRUE, scale = TRUE)
  )

# (Optional) medians of scaled vars
scaled_medians <- LBW_clean %>%
  summarize(
    median_scaled_tmean = median(scaled_tmean, na.rm = TRUE),
    median_scaled_humidity = median(scaled_humidity, na.rm = TRUE),
    median_scaled_precipitation = median(scaled_precipitation, na.rm = TRUE),
    median_PM25_scaled = median(PM25_scaled, na.rm = TRUE),
    median_Education_prop_scaled = median(Education_prop_scaled, na.rm = TRUE),
    median_Wealth_Index_prop_scaled = median(Wealth_Index_prop_scaled, na.rm = TRUE),
    median_Urban_prop_scaled    = median(Urban_prop_scaled, na.rm = TRUE),
    median_Rural_prop_scaled    = median(Rural_prop_scaled, na.rm = TRUE)
  )
print(scaled_medians)

# Save clean dataset
write.csv(LBW_clean, "LBW_clean.csv", row.names = FALSE)
cat("The LBW_clean dataset has been successfully saved as 'LBW_clean.csv'.\n")

# Ordinal encoding and quick checks
LBW_clean <- LBW_clean %>%
  mutate(
    Education   = as.ordered(Education),
    Wealth_Index= as.ordered(Wealth_Index)
  )

library(DescTools)
library(car)
spearman_correlation <- cor(as.numeric(LBW_clean$Education),
                            as.numeric(LBW_clean$Wealth_Index),
                            method = "spearman",
                            use = "complete.obs")
cat("Spearman's rank correlation between Education and Wealth Index:", spearman_correlation, "\n")
cramers_v <- CramerV(LBW_clean$Education, LBW_clean$Wealth_Index)
cat("Cramér's V for the association between Education and Wealth Index:", cramers_v, "\n")
mod_lm <- stats::lm(LBW_count ~ Education + Wealth_Index, data = LBW_clean)
vif_results <- car::vif(mod_lm)
print(vif_results)


###############################################
# Stage 5: models with various parameters (Rubin-pooled)
###############################################

extract_cb_parameters <- function(model, cb_name) {
  coef_full <- fixef(model)$cond
  cb_coef_names <- grep(paste0("^", cb_name), names(coef_full), value = TRUE)
  cb_coefs <- coef_full[cb_coef_names]
  vcov_full <- vcov(model)$cond
  cb_vcov <- vcov_full[cb_coef_names, cb_coef_names, drop = FALSE]
  list(coefs = cb_coefs, vcov = cb_vcov)
}

predict_rr <- function(cb, coefs, vcov, data, cb_name, scaled_median = NULL) {
  tmean_percentiles <- quantile(data$scaled_tmean, probs = c(0.01, 0.10, 0.90, 0.99), na.rm = TRUE)
  cen_med <- median(data$scaled_tmean, na.rm = TRUE)
  pred <- crosspred(cb, coef = coefs, vcov = vcov, at = tmean_percentiles, cen = cen_med)
  rr_results <- data.frame(
    Percentile = names(tmean_percentiles),
    Scaled_Temperature = as.numeric(tmean_percentiles),
    RR = exp(pred$allfit),
    LCI = exp(pred$alllow),
    UCI = exp(pred$allhigh)
  )
  cat("\nRelative Risk Predictions for", cb_name, ":\n")
  print(rr_results)
  return(rr_results)
}

rubin_pool <- function(coef_list, vcov_list) {
  m <- length(coef_list)
  nm <- names(coef_list[[1]])
  stopifnot(all(vapply(coef_list, function(x) identical(names(x), nm), logical(1))))
  q_bar <- Reduce(`+`, coef_list) / m
  W <- Reduce(`+`, vcov_list) / m
  est_mat <- do.call(cbind, coef_list)
  est_centered <- sweep(est_mat, 1, q_bar, "-")
  B <- (est_centered %*% t(est_centered)) / (m - 1)
  Tmat <- W + (1 + 1/m) * B
  list(coefs = q_bar, vcov = Tmat)
}

lag <- 7
degree_exposure <- 2
degree_lag <- 1

# Crossbasis for pooled prediction (built on final LBW_clean)
cb_for_pred_linear <- crossbasis(
  LBW_clean$scaled_tmean,
  lag = lag,
  argvar = list(fun = "lin"),
  arglag = list(fun = "poly", degree = degree_lag),
  group = LBW_clean$province
)

cb_for_pred_poly <- crossbasis(
  LBW_clean$scaled_tmean,
  lag = lag,
  argvar = list(fun = "poly", degree = degree_exposure),
  arglag = list(fun = "poly", degree = degree_lag),
  group = LBW_clean$province
)

fit_one <- function(data_i, basis_type = c("linear","poly"), model_id = 2) {
  basis_type <- match.arg(basis_type)
  data_i <- data_i %>%
    arrange(province, date) %>%
    mutate(
      province = as.factor(province),
      date = as.Date(date),
      Education = as.ordered(Education),
      Wealth_Index = as.ordered(Wealth_Index)
    )
  cb <- switch(
    basis_type,
    "linear" = crossbasis(
      data_i$scaled_tmean,
      lag = lag,
      argvar = list(fun = "lin"),
      arglag = list(fun = "poly", degree = degree_lag),
      group = data_i$province
    ),
    "poly" = crossbasis(
      data_i$scaled_tmean,
      lag = lag,
      argvar = list(fun = "poly", degree = degree_exposure),
      arglag = list(fun = "poly", degree = degree_lag),
      group = data_i$province
    )
  )
  # Add models
  form <- switch(
    as.character(model_id),
    "0" = LBW_count ~ poly(PM25_scaled, 2) + ordered(Education) + ordered(Wealth_Index) +
      poly(Year, 3) + offset(log(WRA_pop)) + (1 | province),
    "1" = LBW_count ~ cb + poly(PM25_scaled, 2) + ordered(Education) + ordered(Wealth_Index) +
      poly(Year, 3) + offset(log(WRA_pop)) + (1 | province),
    "2" = LBW_count ~ cb + poly(PM25_scaled, 2) + ordered(Education) + ordered(Wealth_Index) +
      poly(Year, 3) + offset(log(WRA_pop)) + (1 | province),
    "3" = LBW_count ~ cb + poly(PM25_scaled, 2) + poly(scaled_humidity, 2) +
      ordered(Education) + ordered(Wealth_Index) +
      poly(Year, 3) + offset(log(WRA_pop)) + (1 | province),
    "4" = LBW_count ~ cb + poly(PM25_scaled, 2) + poly(scaled_precipitation, 2) +
      ordered(Education) + ordered(Wealth_Index) +
      poly(Year, 3) + offset(log(WRA_pop)) + (1 | province),
    "5" = LBW_count ~ cb + poly(PM25_scaled, 2) + poly(scaled_humidity, 2) + poly(scaled_precipitation, 2) +
      ordered(Education) + ordered(Wealth_Index) +
      poly(Year, 3) + offset(log(WRA_pop)) + (1 | province)
  )
  fam <- nbinom2()
  mod <- glmmTMB(form, data = data_i, family = fam)
  aic_i <- AIC(mod)
  if (model_id == 0) {
    return(list(cb = NULL, coefs = NULL, vcov = NULL, AIC = aic_i))
  } else {
    par <- extract_cb_parameters(mod, "cb")
    return(list(cb = cb, coefs = par$coefs, vcov = par$vcov, AIC = aic_i))
  }
}

# Fit across imputations and pool
LBW_base <- as.data.frame(LBW_clean)

fit_pool_model <- function(imputed, basis_type = c("linear","poly"), model_id = 2) {
  basis_type <- match.arg(basis_type)
  m <- imputed$m
  fit_list <- vector("list", m)
  for (i in 1:m) {
    imp_i <- complete(imputed, i)
    dat_i <- LBW_base
    stopifnot(nrow(dat_i) == nrow(imp_i))
    shared <- intersect(names(dat_i), names(imp_i))
    for (nm in shared) dat_i[[nm]] <- imp_i[[nm]]
    if ("Urban_prop" %in% names(dat_i) && "Rural_prop" %in% names(dat_i)) {
      dat_i$Rural_prop <- pmin(pmax(100 - dat_i$Urban_prop, 0), 100)
    }
    fit_list[[i]] <- fit_one(dat_i, basis_type = basis_type, model_id = model_id)
  }
  AICs <- vapply(fit_list, `[[`, numeric(1), "AIC")
  AIC_mean <- mean(AICs); AIC_sd <- sd(AICs)
  if (model_id == 0) return(list(pooled = NULL, AIC_mean = AIC_mean, AIC_sd = AIC_sd))
  coef_list <- lapply(fit_list, `[[`, "coefs")
  vcov_list <- lapply(fit_list, `[[`, "vcov")
  pooled <- rubin_pool(coef_list, vcov_list)
  list(pooled = pooled, AIC_mean = AIC_mean, AIC_sd = AIC_sd)
}

# Run models
null_res   <- fit_pool_model(imputed, basis_type = "poly",   model_id = 0)
linear_res <- fit_pool_model(imputed, basis_type = "linear", model_id = 1)
m2_res     <- fit_pool_model(imputed, basis_type = "poly",   model_id = 2)
m3_res     <- fit_pool_model(imputed, basis_type = "poly",   model_id = 3)
m4_res     <- fit_pool_model(imputed, basis_type = "poly",   model_id = 4)
m5_res     <- fit_pool_model(imputed, basis_type = "poly",   model_id = 5)

# RR tables (pooled)
if (!is.null(linear_res$pooled)) {
  linear_params <- list(coefs = linear_res$pooled$coefs, vcov = linear_res$pooled$vcov)
  linear_rr <- predict_rr(cb_for_pred_linear, linear_params$coefs, linear_params$vcov, LBW_clean, "Linear Model (Rubin)")
}
if (!is.null(m2_res$pooled)) {
  model_2_params <- list(coefs = m2_res$pooled$coefs, vcov = m2_res$pooled$vcov)
  model_2_rr <- predict_rr(cb_for_pred_poly, model_2_params$coefs, model_2_params$vcov, LBW_clean, "Model 2 (Rubin)")
}
if (!is.null(m3_res$pooled)) {
  model_3_params <- list(coefs = m3_res$pooled$coefs, vcov = m3_res$pooled$vcov)
  model_3_rr <- predict_rr(cb_for_pred_poly, model_3_params$coefs, model_3_params$vcov, LBW_clean, "Model 3 (Rubin)")
}
if (!is.null(m4_res$pooled)) {
  model_4_params <- list(coefs = m4_res$pooled$coefs, vcov = m4_res$pooled$vcov)
  model_4_rr <- predict_rr(cb_for_pred_poly, model_4_params$coefs, model_4_params$vcov, LBW_clean, "Model 4 (Rubin)")
}
if (!is.null(m5_res$pooled)) {
  model_5_params <- list(coefs = m5_res$pooled$coefs, vcov = m5_res$pooled$vcov)
  model_5_rr <- predict_rr(cb_for_pred_poly, model_5_params$coefs, model_5_params$vcov, LBW_clean, "Model 5 (Rubin)")
}

# Model selection table (all models)
aic_means <- c(null_res$AIC_mean, linear_res$AIC_mean, m2_res$AIC_mean, m3_res$AIC_mean, m4_res$AIC_mean, m5_res$AIC_mean)
names(aic_means) <- c("Null_Model","Linear_Model","Model_2","Model_3","Model_4","Model_5")
aic_sds <- c(null_res$AIC_sd, linear_res$AIC_sd, m2_res$AIC_sd, m3_res$AIC_sd, m4_res$AIC_sd, m5_res$AIC_sd)
min_aic <- min(aic_means); delta_aic <- aic_means - min_aic
weights_all <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
model_summary <- data.frame(Model = names(aic_means), AIC_mean = as.numeric(aic_means),
                            AIC_sd = as.numeric(aic_sds), Delta_AIC = as.numeric(delta_aic),
                            Weight = as.numeric(weights_all), row.names = NULL)
print(model_summary)

# ---- Model-averaged DLNM parameters (poly models only; same basis dimension)
poly_res <- list(m2 = m2_res, m3 = m3_res, m4 = m4_res, m5 = m5_res)
aic_poly <- sapply(poly_res, function(x) x$AIC_mean)
delta_poly <- aic_poly - min(aic_poly)
w_poly <- exp(-0.5 * delta_poly) / sum(exp(-0.5 * delta_poly))

nm <- names(poly_res[[1]]$pooled$coefs)
stopifnot(all(vapply(poly_res, function(x) identical(names(x$pooled$coefs), nm), logical(1))))

# ---- Model-averaged DLNM parameters (poly models only; same basis dimension)
poly_res <- list(m2 = m2_res, m3 = m3_res, m4 = m4_res, m5 = m5_res)

# (optional guard if any model failed)
# poly_res <- Filter(function(x) !is.null(x$pooled), poly_res)

aic_poly  <- sapply(poly_res, function(x) x$AIC_mean)
delta_poly <- aic_poly - min(aic_poly)
w_poly    <- exp(-0.5 * delta_poly) / sum(exp(-0.5 * delta_poly))

## model-averaged coefficients and unconditional variance
stopifnot(all(vapply(poly_res,
                     function(x) identical(names(x$pooled$coefs),
                                           names(poly_res[[1]]$pooled$coefs)),
                     TRUE)))

coefs_list <- lapply(poly_res, function(x) x$pooled$coefs)
vcov_list  <- lapply(poly_res, function(x) x$pooled$vcov)

# weighted mean of coefficients
weighted_coefs <- Reduce(`+`, Map(function(b, w) w * b, coefs_list, w_poly))

# within-model component
V_within <- Reduce(`+`, Map(function(V, w) w * V, vcov_list, w_poly))

# between-model component
V_between <- Reduce(`+`, Map(function(b, w) {
  d <- b - weighted_coefs
  w * (d %*% t(d))
}, coefs_list, w_poly))

weighted_vcov <- V_within + V_between

cat("\nWeighted Coefficients (poly models):\n"); print(weighted_coefs)
cat("\nWeighted Variance-Covariance Matrix (poly models, unconditional):\n"); print(weighted_vcov)
cat("\nWeighted Coefficients (poly models):\n"); print(weighted_coefs)
cat("\nWeighted Variance-Covariance Matrix (poly models):\n"); print(weighted_vcov)


###############################################
# Stage 6: Province-level RR curves (original °C output)
###############################################

###############################################
# Stage 6: Province-level RR curves (original °C output)
###############################################

province_full_rr_results       <- list()
province_percentile_rr_results <- list()
prov_medians                   <- list()

provinces <- unique(LBW_clean$province)

# Scaling parameters for °C back-transformation
mean_tmean <- attr(LBW_clean$scaled_tmean, "scaled:center")
sd_tmean   <- attr(LBW_clean$scaled_tmean, "scaled:scale")

for (prov in provinces) {

  cat("\n--- Processing Province:", prov, "---\n")

  prov_data <- subset(LBW_clean, province == prov)

  median_tmean_prov <- median(prov_data$scaled_tmean, na.rm = TRUE)

  percentiles <- c(0.01, 0.10, 0.90, 0.99)
  tmean_percentiles_scaled <- quantile(prov_data$scaled_tmean,
                                       probs = percentiles,
                                       na.rm = TRUE)

  temp_seq_scaled <- seq(min(prov_data$scaled_tmean, na.rm = TRUE),
                         max(prov_data$scaled_tmean, na.rm = TRUE),
                         length.out = 100)

  temp_seq_original <- temp_seq_scaled * sd_tmean + mean_tmean

  # Province median (scaled + raw) + sanity check
  median_tmean_prov_scaled <- median(prov_data$scaled_tmean, na.rm = TRUE)
  median_tmean_prov_C      <- median(prov_data$tmean_C,      na.rm = TRUE)

  stopifnot(isTRUE(all.equal(
    median_tmean_prov_C,
    median_tmean_prov_scaled * sd_tmean + mean_tmean,
    tolerance = 1e-8
  )))

  # Store medians for export
  prov_medians[[prov]] <- data.frame(
    Province              = prov,
    Median_scaled_for_cen = median_tmean_prov_scaled,
    Median_temperature_C  = median_tmean_prov_C
  )

  # Full RR curve across temperature sequence (centered at province median, scaled)
  pred_full <- crosspred(
    cb_for_pred_poly,
    coef = weighted_coefs,
    vcov = weighted_vcov,
    at   = temp_seq_scaled,
    cen  = median_tmean_prov
  )

  rr_full <- data.frame(
    Province    = prov,
    Temperature = temp_seq_original,
    Percentile  = ecdf(prov_data$tmean_C)(temp_seq_original),
    RR          = exp(pred_full$allfit),
    LCI         = exp(pred_full$alllow),
    UCI         = exp(pred_full$allhigh)
  )

  province_full_rr_results[[prov]] <- rr_full

  # RR at selected percentiles (centered at province median, scaled)
  pred_percentiles <- crosspred(
    cb_for_pred_poly,
    coef = weighted_coefs,
    vcov = weighted_vcov,
    at   = tmean_percentiles_scaled,
    cen  = median_tmean_prov
  )

  tmean_percentiles_original <- tmean_percentiles_scaled * sd_tmean + mean_tmean

  rr_at_percentiles <- data.frame(
    Province   = prov,
    Percentile = names(tmean_percentiles_scaled),
    tmean_C    = tmean_percentiles_original,
    RR         = exp(pred_percentiles$allfit),
    LCI        = exp(pred_percentiles$alllow),
    UCI        = exp(pred_percentiles$allhigh)
  )

  province_percentile_rr_results[[prov]] <- rr_at_percentiles

  cat("\nFull Array Preview:\n"); print(head(rr_full))
  cat("\nSpecific Percentiles Preview:\n"); print(rr_at_percentiles)
}

all_full_rr_results       <- do.call(rbind, province_full_rr_results)
all_percentile_rr_results <- do.call(rbind, province_percentile_rr_results)

write.csv(
  all_full_rr_results,
  file = file.path(dir_outputs, "Province_Specific_Temperature_RR_Full_Array_Scaled.csv"),
  row.names = FALSE
)

write.csv(
  all_percentile_rr_results,
  file = file.path(dir_outputs, "Province_Specific_Temperature_RR_Specific_Percentiles_Scaled.csv"),
  row.names = FALSE
)

cat("Results saved.\n")

prov_medians_df <- do.call(rbind, prov_medians)
# write.csv(prov_medians_df,
#           file = file.path(dir_outputs, "Province_Median_Temperature_Reference.csv"),
#           row.names = FALSE)

# 90th percentile RR per province
rr_90th_percentile <- all_percentile_rr_results %>%
  filter(Percentile == "90%") %>%
  dplyr::select(Province, RR, LCI, UCI)

print(rr_90th_percentile)

# 99th percentile RR per province
rr_99th_percentile <- all_percentile_rr_results %>%
  filter(Percentile == "99%") %>%
  dplyr::select(Province, RR, LCI, UCI)

print(rr_99th_percentile)

write.csv(
  rr_99th_percentile,
  file = file.path(dir_outputs, "Province_RR.csv"),
  row.names = FALSE
)

cat("\nRR at the 99th percentile saved as 'Province_RR.csv'.\n")


#####################
# Stage 7: Province-level RR figure (curve + histogram)
#####################

# Ensure proper column naming
all_full_rr_results <- all_full_rr_results %>% rename(province = Province)

# Compute AlignTemp and percentiles on the curve grid
all_full_rr_results <- all_full_rr_results %>%
  group_by(province) %>%
  mutate(
    AlignTemp = Temperature[which.min(abs(RR - 1) + (UCI - LCI))],
    P1  = quantile(Temperature, 0.01, na.rm = TRUE),
    P10 = quantile(Temperature, 0.10, na.rm = TRUE),
    P90 = quantile(Temperature, 0.90, na.rm = TRUE),
    P99 = quantile(Temperature, 0.99, na.rm = TRUE)
  ) %>%
  ungroup()

# Histogram data from observed baseline temps
hist_data <- LBW_clean %>%
  group_by(province) %>%
  mutate(tmean_C_bin = cut(tmean_C, breaks = seq(min(tmean_C, na.rm = TRUE), max(tmean_C, na.rm = TRUE), by = 1))) %>%
  count(province, tmean_C_bin) %>%
  ungroup() %>%
  mutate(tmean_C_mid = as.numeric(sub("\\((.+),.*", "\\1", as.character(tmean_C_bin))) + 0.5)

combined_plots <- list()

for (prov in unique(all_full_rr_results$province)) {
  rr_data <- all_full_rr_results %>% filter(province == prov)
  hist_data_subset <- hist_data %>% filter(province == prov)
  max_count <- max(hist_data_subset$n, na.rm = TRUE)
  mid_count <- round(max_count / 2)
  y_breaks <- c(0, mid_count, max_count)
  
  risk_plot <- ggplot(rr_data, aes(x = Temperature, y = RR)) +
    geom_line(data = subset(rr_data, Temperature <= AlignTemp), aes(color = "Cold"), size = 1) +
    geom_line(data = subset(rr_data, Temperature >  AlignTemp), aes(color = "Hot"),  size = 1) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), fill = "grey", alpha = 0.2) +
    geom_vline(aes(xintercept = AlignTemp), linetype = "solid", size = 0.8, color = "black") +
    geom_vline(aes(xintercept = P1),  linetype = "dotted", size = 0.5, color = "black") +
    geom_vline(aes(xintercept = P10), linetype = "dotted", size = 0.5, color = "black") +
    geom_vline(aes(xintercept = P90), linetype = "dotted", size = 0.5, color = "black") +
    geom_vline(aes(xintercept = P99), linetype = "dotted", size = 0.5, color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.5, color = "black") +
    labs(title = prov, x = NULL, y = "RR") +
    scale_color_manual(values = c("Cold" = "blue", "Hot" = "red")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
  
  hist_plot <- ggplot(hist_data_subset, aes(x = tmean_C_mid, y = n)) +
    geom_bar(stat = "identity", fill = "lightblue", alpha = 0.6, color = "lightblue") +
    scale_y_continuous(breaks = y_breaks, labels = scales::comma) +
    labs(x = "Temperature (°C)", y = "Months") +
    theme_minimal(base_size = 14) +
    theme(axis.title.y = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 12),
          axis.text.y  = element_text(color = "black"))
  
  combined_plot <- risk_plot / hist_plot + plot_layout(heights = c(3, 1))
  combined_plots[[prov]] <- combined_plot
}

final_combined_plot <- wrap_plots(combined_plots, ncol = 2)
ggsave("Figure1_Final_Province_ER_Curves_Scaled.png", final_combined_plot, width = 16, height = 12)
print(final_combined_plot)



# Small helper: build lines & ribbons on the ORIGINAL °C scale
generate_plot_data <- function(cb, coefs, vcov, data, model_name, color_hex) {
  # scaling params from the scaled variable
  mean_tmean <- attr(data$scaled_tmean, "scaled:center")
  sd_tmean   <- attr(data$scaled_tmean, "scaled:scale")
  stopifnot(is.finite(mean_tmean), is.finite(sd_tmean))
  
  # sequence in original °C, then transform to scaled
  temp_seq <- seq(min(data$tmean_C, na.rm = TRUE),
                  max(data$tmean_C, na.rm = TRUE),
                  length.out = 200)
  temp_seq_scaled <- (temp_seq - mean_tmean) / sd_tmean
  
  # center at median of the scaled variable
  cen_scaled <- median(data$scaled_tmean, na.rm = TRUE)
  
  pr <- dlnm::crosspred(cb, coef = coefs, vcov = vcov,
                        at = temp_seq_scaled, cen = cen_scaled)
  
  tibble::tibble(
    Temperature = temp_seq,                     # original °C
    RR  = as.numeric(exp(pr$allfit)),
    LCI = as.numeric(exp(pr$alllow)),
    UCI = as.numeric(exp(pr$allhigh)),
    Model = model_name,
    Color = color_hex
  )
}

# Collect available models safely (skip NULL ones)
avail <- list()

if (exists("linear_params") && !is.null(linear_params)) {
  avail[["Linear Model"]] <- list(cb = cb_for_pred_linear,
                                  coefs = linear_params$coefs,
                                  vcov = linear_params$vcov,
                                  col = "#1f77b4")
}
if (exists("model_2_params") && !is.null(model_2_params)) {
  avail[["Model 2"]] <- list(cb = cb_for_pred_poly,
                             coefs = model_2_params$coefs,
                             vcov = model_2_params$vcov,
                             col = "#d62728")
}
if (exists("model_3_params") && !is.null(model_3_params)) {
  avail[["Model 3"]] <- list(cb = cb_for_pred_poly,
                             coefs = model_3_params$coefs,
                             vcov = model_3_params$vcov,
                             col = "#2ca02c")
}
if (exists("model_4_params") && !is.null(model_4_params)) {
  avail[["Model 4"]] <- list(cb = cb_for_pred_poly,
                             coefs = model_4_params$coefs,
                             vcov = model_4_params$vcov,
                             col = "#9467bd")
}
if (exists("model_5_params") && !is.null(model_5_params)) {
  avail[["Model 5"]] <- list(cb = cb_for_pred_poly,
                             coefs = model_5_params$coefs,
                             vcov = model_5_params$vcov,
                             col = "#ff7f0e")
}

if (length(avail) == 0) stop("No available model parameters to plot.")

# Build plotting data
plot_data <- dplyr::bind_rows(lapply(names(avail), function(nm) {
  spec <- avail[[nm]]
  generate_plot_data(spec$cb, spec$coefs, spec$vcov,
                     LBW_clean, nm, spec$col)
}))

# Named color vectors for ggplot
cols <- stats::setNames(unique(plot_data$Color), unique(plot_data$Model))

combined_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Temperature, y = RR, color = Model, fill = Model)) +
  ggplot2::geom_line(size = 1.1) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = LCI, ymax = UCI), alpha = 0.20, color = NA) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
  ggplot2::labs(
    title = "Combined Exposure–Response Curves (original °C)",
    x = "Temperature (°C)",
    y = "Relative Risk (RR)"
  ) +
  ggplot2::scale_color_manual(values = cols) +
  ggplot2::scale_fill_manual(values = cols) +
  ggplot2::theme_classic(base_size = 12)

ggplot2::ggsave("Combined_Exposure_Response_Original_Scale_2025.png", combined_plot, width = 10, height = 8, dpi = 300)
print(combined_plot)


## ===============================
## Stage 8: AF estimation & projections
## ===============================


# --- Safety checks
stopifnot(exists("LBW_clean"), exists("all_full_rr_results"),
          exists("weighted_coefs"), exists("weighted_vcov"))


# --- READ climate projection files (make sure these paths match your files) ---
canonise_prov <- function(x) {
  xu <- toupper(trimws(as.character(x)))
  dplyr::recode(
    xu,
    "KHYBER PAKHTUNKHWA" = "KPK",
    "GILGIT BALTISTAN"   = "GB",
    "AZAD KASHMIR"       = "AJK",
    "BALOCHISTAN"        = "Baluchistan",
    "PUNJAB"             = "Punjab",
    "SINDH"              = "Sindh",
    "ISLAMABAD CAPITAL TERRITORY" = "ICT",
    "CAPITAL TERRITORY"  = "ICT",
    "ISLAMABAD"          = "ICT",
    .default = xu
  )
}

# Helper to read and standardise a projection CSV (expects province/date/tmean; humidity is optional)
read_proj <- function(path) {
  df <- fread(path)
  # Try to find province column and standardise the name
  nm <- names(df)
  prov_col <- if ("province" %in% nm) "province" else if ("ADM1_EN" %in% nm) "ADM1_EN" else if ("ADM1_EN.x" %in% nm) "ADM1_EN.x" else stop("No province column in: ", path)
  setnames(df, prov_col, "province")
  # Date
  date_col <- if ("date" %in% nm) "date" else if ("date.y" %in% nm) "date.y" else grep("^date", nm, value = TRUE)[1]
  if (is.na(date_col)) stop("No date column in: ", path)
  setnames(df, date_col, "date")
  df[, date := as.Date(date)]
  # Temperature
  if (!"tmean_C" %in% names(df)) {
    if ("mean_C" %in% names(df)) setnames(df, "mean_C", "tmean_C") else stop("Need tmean_C/mean_C in: ", path)
  }
  # Standardise province values to study naming
  df[, province := factor(dplyr::recode(province,
                                        "Balochistan"="Baluchistan",
                                        "Azad Kashmir"="AJK",
                                        "Khyber Pakhtunkhwa"="KPK",
                                        "Gilgit Baltistan"="GB",
                                        .default = as.character(province)
  ))]
  df[]
}

# ==== EDIT THESE PATHS ====
path_RCP45_48 <- "prov_RCP4.5_48.csv"
path_RCP85_48 <- "prov_RCP8.5_48.csv"
path_RCP45_68 <- "prov_RCP4.5_68.csv"
path_RCP85_68 <- "prov_RCP8.5_68.csv"

prov_RCP4.5_48 <- read_proj(path_RCP45_48)
prov_RCP8.5_48 <- read_proj(path_RCP85_48)
prov_RCP4.5_68 <- read_proj(path_RCP45_68)
prov_RCP8.5_68 <- read_proj(path_RCP85_68)

# quick sanity prints
lapply(list(prov_RCP4.5_48, prov_RCP8.5_48, prov_RCP4.5_68, prov_RCP8.5_68), function(x) {
  cat(unique(x$province), "\n"); print(summary(x$tmean_C))
})


# --- Crossbasis matching your main spec
cb_for_pred_poly <- crossbasis(
  LBW_clean$scaled_tmean,
  lag   = 7,
  argvar = list(fun = "poly", degree = 2),
  arglag = list(fun = "poly", degree = 1),
  group = LBW_clean$province
)

# --- Clean RR grid and standardized names (kept as-is; used for province filtering later)
if ("Province" %in% names(all_full_rr_results) && !"province" %in% names(all_full_rr_results)) {
  all_full_rr_results <- all_full_rr_results %>% dplyr::rename(province = Province)
}
stopifnot(all(c("province","Temperature","RR","LCI","UCI") %in% names(all_full_rr_results)))

rr_clean <- all_full_rr_results %>%
  mutate(
    Temperature = as.numeric(Temperature),
    RR  = as.numeric(RR),
    LCI = as.numeric(LCI),
    UCI = as.numeric(UCI)
  ) %>%
  filter(is.finite(Temperature) & is.finite(RR) & is.finite(LCI) & is.finite(UCI)) %>%
  group_by(province) %>%
  arrange(Temperature, .by_group = TRUE) %>%
  ungroup()

# --- Province ALIGN temp for AF (UPDATED: use BASELINE MEDIAN °C)
align_tbl <- LBW_clean %>%
  group_by(province) %>%
  summarise(AlignTemp = median(tmean_C, na.rm = TRUE), .groups = "drop")

# province-specific centers in scaled space (unchanged)
prov_centers <- LBW_clean %>%
  group_by(province) %>%
  summarise(cen_scaled = median(scaled_tmean, na.rm = TRUE), .groups = "drop")

# Baseline scaling params
mean_s <- attr(LBW_clean$scaled_tmean, "scaled:center")
sd_s   <- attr(LBW_clean$scaled_tmean, "scaled:scale")

# --- Helpers (robust)
safe_crosspred <- function(cb, coefs, vcov, at, cen) {
  at <- at[is.finite(at) & !is.na(at)]
  if (length(at) == 0) return(NULL)
  if (length(cen) != 1 || !is.finite(cen) || is.na(cen)) {
    stop("safe_crosspred(): 'cen' must be a single finite number.")
  }
  suppressWarnings({
    cp <- try(crosspred(cb, coef = coefs, vcov = vcov, at = at, cen = cen), silent = TRUE)
  })
  if (inherits(cp, "try-error")) stop("safe_crosspred(): crosspred failed.")
  cp
}

# MC simulator in log space using crosspred()’s SEs
simulate_af_log <- function(mu_log, se_log, pexp, n = 1000) {
  mu_log <- as.numeric(mu_log)
  se_log <- as.numeric(se_log)
  se_log[!is.finite(se_log) | se_log <= 0] <- 1e-6
  pexp <- ifelse(is.finite(pexp), pexp, 0)
  
  sims <- replicate(n, {
    rr_draw <- exp(rnorm(length(mu_log), mean = mu_log, sd = se_log))
    mean(((rr_draw - 1) / rr_draw) * pexp * 100, na.rm = TRUE)
  })
  sims <- sims[is.finite(sims)]
  if (!length(sims)) return(c(AF = NA_real_, LCI = NA_real_, UCI = NA_real_))
  c(AF = mean(sims), LCI = unname(quantile(sims, 0.025)), UCI = unname(quantile(sims, 0.975)))
}


af_from_cp <- function(cp, prop, scale100 = TRUE) {
  rr <- exp(cp$allfit)
  af_day <- (rr - 1) / rr
  out <- prop * mean(af_day, na.rm = TRUE)
  if (scale100) out <- 100 * out
  out
}
# Core AF calculator (returns a list with heat/cold named vectors)
compute_AF_for_dataset <- function(df, province_name, cb_pred, coefs, V,
                                   align_temp, cen_scaled, n_mc = 1000) {
  
  dprov <- df %>% dplyr::filter(province == province_name)
  if (nrow(dprov) == 0) {
    return(list(heat = c(AF = NA, LCI = NA, UCI = NA),
                cold = c(AF = NA, LCI = NA, UCI = NA)))
  }
  
  prop_heat <- mean(dprov$tmean_C >  align_temp, na.rm = TRUE)
  prop_cold <- mean(dprov$tmean_C <= align_temp, na.rm = TRUE)
  
  at_hot_scaled  <- ((dprov$tmean_C[dprov$tmean_C >  align_temp] - mean_s) / sd_s)
  at_cold_scaled <- ((dprov$tmean_C[dprov$tmean_C <= align_temp] - mean_s) / sd_s)
  
  if (length(cen_scaled) != 1 || !is.finite(cen_scaled) || is.na(cen_scaled)) {
    warning("compute_AF_for_dataset(): bad 'cen_scaled' for ", province_name,
            ". Using overall median of baseline scaled_tmean.")
    cen_scaled <- median(LBW_clean$scaled_tmean, na.rm = TRUE)
  }
  
  use_mvn <- TRUE   # toggle here
  
  # HOT
  if (length(at_hot_scaled) > 0) {
    if (use_mvn) {
      draws <- rmvnorm(n_mc, mean = coefs, sigma = V)
      sims <- apply(draws, 1, function(b) {
        cp <- dlnm::crosspred(cb_pred, coef = b, vcov = V, at = at_hot_scaled, cen = cen_scaled)
        af_from_cp(cp, prop_heat)
      })
      heat_res <- c(
        AF  = mean(sims, na.rm = TRUE),
        LCI = unname(quantile(sims, 0.025, na.rm = TRUE, names = FALSE)),
        UCI = unname(quantile(sims, 0.975, na.rm = TRUE, names = FALSE))
      )
      
    } else {
      pr_hot   <- safe_crosspred(cb_pred, coefs, V, at = at_hot_scaled, cen = cen_scaled)
      heat_res <- simulate_af_log(pr_hot$allfit, pr_hot$allse, prop_heat, n_mc)
    }
  } else heat_res <- c(AF = 0, LCI = 0, UCI = 0)
  
  # COLD (same pattern)
  if (length(at_cold_scaled) > 0) {
    if (use_mvn) {
      draws <- rmvnorm(n_mc, mean = coefs, sigma = V)
      sims <- apply(draws, 1, function(b) {
        cp <- dlnm::crosspred(cb_pred, coef = b, vcov = V, at = at_cold_scaled, cen = cen_scaled)
        af_from_cp(cp, prop_cold)
      })
      cold_res <- c(
        AF  = mean(sims, na.rm = TRUE),
        LCI = unname(quantile(sims, 0.025, na.rm = TRUE, names = FALSE)),
        UCI = unname(quantile(sims, 0.975, na.rm = TRUE, names = FALSE))
      )
    } else {
      pr_cold  := safe_crosspred(cb_pred, coefs, V, at = at_cold_scaled, cen = cen_scaled)
      cold_res <- simulate_af_log(pr_cold$allfit, pr_cold$allse, prop_cold, n_mc)
    }
  } else cold_res <- c(AF = 0, LCI = 0, UCI = 0)
  return(list(
    heat = heat_res,
    cold = cold_res
  ))
}  


# --- Baseline AF (now using MEDIAN °C threshold)
baseline_AF <- bind_rows(lapply(unique(LBW_clean$province), function(prov) {
  align <- align_tbl$AlignTemp[align_tbl$province == prov]
  cen_s <- prov_centers$cen_scaled[prov_centers$province == prov]
  
  res <- compute_AF_for_dataset(
    df         = LBW_clean %>% dplyr::select(province, date, tmean_C),
    province_name = prov,
    cb_pred    = cb_for_pred_poly,
    coefs      = weighted_coefs,
    V          = weighted_vcov,
    align_temp = align,
    cen_scaled = cen_s,
    n_mc       = 1000
  )
  if (is.null(res)) return(NULL)
  tibble(
    Province = prov,
    AF_Heat = res$heat["AF"],  LCI_Heat = res$heat["LCI"],  UCI_Heat = res$heat["UCI"],
    AF_Cold = res$cold["AF"],  LCI_Cold = res$cold["LCI"],  UCI_Cold = res$cold["UCI"],
    Scenario = "NA",
    Period   = "Baseline"
  )
}))

# --- Projections: standardize to (province, date, tmean_C)
prep_proj <- function(df) {
  df <- as.data.frame(df)
  nn <- names(df)
  
  if ("province" %in% nn) {
    # ok
  } else if ("ADM1_EN" %in% nn) {
    names(df)[names(df) == "ADM1_EN"] <- "province"
  } else if ("ADM1_EN.x" %in% nn) {
    names(df)[names(df) == "ADM1_EN.x"] <- "province"
  } else {
    stop("prep_proj(): can't find province/ADM1_EN")
  }
  
  if ("date" %in% nn) {
    # ok
  } else if ("date.y" %in% nn) {
    names(df)[names(df) == "date.y"] <- "date"
  } else {
    cand <- grep("^date", nn, value = TRUE)
    if (length(cand) >= 1) names(df)[names(df) == cand[1]] <- "date" else
      stop("prep_proj(): can't find a date column")
  }
  
  if (!inherits(df$date, "Date")) df$date <- as.Date(df$date)
  
  if (!("tmean_C" %in% names(df))) stop("prep_proj(): need 'tmean_C'")
  
  df[, c("province","date","tmean_C"), drop = FALSE]
}

proj_list <- list(
  RCP4.5_48 = prep_proj(prov_RCP4.5_48),
  RCP8.5_48 = prep_proj(prov_RCP8.5_48),
  RCP4.5_68 = prep_proj(prov_RCP4.5_68),
  RCP8.5_68 = prep_proj(prov_RCP8.5_68)
)

# --- Scenario AF (median threshold carried through)
proj_AF <- bind_rows(lapply(names(proj_list), function(nm) {
  df <- proj_list[[nm]] %>%
    dplyr::filter(province %in% unique(rr_clean$province))
  
  bind_rows(lapply(unique(df$province), function(prov) {
    cen_s <- prov_centers$cen_scaled[prov_centers$province == prov]
    if (!length(cen_s) || !is.finite(cen_s)) cen_s <- median(LBW_clean$scaled_tmean, na.rm = TRUE)
    
    align <- align_tbl$AlignTemp[align_tbl$province == prov]
    res <- compute_AF_for_dataset(
      df         = df,
      province_name = prov,
      cb_pred    = cb_for_pred_poly,
      coefs      = weighted_coefs,
      V          = weighted_vcov,
      align_temp = align,
      cen_scaled = cen_s,
      n_mc       = 1000
    )
    if (is.null(res)) return(NULL)
    tibble(
      Province = prov,
      AF_Heat = res$heat["AF"],  LCI_Heat = res$heat["LCI"],  UCI_Heat = res$heat["UCI"],
      AF_Cold = res$cold["AF"],  LCI_Cold = res$cold["LCI"],  UCI_Cold = res$cold["UCI"],
      Scenario = sub("_.*","", nm),
      Period   = ifelse(grepl("48", nm), "2048–2057", "2068–2077")
    )
  }))
}))

# --- Combine
all_combined_af_results <- bind_rows(baseline_AF, proj_AF) %>%
  mutate(
    Period   = factor(Period,   levels = c("Baseline","2048–2057","2068–2077")),
    Scenario = factor(Scenario, levels = c("NA","RCP4.5","RCP8.5"))
  )


# If you want to exclude AJK:
# all_combined_af_results <- all_combined_af_results %>% filter(Province != "AJK")

cat("AF calculation complete.\n")
print(head(all_combined_af_results))


# 0) Output folder
out_dir <- "outputs_af"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1) Tidy + round for export
round_cols <- c("AF_Heat","LCI_Heat","UCI_Heat",
                "AF_Cold","LCI_Cold","UCI_Cold")

af_all <- all_combined_af_results %>%
  mutate(
    across(all_of(round_cols), ~ round(as.numeric(.x), 3)),
    Province = as.character(Province),
    Period   = as.character(Period),
    Scenario = as.character(Scenario)
  ) %>%
  relocate(Province, Period, Scenario)

# 2) Full table (heat + cold)
write.csv(af_all,
          file.path(out_dir, "Combined_AF_Results.csv"),
          row.names = FALSE)

# 3) Heat-only table (like before)
af_heat_only <- af_all %>%
  dplyr::select(Province, Period, Scenario,
                AF_Heat, LCI_Heat, UCI_Heat)

write.csv(af_heat_only,
          file.path(out_dir, "Final_AF_Combined_Heat_Only.csv"),
          row.names = FALSE)

# 4) Optional: burden vs prevented fraction for heat
af_heat_bp <- af_all %>%
  mutate(
    AF_Heat_Burden = pmax(AF_Heat, 0),   # set negative AFs to 0 for burden
    PF_Heat        = pmax(-AF_Heat, 0)   # prevented fraction from heat
  ) %>%
  dplyr::select(Province, Period, Scenario, AF_Heat_Burden, PF_Heat)

write.csv(af_heat_bp,
          file.path(out_dir, "AF_Heat_Burden_and_Prevented.csv"),
          row.names = FALSE)

# 5) Optional: versions excluding AJK
af_all_noAJK <- af_all %>% filter(Province != "AJK")
write.csv(af_all_noAJK,
          file.path(out_dir, "Combined_AF_Results_noAJK.csv"),
          row.names = FALSE)

af_heat_only_noAJK <- af_heat_only %>% filter(Province != "AJK")
write.csv(af_heat_only_noAJK,
          file.path(out_dir, "Final_AF_Combined_Heat_Only_noAJK.csv"),
          row.names = FALSE)

cat("AF CSVs saved in: ", normalizePath(out_dir), "\n")

# --- Sensitivity: use MIN-RISK/MMT (from RR curves) as the heat/cold threshold ---

# 1) Build MMT/min-risk align table from the RR grid
align_tbl_mmt <- rr_clean %>%
  dplyr::group_by(province) %>%
  dplyr::summarise(
    AlignTemp = Temperature[which.min(abs(RR - 1) + (UCI - LCI))],
    .groups = "drop"
  )

# 2) Helper (unchanged)
run_af_with_align <- function(align_tbl_in) {
  # Baseline
  base <- dplyr::bind_rows(lapply(unique(LBW_clean$province), function(prov) {
    align <- align_tbl_in$AlignTemp[align_tbl_in$province == prov]
    cen_s <- prov_centers$cen_scaled[prov_centers$province == prov]
    res   <- compute_AF_for_dataset(
      df         = dplyr::select(LBW_clean, province, date, tmean_C),
      province_name = prov,
      cb_pred    = cb_for_pred_poly,
      coefs      = weighted_coefs, V = weighted_vcov,
      align_temp = align, cen_scaled = cen_s, n_mc = 1000
    )
    if (is.null(res)) return(NULL)
    tibble::tibble(
      Province = prov,
      AF_Heat = res$heat["AF"],  LCI_Heat = res$heat["LCI"],  UCI_Heat = res$heat["UCI"],
      AF_Cold = res$cold["AF"],  LCI_Cold = res$cold["LCI"],  UCI_Cold = res$cold["UCI"],
      Scenario = "NA", Period = "Baseline"
    )
  }))
  
  # Projections
  proj <- dplyr::bind_rows(lapply(names(proj_list), function(nm) {
    dfp <- proj_list[[nm]] %>% dplyr::filter(province %in% unique(LBW_clean$province))
    dplyr::bind_rows(lapply(unique(dfp$province), function(prov) {
      align <- align_tbl_in$AlignTemp[align_tbl_in$province == prov]
      cen_s <- prov_centers$cen_scaled[prov_centers$province == prov]
      res   <- compute_AF_for_dataset(
        df = dfp, province_name = prov,
        cb_pred = cb_for_pred_poly, coefs = weighted_coefs, V = weighted_vcov,
        align_temp = align, cen_scaled = cen_s, n_mc = 1000
      )
      if (is.null(res)) return(NULL)
      tibble::tibble(
        Province = prov,
        AF_Heat = res$heat["AF"],  LCI_Heat = res$heat["LCI"],  UCI_Heat = res$heat["UCI"],
        AF_Cold = res$cold["AF"],  LCI_Cold = res$cold["LCI"],  UCI_Cold = res$cold["UCI"],
        Scenario = sub("_.*","", nm),
        Period   = ifelse(grepl("48", nm), "2048–2057", "2068–2077")
      )
    }))
  }))
  
  dplyr::bind_rows(base, proj)
}

# 3) Run MMT-threshold AF and compare to main (median-based)
af_mmt_thresh <- run_af_with_align(align_tbl_mmt)

af_compare <- dplyr::left_join(
  all_combined_af_results %>%
    dplyr::select(Province, Period, Scenario, AF_Heat_MED = AF_Heat, AF_Cold_MED = AF_Cold),
  af_mmt_thresh %>%
    dplyr::select(Province, Period, Scenario, AF_Heat_MMT = AF_Heat, AF_Cold_MMT = AF_Cold),
  by = c("Province","Period","Scenario")
) %>%
  dplyr::mutate(
    Diff_Heat = AF_Heat_MED - AF_Heat_MMT,
    Diff_Cold = AF_Cold_MED - AF_Cold_MMT
  )

print(af_compare)


# --- 1) Prep & relabel for plotting
af_df <- all_combined_af_results %>%
  mutate(
    Scenario = as.character(Scenario),
    Scenario = ifelse(is.na(Scenario) | Scenario == "NA", "Baseline", Scenario),
    Scenario = dplyr::recode(Scenario,
                             "RCP4.5" = "SSP2-4.5",
                             "RCP8.5" = "SSP5-8.5"),
    Scenario = factor(Scenario, levels = c("Baseline","SSP2-4.5","SSP5-8.5")),
    Period   = factor(Period,   levels = c("Baseline","2048–2057","2068–2077"))
  )


# Heat-only table ready for plotting
af_heat <- af_df %>%
  transmute(
    Province, Period, Scenario,
    AF = AF_Heat, LCI = LCI_Heat, UCI = UCI_Heat
  ) %>%
  # Build a tidy y-axis label that orders exactly as the paper:
  mutate(
    Period_Scenario = case_when(
      Period == "Baseline" ~ "Baseline",
      Period == "2048–2057" & Scenario == "SSP2-4.5" ~ "2048–2057 · SSP2-4.5",
      Period == "2048–2057" & Scenario == "SSP5-8.5" ~ "2048–2057 · SSP5-8.5",
      Period == "2068–2077" & Scenario == "SSP2-4.5" ~ "2068–2077 · SSP2-4.5",
      Period == "2068–2077" & Scenario == "SSP5-8.5" ~ "2068–2077 · SSP5-8.5",
      TRUE ~ NA_character_
    ),
    Period_Scenario = factor(
      Period_Scenario,
      levels = c("2068–2077 · SSP5-8.5","2068–2077 · SSP2-4.5",
                 "2048–2057 · SSP5-8.5","2048–2057 · SSP2-4.5",
                 "Baseline")
    ),
    Label = sprintf("%.2f (%.2f, %.2f)", AF, LCI, UCI)
  )

# --- 2) Single-province plot helper (auto x-limits; labels to the right)
plot_af_province <- function(prov, df = af_heat, which = c("Heat")){
  d <- df %>% filter(Province == prov)
  if (nrow(d) == 0) stop("No rows for province: ", prov)
  
  rng <- range(c(d$LCI, d$UCI, 0), na.rm = TRUE)
  pad <- diff(rng) * 0.10
  xr  <- c(rng[1] - pad, rng[2] + pad)
  
  ggplot(d, aes(x = AF, y = Period_Scenario, color = Scenario)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.20, size = 0.5) +
    geom_point(size = 2) +
    geom_text(aes(x = xr[2] - pad*0.25, label = Label),
              hjust = 1, size = 3, color = "black") +
    scale_color_manual(values = c("Baseline" = "black",
                                  "SSP2-4.5" = "#2c7bb6",
                                  "SSP5-8.5" = "#d7191c")) +
    coord_cartesian(xlim = xr) +
    labs(
      title = paste0(prov, " — Heat-related AF (%), baseline-anchored"),
      x = "AF (%)", y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank())
}

# --- 3) Make the four-province panel (as before)
p_baluchistan <- plot_af_province("Baluchistan")
p_kpk         <- plot_af_province("KPK")
p_punjab      <- plot_af_province("Punjab")
p_sindh       <- plot_af_province("Sindh")

combined_plot <- (p_baluchistan / p_kpk / p_punjab / p_sindh) + plot_layout(heights = rep(1,4))
ggsave("Figure2_proj_AF_plot_updated.png", combined_plot, width = 10, height = 12, dpi = 300)
print(combined_plot)

# --- (Optional) Faceted overview of all provinces at once
p_all <- ggplot(af_heat, aes(x = AF, y = Period_Scenario, color = Scenario)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.20, size = 0.4) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c("Baseline" = "black",
                                "SSP2-4.5" = "#2c7bb6",
                                "SSP5-8.5" = "#d7191c")) +
  labs(title = "Heat-related AF (%), baseline-anchored",
       x = "AF (%)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank()) +
  facet_wrap(~ Province, ncol = 2, scales = "free_x")

# ggsave("Figure2_proj_AF_all_provinces.png", p_all, width = 12, height = 10, dpi = 300)
# print(p_all)


## ===============================
## Stage 9: Heat Vulnerability Index (HVI)
## ===============================

############################################################
## Heat Vulnerability Index (HVI)
## - Uses province-level RR at 99th percentile as susceptibility
############################################################

## 1) Load inputs

# District-level exposures & sociodemographics
path_dist_births <- file.path(dir_data_raw, "dist_births_mor_mpi.csv")
stopifnot(file.exists(path_dist_births))

dist_births <- data.table::fread(path_dist_births)

# Province-level RR @ 99th percentile (from your main DLNM stage)
# If the object doesn't already exist in memory, load the CSV we saved earlier.
if (!exists("rr_99th_percentile")) {
  if (file.exists("Province_RR_99th_Percentile_Celsius.csv")) {
    rr_99th_percentile <- read.csv("Province_RR_99th_Percentile_Celsius.csv")
  } else stop("rr_99th_percentile not in memory and 'Province_RR_99th_Percentile_Celsius.csv' not found.")
}


# Keep only the columns we need from district data
columns_to_keep <- c(
  "ADM2_EN", "ADM1_EN.x",               # district & province names
  "avg_pm", "avg_tmean",
  "avg_prec",
  "SUM",                                # births (for later summaries; not used in index)
  "IHME_LMICS", "mpi"                   # sociodemographic/susceptibility
)
dist_births <- dist_births %>%
  dplyr::select(all_of(columns_to_keep))

## 2) Standardise province names on both sides before join

# District data -> add standardised province
dist_births_std <- dist_births %>%
  mutate(prov_std = canonise_prov(ADM1_EN.x))

RR_std <- rr_99th_percentile %>%
  mutate(Province_up = toupper(Province),
         prov_std    = canonise_prov(Province_up)) %>%
  dplyr::select(prov_std, RR, LCI, UCI)

merged_data <- inner_join(dist_births_std, RR_std, by = "prov_std")

if (nrow(merged_data) == 0) {
  stop("Join produced 0 rows. Check province names in dist_births and rr_99th_percentile.")
}

## 3) Construct HVI (no double-counting; propagate RR uncertainty)

# Choose HVI components:
# - Use tmean (avg_tmean) as the heat hazard
# - Add PM (avg_pm) and two socio metrics (IHME_LMICS, mpi)
hvi_vars <- c("avg_tmean", "avg_pm", "IHME_LMICS", "mpi")
missing_hvi <- setdiff(hvi_vars, names(merged_data))
if (length(missing_hvi)) stop("Missing HVI columns: ", paste(missing_hvi, collapse = ", "))

merged_data <- merged_data %>%
  mutate(across(all_of(hvi_vars), ~ as.numeric(scale(.x)), .names = "z_{col}"))

# Log-normal RR uncertainty (unchanged)
set.seed(123)
num_samples <- 1000
mu_log <- log(merged_data$RR)
sd_log <- (log(merged_data$UCI) - log(merged_data$LCI)) / (2 * qnorm(0.975))
sd_log[!is.finite(sd_log) | sd_log <= 0] <- 1e-6
rr_draws <- replicate(num_samples, exp(rnorm(nrow(merged_data), mu_log, sd_log)))

# Equal weights
w <- c(
  z_avg_tmean  = 0.25,
  z_avg_pm     = 0.25,
  z_IHME_LMICS = 0.25,
  z_mpi        = 0.25
)

# Build base_score
merged_data <- as.data.frame(merged_data)
for (nm in names(w)) {
  if (!nm %in% names(merged_data)) stop("Column missing for HVI: ", nm)
  if (!is.numeric(merged_data[[nm]])) merged_data[[nm]] <- as.numeric(merged_data[[nm]])
}
Z <- as.matrix(merged_data[, names(w), drop = FALSE]); storage.mode(Z) <- "double"
w_vec <- matrix(as.numeric(w), ncol = 1); rownames(w_vec) <- names(w)
Z <- Z[, rownames(w_vec), drop = FALSE]
base_score <- as.numeric(Z %*% w_vec)

# Combine with RR draws
HVI_sims <- sweep(rr_draws, 1, base_score, `*`)
merged_data$mean_HVI <- rowMeans(HVI_sims)
merged_data$HVI_LCI  <- apply(HVI_sims, 1, quantile, 0.025)
merged_data$HVI_UCI  <- apply(HVI_sims, 1, quantile, 0.975)
merged_data$HVI_quintile <- dplyr::ntile(merged_data$mean_HVI, 5)

hvi_data <- merged_data %>%
  dplyr::select(ADM2_EN, mean_HVI, HVI_LCI, HVI_UCI, HVI_quintile)

write.csv(hvi_data, "HVI_by_district.csv", row.names = FALSE)
message("Saved: HVI_by_district.csv")


## 4) Map: merge with shapefiles and plot

# District shapefile (ADM2)
adm2_path <- file.path(dir_data_raw, "shapefiles", "ADM2", "pak_admbnda_adm2_wfp_20220909.shp")
districts_shp <- sf::st_read(adm2_path, quiet = TRUE)

# Merge HVI with districts
district_data <- merge(districts_shp, hvi_data, by.x = "ADM2_EN", by.y = "ADM2_EN", all.x = TRUE)

# Also merge back the original dist_births (for panel maps)
final_district_data <- merge(district_data, dist_births, by = "ADM2_EN", all.x = TRUE)

# Drop non-essential columns for clarity (optional)
cols_drop <- c("Shape_Leng", "Shape_Area", "ADM2_PCODE", "ADM2_REF",
               "ADM2ALT1EN", "ADM2ALT2EN", "ADM1_EN", "ADM1_PCODE",
               "ADM0_EN", "ADM0_PCODE", "date", "validOn", "validTo")
cols_drop <- intersect(cols_drop, names(final_district_data))
final_district_data <- final_district_data %>% dplyr::select(-all_of(cols_drop))

# Province boundaries (ADM1) for overlay and labels
adm1_path <- file.path(dir_data_raw, "shapefiles", "ADM1", "pak_admbnda_adm1_wfp_20220909.shp")
province_boundaries <- sf::st_read(adm1_path, quiet = TRUE)


# Province label positions (centroids)
province_centroids <- province_boundaries %>%
  mutate(centroid = st_centroid(geometry)) %>%
  st_as_sf() %>%
  mutate(
    lon = st_coordinates(centroid)[,1],
    lat = st_coordinates(centroid)[,2]
  )

# Main HVI map
hvi_plot <- ggplot(data = final_district_data) +
  geom_sf(aes(fill = mean_HVI), color = NA) +
  geom_sf(data = province_boundaries, fill = NA, color = "black", size = 0.7) +
  geom_text(
    data = province_centroids %>% filter(!ADM1_EN %in% c("Azad Kashmir", "Islamabad")),
    aes(
      x = ifelse(ADM1_EN == "Khyber Pakhtunkhwa", lon - 1.0, lon),
      y = ifelse(ADM1_EN == "Khyber Pakhtunkhwa", lat - 0.5, lat),
      label = ADM1_EN
    ),
    size = 4.5, fontface = "bold", color = "white"
  ) +
  scale_fill_viridis(option = "turbo", na.value = "gray80", name = "Heat Vulnerability Index") +
  theme_light() +
  theme(
    panel.background = element_rect(fill = "gray40", color = NA),
    legend.position  = c(0.82, 0.18),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.height = unit(1.0, "cm"),
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.grid.minor = element_blank()
  )

ggsave("HVI_Map_Final.png", hvi_plot, width = 12, height = 8, dpi = 300)
print(hvi_plot)

## 5) Optional: panel maps for environmental context (unchanged)

plot_pm <- ggplot(final_district_data) +
  geom_sf(aes(fill = avg_pm)) +
  scale_fill_viridis(option = "plasma", na.value = "gray80") +
  labs(title = NULL, fill = "PM2.5") +
  theme_minimal() + theme(legend.position = "right", legend.key.height = unit(0.8, "cm"))

plot_precip <- ggplot(final_district_data) +
  geom_sf(aes(fill = avg_prec)) +
  scale_fill_viridis(option = "viridis", na.value = "gray80") +
  labs(title = NULL, fill = "Precipitation") +
  theme_minimal() + theme(legend.position = "right", legend.key.height = unit(0.8, "cm"))

plot_births <- ggplot(final_district_data) +
  geom_sf(aes(fill = SUM)) +
  scale_fill_viridis(option = "mako", na.value = "gray80") +
  labs(title = NULL, fill = "Births") +
  theme_minimal() + theme(legend.position = "right", legend.key.height = unit(0.8, "cm"))

plot_mpi <- ggplot(final_district_data) +
  geom_sf(aes(fill = mpi)) +
  scale_fill_viridis(option = "inferno", na.value = "gray80") +
  labs(title = NULL, fill = "MPI") +
  theme_minimal() + theme(legend.position = "right", legend.key.height = unit(0.8, "cm"))

plot_deaths <- ggplot(final_district_data) +
  geom_sf(aes(fill = IHME_LMICS)) +
  scale_fill_viridis(option = "cividis", na.value = "gray80") +
  labs(title = NULL, fill = "Deaths") +
  theme_minimal() + theme(legend.position = "right", legend.key.height = unit(0.8, "cm"))


plot_tmean <- ggplot(final_district_data) +
  geom_sf(aes(fill = avg_tmean)) +
  scale_fill_viridis(option = "turbo", na.value = "gray80") +
  labs(title = NULL, fill = "Mean Temperature (°C)") +
  theme_minimal() + theme(legend.position = "right", legend.key.height = unit(0.8, "cm"))


bottom_panel <- (plot_pm + plot_precip + plot_births) /
  (plot_tmean + plot_mpi + plot_deaths) +
  plot_annotation(tag_levels = "a")

ggsave("Paneled_Environmental_Maps.png", bottom_panel, width = 18, height = 12, dpi = 300)

#projected HVI
## ------------------------------------------------------------
## Projected HVI using province-level RR from projections
## Paste this AFTER your AF/RR pipeline and BEFORE HVI mapping.
## Requires:
##   - cb_for_pred_poly, weighted_coefs, weighted_vcov
##   - prov_centers (province, cen_scaled), mean_s, sd_s
##   - projection frames: prov_RCP4.5_48, prov_RCP8.5_48, prov_RCP4.5_68, prov_RCP8.5_68
##   - district frame: dist_births with ADM1_EN.x + HVI components (avg_tmean, avg_pm, IHME_LMICS, mpi, SUM)
## ------------------------------------------------------------

## --- (A) Helpers -------------------------------------------------------------

# Standardize a projection frame to (province, date, tmean_C)
prep_proj_t <- function(df) {
  df <- as.data.frame(df); nm <- names(df)
  if (!"province" %in% nm) {
    if ("ADM1_EN" %in% nm) names(df)[nm == "ADM1_EN"] <- "province" else
      if ("ADM1_EN.x" %in% nm) names(df)[nm == "ADM1_EN.x"] <- "province" else
        stop("prep_proj_t(): need province/ADM1_EN/ADM1_EN.x")
  }
  if (!"date" %in% nm) {
    if ("date.y" %in% nm) names(df)[nm == "date.y"] <- "date" else {
      cand <- grep("^date", nm, value = TRUE)
      if (length(cand) >= 1) names(df)[nm == cand[1]] <- "date" else
        stop("prep_proj_t(): need a date column")
    }
  }
  if (!inherits(df$date, "Date")) df$date <- as.Date(df$date)
  if (!"tmean_C" %in% nm) {
    if ("mean_C" %in% nm) names(df)[nm == "mean_C"] <- "tmean_C" else
      stop("prep_proj_t(): need tmean_C or mean_C")
  }
  df[, c("province","date","tmean_C"), drop = FALSE]
}

# Safe crosspred at one 'at' point
safe_one_crosspred <- function(cb, coefs, vcov, at_scaled, cen_scaled) {
  stopifnot(length(at_scaled) == 1, length(cen_scaled) == 1)
  cp <- try(crosspred(cb, coef = coefs, vcov = vcov, at = at_scaled, cen = cen_scaled), silent = TRUE)
  if (inherits(cp, "try-error")) return(NULL)
  list(mu_log = as.numeric(cp$allfit), se_log = as.numeric(cp$allse))
}

# Compute province RR at projected P99 temperature for a scenario/period
proj_rr_at_p99 <- function(df_proj, scenario_lab, period_lab) {
  df_proj %>%
    group_by(province) %>%
    summarise(t99 = quantile(tmean_C, 0.99, na.rm = TRUE), .groups = "drop") %>%
    rowwise() %>%
    mutate(
      at_scaled = (t99 - mean_s) / sd_s,
      cen_scaled = {
        cs <- prov_centers$cen_scaled[prov_centers$province == province]
        if (length(cs) != 1 || !is.finite(cs)) median(LBW_clean$scaled_tmean, na.rm = TRUE) else cs
      },
      pred = list(safe_one_crosspred(cb_for_pred_poly, weighted_coefs, weighted_vcov, at_scaled, cen_scaled))
    ) %>%
    ungroup() %>%
    mutate(
      mu_log = map_dbl(pred, ~ if (is.null(.x)) NA_real_ else .x$mu_log),
      se_log = map_dbl(pred, ~ if (is.null(.x)) NA_real_ else .x$se_log),
      RR     = exp(mu_log),
      LCI    = exp(mu_log - 1.96 * se_log),
      UCI    = exp(mu_log + 1.96 * se_log),
      Scenario = scenario_lab,
      Period   = period_lab
    ) %>%
    dplyr::select(province, Scenario, Period, t99, mu_log, se_log, RR, LCI, UCI)
}

# Lognormal RR draws
rr_draws_by_prov <- function(mu_log, se_log, n = 1000) {
  se_log[!is.finite(se_log) | se_log <= 0] <- 1e-6
  matrix(exp(rnorm(n = n, mean = mu_log, sd = se_log)), nrow = 1)
}

label_period <- function(nm) if (grepl("48", nm)) "2048–2057" else "2068–2077"
label_scn    <- function(nm) sub("_.*", "", nm)              # "RCP4.5" / "RCP8.5"
label_ssp    <- function(x) dplyr::recode(x, "RCP4.5"="SSP2-4.5", "RCP8.5"="SSP5-8.5", .default = x)

## --- (B) Inputs needed here --------------------------------------------------
# Expect these four objects (lists/data.frames of province-date climate projections):
#   prov_RCP4.5_48, prov_RCP8.5_48, prov_RCP4.5_68, prov_RCP8.5_68
# Each must have province + date + tmean_C
stopifnot(exists("prov_RCP4.5_48"), exists("prov_RCP8.5_48"),
          exists("prov_RCP4.5_68"), exists("prov_RCP8.5_68"))

proj_raw <- list(
  RCP4.5_48 = prov_RCP4.5_48,
  RCP8.5_48 = prov_RCP8.5_48,
  RCP4.5_68 = prov_RCP4.5_68,
  RCP8.5_68 = prov_RCP8.5_68
)

# Ensure district data present
stopifnot(exists("dist_births"))
stopifnot(all(c("ADM1_EN.x","ADM2_EN","avg_tmean","avg_pm","IHME_LMICS","mpi","SUM") %in% names(dist_births)))

## --- (C) Province-level projected RR tables ---------------------------------

proj_list_t <- lapply(proj_raw, prep_proj_t)

proj_t_tbl <- bind_rows(lapply(names(proj_list_t), function(nm) {
  dfp <- proj_list_t[[nm]] %>% dplyr::filter(province %in% unique(prov_centers$province))
  dfp %>%
    group_by(province) %>%
    summarise(T_mean_proj = mean(tmean_C, na.rm = TRUE), .groups = "drop") %>%
    mutate(Scenario = label_ssp(label_scn(nm)), Period = label_period(nm))
}))

proj_t_tbl_std <- proj_t_tbl %>%
  mutate(prov_std = canonise_prov(province)) %>%
  dplyr::select(prov_std, Scenario, Period, T_mean_proj)
# Compute RR tables
proj_rr_tbl <- bind_rows(lapply(names(proj_list_t), function(nm) {
  dfp <- proj_list_t[[nm]] %>% filter(province %in% unique(prov_centers$province))
  proj_rr_at_p99(dfp, scenario_lab = label_ssp(label_scn(nm)), period_lab = label_period(nm))
}))


# Canonicalise province names for joins
proj_rr_tbl_std <- proj_rr_tbl %>%
  mutate(prov_std = canonise_prov(province))

# proj_hi_tbl_std <- proj_hi_tbl %>%
#   mutate(prov_std = canonise_prov(province)) %>%
#   dplyr::select(prov_std, Scenario, Period, HI_mean_proj)
stopifnot(
  exists("cb_for_pred_poly"),
  exists("weighted_coefs"),
  exists("weighted_vcov"),
  exists("prov_centers"),
  exists("mean_s"), is.finite(mean_s),
  exists("sd_s"),   is.finite(sd_s), sd_s > 0
)

stopifnot(exists("LBW_clean"), "scaled_tmean" %in% names(LBW_clean))

## --- (E) Baseline tmean reference (for z-scoring projected tmean) ------------------

baseline_T_mean <- mean(dist_births$avg_tmean, na.rm = TRUE)
baseline_T_sd   <- sd(dist_births$avg_tmean,   na.rm = TRUE)
if (!is.finite(baseline_T_sd) || baseline_T_sd == 0) {
  stop("Baseline T SD is zero/NA. Check dist_births$avg_tmean.")
}

# Static z for non-climate components
dist_scaled_static <- dist_births %>%
  mutate(
    z_avg_pm     = as.numeric(scale(avg_pm)),
    z_IHME_LMICS = as.numeric(scale(IHME_LMICS)),
    z_mpi        = as.numeric(scale(mpi))
  ) %>%
  dplyr::select(ADM2_EN, z_avg_pm, z_IHME_LMICS, z_mpi)

# # Baseline distribution from district data
# baseline_HI_mean  <- mean(dist_births$avg_heat_index, na.rm = TRUE)
# baseline_HI_sd    <- sd(dist_births$avg_heat_index,   na.rm = TRUE)
# if (!is.finite(baseline_HI_sd) || baseline_HI_sd == 0) {
#   stop("Baseline HI SD is zero/NA. Check dist_births$avg_heat_index.")
# }

# Districts with province key
dist_births_std <- dist_births %>%
  mutate(prov_std = canonise_prov(ADM1_EN.x)) %>%
  dplyr::select(ADM2_EN, prov_std, avg_pm, IHME_LMICS, mpi, SUM)

## --- (F) Build projected district HVI per Scenario × Period ------------------

# Equal weights (same as baseline)
w_pm   <- 0.25; w_IHME <- 0.25; w_mpi <- 0.25; w_T <- 0.25


set.seed(123); nsim <- 1000

HVI_proj_list <- proj_rr_tbl_std %>%
  split(list(.$Scenario, .$Period), drop = TRUE) %>%
  purrr::imap(function(df_rr, key) {
    scn <- unique(df_rr$Scenario); per <- unique(df_rr$Period)
    
    df_t <- proj_t_tbl_std %>% dplyr::filter(Scenario == scn, Period == per)
    df_rrt <- df_rr %>% inner_join(df_t, by = c("prov_std","Scenario","Period"))
    
    rr_sims <- lapply(seq_len(nrow(df_rrt)), function(i) {
      rr_draws_by_prov(df_rrt$mu_log[i], df_rrt$se_log[i], n = nsim)
    })
    names(rr_sims) <- df_rrt$prov_std
    
    dist_base <- dist_births %>%
      mutate(prov_std = canonise_prov(ADM1_EN.x)) %>%
      dplyr::select(ADM2_EN, prov_std, avg_pm, IHME_LMICS, mpi, SUM) %>%
      left_join(dist_scaled_static, by = "ADM2_EN") %>%
      left_join(df_t, by = "prov_std") %>%
      mutate(z_T_proj = (T_mean_proj - baseline_T_mean) / baseline_T_sd)
    
    #  w_T <- 0.25; w_pm <- 0.25; w_IHME <- 0.25; w_mpi <- 0.25
    
    dist_out <- dist_base %>%
      rowwise() %>%
      mutate(
        rr_vec = list(as.numeric(rr_sims[[prov_std]])),
        base_score_proj = ifelse(is.na(z_T_proj), NA_real_,
                                 w_T * z_T_proj + w_pm * z_avg_pm + w_IHME * z_IHME_LMICS + w_mpi * z_mpi),
        HVI_mean = ifelse(length(rr_vec) == 0 || is.na(base_score_proj), NA_real_,
                          base_score_proj * mean(rr_vec, na.rm = TRUE)),
        HVI_LCI  = ifelse(length(rr_vec) == 0 || is.na(base_score_proj), NA_real_,
                          base_score_proj * quantile(rr_vec, 0.025, na.rm = TRUE)),
        HVI_UCI  = ifelse(length(rr_vec) == 0 || is.na(base_score_proj), NA_real_,
                          base_score_proj * quantile(rr_vec, 0.975, na.rm = TRUE))
      ) %>%
      ungroup() %>%
      mutate(Scenario = scn, Period = per) %>%
      dplyr::select(ADM2_EN, prov_std, Scenario, Period, base_score_proj, HVI_mean, HVI_LCI, HVI_UCI, SUM)
    
    dist_out
  })

HVI_proj_districts <- bind_rows(HVI_proj_list)
readr::write_csv(HVI_proj_districts, "Projected_HVI_by_District.csv")

HVI_proj_province <- HVI_proj_districts %>%
  group_by(prov_std, Scenario, Period) %>%
  summarise(
    HVI_mean_pw = weighted.mean(HVI_mean, w = SUM, na.rm = TRUE),
    HVI_LCI_pw  = weighted.mean(HVI_LCI,  w = SUM, na.rm = TRUE),
    HVI_UCI_pw  = weighted.mean(HVI_UCI,  w = SUM, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(HVI_proj_province, "Projected_HVI_by_Province.csv")

cat("✓ Projected HVI tables written:\n  - Projected_HVI_by_District.csv\n  - Projected_HVI_by_Province.csv\n")


proj_hvi <- readr::read_csv(file.path(dir_data_proc, "Projected_HVI_by_District.csv"))

districts_shp <- st_read(adm2_path, quiet = TRUE)  # same as earlier
proj_map <- merge(districts_shp, proj_hvi, by.x = "ADM2_EN", by.y = "ADM2_EN", all.x = TRUE)

ggplot(proj_map) +
  geom_sf(aes(fill = HVI_mean), color = NA) +
  scale_fill_viridis(option = "turbo", na.value = "gray85") +
  facet_grid(Scenario ~ Period) +
  labs(fill = "Projected HVI") +
  theme_minimal()


# --- read baseline & projected ---
base_hvi <- readr::read_csv(
  file.path(dir_data_proc, "HVI_by_district.csv"),
  show_col_types = FALSE
) %>%
  dplyr::select(ADM2_EN, base_HVI = mean_HVI)

proj_hvi <- readr::read_csv(
  file.path(dir_data_proc, "Projected_HVI_by_District.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Scenario = forcats::fct_relevel(Scenario, "SSP2-4.5", "SSP5-8.5"),
    Period   = forcats::fct_relevel(Period, "2048–2057", "2068–2077")
  )

# Join baseline
proj_vs_base <- proj_hvi %>%
  left_join(base_hvi, by = "ADM2_EN") %>%
  mutate(
    dHVI      = HVI_mean - base_HVI,                  # absolute change
    pct_change = 100 * (HVI_mean - base_HVI) / abs(base_HVI)
  )

# Shapes
adm2 <- st_read(adm2_path, quiet = TRUE)
adm1 <- st_read(adm1_path, quiet = TRUE)

# Merge to shapes
map_delta <- adm2 %>% left_join(proj_vs_base, by = c("ADM2_EN" = "ADM2_EN"))

# (i) Absolute change map (diverging around 0)
lims_d <- max(abs(range(map_delta$dHVI, na.rm = TRUE)))
# Reusable white theme
theme_white_facets <- theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key       = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  )

# (i) Absolute change map (Δ HVI) with white background
p_dHVI <- ggplot(map_delta) +
  geom_sf(aes(fill = dHVI), color = NA) +
  geom_sf(data = adm1, fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_viridis(option = "plasma", na.value = "grey85",
                     limits = c(-lims_d, lims_d), name = "Δ HVI") +
  facet_grid(Scenario ~ Period) +
  theme_white_facets

ggsave(
  filename = file.path(dir_figures, "Projected_vs_Baseline_DeltaHVI.png"),
  plot     = p_dHVI,
  width    = 12,
  height   = 8,
  dpi      = 300
)


# (ii) Percent change map
lims_p <- max(abs(range(map_delta$pct_change, na.rm = TRUE)))
# (ii) Percent change map with white background
p_pct <- ggplot(map_delta) +
  geom_sf(aes(fill = pct_change), color = NA) +
  geom_sf(data = adm1, fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_viridis(option = "magma", na.value = "grey85",
                     limits = c(-lims_p, lims_p), name = "% change") +
  facet_grid(Scenario ~ Period) +
  theme_white_facets

ggsave(
  filename = file.path(dir_figures, "Projected_vs_Baseline_PercentChange.png"),
  plot     = p_pct,
  width    = 12,
  height   = 8,
  dpi      = 300
)
p_dHVI; p_pct



#############################
# Stage 10: New Subgroup analysis via single interaction model + linear combos
#############################

# --- Settings (match main model) ---
lag <- 7
degree_exposure <- 2
degree_lag <- 1

# Ensure scaled_tmean exists (as in your main pipeline)
stopifnot("scaled_tmean" %in% names(LBW_clean))

# Define subgroup flags used in Table 2
LBW_clean <- LBW_clean %>%
  mutate(
    PM25_Group = factor(ifelse(PM25 < 25, "Fair", "Poor/Hazardous"),
                        levels = c("Fair","Poor/Hazardous")),
    Education_Group = factor(case_when(
      Education %in% c("No education","Primary") ~ "Low Education",
      Education %in% c("Secondary","Higher")     ~ "Secondary+",
      TRUE ~ NA_character_), levels = c("Low Education","Secondary+")),
    Wealth_Group = factor(case_when(
      Wealth_Index %in% c("Poorest","Poorer","Middle") ~ "Low Wealth",
      Wealth_Index %in% c("Richer","Richest")          ~ "High Wealth",
      TRUE ~ NA_character_), levels = c("Low Wealth","High Wealth")),
    Area_Group = factor(Area, levels = c("Urban","Rural"))
  )

# Helper to fit an interaction model for one domain and return RR at 90/99 + p-het
fit_interaction_domain <- function(data, domain_var, labels_map,
                                   cb_name = "cb_all",
                                   pct_vec = c(`90%`=0.90, `99%`=0.99)) {
  dat <- data %>% dplyr::filter(!is.na(.data[[domain_var]])) %>% as.data.frame()
  
  # housekeeping
  if (!"PM25_scaled" %in% names(dat)) dat$PM25_scaled <- as.numeric(scale(dat$PM25))
  dat$Education_Group <- droplevels(factor(dat$Education_Group))
  dat$province        <- droplevels(factor(dat$province))
  
  # Use treatment contrasts (avoid polynomial contrasts for ordered factors)
  old_contr <- options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment"))
  on.exit(options(old_contr), add = TRUE)
  
  # ---- crossbasis (unchanged) ----
  cb <- crossbasis(dat$scaled_tmean, lag = lag,
                   argvar = list(fun="poly", degree=degree_exposure),
                   arglag = list(fun="poly", degree=degree_lag),
                   group = dat$province)
  dat[[cb_name]] <- cb
  cb_cols <- colnames(cb)
  main_expected <- paste0(cb_name, cb_cols)
  
  # two levels for the domain
  levs <- levels(dat[[domain_var]]); stopifnot(length(levs) == 2)
  base_lev <- levs[1]; alt_lev <- levs[2]
  
  # percentiles & center (define BEFORE using them)
  get_group_pcts <- function(g) {
    x <- dat %>% dplyr::filter(.data[[domain_var]] == g) %>% dplyr::pull(scaled_tmean)
    stats::quantile(x, probs = unname(pct_vec), na.rm = TRUE)
  }
  at_base <- get_group_pcts(base_lev)
  at_alt  <- get_group_pcts(alt_lev)
  cen_all <- stats::median(dat$scaled_tmean, na.rm = TRUE)
  
  # ----- build covariate string: include Education_Group unless it IS the domain -----
  covars <- c("poly(PM25_scaled, 2)", "poly(Year, 3)")
  if (domain_var != "Education_Group" && "Education_Group" %in% names(dat)) {
    # only add if both levels present; otherwise skip to avoid singularities
    if (nlevels(dat$Education_Group) >= 2) covars <- c(covars, "Education_Group")
  }
  covar_str <- paste(covars, collapse = " + ")
  
  # ---- model (no ordered(Education)) ----
  form <- as.formula(paste0(
    "LBW_count ~ ", cb_name, "*", domain_var,
    " + ", covar_str,
    " + offset(log(WRA_pop)) + (1 | province)"
  ))
  fit <- glmmTMB(form, data = dat, family = nbinom2(), na.action = na.exclude)
  
  beta <- fixef(fit)$cond; V <- vcov(fit)$cond; bn <- names(beta)
  
  # prediction-row helper
  pred_row_cb <- function(cb, at, cen) {
    q <- ncol(cb); Iq <- diag(q)
    vapply(seq_len(q), function(j) {
      cp <- dlnm::crosspred(cb, coef = Iq[j,], vcov = diag(q), at = at, cen = cen)
      as.numeric(cp$allfit)
    }, numeric(1))
  }
  
  # align main cb terms
  main_idx <- match(main_expected, bn)
  if (any(is.na(main_idx))) stop("Missing CB coef: ", main_expected[which(is.na(main_idx))[1]])
  beta_base <- beta[main_idx]
  V_base    <- V[main_idx, main_idx, drop = FALSE]
  
  # align interaction terms  cb_col:DomainVarAltLevel
  find_int_for_main <- function(main_nm) {
    hits <- which(startsWith(bn, paste0(main_nm, ":")) &
                    grepl(paste0(":", domain_var, alt_lev), bn))
    if (length(hits) == 1) hits else NA_integer_
  }
  int_idx_aligned <- vapply(main_expected, find_int_for_main, integer(1))
  beta_alt <- beta_base
  add_ok   <- !is.na(int_idx_aligned)
  beta_alt[add_ok] <- beta_alt[add_ok] + beta[int_idx_aligned[add_ok]]
  
  V_alt <- V_base
  if (any(add_ok)) {
    mm <- main_idx[add_ok]; ii <- int_idx_aligned[add_ok]
    V_alt[add_ok, add_ok] <-
      V[mm, mm, drop = FALSE] + V[ii, ii, drop = FALSE] +
      V[mm, ii, drop = FALSE] + V[ii, mm, drop = FALSE]
  }
  
  # predictions
  pr_base <- crosspred(cb, coef = beta_base, vcov = V_base, at = unname(at_base), cen = cen_all)
  pr_alt  <- crosspred(cb, coef = beta_alt,  vcov = V_alt,  at = unname(at_alt),  cen = cen_all)
  
  tab_base <- tibble::tibble(
    Subgroup   = paste0(domain_var, "=", base_lev),
    Percentile = names(pct_vec),
    logRR      = as.numeric(pr_base$allfit),
    SE         = as.numeric(pr_base$allse),
    RR         = exp(logRR),
    LCI        = exp(logRR - 1.96*SE),
    UCI        = exp(logRR + 1.96*SE)
  )
  tab_alt <- tibble::tibble(
    Subgroup   = paste0(domain_var, "=", alt_lev),
    Percentile = names(pct_vec),
    logRR      = as.numeric(pr_alt$allfit),
    SE         = as.numeric(pr_alt$allse),
    RR         = exp(logRR),
    LCI        = exp(logRR - 1.96*SE),
    UCI        = exp(logRR + 1.96*SE)
  )
  preds <- dplyr::bind_rows(tab_base, tab_alt)
  
  # covariance-aware heterogeneity (GLHT)
  make_glht_p <- function(label) {
    Lb <- pred_row_cb(cb, unname(at_base[label]), cen_all)
    La <- pred_row_cb(cb, unname(at_alt[label]),  cen_all)
    
    K <- matrix(0, nrow = 1, ncol = length(bn)); colnames(K) <- bn
    K[1, main_idx] <- (La - Lb)
    ok <- !is.na(int_idx_aligned)
    K[1, int_idx_aligned[ok]] <- K[1, int_idx_aligned[ok]] + La[ok]
    
    ht <- multcomp::glht(fit, linfct = K)
    s  <- summary(ht)   # <-- FIX: use the generic summary(), not multcomp::summary
    
    tibble::tibble(
      Percentile = label,
      logRR_diff = drop(K %*% beta),
      SE         = s$test$sigma,
      z          = s$test$tstat,
      p_heterogeneity = s$test$pvalues
    )
  }
  
  het_covaware <- dplyr::bind_rows(make_glht_p("90%"), make_glht_p("99%")) %>%
    mutate(Domain = domain_var,
           Contrast = paste0(base_lev, " vs ", alt_lev))
  
  list(preds = preds, het = het_covaware)
}

# Run per-domain interaction models
res_pm25   <- fit_interaction_domain(LBW_clean, "PM25_Group",
                                     labels_map = c("Fair"="fair (< 25 μg m-³)",
                                                    "Poor/Hazardous"="poor–hazardous (≥ 25 μg m-³)"))
res_edu    <- fit_interaction_domain(LBW_clean, "Education_Group",
                                     labels_map = c("Low Education"="primary or no education",
                                                    "Secondary+"="secondary or above"))
res_wealth <- fit_interaction_domain(LBW_clean, "Wealth_Group",
                                     labels_map = c("Low Wealth"="low","High Wealth"="high"))
res_area   <- fit_interaction_domain(LBW_clean, "Area_Group",
                                     labels_map = c("Urban"="urban","Rural"="rural"))

# Combine outputs
pred_by_sub <- dplyr::bind_rows(res_pm25$preds, res_edu$preds, res_wealth$preds, res_area$preds)
het_tbl     <- dplyr::bind_rows(res_pm25$het,   res_edu$het,   res_wealth$het,   res_area$het)



# --- pretty labels as a lookup table ---
label_pretty <- c(
  "PM25_Group=Fair"               = "fair (< 25 μg m-³)",
  "PM25_Group=Poor/Hazardous"     = "poor–hazardous (≥ 25 μg m-³)",
  "Education_Group=Low Education" = "primary or no education",
  "Education_Group=Secondary+"    = "secondary or above",
  "Wealth_Group=Low Wealth"       = "low",
  "Wealth_Group=High Wealth"      = "high",
  "Area_Group=Urban"              = "urban",
  "Area_Group=Rural"              = "rural"
)

lab_df <- enframe(label_pretty, name = "Subgroup", value = "Subgroup_pretty")

# --- RR table in the same shape as your manuscript table ---
rr_table <- pred_by_sub %>%
  filter(Percentile %in% c("90%","99%")) %>%
  left_join(lab_df, by = "Subgroup") %>%
  mutate(
    Subgroup = coalesce(Subgroup_pretty, Subgroup),
    RR_CI    = sprintf("%.2f (%.2f–%.2f)", RR, LCI, UCI)
  ) %>%
  dplyr::select(Subgroup, Percentile, RR_CI) %>%
  pivot_wider(names_from = Percentile, values_from = RR_CI) %>%
  rename(`risk estimates (90th Percentile)` = `90%`,
         `risk estimates (99th Percentile)` = `99%`)

print(rr_table)

# --- Heterogeneity p-values (already built as het_tbl) in a wide view ---
het_table_wide <- het_tbl %>%
  mutate(PercentileCol = ifelse(Percentile == "90%", "p-het (90th)", "p-het (99th)")) %>%
  dplyr::select(Domain, Contrast, PercentileCol, p_heterogeneity) %>%
  pivot_wider(names_from = PercentileCol, values_from = p_heterogeneity)

print(het_table_wide)


#load packages (if needed)
library(dplyr)
library(tidyr)
library(stringr)

# 1) Add Domain parsed from "Something=Level"
pred_by_sub_fix <- pred_by_sub %>%
  mutate(
    Domain  = sub("=.+$", "", Subgroup),      # e.g., "PM25_Group"
    Level   = sub("^.+?=", "", Subgroup)      # e.g., "Fair"
  )

# 2) Pretty labels
pretty_sub <- c(
  "PM25_Group=Fair"               = "Fair (<25 μg/m³)",
  "PM25_Group=Poor/Hazardous"     = "Poor–Hazardous (≥25 μg/m³)",
  "Education_Group=Low Education" = "Low education",
  "Education_Group=Secondary+"    = "Secondary+",
  "Wealth_Group=Low Wealth"       = "Low wealth",
  "Wealth_Group=High Wealth"      = "High wealth",
  "Area_Group=Urban"              = "Urban",
  "Area_Group=Rural"              = "Rural"
)

# ---------- Table 1: compact (one row per domain) ----------
pred_wide <- pred_by_sub_fix %>%
  mutate(
    Pretty = dplyr::coalesce(pretty_sub[Subgroup], Subgroup),
    RR_CI  = sprintf("%.2f (%.2f–%.2f)", RR, LCI, UCI)
  ) %>%
  dplyr::select(Domain, Subgroup = Pretty, Percentile, RR_CI) %>%
  tidyr::pivot_wider(names_from = c(Subgroup, Percentile),
                     values_from = RR_CI) %>%
  relocate(matches("90%$"), .before = matches("99%$"))

het_wide <- het_tbl %>%
  group_by(Domain) %>%
  summarise(
    `p-het (90%)` = p_heterogeneity[Percentile == "90%"][1],
    `p-het (99%)` = p_heterogeneity[Percentile == "99%"][1],
    .groups = "drop"
  )

table_compact <- pred_wide %>%
  left_join(het_wide, by = "Domain") %>%
  mutate(
    Domain = dplyr::recode(
      Domain,
      "PM25_Group"      = "PM₂.₅",
      "Education_Group" = "Education",
      "Wealth_Group"    = "Wealth",
      "Area_Group"      = "Area",
      .default = Domain
    )
  )


table_compact
# write.csv(table_compact, "Table_Subgroups_RR_and_phet_compact.csv", row.names = FALSE)

# ---------- Table 2: long (two rows per domain) ----------
rr_rows <- pred_by_sub_fix %>%
  mutate(
    Pretty = dplyr::coalesce(pretty_sub[Subgroup], Subgroup),
    RR90 = ifelse(Percentile == "90%",
                  sprintf("%.2f (%.2f–%.2f)", RR, LCI, UCI), NA_character_),
    RR99 = ifelse(Percentile == "99%",
                  sprintf("%.2f (%.2f–%.2f)", RR, LCI, UCI), NA_character_)
  ) %>%
  group_by(Domain, Pretty) %>%
  summarise(
    `RR (90%)` = dplyr::first(na.omit(RR90), default = NA_character_),
    `RR (99%)` = dplyr::first(na.omit(RR99), default = NA_character_),
    .groups = "drop"
  )

phet_cols <- het_tbl %>%
  dplyr::select(Domain, Percentile, p_heterogeneity) %>%
  mutate(Percentile = ifelse(Percentile == "90%", "p-het (90%)", "p-het (99%)")) %>%
  tidyr::pivot_wider(names_from = Percentile, values_from = p_heterogeneity)

table_long <- rr_rows %>%
  left_join(phet_cols, by = "Domain") %>%
  mutate(
    Domain = dplyr::recode(
      Domain,
      "PM25_Group"      = "PM₂.₅",
      "Education_Group" = "Education",
      "Wealth_Group"    = "Wealth",
      "Area_Group"      = "Area",
      .default = Domain
    )
  ) %>%
  arrange(factor(Domain, levels = c("PM₂.₅","Wealth","Education","Area")),
          Pretty)

table_long
#write.csv(table_long, "Table_Subgroups_RR_and_phet_long.csv", row.names = FALSE)

save.image()
 
