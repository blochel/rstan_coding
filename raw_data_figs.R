# data figures
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# Standardize water data column name to subregion
water_standardized <- water |>
  rename(subregion = region)

# raw water dynamics for each water year ----------------------------------
water_long <- water_standardized |>
  pivot_longer(
    cols = c(init_depth, breed_season_depth, recession,
             pre_recession, post_recession, dry_days, reversals),
    names_to = "metric",
    values_to = "value"
  )

ggplot(water_long, aes(x = year, y = value, color = subregion, group = subregion)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    x = "Water Year",
    y = "Value",
    color = "Subregion",
    title = "Raw Water Dynamics by Water Year"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# each hydro year initiation time -----------------------------------------
colonies_url <- "https://raw.githubusercontent.com/weecology/EvergladesWadingBird/main/SiteandMethods/colonies.csv"
colonies <- read_csv(colonies_url, show_col_types = FALSE)

bird_subregion_summary <- bird_data |>
  left_join(colonies, by = "colony") |>
  group_by(water_year, subregion) |>
  summarise(y = min(julian_wy, na.rm = TRUE), .groups = "drop") |>
  rename(year = water_year)

bird_average_summary <- bird_data |>
  group_by(water_year) |>
  summarise(y = mean(julian_wy, na.rm = TRUE), .groups = "drop") |>
  rename(year = water_year)

ggplot() +
  geom_smooth(data = bird_average_summary, aes(x = year, y = y)) +
  geom_point(data = bird_average_summary, aes(x = year, y = y),
             colour = "black", size = 4) +
  geom_line(data = bird_subregion_summary,
            aes(x = year, y = y, color = subregion, group = subregion),
            alpha = 0.7) +
  geom_point(data = bird_subregion_summary,
             aes(x = year, y = y, color = subregion, group = subregion),
             size = 1.5) +
  labs(
    x = "Water Year",
    y = "First Nest Initiation Date (Day of Water Year)",
    color = "Subregion",
    title = "Wood Stork Nest Initiation by Subregion Across Water Years\n[black points & smooth average initiation]"
  ) +
  theme_minimal()

# both in one ------------------------------------------------------------
target_subregions <- c("2a", "2b", "3an", "3as", "3ase", "3b")

# Combine observed initiation subregions with explicit target subregions
selected_subregions <- unique(c(
  bird_subregion_summary |> filter(!is.na(y)) |> pull(subregion),
  target_subregions
))

# Filter and format nest initiation by subregion
bird_initiation_long <- bird_subregion_summary |>
  filter(subregion %in% selected_subregions) |>
  rename(value = y) |>
  mutate(metric = "nest_initiation")

# Filter water dynamics by subregion and reshape
water_long_filtered <- water_standardized |>
  filter(subregion %in% selected_subregions) |>
  pivot_longer(
    cols = c(init_depth, breed_season_depth, recession,
             pre_recession, post_recession, dry_days, reversals),
    names_to = "metric",
    values_to = "value"
  )

# Combine datasets
combined_data_filtered <- bind_rows(water_long_filtered, bird_initiation_long)

# Create a consistent color palette across all selected subregions
unique_subregions <- sort(unique(combined_data_filtered$subregion))
subregion_palette <- setNames(
  scales::hue_pal()(length(unique_subregions)),
  unique_subregions
)

# Plot raw water dynamics and nest initiation by subregion
ggplot(combined_data_filtered, aes(x = year, y = value, color = subregion, group = subregion)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_color_manual(values = subregion_palette) +
  labs(
    x = "Water Year",
    y = "Value",
    color = "Subregion",
    title = "Raw Water Dynamics and Nest Initiation by Subregion"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")