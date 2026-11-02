#data figures


library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)


# raw water dynamics for each water year ----------------------------------


water_long <- water |>
  pivot_longer(
    cols = c(init_depth, breed_season_depth, recession, 
             pre_recession, post_recession, dry_days, reversals),
    names_to = "metric",
    values_to = "value"
  )

ggplot(water_long, aes(x = year, y = value, color = region, group = region)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    x = "Water Year",
    y = "Value",
    color = "Region",
    title = "Raw Water Dynamics by Water Year"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


# each hydro year initation time  -----------------------------------------


colonies_url <- "https://raw.githubusercontent.com/weecology/EvergladesWadingBird/main/SiteandMethods/colonies.csv"
colonies <- read_csv(colonies_url)

bird_region_summary <- bird_data |>
  left_join(colonies, by = "colony") |>
  group_by(water_year, region) |>
  summarise(y = min(julian_wy, na.rm = TRUE), .groups = "drop") |>
  rename(year = water_year)

bird_average_summary <- bird_data |>
  left_join(colonies, by = "colony") |>
  group_by(water_year) |>
  summarise(y = mean(julian_wy, na.rm = TRUE), .groups = "drop") |>
  rename(year = water_year)

ggplot() +
  geom_smooth(data = bird_average_summary,
              aes(x = year, y = y))+
  # geom_line(data = bird_average_summary, colour = 'grey',
  #           aes(x = year, y = y), linewidth = 1.5) +
  geom_point(data = bird_average_summary, colour = 'black',
             aes(x = year, y = y), size = 4) +
  geom_line(data = bird_region_summary, 
            aes(x = year, y = y, color = region, group = region), alpha = 0.7) +
  geom_point(data = bird_region_summary, 
             aes(x = year, y = y, color = region, group = region), size = 1.5) +
  labs(
    x = "Water Year",
    y = "First Nest Initiation Date (Day of Water Year)",
    color = "Region",
    title = "Wood Stork Nest Initiation by Region Across Water Years
    [black points & smooth average initation]"
  ) +
  
  theme_minimal()



# both in one  ------------------------------------------------------------



# 1. Format regional initiation to match water_long columns
bird_initiation_long <- bird_region_summary |>
  rename(value = y) |>
  mutate(metric = "nest_initiation")

# 2. Reshape water dynamics into long format
water_long <- water |>
  pivot_longer(
    cols = c(init_depth, breed_season_depth, recession,
             pre_recession, post_recession, dry_days, reversals),
    names_to = "metric",
    values_to = "value"
  )

# 3. Combine both datasets
combined_data <- bind_rows(water_long, bird_initiation_long)

# 4. Generate a shared color palette across all unique regions
all_regions <- sort(unique(combined_data$region))
region_palette <- setNames(
  scales::hue_pal()(length(all_regions)),
  all_regions
)

# 5. Plot water metrics and regional nest initiation together
ggplot(combined_data, aes(x = year, y = value, color = region, group = region)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_color_manual(values = region_palette) +
  labs(
    x = "Water Year",
    y = "Value",
    color = "Region",
    title = "Raw Water Dynamics and Nest Initiation by Water Year"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
