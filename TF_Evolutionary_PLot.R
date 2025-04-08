##############################################################################
# 0) Load packages
##############################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

##############################################################################
# 1) Define species + groups
##############################################################################
# First, define your original species_groups
species_groups <- data.frame(
  species = c("Cnc", "Cre", "Csc", "Mpu", "Olu", "Orc", "Ota",
              "Mpo", 
              "Ppa", 
              "Smo",
              "Pab", "Pme", "Psi",
              "Osi", "Zma",
              "Ath", "Ptr", "Vvi"),
  group = c(
    rep("Green Algae", 7), 
    "Liverwort", 
    "Moss", 
    "Lycophyte",  # <-- We'll rename this shortly
    rep("Gymnosperm", 3),
    rep("Monocots", 2),
    rep("Eudicots", 3)
  ),
  stringsAsFactors = FALSE
)

# Swap the role: rename group="Lycophyte" to group="Spike Moss"
species_groups <- species_groups %>%
  mutate(
    group = ifelse(group == "Lycophyte", "Spike Moss", group)
  )

total_width <- 100
n_sp       <- nrow(species_groups)

species_positions <- species_groups %>%
  mutate(
    x_position = seq(
      from  = total_width/(2*n_sp),
      by    = total_width/n_sp,
      length.out = n_sp
    )
  )

species_width <- total_width / n_sp

group_boundaries <- species_positions %>%
  group_by(group) %>%
  summarize(
    start_x = min(x_position) - species_width/2,
    end_x   = max(x_position) + species_width/2,
    center_x = (start_x + end_x)/2,
    .groups  = "drop"
  )

boundary_lines <- group_boundaries %>%
  arrange(end_x) %>%
  mutate(boundary_x = end_x)

##############################################################################
# 2) Super‐groups + ordering
##############################################################################
# Also swap the lines so that the group "Spike Moss" belongs to super_group "Lycophytes"
taxonomic_groups <- data.frame(
  group = c("Green Algae", "Liverwort", "Moss", "Spike Moss", 
            "Gymnosperm", "Monocots", "Eudicots"),
  super_group = c("Green Algae", "Bryophytes", "Bryophytes", 
                  "Lycophytes",  # <-- was "Spike Moss" before
                  "Gymnosperms", "Angiosperms", "Angiosperms"),
  super_group_order = c(1, 2, 2, 3, 4, 5, 5),
  stringsAsFactors = FALSE
)

species_positions <- species_positions %>%
  left_join(taxonomic_groups, by = "group")

group_boundaries <- group_boundaries %>%
  left_join(taxonomic_groups, by = "group")

# Summaries for super‐group labeling
super_group_boundaries <- group_boundaries %>%
  group_by(super_group) %>%
  summarize(
    start_x = min(start_x),
    end_x   = max(end_x),
    center_x = (start_x + end_x)/2,
    super_group_order = first(super_group_order),
    .groups = "drop"
  ) %>%
  arrange(super_group_order) %>%
  # SHIFT "Angiosperms" label if you like (adjust numeric offset as needed)
  mutate(
    center_x = ifelse(super_group == "Angiosperms", center_x - 1, center_x)
  )

##############################################################################
# 3) Presence data + pivot
##############################################################################
presence_mat_long <- presence_mat %>%
  pivot_longer(
    cols = -Family,
    names_to = "species",
    values_to = "presence"
  ) %>%
  left_join(species_positions, by = "species")

tf_origin <- presence_mat_long %>%
  filter(presence == 1) %>%
  group_by(Family) %>%
  summarize(
    first_super_group = super_group[which.min(
      match(super_group, unique(super_group_boundaries$super_group))
    )],
    first_group = group[which.min(
      match(group, unique(group_boundaries$group))
    )],
    .groups = "drop"
  )

##############################################################################
# 4) Factor levels for bottom->top
##############################################################################
# Also rename "Spike Moss" -> "Lycophytes" in your super_group_levels
super_group_levels <- c("Green Algae", "Bryophytes", "Gymnosperms",
                        "Angiosperms", "Lycophytes")

tf_origin <- tf_origin %>%
  mutate(
    first_super_group = factor(first_super_group, levels = super_group_levels)
  )

ordered_families <- tf_origin %>%
  arrange(first_super_group, first_group, Family) %>%
  pull(Family)

family_positions <- data.frame(
  Family     = ordered_families,
  y_position = seq_along(ordered_families)
) %>%
  left_join(tf_origin, by = "Family")

##############################################################################
# 5) Continuous presence spans
##############################################################################
full_grid <- expand.grid(
  Family  = ordered_families,
  species = species_positions$species,
  stringsAsFactors = FALSE
) %>%
  left_join(species_positions, by = "species") %>%
  left_join(
    presence_mat_long %>% select(Family, species, presence),
    by = c("Family","species")
  ) %>%
  mutate(presence = ifelse(is.na(presence), 0, presence)) %>%
  left_join(family_positions, by = "Family")

span_data <- full_grid %>%
  arrange(Family, x_position) %>%
  group_by(Family) %>%
  mutate(
    change = presence != lag(presence, default = presence[1]),
    run_id = cumsum(change)
  ) %>%
  group_by(Family, run_id, presence) %>%
  summarize(
    start_x = min(x_position) - species_width/2,
    end_x   = max(x_position) + species_width/2,
    y_pos   = first(y_position),
    first_super_group = first(first_super_group),
    .groups = "drop"
  ) %>%
  filter(presence == 1)

##############################################################################
# 6) Darker color palette for stronger outlines
##############################################################################
# Swap the "Spike Moss" entry with "Lycophytes"
origin_colors <- c(
  "Green Algae"  = "#339999",  
  "Bryophytes"   = "#FFB266",  
  "Gymnosperms"  = "#FF704D",  
  "Angiosperms"  = "#FF3333",
  "Lycophytes"   = "#B266B2"   # now called "Lycophytes"
)

x_max <- total_width
y_min <- 0
y_max <- length(ordered_families) + 1

##############################################################################
# 7) LINEAGE-BASED label positions
#    "Lycophytes" => center under Osi
#    "Bryophytes" => center under Psi
#    else => x_max/2
##############################################################################
osi_x <- species_positions %>%
  filter(species == "Osi") %>%
  pull(x_position)

psi_x <- species_positions %>%
  filter(species == "Psi") %>%
  pull(x_position)

# Adjust case_when to match new super_group name "Lycophytes" instead of "Spike Moss"
family_positions <- family_positions %>%
  mutate(
    label_x = case_when(
      first_super_group == "Lycophytes"   ~ osi_x,      # formerly "Spike Moss"
      first_super_group == "Bryophytes"   ~ psi_x,
      TRUE                                ~ x_max / 2
    )
  )

##############################################################################
# 8) Main plot
##############################################################################
p_main <- ggplot() +
  # No background rect
  geom_vline(
    data = boundary_lines,
    aes(xintercept = boundary_x),
    color = "gray40",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = seq(y_min + 0.5, y_max - 0.5, by = 1),
    color = "gray90"
  ) +
  geom_rect(
    data = span_data,
    aes(
      xmin = start_x,
      xmax = end_x,
      ymin = y_pos - 0.4,
      ymax = y_pos + 0.4,
      color = first_super_group
    ),
    fill = "white",
    linewidth = 0.7
  ) +
  # Use our new label_x for gene family text
  geom_text(
    data = family_positions,
    aes(
      x = label_x,
      y = y_position,
      label = Family
    ),
    size = 2.5,
    hjust = 0.5
  ) +
  # Only 3 keys in the legend: "Green Algae", "Bryophytes", "Lycophytes"
  scale_color_manual(
    values = origin_colors,
    breaks = c("Green Algae", "Bryophytes", "Lycophytes"),
    name   = "Lineages",
    guide  = guide_legend(title.position = "top", nrow=1)
  ) +
  scale_x_continuous(
    breaks = species_positions$x_position,
    labels = species_positions$species,
    limits = c(0, x_max),
    expand = c(0,0),
    position = "top"
  ) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    expand = c(0,0)
  ) +
  theme_minimal() +
  theme(
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    panel.grid    = element_blank(),
    panel.border  = element_rect(color = "black", fill = NA),
    plot.margin   = margin(t = 0, r = 10, b = 0, l = 10),
    legend.position = "bottom"
  ) +
  labs(x = NULL, y = NULL)

##############################################################################
# 9) Header plot
##############################################################################
header_plot <- ggplot() +
  geom_vline(
    data = boundary_lines,
    aes(xintercept = boundary_x),
    color = "gray40",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  # Super-group text (with "Lycophytes" included)
  geom_text(
    data = super_group_boundaries,
    aes(x = center_x, y = 0.7, label = super_group),
    size = 3.5,
    fontface = "bold"
  ) +
  # Sub-group labels, skip duplicates + "Gymnosperm"
  geom_text(
    data = group_boundaries %>%
      filter(group != super_group & group != "Gymnosperm"),
    aes(x = center_x, y = 0.3, label = group),
    size = 3
  ) +
  scale_x_continuous(limits = c(0, x_max), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1)) +
  theme_void() +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 0, l = 10)
  )

##############################################################################
# 10) Combine header + main
##############################################################################
final_plot <- cowplot::plot_grid(
  header_plot,
  p_main,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(0.1, 0.9)
)

ggsave("tf_distribution_lycophytes.png",
       final_plot, width = 12, height = 20, dpi = 300)
print(final_plot)
