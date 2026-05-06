library(readr)
library(dplyr)
library(tidyr)
library(janitor)
library(stringr)
library(ggplot2)
library(hms)

# Read CSV
raw <- read_csv("data/ethanol2_3_5Bleach50_100_200ppm_PQ50_100_150uM_D74D105_ex3.csv")

# Extract data
time_row <- 26
dat <- raw[(time_row +1):nrow(raw), 2:93]
colnames(dat) <- as.character(raw[time_row, 2:93])
dat <- dat[complete.cases(dat[[1]]),]

# Rename wells
colnames(dat)[1:92] <- c(
  "time", "Temperature",
  "A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12",
  "B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12",
  "C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12",
  "D4","D5","D6","D7","D8","D9","D10","D11","D12",
  "E4","E5","E6","E7","E8","E9","E10","E11","E12",
  "F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12",
  "G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12",
  "H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12"
)

# Convert wells to numeric
well_cols <- 3:ncol(dat)
dat[well_cols] <- lapply(dat[well_cols], function (x) as.numeric(as.character(x)))

# Long format
dat_long <- dat %>%
  pivot_longer(
    cols = 3:ncol(dat),
    names_to = "Well",
    values_to = "OD"
  )

# BLEACH
dat_long_bleach <- dat_long |>
  mutate(
    Condition = case_when(
      Well %in% c("A1","A2","A3") ~ "A909_Untreated",
      Well %in% c("A4","A5","A6") ~ "A909_50ppm",
      Well %in% c("A7","A8","A9") ~ "A909_100ppm",
      Well %in% c("A10","A11","A12") ~ "A909_200ppm",
      Well %in% c("B1","B2","B3") ~ "D74_Untreated",
      Well %in% c("B4","B5","B6") ~ "D74_50ppm",
      Well %in% c("B7","B8","B9") ~ "D74_100ppm",
      Well %in% c("B10","B11","B12") ~ "D74_200ppm",
      Well %in% c("C1","C2","C3") ~ "D105_Untreated",
      Well %in% c("C4","C5","C6") ~ "D105_50ppm",
      Well %in% c("C7","C8","C9") ~ "D105_100ppm",
      Well %in% c("C10","C11","C12") ~ "D105_200",
      TRUE ~ NA_character_
    ),
    Replicate = case_when(
      Well %in% c("A1","A4","A7","A10") ~ "R1",
      Well %in% c("A2","A5","A8","A11") ~ "R2",
      Well %in% c("A3","A6","A9","A12") ~ "R3",
      Well %in% c("B1","B4","B7","B10") ~ "R4",
      Well %in% c("B2","B5","B8","B11") ~ "R5",
      Well %in% c("B3","B6","B9","B12") ~ "R6",
      Well %in% c("C1","C4","C7","C10") ~ "R7",
      Well %in% c("C2","C5","C8","C11") ~ "R8",
      Well %in% c("C3","C6","C9","C12") ~ "R9",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(Condition))

# Convert time (BLEACH)
dat_long_bleach <- dat_long_bleach |>
  mutate(time_hours = as.numeric(hms::as_hms(time)) / 3600) |>
  filter(!is.na(time_hours))

# Summary (BLEACH)
summary_data_bleach <- dat_long_bleach |>
  group_by(time_hours, Condition) |>
  summarise(
    mean_OD = mean(OD, na.rm = TRUE),
    sd_OD   = sd(OD, na.rm = TRUE),
    n       = sum(!is.na(OD)),
    se_OD   = sd_OD / sqrt(n),
    .groups = "drop"
  )

# Plot 1 (BLEACH)
ggplot(summary_data_bleach, aes(x = time_hours, y = mean_OD, color = Condition)) +
  geom_line(linewidth = 0.7, alpha = 0.5) +
  geom_point(size = 1.3, alpha = 0.8) +
  geom_smooth(se = FALSE, linewidth = 1.3, span = 0.25) +
  geom_ribbon(
    aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD, fill = Condition),
    alpha = 0.04,
    color = NA
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "GBS Growth Curves For Bleach",
    x = "Time (hours)",
    y = "OD600"
  )
ggsave("figure/Bleach_Treatment_Growth_Curve.png", width = 6, height = 4, dpi = 1200)

# ETHANOL ANALYSIS
dat_long_ethanol <- dat_long %>%
  mutate(
    Condition = case_when(
      Well %in% c("A1","A2","A3") ~ "A909_Untreated",
      Well %in% c("D4","D5","D6") ~ "A909_2%",
      Well %in% c("D7","D8","D9") ~ "A909_3%",
      Well %in% c("D10","D11","D12") ~ "A909_5%",
      
      Well %in% c("B1","B2","B3") ~ "D74_Untreated",
      Well %in% c("E4","E5","E6") ~ "D74_2%",
      Well %in% c("E7","E8","E9") ~ "D74_3%",
      Well %in% c("E10","E11","E12") ~ "D74_5%",
      
      Well %in% c("C1","C2","C3") ~ "D105_Untreated",
      Well %in% c("F1","F2","F3") ~ "D105_2%",
      Well %in% c("F4","F5","F6") ~ "D105_3%",
      Well %in% c("F7","F8","F9") ~ "D105_5%",
      
      TRUE ~ NA_character_
    ),
    Replicate = case_when(
      Well %in% c("A1","D4","D7","D10") ~ "R1",
      Well %in% c("B1","E4","E7","E10") ~ "R2",
      Well %in% c("C1","F1","F4","F7") ~ "R3",
      Well %in% c("A2","D5","D8","D6") ~ "R4",
      Well %in% c("B2","E5","E8","E11") ~ "R5",
      Well %in% c("C2","F2","F5","F8") ~ "R6",
      Well %in% c("A3","D6","D9","D12") ~ "R7",
      Well %in% c("B3","E6","E9","E12") ~ "R8",
      Well %in% c("C3","F3","F6","F9") ~ "R9",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(Condition))

# Convert time (Ethanol)
dat_long_ethanol <- dat_long_ethanol |>
  mutate(time_hours = as.numeric(hms::as_hms(time)) / 3600) |>
  filter(!is.na(time_hours))

# summary
summary_ethanol <- dat_long_ethanol %>%
  group_by(time_hours, Condition) %>%
  summarise(
    mean_OD = mean(OD, na.rm = TRUE),
    sd_OD   = sd(OD, na.rm = TRUE),
    n       = sum(!is.na(OD)),
    se_OD   = sd_OD / sqrt(n),
    .groups = "drop"
  )

# Plot 2 (ETHANOL)
ggplot(summary_ethanol, aes(x = time_hours, y = mean_OD, color = Condition)) +
  geom_line(linewidth = 0.7, alpha = 0.5) +
  geom_point(size = 1.3, alpha = 0.8) +
  geom_smooth(se = FALSE, linewidth = 1.3, span = 0.25) +
  geom_ribbon(
    aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD, fill = Condition),
    alpha = 0.04,
    color = NA
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "GBS Growth Curves For Ethanol",
    x = "Time (hours)",
    y = "OD600"
  )

ggsave("figure/Ethanol_Treatment_Growth_Curve.png", width = 6, height = 4, dpi = 1200)

#PQ ANALYSIS
dat_long_pq <- dat_long %>%
  mutate(
    Condition = case_when(
      Well %in% c("A1","A2","A3") ~ "A909_Untreated",
      Well %in% c("F10","F11","F12") ~ "A909_50PQ",
      Well %in% c("G1","G2","G3") ~ "A909_100PQ",
      Well %in% c("G4","G5","G6") ~ "A909_150PQ",
      
      Well %in% c("B1","B2","B3") ~ "D74_Untreated",
      Well %in% c("G7","G8","G9") ~ "D74_50PQ",
      Well %in% c("G10","G11","G12") ~ "D74_100PQ",
      Well %in% c("H1","H2","H3") ~ "D74_150PQ",
      
      Well %in% c("C1","C2","C3") ~ "D105_Untreated",
      Well %in% c("H4","H5","H6") ~ "D105_50PQ",
      Well %in% c("H7","H8","H9") ~ "D105_100PQ",
      Well %in% c("H10","H11","H12") ~ "D105_150PQ",
      
      TRUE ~ NA_character_
    ),
    Replicate = case_when(
      Well %in% c("A1","F10","G1","G4") ~ "R1",
      Well %in% c("B1","F11","G2","G5") ~ "R2",
      Well %in% c("C1","F12","G3","G6") ~ "R3",
      Well %in% c("A2","G7","G10","H1") ~ "R4",
      Well %in% c("B2","G8","G11","H2") ~ "R5",
      Well %in% c("C2","G9","G12","H3") ~ "R6",
      Well %in% c("A3","H4","H5","H6") ~ "R7",
      Well %in% c("B3","H7","H8","H9") ~ "R8",
      Well %in% c("C3","H10","H11","H12") ~ "R9",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(Condition))

# Convert time (Ethanol)
dat_long_pq <- dat_long_pq |>
  mutate(time_hours = as.numeric(hms::as_hms(time)) / 3600) |>
  filter(!is.na(time_hours))

# Summary
summary_pq <- dat_long_pq %>%
  group_by(time_hours, Condition) %>%
  summarise(
    mean_OD = mean(OD, na.rm = TRUE),
    sd_OD   = sd(OD, na.rm = TRUE),
    n       = sum(!is.na(OD)),
    se_OD   = sd_OD / sqrt(n),
    .groups = "drop"
  )

# Plot 3 (PQ)
ggplot(summary_pq, aes(x = time_hours, y = mean_OD, color = Condition)) +
  geom_line(linewidth = 0.7, alpha = 0.5) +
  geom_point(size = 1.3, alpha = 0.8) +
  geom_smooth(se = FALSE, linewidth = 1.3, span = 0.25) +
  geom_ribbon(
    aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD, fill = Condition),
    alpha = 0.04,
    color = NA
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "GBS Growth Curves For PQ",
    x = "Time (hours)",
    y = "OD600"
  )
ggsave("figure/PQ_Treatment_Growth_Curve.png", width = 6, height = 4, dpi = 1200)

# Add a Treatment label
summary_data_bleach <- summary_data_bleach %>%
  mutate(Treatment = "Bleach")

summary_ethanol <- summary_ethanol %>%
  mutate(Treatment = "Ethanol")

summary_pq <- summary_pq %>%
  mutate(Treatment = "PQ")

#Combine
combined_summary <- bind_rows(
  summary_data_bleach,
  summary_ethanol,
  summary_pq
)
# Plot 4 (Combined)
ggplot(combined_summary, aes(x = time_hours, y = mean_OD, color = Condition)) +
  geom_line(linewidth = 0.7, alpha = 0.5) +
  geom_point(size = 1.3, alpha = 0.8) +
  geom_smooth(
    se = FALSE,
    linewidth = 1.3,
    span = 0.25   # smaller = tighter curve (adjust 0.2–0.4)
  ) +
  geom_ribbon(
    aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD, fill = Condition),
    alpha = 0.04,
    color = NA
  ) +
  
  facet_wrap(~Treatment, scales = "free_y") +
  
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  
  labs(
    title = "GBS Growth Curves Across Treatments",
    x = "Time (hours)",
    y = "OD600"
  )
ggsave("figure/Combined.png", width = 18, height = 4, dpi = 1200)