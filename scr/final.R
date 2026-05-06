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
dat_long_bleach <- dat %>%
  pivot_longer(
    cols = 3:ncol(dat),
    names_to = "Well",
    values_to = "OD"
  )

# BLEACH
dat_long_bleach <- dat_long_bleach %>%
  mutate(
    Condition = case_when(
      Well %in% c("A1","A2","A3") ~ "A909_un",
      Well %in% c("A4","A5","A6") ~ "A909_50ppm",
      Well %in% c("A7","A8","A9") ~ "A909_100ppm",
      Well %in% c("A10","A11","A12") ~ "A909_200ppm",
      Well %in% c("B1","B2","B3") ~ "A909_D75_un",
      Well %in% c("B4","B5","B6") ~ "A909_D75_50ppm",
      Well %in% c("B7","B8","B9") ~ "A909_D75_100ppm",
      Well %in% c("B10","B11","B12") ~ "A909_D75_200ppm",
      Well %in% c("C1","C2","C3") ~ "A909_D105_un",
      Well %in% c("C4","C5","C6") ~ "A909_D105_50ppm",
      Well %in% c("C7","C8","C9") ~ "A909_D105_100ppm",
      Well %in% c("C10","C11","C12") ~ "A909_D105_200",
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
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD, fill = Condition),
              alpha = 0.09, color = NA) +
  theme_classic() +
  labs(
    title = "GBS Growth Curve (Bleach)",
    x = "Time (hours)",
    y = "OD600"
  )

ggsave("figure/Bleach_Treatment_Growth_Curve.png", width = 6, height = 4, dpi = 1200)

# ETHANOL ANALYSIS
dat_long_ethanol <- dat_long %>%
  mutate(
    Condition = case_when(
      Well %in% c("A1","A2","A3") ~ "A909_2%",
      Well %in% c("A4","A5","A6") ~ "A909_3%",
      Well %in% c("A7","A8","A9") ~ "A909_5%",
      
      Well %in% c("B1","B2","B3") ~ "A909_D74_2%",
      Well %in% c("B4","B5","B6") ~ "A909_D74_3%",
      Well %in% c("B7","B8","B9") ~ "A909_D74_5%",
      
      Well %in% c("C1","C2","C3") ~ "A909_D105_2%",
      Well %in% c("C4","C5","C6") ~ "A909_D105_3%",
      Well %in% c("C7","C8","C9") ~ "A909_D105_5%",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Condition)) %>%
  mutate(time_hours = as.numeric(hms::as_hms(time)) / 3600) %>%
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
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD, fill = Condition),
              alpha = 0.2, color = NA) +
  theme_classic() +
  labs(
    title = "GBS Ethanol Growth Curve",
    x = "Time (hours)",
    y = "OD600"
  )

ggsave("figure/Ethanol_Treatment_Growth_Curve.png", width = 6, height = 4, dpi = 1200)

#PQ ANALYSIS
dat_long_pq <- dat_long %>%
  mutate(
    Condition = case_when(
      # F-row equivalent of ethanol-style grouping
      Well %in% c("A1","A2","A3") ~ "A909_50PQ",
      Well %in% c("A4","A5","A6") ~ "A909_100PQ",
      Well %in% c("A7","A8","A9") ~ "A909_150PQ",
      
      Well %in% c("B1","B2","B3") ~ "A909_D74_50PQ",
      Well %in% c("B4","B5","B6") ~ "A909_D74_100PQ",
      Well %in% c("B7","B8","B9") ~ "A909_D74_150PQ",
      
      Well %in% c("C1","C2","C3") ~ "A909_D105_50PQ",
      Well %in% c("C4","C5","C6") ~ "A909_D105_100PQ",
      Well %in% c("C7","C8","C9") ~ "A909_D105_150PQ",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Condition)) %>%
  mutate(time_hours = as.numeric(hms::as_hms(time)) / 3600) %>%
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
  geom_line(size = 1.2) +
  geom_ribbon(
    aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD, fill = Condition),
    alpha = 0.2, color = NA
  ) +
  theme_classic() +
  labs(
    title = "GBS PQ Growth Curve",
    x = "Time (hours)",
    y = "OD600"
  )
ggsave("figure/PQ_Treatment_Growth_Curve.png", width = 6, height = 4, dpi = 1200)
