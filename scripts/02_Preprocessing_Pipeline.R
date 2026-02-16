# ------------------------------------------------------------
# 02_Preprocessing_Pipeline.R
#
# Cardiac feedback processing – preprocessing pipeline
#
# This script implements the preprocessing steps for
# feedback-locked cardiac responses, including:
#   - marker cleaning and trial segmentation
#   - IBI computation and alignment to feedback onset
#   - exclusion of invalid trials
#
# Due to data protection regulations, raw data files,
# subject identifiers, and exclusion lists are not included.
# ------------------------------------------------------------


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

rm(list = ls())

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# ------------------------------------------------------------
# 2. Paths & Global Parameters
# ------------------------------------------------------------

data_dir   <- "data"
output_dir <- "outputs"

if (!dir.exists(output_dir)) dir.create(output_dir)

# Length of one task segment (-4s to +5s)
segment_length  <- 9000

# Expected relative position of feedback trigger
trigger_rel_pos <- 4001

# ------------------------------------------------------------
# 3. Import Marker Files
# ------------------------------------------------------------

marker_ibi <- list.files(
  path = data_dir,
  pattern = "_2.Markers",
  full.names = TRUE
) %>%
  set_names() %>%
  map_df(fread, .id = "ID")

marker_ibi <- marker_ibi %>%
  mutate(
    ID = basename(ID),
    ID = gsub("_2\\.Markers$", "", marker_ibi$ID) 
  )

# Fix known ID inconsistencies
marker_ibi <- marker_ibi %>%
  mutate(
    ID = case_when(
      ID == "Doors_FOR11141x_AllTriggers"       ~ "Doors_FOR11141_AllTriggers",
      ID == "DOORS_FOR14166_AllTriggers"        ~ "Doors_FOR14166_AllTriggers",
      ID == "EmoEEG_Doors_M1003_AllTriggers"    ~ "Doors_Emo_M1003_AllTriggers",
      ID == "EmoEEG_Reversal_M1003_AllTriggers" ~ "Reversal_Emo_M1003_AllTriggers",
      TRUE ~ ID
    )
  )

# Backup raw version
marker_ibi_raw <- marker_ibi

fwrite(
  marker_ibi_raw,
  file = file.path(output_dir, paste0("marker_ibi_raw_", Sys.Date(), ".csv"))
)

# ------------------------------------------------------------
# 4. Basic Marker Cleanup
# ------------------------------------------------------------

marker_ibi <- marker_ibi %>%
  mutate(
    Type        = str_remove_all(Type, " "),
    Description = str_remove_all(Description, " ")
  )

# ------------------------------------------------------------
# 5. Trial Definition
# ------------------------------------------------------------
# Trials are defined exclusively by NewSegment markers.
# No fixed time-based segmentation is used in this analysis.

setDT(marker_ibi)
marker_ibi[, trialnr := cumsum(Type == "NewSegment"), by = ID]

marker_ibi <- marker_ibi %>%
  mutate(
    task = case_when(
      str_detect(ID, "Doors")    ~ "Doors",
      str_detect(ID, "Reversal") ~ "Reversal",
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# 6. Remove Practice Trials (Reversal Only)
# ------------------------------------------------------------

practice_trials <- marker_ibi %>%
  filter(task == "Reversal", Description == "S111") %>%
  pull(trialnr)

marker_ibi <- marker_ibi %>%
  filter(!(task == "Reversal" & trialnr %in% practice_trials))

# ------------------------------------------------------------
# 7. Clean Subject IDs & Merge Demographics
# ------------------------------------------------------------

marker_ibi <- marker_ibi %>%
  mutate(
    CleanID = str_extract(ID, "(?<=Doors_|Reversal_).*?(?=_AllTriggers)")
  )

# Correct known mismatches
marker_ibi <- marker_ibi %>%
  mutate(
    CleanID = case_when(
      CleanID == "FOR11907" ~ "FOR11002",
      CleanID == "FOR11908" ~ "FOR11003",
      CleanID == "FOR11910" ~ "FOR11005",
      TRUE ~ CleanID
    )
  )

demographics <- fread(file.path(data_dir, "Demographics.csv"))

marker_ibi <- marker_ibi %>%
  left_join(
    demographics %>% select(ID, Group, Age, Gender, BMI),
    by = c("CleanID" = "ID")
  )

# ------------------------------------------------------------
# 8. R-Peaks & IBI Computation
# ------------------------------------------------------------

r_peaks <- marker_ibi %>%
  filter(Description == "R") %>%
  arrange(ID, trialnr, Position) %>%
  group_by(ID, trialnr) %>%
  mutate(ibi_ms = Position - lag(Position)) %>%
  ungroup()

# ------------------------------------------------------------
# 9. Feedback Trigger Detection
# ------------------------------------------------------------

doors_triggers <- marker_ibi %>%
  filter(
    Description %in% c("S44", "S55"),
    Position %% segment_length == trigger_rel_pos
  )

reversal_triggers <- marker_ibi %>%
  filter(
    Description %in% c("S100", "S200"),
    Position %% segment_length == trigger_rel_pos
  )

triggers <- bind_rows(doors_triggers, reversal_triggers) %>%
  transmute(
    ID,
    trialnr,
    trigger     = Description,
    trigger_pos = Position,
    feedback = case_when(
      trigger %in% c("S44", "S100") ~ "pos",
      trigger %in% c("S55", "S200") ~ "neg"
    )
  )

# ------------------------------------------------------------
# 10. IBI0 Identification & IBI Change
# ------------------------------------------------------------

r_peaks <- r_peaks %>%
  group_by(ID, trialnr) %>%
  arrange(Position) %>%
  mutate(
    ibi_start = lag(Position),
    ibi_end   = Position
  ) %>%
  ungroup() %>%
  left_join(triggers, by = c("ID", "trialnr")) %>%
  mutate(
    ibi_0_flag = trigger_pos > ibi_start & trigger_pos <= ibi_end
  )

r_peaks <- r_peaks %>%
  group_by(ID, trialnr) %>%
  mutate(
    ibi0_index = which(ibi_0_flag),
    ibi_name   = row_number() - ibi0_index
  ) %>%
  ungroup()

r_peaks <- r_peaks %>%
  filter(between(ibi_name, -2, 2)) %>%
  group_by(ID, trialnr) %>%
  mutate(
    ibi_minus2 = ibi_ms[ibi_name == -2],
    ibi_change = ifelse(ibi_name == -2, NA, ibi_ms - ibi_minus2)
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 11. Master Table
# ------------------------------------------------------------

marker_info <- marker_ibi %>%
  select(CleanID, task, trialnr, Group, Age, Gender, BMI) %>%
  distinct()

r_peaks_info <- r_peaks %>%
  select(CleanID, task, trialnr, ibi_name, ibi_change)

triggers_info <- triggers %>%
  left_join(
    marker_ibi %>% select(ID, CleanID, task) %>% distinct(),
    by = "ID"
  ) %>%
  select(CleanID, task, trialnr, feedback)

final_table <- marker_info %>%
  left_join(r_peaks_info,  by = c("CleanID", "task", "trialnr")) %>%
  left_join(triggers_info, by = c("CleanID", "task", "trialnr"))

# ------------------------------------------------------------
# 12. Save Outputs
# ------------------------------------------------------------

fwrite(
  final_table,
  file = file.path(output_dir, paste0("final_table_", Sys.Date(), ".csv"))
)

# ------------------------------------------------------------
# 12. Trial Exclusion (Documentation Only)
# ------------------------------------------------------------

# A predefined list of trials with technical artifacts or
# implausible cardiac intervals was excluded prior to the
# statistical analyses.
#
# The exclusion list is not included in this public repository
# due to data protection regulations.

