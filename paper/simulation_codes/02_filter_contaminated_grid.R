library(readr)
library(dplyr)

scenario <- read_csv("results/tables/scenario_grid_rrhfd_main.csv", show_col_types = FALSE)
job <- read_csv("results/tables/job_grid_rrhfd_main.csv", show_col_types = FALSE)

scenario_cont <- scenario %>%
  filter(study == "contaminated_rolling")

job_cont <- job %>%
  filter(study == "contaminated_rolling") %>%
  arrange(scenario_id, block_id) %>%
  mutate(job_id = row_number()) %>%
  select(job_id, everything())

write_csv(scenario_cont, "results/tables/contaminated/scenario_grid_rrhfd_contaminated.csv")
write_csv(job_cont, "results/tables/contaminated/job_grid_rrhfd_contaminated.csv")

cat("Contaminated scenarios:", nrow(scenario_cont), "\n")
cat("Contaminated jobs:", nrow(job_cont), "\n")
