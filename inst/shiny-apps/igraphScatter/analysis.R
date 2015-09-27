# Import data and prepare objects for shiny app
# Data already cleaned. Lowly expressed or highly expressed genes were removed.
# save(anno_single, molecules_single_collision, molecules_cv,
#       file = "../data/igraphScatter-demo.rda")

#' @export
demo_data <- function(molecules_cv,
                          anno_single, molecules_single_collision,
                          per_person) {

  ii_person <- match(per_person, unique(anno_single$individual))
  # This object contains coefficieint of variation and mean
  # molecules_cv
  df_cvplot <- data.frame(mean = log10(molecules_cv[[ii_person]]$mean),
                          cv = molecules_cv[[ii_person]]$cv)

  # This object contains per gene count data across cells
  per_person <- unique(anno_single$individual)[ii_person]
  df_geneplot <- data.frame(
    molecules_single_collision[, anno_single$individual == per_person])
  df_geneplot <- as.matrix(df_geneplot)

  list(molecules_cv = molecules_cv,
       anno_single = anno_single,
       molecules_single_collision = molecules_single_collision,
       per_person = per_peson)
}


