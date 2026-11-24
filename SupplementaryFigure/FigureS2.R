# Supplementary Figure S2. Prediction in the sample with both clock layers.
# Rscript --vanilla SupplementaryFigure/FigureS2.R [output_directory]
# Use Unicode paths even when launched from a non-UTF-8 Windows shell.
if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
script_dir <- dirname(normalizePath(script_path, mustWork=TRUE))
# FIGURE_DATA_DIR points to a local directory containing the plotting inputs.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 1L) stop("Supply at most one output directory.", call.=FALSE)
input_root <- Sys.getenv("FIGURE_DATA_DIR")
if (!nzchar(input_root) || !dir.exists(path.expand(input_root))) {
  stop("Set FIGURE_DATA_DIR to your local plotting input directory.", call.=FALSE)
}
input_dir <- file.path(path.expand(input_root), "SupplementaryFigure")
if (!dir.exists(input_dir)) stop("Missing input directory: ", input_dir, call.=FALSE)
output_dir <- if (length(args)) path.expand(args[[1L]]) else file.path(dirname(script_dir), "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

source(file.path(script_dir, "plot_s2_s3.R"), encoding="UTF-8")
data_dir <- file.path(dirname(input_dir), "Figure6")
cross <- load_cross(data_dir, file.path(input_dir, "FigureS2"))
export_plot(file.path(output_dir, "FigureS2_two_layers"), 190.5, function() draw_cross(cross))
