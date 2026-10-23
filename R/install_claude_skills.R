#' Install Qploidy2 Claude AI Skills
#'
#' Downloads the Qploidy2 Claude Code skill files from GitHub and places them
#' in the \code{.claude/agents/} directory of your current R working directory.
#' Once installed, the skills are automatically available to Claude Code when
#' you open that directory as a project.
#'
#' Two skills are installed:
#' \itemize{
#'   \item \strong{qploidy2-input-prep}: guides input preparation from VCF,
#'     Axiom, or Illumina array files.
#'   \item \strong{qploidy2-analysis}: covers the full analysis pipeline —
#'     standardization, BAF model selection, HMM copy number estimation,
#'     filtering, visualization, and export.
#' }
#'
#' @param dest_dir Destination directory for the skill files. Defaults to
#'   \code{.claude/agents/} inside the current working directory.
#' @param overwrite Logical. If \code{TRUE}, overwrite existing skill files.
#'   Default is \code{FALSE}.
#'
#' @return Invisibly returns the paths of the installed skill files.
#'
#' @examples
#' \dontrun{
#' install_claude_skills()
#' }
#'
#' @export
install_claude_skills <- function(dest_dir = file.path(getwd(), ".claude", "agents"),
                                  overwrite = FALSE) {
  base_url <- paste0(
    "https://raw.githubusercontent.com/",
    "Breeding-Insight/Qploidy2/main/.claude/agents/"
  )

  skills <- c("qploidy2-input-prep.md", "qploidy2-analysis.md")

  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
    message("Created directory: ", dest_dir)
  }

  installed <- character(0)

  for (skill in skills) {
    dest_file <- file.path(dest_dir, skill)

    if (file.exists(dest_file) && !overwrite) {
      message("Skipping ", skill, " (already exists). Use overwrite = TRUE to replace.")
      next
    }

    url <- paste0(base_url, skill)
    tryCatch({
      utils::download.file(url, destfile = dest_file, quiet = TRUE)
      message("Installed: ", dest_file)
      installed <- c(installed, dest_file)
    }, error = function(e) {
      warning("Failed to download ", skill, ": ", conditionMessage(e))
    })
  }

  if (length(installed) > 0) {
    message(
      "\nSkills installed successfully!\n",
      "Open this directory in Claude Code and ask Claude about the Qploidy2 workflow.\n",
      "Or copy the skill files to any other project's .claude/agents/ folder."
    )
  }

  invisible(installed)
}
