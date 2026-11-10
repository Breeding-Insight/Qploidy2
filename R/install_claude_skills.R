#' Install Qploidy2 Claude AI Skills
#'
#' Downloads the Qploidy2 Claude Code skill files from GitHub and places them
#' in the \code{.claude/skills/} directory of your current R working directory,
#' one folder per skill (\code{.claude/skills/<skill>/SKILL.md}). Once
#' installed, the skills are automatically available to Claude Code when you
#' open that directory as a project.
#'
#' Two skills are installed:
#' \itemize{
#'   \item \strong{qploidy2-input-prep}: guides input preparation from VCF,
#'     Axiom, or Illumina array files.
#'   \item \strong{qploidy2-analysis}: covers the full analysis pipeline --
#'     standardization, BAF model selection, HMM copy number estimation,
#'     filtering, visualization, and export.
#' }
#'
#' @param dest_dir Destination directory for the skills. Defaults to
#'   \code{.claude/skills/} inside the current working directory. Each skill is
#'   written to \code{<dest_dir>/<skill>/SKILL.md}.
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
install_claude_skills <- function(dest_dir = file.path(getwd(), ".claude", "skills"),
                                  overwrite = FALSE) {
  base_url <- paste0(
    "https://raw.githubusercontent.com/",
    "Breeding-Insight/Qploidy2/main/.claude/skills/"
  )

  skills <- c("qploidy2-input-prep", "qploidy2-analysis")

  installed <- character(0)

  for (skill in skills) {
    skill_dir <- file.path(dest_dir, skill)
    dest_file <- file.path(skill_dir, "SKILL.md")

    if (file.exists(dest_file) && !overwrite) {
      message("Skipping ", skill, " (already exists). Use overwrite = TRUE to replace.")
      next
    }

    if (!dir.exists(skill_dir)) dir.create(skill_dir, recursive = TRUE)

    url <- paste0(base_url, skill, "/SKILL.md")
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
      "Or copy the skill folders to any other project's .claude/skills/ folder."
    )
  }

  invisible(installed)
}
