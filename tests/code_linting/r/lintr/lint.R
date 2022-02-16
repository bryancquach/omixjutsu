library(optparse, quietly = T)
options(stringsAsFactors = F)

p <- OptionParser(
  usage = "\n%prog [options] code_file",
  description = "\nR Code Linter\n\nA wrapper script for the lintr R package.\n\nChecks an R code file for style consistency with a specified style guide. By default it uses a modified version of Google's R style guide.",
  prog = "Rscript lint.R"
)
p <- add_option(
  object = p, 
  opt_str = c("-c", "--config"), 
  default = NULL, 
  type = "character",
  help = "A .lintr DCF file that specifies the linters to use."
)
argv <- parse_args(p, positional_arguments = 1)
config_file <- argv$options$config
code_file <- argv$args
if (is.null(config_file)) {
  config_file <- normalizePath("tests/code_linting/r/lintr/.lintr")
} else {
  config_file <- normalizePath(config_file)
}
cat("\nUsing .lintr file", config_file, "...\n")
if (!file.exists(config_file)) {
  warnings(
    paste0(
      "Default .lintr file could not be found at ", 
      config_file, 
      ". Reverting to default lintr package settings"
    )
  )
}
options(lintr.linter_file = config_file)
cat("\nLinting file", code_file, "...\n\n")
lintr::lint(code_file)
