options(repos = "https://cloud.r-project.org")

needed <- c(
  "caret","e1071","ggplot2","lattice","pROC","nnet","foreach",
  "ModelMetrics","recipes","lubridate","Matrix","plyr","stringr",
  "survival","tibble","withr","nlme","randomForest"
)

to_install <- setdiff(needed, rownames(installed.packages()))
if (length(to_install)) {
  install.packages(to_install, dependencies = c("Depends","Imports"))
}
