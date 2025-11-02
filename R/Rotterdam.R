library(ggplot2)
library(ggpubr)
library(devtools)
install_github("grichards230/smoothC")
library(smoothC)
library(survival)

compare_models(
  c("age","meno","hormon"),
  c("nodes","pgr"),
  "rtime", "recur", rotterdam
)

ggsave(filename="Rotterdam.pdf",height=7,width=14)



