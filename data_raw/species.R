#different species datasets

library(usethis)
library(readr)

human = read.csv("./data_raw/human.txt", sep = "\t")
usethis::use_data(human, overwrite = TRUE)

c.elegans = read.csv("./data_raw/c.elegans.txt", sep = "\t")
usethis::use_data(c.elegans, overwrite = TRUE)

drosophila = read.csv("./data_raw/drosophila.txt", sep = "\t")
usethis::use_data(drosophila, overwrite = TRUE)

mouse = read.csv("./data_raw/mouse.txt", sep = "\t")
usethis::use_data(mouse, overwrite = TRUE)

rat = read.csv("./data_raw/rat.txt", sep = "\t")
usethis::use_data(rat, overwrite = TRUE)

zebrafish = read.csv("./data_raw/zebrafish.txt", sep = "\t")
usethis::use_data(zebrafish, overwrite = TRUE)

