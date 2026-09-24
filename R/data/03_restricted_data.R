library(tidyverse)

c('./data/gs_data.csv','./data/wide_data.csv','./data/narrow_data.csv') %>%
  walk({
    function(datafile) {
      outfile <- sub("\\.csv$", "_restricted.csv", datafile)
      read_csv(datafile) %>%
        filter(year > 2022) %>%
        write_csv(outfile)
    }
  })
