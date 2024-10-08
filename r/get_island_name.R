get_island_name <- function(island_num) {

  switch(paste(island_num, collapse = ","),
         "20" = "Nesøy",
         "22" = "Myken",
         "23" = "Træna",
         "24" = "Selvær",
         "25" = "Sanna",
         "22,23,24,34,35" = "Non-farm islands",
         "26" = "Gjerøy",
         "27" = "Hestmannøy",
         "28" = "Indre Kvarøy",
         "30" = "Selsøyvik",
         "33" = "Lurøy-Onøy",
         "331" = "Lurøy",
         "332" = "Onøy",
         "34" = "Lovund",
         "35" = "Sleneset",
         "38" = "Aldra",
         "20,26,27,28,33,38" = "Farm islands",
         "20,22,23,24,26,27,28,33,34,35,38" = "Helgeland system", # TODO: change to Helgeland System
         "60" = "Leka",
         "61" = "Vega",
         "67" = "Lauvøya-Selnes-Flenstad",
         "68" = "Rånes",
         "60,61,63,67,68" = "Southern system",
         "77" = "various1",
         "88" = "various2",
         NA)
}
