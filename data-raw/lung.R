library("GEOquery")

# Get the GSE object
gse <- getGEO("GSE102287", GSEMatrix=FALSE)

# Get the GSM list
gsmlist <- GSMList(gse)

# Initialize the dataframe
df <- data.frame()
df2 <- data.frame()

# Loop through each item in gsmlist
for (i in 1:length(gsmlist)) {
  mygsm <- gsmlist[[i]]  # Get the current sample from the list

  # Extract characteristics of each sample
  data <- Meta(mygsm)$characteristics_ch1
  series_id <- Meta(mygsm)$series_id

  # Convert the data into a dictionary
  dict <- setNames(
    sapply(data, function(x) strsplit(x, ": ")[[1]][2]),  # Extract the value
    sapply(data, function(x) strsplit(x, ": ")[[1]][1])   # Extract the key
  )

  # Convert the dictionary to a row
  new_row <- as.data.frame(t(dict), stringsAsFactors = FALSE)

  # Add the new row to the dataframe
  if (series_id[1] == "GSE101929")
    # mRNA data
    df <- rbind(df, new_row)
  else
    # miRNA data
    df2 <- rbind(df2, new_row)
}

# mRNA data
lung <- df[df$race == 'AA' & !duplicated(df$individual), ]
lung$age <- as.numeric(lung$age)
lung$time <- as.numeric(lung$`survival (days)`)
lung$status <- ifelse(lung$`death due to lung cancer (all years)` == "Alive", 0, 1)

# miRNA data
lung2 <- df2[df2$race == 'AA' & !duplicated(df2$patient) & df2$`survival (months)` != 'Unknown', ]
lung2$age <- as.numeric(lung2$age)
lung2$time <- as.numeric(lung2$`survival (months)`) * 30
lung2$status <- 0

lung <- rbind(lung[, c(2, 5, 12, 13)], lung2[, c(6, 3, 10, 11)])
rownames(lung) <- NULL
usethis::use_data(lung, overwrite = TRUE)
