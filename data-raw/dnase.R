## code to prepare `dnase` dataset goes here

library(survival)

# Extract the first row for each subject
first <- subset(rhDNase, !duplicated(id))
dnase <- tmerge(first, first, id = id, tstop = as.numeric(end.dt - entry.dt))

# Subject whose fu ended during the 6 day window after finishing their
# antibiotics are not counted as new infections
temp.end <- with(rhDNase, pmin(ivstop + 6, end.dt - entry.dt))

# Create the event indicator
dnase <- tmerge(dnase, rhDNase,
  id = id,
  infect = event(ivstart),
  end = event(temp.end)
)

# Toss out the non-at-risk intervals, and extra variables
dnase <- subset(dnase, (infect == 1 | end == 0), c(id:trt, fev:infect))

# Remove duplicated subjects
dnase <- subset(dnase, !duplicated(id))

# Set the survival time
dnase$time <- dnase$tstop - dnase$tstart
dnase <- subset(dnase, select = c(trt, fev, infect, time))

usethis::use_data(dnase, overwrite = TRUE)
