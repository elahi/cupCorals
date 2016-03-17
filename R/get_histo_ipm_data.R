
# Data from ipmData.csv:
dat <- read.csv("./data/bael_ipmData.csv", header=TRUE, na.strings="NA")

# Select the relevant columns
d <- dat[, which(names(dat) %in% 
                   c("quad", "date", "date.no", "coral.id", "area", 
                     "feret", "code", "sizeOK", "surv", "growth", "recruit"))]

dHIST <- droplevels(d[d$code != "angle" & d$code != "algae" &
                        d$code != "nv" & d$code != "dead", ])
dHIST <- droplevels(dHIST[complete.cases(dHIST$area), ]) # drop NAs

# include only 2007 and 2010: years for IPM
hist07 <- droplevels(dHIST[dHIST$date.no == 39426, ])
hist10 <- droplevels(dHIST[dHIST$date.no == 40515, ])
hist0710 <- rbind(hist07, hist10)