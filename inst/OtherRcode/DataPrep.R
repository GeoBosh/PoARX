#--- Call libraries -----------------------------------------------------------#

library("stringr")
library("XML")
library("RCurl")
library("lubridate")
library("cranlogs")
library("PoARX")

#--- Source code --------------------------------------------------------------#

getwebpage <- function(urlAddress) {
  res <- as.character(getURL(urlAddress))
}

extractPackageDates <- function(pkgName){
  url <- paste0("https://cran.r-project.org/src/contrib/Archive/", pkgName, "/")
  webpage <- getwebpage(url)
  marks <- unlist(gregexpr(paste0(">", pkgName, "_"), webpage))+1
  blurbLength <- nchar("</a></td><td align=\"right\">")
  dateLength <- nchar("2017-11-14")
  dates <- lapply(marks, function(x){
    archiveEnd <- unlist(gregexpr("tar.gz", webpage))
    extEnd <- archiveEnd[archiveEnd > x][1] + nchar("tar.gz")
    st <- extEnd + blurbLength
    en <- st + dateLength - 1
    substr(webpage, st, en)
  })
  dates <- unlist(dates)
  dates <- tryCatch(ymd(dates),
                    warning = function(cond){
                      message(cond)
                      return(dates)
                    })
  return(sort(dates))
}

downloadCounts <- function(pkgName, from, to){
  downloadsData <- cran_downloads(package = pkgName, from = from, to = to)
  releaseDates <- extractPackageDates(pkgName)
  downloadsData <- downloadsData[downloadsData$date >= releaseDates[1],]
  downloadsData$update <- as.numeric(downloadsData$date %in% releaseDates)
  return(downloadsData)
}

mycumsum <- function(x){
  k <- c(which(x == 0), length(x)+1)
  y <- lapply(seq_len(length(k)-1), function(i){
    cumsum(x[k[i]:(k[i+1]-1)])
  })
  y <- unlist(y)
  return(y)
}

#--- Countr data --------------------------------------------------------------#

countrData <- downloadCounts("Countr", from = "2016-03-23", to = "2017-12-12")

### Add 6 'lagged' update variables
countrData$updateLag1 <- c(rep(0,1), countrData$update[-c(600)])
countrData$updateLag2 <- c(rep(0,2), countrData$update[-c(599:600)])
countrData$updateLag3 <- c(rep(0,3), countrData$update[-c(598:600)])
countrData$updateLag4 <- c(rep(0,4), countrData$update[-c(597:600)])
countrData$updateLag5 <- c(rep(0,5), countrData$update[-c(596:600)])
countrData$updateLag6 <- c(rep(0,6), countrData$update[-c(595:600)])

### And create an indicator variable of "update occurred this week"
countrData$recentUpdate <- rowSums(countrData[,grepl("update", names(countrData))])

countrData <- countrData[,c(1:3, 11)]

rownames(countrData) <- NULL

#--- Flexsurv data ------------------------------------------------------------#

flexsurvData <- downloadCounts("flexsurv", from = "2016-03-23", to = "2017-12-12")

### Add 6 'lagged' update variables
flexsurvData$updateLag1 <- c(rep(0,1), flexsurvData$update[-c(600)])
flexsurvData$updateLag2 <- c(rep(0,2), flexsurvData$update[-c(599:600)])
flexsurvData$updateLag3 <- c(rep(0,3), flexsurvData$update[-c(598:600)])
flexsurvData$updateLag4 <- c(rep(0,4), flexsurvData$update[-c(597:600)])
flexsurvData$updateLag5 <- c(rep(0,5), flexsurvData$update[-c(596:600)])
flexsurvData$updateLag6 <- c(rep(0,6), flexsurvData$update[-c(595:600)])

### And create an indicator variable of "update occurred this week"
flexsurvData$recentUpdate <- rowSums(flexsurvData[,grepl("update", names(flexsurvData))])

flexsurvData <- flexsurvData[,c(1:3, 11)]

rownames(flexsurvData) <- NULL

#--- Save data ----------------------------------------------------------------#

save(countrData, file = "CountrData.RData")
save(flexsurvData, file = "flexsurvData.RData")

#--- Join data ----------------------------------------------------------------#

countData <- cbind(countrData[,c(1,2,4)], flexsurvData[,c(2,4)])
names(countData) <- c("date", "CountrDownloads", "CountrUpdate",
                      "flexsurvDownloads", "flexsurvUpdate")
save(countData, file = "countData.RData")

