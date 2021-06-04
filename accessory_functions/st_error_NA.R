# Standard error without NAs
st_error_NA <- function(x) sd(x,na.rm=TRUE)/sqrt(length(x[!is.na(x)]))
