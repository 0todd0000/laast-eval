

read.data <- function(file.name){
    a         <- read.csv(file.name, header=FALSE)
    a         <- as.matrix( a )
    group     <- a[,1]
    y         <- a[,2:ncol(a)]
    u         <- unique(group)
    y1        <- y[group==u[1] , ]   # first group
    y2        <- y[group==u[2] , ]   # second group
    return( list(y1, y2) )
}