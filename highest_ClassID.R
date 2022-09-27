# Creates a vector for each ASV labeling its highest classification with a unique ASV number
highest_ClassID <- function (x) {
  y <- ifelse(is.na(x$Species),
              ifelse(is.na(x$Genus),
                     ifelse(is.na(x$Family),
                            ifelse(is.na(x$Order),
                                   ifelse(is.na(x$Class),
                                          ifelse(is.na(x$Phylum),
                                                 ifelse(is.na(x$Kingdom),
                                                        paste0("Unknown|",c(1:dim(x)[1])),
                                                        paste0(x$Kingdom, "|k-",c(1:dim(x)[1])) ),
                                                 paste0(x$Phylum, "|p-",c(1:dim(x)[1])) ),
                                          paste0(x$Class, "|c-",c(1:dim(x)[1])) ),
                                   paste0(x$Order, "|o-",c(1:dim(x)[1])) ),
                            paste0(x$Family, "|f-",c(1:dim(x)[1])) ),
                     paste0(x$Genus, "|g-",c(1:dim(x)[1])) ),
              paste0(x$Genus, " ", x$Species, "|s-",c(1:dim(x)[1])) )
  return(y)
}