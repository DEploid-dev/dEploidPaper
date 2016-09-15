fun.divide.to.seg <- function(hapLength, numSeg = 51){
    return(floor(seq(1, hapLength, length.out =numSeg)))
}
