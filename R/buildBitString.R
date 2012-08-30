buildBitString <- function(bStrings){
    # example:
    # buildBitString(list(c(1,1,0,1,0,0), c(0,1,0), c(0,1)))
    #$bs
    #[1] 1 1 0 1 1 1
    #$bsTimes
    #[1] 1 1 0 1 2 3

    # build bitstrings from different steady states
    if (is.list(bStrings)==FALSE){
        if (is.vector(bStrings, mode="numeric")==TRUE){
            bStrings = list(bStrings)
        } else{ 
            stop("CellNOpt Error: argument bString must be a list of bitstring or a single bitstring. Each bitstring must be a vector. See manual for details.")

        }
    }
    N = length(bStrings)

    if (N < 1 ){
        stop("CellNOpt Error: bStrings argument must be a list of length greater or equal to 1")
    }

    if (N>=1){
        bs = bStrings[[1]]
        bsTimes = bStrings[[1]]
    }

    if (N>1){
        for (i in 2:N){
            bs[which(bs==0)] <- bStrings[[i]]
            bsTimes[which(bsTimes==0)]<-bStrings[[i]]*(i)
        }
    }

    return(list(bs=bs, bsTimes=bsTimes))
}
