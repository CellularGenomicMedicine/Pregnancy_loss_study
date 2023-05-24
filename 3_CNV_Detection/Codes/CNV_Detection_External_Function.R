distcompseg <- function(M1Chr,M2Chr){

Inds <- colnames(M1Chr)[-c(1:3)]

A <- M1Chr[,"Position"]
B <- M2Chr[,"Position"]

Int <- 10

d = rep(0,length(A))

dMatChr <- cbind(M1Chr[,c("Name","Chr","Position")],matrix(NA,nrow(M1Chr),ncol(M1Chr)-3))
colnames(dMatChr) <- colnames(M1Chr)

for(i in 1:length(A)){
#print(i)
        Ph <- abs(B-A[i])
        sPh <- sort(Ph)
        d[i] <- which(Ph==min(Ph))
        
        for(ind in Inds){
				
		if(d[i]<=((Int/2)+1)){
					
			dMatChr[i,ind] <- median(abs(M1Chr[i,ind]-M2Chr[1:(d[i]+(Int/2)),ind]))	
				
		} else if(d[i]>=(nrow(M2Chr)-(Int/2))){
                
	                dMatChr[i,ind] <- median(abs(M1Chr[i,ind]-M2Chr[(d[i]-(Int/2)):nrow(M2Chr),ind]))
				
		} else {
                
	                dMatChr[i,ind] <- median(abs(M1Chr[i,ind]-M2Chr[(d[i]-(Int/2)):(d[i]+(Int/2)),ind]))
				
		}

        }#end ind
        #M2adj[M2adj$Chr==chr,]<-  M2adjChr

}#end i


dMatChr
}#end chr loop

