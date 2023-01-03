#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%        					Function: CallRate         		    		   %%%%%%%%%%%%%#
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(->): 
# 
#
#(<-): 
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
callrate <- function(Gtype){
	
	CallRate <- (sum(Gtype!="NC" )/length(Gtype))*100
	CallRateAA <- (sum(Gtype=="AA")/length(Gtype))*100
	CallRateBB <- (sum(Gtype=="BB")/length(Gtype))*100
	CallRateHet <- (sum(Gtype=="AB")/length(Gtype))*100
	CallRateHom <- (sum(Gtype=="AA"|Gtype=="BB")/length(Gtype))*100
	CallRateNoCall <- (sum(Gtype=="NC")/length(Gtype))*100
	#Test <- CallRateAA + CallRateBB + CallRateHet + CallRateNoCall
	#Test1 <- CallRateHet + CallRateHom + CallRateNoCall
	CallRates <- rbind(CallRate,CallRateAA, CallRateBB,CallRateHet, CallRateHom,CallRateNoCall)
	CallRates

}#end function
