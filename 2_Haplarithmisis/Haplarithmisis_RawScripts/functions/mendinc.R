#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%   			Function: Mendelian Inconsistancay Rate             	   %%%%%%%%%%%%%#
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
 
 mendinc <- function(Father,Mother,Child){
 	
 	Mdc <- sum((Child=="BB" & (Father=="AB"|Father=="BA") & Mother=="AA") | 
 				 (Child=="AA" & (Father=="AB"|Father=="BA") & Mother=="BB") | 
 	             (Child=="BB" & (Mother=="AB"|Mother=="BA") & Father=="AA") | 
 	             (Child=="AA" & (Mother=="AB"|Mother=="BA") & Father=="BB") | 
 	             (Father=="AA" & (Child=="AA"|Child=="BB")  & Mother=="BB") | 
 	             (Father=="BB" & (Child=="AA"|Child=="BB")  & Mother=="AA") | 
 	             (Child=="AA" & (Father=="NC"|Father=="NoCall") & Mother=="BB") | 
 	             (Child=="BB" & (Father=="NC"|Father=="NoCall") & Mother=="AA") | 
 	             (Child=="AA" & (Mother=="NC"|Mother=="NoCall") & Father=="BB") | 
 	             (Child=="BB" & (Mother=="NC"|Mother=="NoCall") & Father=="AA") | 
 	             (Father=="AA"& Mother=="AA" & Child!="AA" & Child!="NC"&Child!="NoCall") | 
 	             (Father=="BB"& Mother=="BB" & Child!="BB" & Child!="NC"&Child!="NoCall") ) 
 	
 	MdcRate <- (Mdc/sum(Child!="NC"))*100
 	MdcRate
 	
 }#end function