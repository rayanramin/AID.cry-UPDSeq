# NDC2
Normalized Differential Coverage (NDC) 2

This is the new and improved function to get normalized differential coverage of two samples. 

The old version of NDC is available by choosing NDC_Version = 1.  
NDC_Version = 2 is what we have used for our analysis in the AID.cry publication.   
NDC_Version = 3 is similar to NDC2 but more useful for the genome has many regions of zero coverage.

Example, running NDC2 with default window sizes for a circular genome:

NDC2 <- NDC(x=treatment$cov , y= control$cov ,  NDC_Version = 2 , circ = T)

