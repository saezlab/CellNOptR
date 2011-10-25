plotCNOlistPDF <-
function(CNOlist,fileName){
	pdf(file=fileName,width=14,height=7)
	plotCNOlist(CNOlist)
	dev.off()
	}

