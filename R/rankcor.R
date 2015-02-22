# Evaluates the significance of rank correlations for a pair of rankings
rankor<-function(x, y, index = "spearman", approx="exact", tiex="woodbury", CC=FALSE, type="two-sided", print=TRUE, sizer=1000000, repl=1000){
	 if (!is.numeric(x)) {stop("Non-numeric argument to mathematical function")}
	 if (!is.matrix(x) & !is.numeric(y)) {stop("Non-numeric argument to mathematical function")}
	 nx<-length(x);ny<-length(y);n<-max(nx,ny)
	 if (!(nx==ny)) {stop("not all arguments have the same length")}
	 ain<-c("spearman","kendall","gini","r4")
	 apx<-c("gaussian","student","vggfr","exact")
	 alter<- c("two-sided", "less", "greater")
	 tos<-c("woodbury","gh","wgh","midrank","dubois")
	 index<-tolower(index); index<-match.arg(index, ain, several.ok = TRUE) 
	 approx<-tolower(approx);approx<-match.arg(approx, apx, several.ok = TRUE)
	 tiex<-tolower(tiex);tiex<-match.arg(tiex, tos, several.ok = TRUE)
	 type=tolower(type);type<-match.arg(type, alter, several.ok = TRUE)
   	 r<-comprank(x,y,index,tiex,sizer,repl,print)
   	 a<-ranktes(r, n, index, approx, CC, type, print)
	return(a)}