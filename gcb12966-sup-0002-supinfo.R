#######################################################################################################
#######################################################################################################

zhiv<-function(x,prior="non-uniform",prior.method = "within-pop" ,global="F",format = "default"){
if (format=="aflpdat"){x <- cbind(x[,2],x[,1],x[,-c(1:2)])}
if (format=="default") {x <- x}
x <- droplevels(x)
nsamp<-dim(x)[1]							## number of individuals
nloci<-dim(x)[2]-2						## number of loci
nsamp2<-apply(x[,3:(nloci+2)],2,function(y){
				yy <- c(y[which(y==0)],y[which(y==1)])
				length(yy)})


if (prior.method == "among-pop"){

if(global=="T"){							## if calculating overall allele frequencies
			popnames<-"Global"
			npop<-1
			xx<-x[,3:(nloci+2)]
			nn<-apply(xx,2,function(y){
				yy <- c(y[which(y==0)],y[which(y==1)])
				length(yy)})
			# names(poplength)<-"global"
			nhom1<-apply(xx,2,function(y)sum(y[which(y==1)]))}

if(global=="F"){
			popnames<-names(table(x[,1])[table(x[,1])>0])
			npop<-length(table(x[,1]))
			xx<-split(x[,3:(nloci+2)],x[,1],drop=T)
			poplength<-lapply(xx,function(z){apply(z,2,function(y){
				yy <- c(y[which(y==0)],y[which(y==1)])
				length(yy)})})
			nn<-unlist(poplength)
			nhom1<-lapply(xx,function(z){apply(z,2,function(y)sum(y[which(y==1)]))})
			nhom1 <- unlist(nhom1)
			}
			
if(prior=="non-uniform"){
			
			if(npop==1){mm<-(nn-nhom1)
					Ri<-((nn-nhom1)/nn)		## fraction of null homozygotes
					if(length(unique(Ri))==1)mm[1]<-0.5
					Ri <- mm/nn	## recalculate fraction of null homozygotes
					wtRi<-Ri  *  (1/nloci)
					Rbar<-sum(wtRi)}


			if(npop>1){
					mm<-(nn-nhom1)
					dim(mm)<-c(nloci,npop)
					Ri<-(nn-nhom1)/nn
					dim(Ri)<-c(nloci,npop)
					all_same<-which(apply(Ri,1,function(x)length(unique(x))==1))
					if(length(all_same)!=0)mm[all_same,1]<-mm[all_same,1]+0.5
					dim(nn)<-c(nloci,npop)
					Ri<-mm/nn
					nsamp3<-rep(nsamp2, npop); dim(nsamp3)<-c(nloci,npop)
					wtRi<-Ri  *  nn/nsamp3
					Rbar<-apply(wtRi,1,sum)
					}

			sqRi<-Ri^2
			if(npop>1){
					wtsqRi<-sqRi *  nn/nsamp3
					E.sqRi<-apply(wtsqRi,1,sum)}
			if(npop==1){
					wtsqRi<-sqRi *  1/nloci
					E.sqRi<-sum(wtsqRi)}

			sigR<-E.sqRi-(Rbar^2)

			ahat<-Rbar*(((Rbar*(1-Rbar))/sigR)-1)
			bhat<-(1-Rbar)*(((Rbar*(1-Rbar))/sigR)-1)
}



if (prior=="uniform"){
			ahat<-rep(0.5,nloci)
			bhat<-rep(1,nloci)}

m<-as.vector(nn)-nhom1
ahat<-rep(ahat,npop)
bhat<-rep(bhat,npop)

qhat<-exp(lgamma(m+ahat+0.5)+lgamma(as.vector(nn)+ahat+bhat)-
lgamma(m+ahat)-lgamma(as.vector(nn)+ahat+bhat+0.5))
sqSE<- exp(lgamma(m+ahat+1)+lgamma(as.vector(nn)+ahat+bhat)-
lgamma(m+ahat)-lgamma(as.vector(nn)+ahat+bhat+1))-(qhat^2)

dim(qhat)<-c(nloci,npop)
colnames(qhat)<-paste(popnames,"_q",sep="")
qhat<-t(qhat)
colnames(qhat)<-names(x)[3:(nloci+2)]

dim(sqSE)<-c(nloci,npop)
colnames(sqSE)<-paste(popnames,"_SE",sep="")
sqSE<-t(sqSE)
colnames(sqSE)<-names(x)[3:(nloci+2)]

res<-rbind(qhat,sqSE)
}

#######################################################################################

if (prior.method == "within-pop"){

if(global=="T"){							## if calculating overall allele frequencies
			popnames<-"Global"
			npop<-1
			xx<-x[,3:(nloci+2)]
			nn<-apply(xx,2,function(y){
				yy <- c(y[which(y==0)],y[which(y==1)])
				length(yy)})
			# names(poplength)<-"global"
			nhom1<-apply(xx,2,function(y)sum(y[which(y==1)]))}

if(global=="F"){
			popnames<-names(table(x[,1])[table(x[,1])>0])
			npop<-length(table(x[,1]))
			xx<-split(x[,3:(nloci+2)],x[,1],drop=T)
			poplength<-lapply(xx,function(z){apply(z,2,function(y){
				yy <- c(y[which(y==0)],y[which(y==1)])
				length(yy)})})
			nn<-unlist(poplength)
			nhom1<-lapply(xx,function(z){apply(z,2,function(y)sum(y[which(y==1)]))})
			nhom1 <- unlist(nhom1)}

if(prior=="non-uniform"){
			

			if(npop==1){
					mm<-(nn-nhom1)
					Ri<-((nn-nhom1)/nn)		## fraction of null homozygotes
					if(length(unique(Ri))==1)mm[1]<-0.5
					Ri <- mm/nn
					wtRi<-Ri  *  (1/nloci)
					Rbar<-sum(wtRi)
					}

			if(npop>1){
					mm<-(nn-nhom1)
					dim(mm)<-c(nloci,npop)
					Ri<-mm/nn
					dim(Ri)<-c(nloci,npop)
					all_same<-which(apply(Ri,2,function(x)length(unique(x))==1))
					if(length(all_same)!=0)mm[1,all_same]<-mm[1,all_same]+0.5
					dim(nn)<-c(nloci,npop)
					Ri<-mm/nn
					wtRi<-Ri  *  1/nloci
					dim(wtRi)<- c(nloci,npop)
					Rbar<-apply(wtRi,2,sum)
					dim(Ri)<-c(nloci,npop)}

			sqRi<-Ri^2
			if(npop>1){
					wtsqRi<-sqRi *  1/nloci
					dim(wtsqRi)<-c(nloci,npop)
					E.sqRi<-apply(wtsqRi,2,sum)}
			if(npop==1){
					wtsqRi<-sqRi *  1/nloci
					E.sqRi<-sum(wtsqRi)}

			sigR<-E.sqRi-(Rbar^2)

			ahat<-Rbar*(((Rbar*(1-Rbar))/sigR)-1)
			bhat<-(1-Rbar)*(((Rbar*(1-Rbar))/sigR)-1)
}



if (prior=="uniform"){
			ahat<-rep(0.5,npop)
			bhat<-rep(1,npop)}

m<-as.vector(nn)-unlist(nhom1)
ahat<-rep(ahat,each=nloci)
bhat<-rep(bhat,each=nloci)

qhat<-exp(lgamma(m+ahat+0.5)+lgamma(as.vector(nn)+ahat+bhat)-
lgamma(m+ahat)-lgamma(as.vector(nn)+ahat+bhat+0.5))
sqSE<- exp(lgamma(m+ahat+1)+lgamma(as.vector(nn)+ahat+bhat)-
lgamma(m+ahat)-lgamma(as.vector(nn)+ahat+bhat+1))-(qhat^2)

dim(qhat)<-c(nloci,npop)
colnames(qhat)<-paste(popnames,"_q",sep="")
qhat<-t(qhat)
colnames(qhat)<-names(x)[3:(nloci+2)]

dim(sqSE)<-c(nloci,npop)
colnames(sqSE)<-paste(popnames,"_SE",sep="")
sqSE<-t(sqSE)
colnames(sqSE)<-names(x)[3:(nloci+2)]

res<-rbind(qhat,sqSE)
}


return (res)
}


#######################################################################################################
#######################################################################################################
## function for computing Cockerham and Weir's 1993 beta statistic of population differentiation
#  Cockerham CC, Weir BS (1993) Estimation of gene flow from F-statistics. Evolution 47, 855–863.

betastat<-function(x,format = "default"){
if (format=="aflpdat"){x <- cbind(x[,2],x[,1],x[,-c(1:2)])}
if (format=="default") {x <- x}
x<-droplevels(x)
poplength<-table(x[,1])
nloci<-dim(x)[2]-2
npop<-length(table(x[,1])[table(x[,1])>0])
alfq <- zhiv(x)
alfq<-rbind(alfq[1:npop,],1-alfq[1:npop,])
alfq<-cbind(alfq,"Allele"=rep(0:1,each=npop))
M<-(dim(x)[1])/npop
X<-aggregate(alfq[,1:nloci],by=list(alfq[,nloci+1]),FUN=function(x)sum(x^2))
Y<-aggregate(alfq[,1:nloci],by=list(alfq[,nloci+1]),FUN=function(x)(sum(x))^2)
X<-as.data.frame(sapply(X[2:(nloci+1)],sum))
Y<-as.data.frame(sapply(Y[2:(nloci+1)],sum))
F_0<-((2*M*X)-npop)/(((2*M)-1)*npop)
F_1<-(Y-X)/(npop*(npop-1))
beta<-(F_0-F_1)/(1-F_1)
names(beta)<-"Beta"
return(beta)}


#######################################################################################################
#######################################################################################################
## function for conducting two-sided (or one-sided) permutation test of significance for beta
# Argument x is the genotypes data, as above
# Argument nperm is the number of random permutations of the individuals among populations
# Argument locus.summary is the summary statistic used to amalgamate locus-specific estimates for beta. Values can be "mean" or "median"
# Argument type specifies the type of p-value to be output, either "one-tailed" or "two-tailed"

betatest<-function(x,nperm = 1000,locus.summary = "mean",type = "two-tailed",format = "default"){
	if (format=="aflpdat"){x <- cbind(x[,2],x[,1],x[,-c(1:2)])}
	if (format=="default") {x <- x}
	res<-vector("numeric",nperm)
	nind<-dim(x)[1]
	pop_ind<-x[,1:2]
	for(i in 1:(nperm-1)){
		s.ord<-sample(1:nind,nind)
		xx<-x[s.ord,]
		xx[,1:2]<-pop_ind
		betax<-betastat(xx)
		if (locus.summary == "mean") {betaxm<-mean(betax[,1])}
		if (locus.summary == "median") {betaxm<-median(betax[,1])}
		res[i]<-betaxm
		if(i/100 == floor(i/100))cat(paste(i,"\n",sep=""))
		}
	Beta_obs<-betastat(x)
	if (locus.summary == "mean") {Beta_obs <- mean(Beta_obs[,1])}
	if (locus.summary == "median") {Beta_obs <- median(Beta_obs[,1])}
	res[nperm] <- Beta_obs
	
	
	names(res)[nperm]<- "Obs"
	obs.rank <- rep(0,2)
	obs.rank[1] <- rank(res,ties.method="average")["Obs"]
	obs.rank[2] <- nperm-rank(res,ties.method="average")["Obs"]+1
	if(type == "two-tailed") {pval <- (min(obs.rank)*2)/nperm}
	if(type == "one-tailed") {pval <- min(obs.rank)/nperm}
	
	cat(paste("\n","\n","Observed beta statistic:","\n",Beta_obs,"\n"))
	cat(paste("\n",type,"p-value:","\n",pval,"\n","\n"))
	cat(paste("0.0005, 0.005, 0.025, 0.975, 0.995 and 0.9995 quantiles of the null distribution for beta:","\n"))
	cat(quantile(res,probs=c(0.0005,0.005,0.025,0.975,0.995,0.9995)))
		return(res)
	}
	
	
# test <- betatest(Fo.trt,100,locus.summary = "mean",type = "two-tailed")

#######################################################################################################
#######################################################################################################
# Function for creating pairwise distance matrix using Cockerham and Weir's 1993 beta statistic of population differentiation

beta.dist.matrix <- function(x,format = "default"){
if (format=="aflpdat"){x <- cbind(x[,2],x[,1],x[,-c(1:2)])}
if (format=="default") {x <- x}
yy <- split(x,x[,1],drop=T)

res <- matrix(NA,length(yy),length(yy))
dimnames(res)[1] <- list(names(yy))
dimnames(res)[2] <- list(names(yy))

for(i in 1:length(yy)){
	for(j in 1:length(yy)){
		xx <- rbind(yy[[i]],yy[[j]])
		res[i,j]<- mean(betastat(xx)[,1])	
	}
}
return (res)
}

#######################################################################################################
#######################################################################################################

FST <- function(x,diversity.between = F,format = "default"){

if (format=="aflpdat"){x <- cbind(x[,2],x[,1],x[,-c(1:2)])}
if (format=="default") {x <- x}

data <- zhiv(x)
npop<-dim(data)[1]/2
nloc<-dim(data)[2]
popnames<-unlist(strsplit(row.names(data)[1:npop],"_q"))
xx<-split(as.data.frame(data),rep(1:npop,2))
xx<-lapply(xx,as.data.frame)
Hjt<-lapply(xx,function(x)2*x[1,]*(1-x[1,])+ 2*x[2,])
Hjt <- unsplit(Hjt,1:npop)

Hjk <- matrix(NA,npop,npop)
for(i in 2:npop){
for(j in 1:npop){
if (j<i){

Hjk[i,j] <- mean(unlist(data[i,] + data[j,] - 2 * data[i,] * data[j,] - ((Hjt[i,] + Hjt[j,])/2)))
}
}}

nvb <- floor(npop/2)
vb.1 <- var (diag(Hjk[1:nvb*2,1:nvb*2-1]))

vb.2.1 <- rep(1:(nvb/2)*4-3,each=2)+c(0,1)
vb.2.2 <- rep(1:(nvb/2)*4-1,each=2)+c(0,1)
vb.2 <- var(diag(Hjk[vb.2.2,vb.2.1]))
vb <- mean(c(vb.1,vb.2))

Hjk.split <- split(t(Hjk),1:npop)
Hjk.split <- Hjk.split[1:(npop-2)]
Hjk.split <- lapply(Hjk.split,function(x)x[-which(is.na(x))])

cb.res <- vector("numeric",10)

for(i in 1:5){
Hjk.split <- lapply(Hjk.split,function(x)sample(x,length(x)))
cb.1 <- unlist(lapply(Hjk.split,function(x)x[1]))
cb.1 <- rep.int(cb.1,c((npop-2):1))
cb.2 <- unlist(lapply(Hjk.split,function(x)x[2:(length(x))]))
cb <- cbind(cb.1,cb.2)
cb.res[i] <- cov(cb)[2,1]
}

Hjk.split2 <- split(Hjk,1:npop)
Hjk.split2 <- Hjk.split2[3:npop]
Hjk.split2 <- lapply(Hjk.split2,function(x)x[-which(is.na(x))])

for(i in 6:10){
Hjk.split2 <- lapply(Hjk.split2,function(x)sample(x,length(x)))
cb.1 <- unlist(lapply(Hjk.split2,function(x)x[1]))
cb.1 <- rep.int(cb.1,(1:(npop-2)))
cb.2 <- unlist(lapply(Hjk.split2,function(x)x[2:(length(x))]))
cb <- cbind(cb.1,cb.2)
cb.res[i] <- cov(cb)[2,1]
}
cb <- mean(cb.res)

var.hb <- (2*(vb+2*(npop-2)*cb))/(npop*(npop-1))

shjk <- rowSums(Hjk,na.rm=T) + colSums(Hjk,na.rm=T)
Hj <- apply(Hjt,1,mean)



Hjk2 <- as.vector(Hjk)
Hjk2 <- Hjk2[-which(is.na(Hjk2))]

Hb <- (2/(npop*(npop-1)))*sum(Hjk2)

cov.hb.hw <- (1/npop)*( (1/(npop*(npop-1))) * sum (Hj*shjk) - (mean(Hj) * Hb))
var.hw <- (1/(npop*(npop-1))) * sum((Hj-mean(Hj))^2)

Hw <- mean(apply(Hjt,1,mean))
Ht <- Hb + Hw
FST <- (Hb/Ht) * (1 + ((Hb*var.hw - Hw*var.hb + (Hb-Hw)*cov.hb.hw) / (Hb*(Ht^2))))^-1

if(diversity.between == T) {return (list(Hjk,"Hj" = apply(Hjt,1,mean)))} else {return(FST)}
}



#######################################################################################################
#######################################################################################################



Hj<-function(x,format = "default"){
	
	if (format=="aflpdat"){x <- cbind(x[,2],x[,1],x[,-c(1:2)])}
	if (format=="default") {x <- x}

	
	x <- zhiv(x)
	npop<-dim(x)[1]/2
	popnames<-row.names(x)[1:npop]
	popnames<-strsplit(popnames,"_")
	popnames<-unlist(lapply(popnames,function(x)x[1]))
	split.index<-rep(1:npop,2)
	xx<-split(x,split.index)
	nloc<-dim(x)[2]
	xx<-lapply(xx,function(z)as.data.frame(matrix(z,2,nloc)))
	Hjtemp<-lapply(xx,function(z)2*z[1,]*(1-z[1,])+ 2*z[2,])
	VHj<-unlist(lapply(Hjtemp,function(x)var(unlist(x[1,]))/nloc))
	SEHj<-sqrt(VHj)
	Hjt<-unlist(lapply(Hjtemp,function(z) mean(unlist(z))))
	Hj.out<-cbind("Hj"=Hjt,"SE_Hj"=SEHj,"Var_Hj"=VHj)
	row.names(Hj.out)<-popnames
	return(Hj.out)
}

#######################################################################################################
#######################################################################################################


Hj_samp <- function(x,size,nsamp,format = "default"){
	
	if (format=="aflpdat"){x <- cbind(x[,2],x[,1],x[,-c(1:2)])}
	if (format=="default") {x <- x}

	x <- droplevels(x)
	popnames <- names(table(x[,1]))
	npop <- length(unique(x[,1]))
	x <- split(x,x[,1],drop=T)
	res <- matrix(NA,npop,4)
	res[,1]<- popnames
	xx <- lapply(x,function(y){
		res2 <- matrix(NA,nsamp,3)
		nind <- dim(y)[1]
		for(i in 1:nsamp){
			yy <- y[sample(1:nind,size=size),]
			res2[i,] <- Hj(yy)  }
		zz <- apply(res2,2,mean)
		zz
	})
	xx <- unlist (xx)
	xx <- matrix(xx,npop,3,byrow=T)
	res[,2:4] <- xx
	res <- data.frame(res)
	names(res) <- c("Population", "Hj", "SE_Hj","Var_Hj")
	return(res)
}

#######################################################################################################
#######################################################################################################

