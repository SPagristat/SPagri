#'Provides RBD analysis and Path analysis
#' @export
#' @param DATA  Experimental data from RBD
GENVA_RBD<- function(DATA){
  N <- nrow(DATA)
  totdf =N-1
  p= ncol(DATA)
  repsum=aggregate(DATA[,3:p], by=list(DATA[,2]), FUN=sum)
  rownames(repsum) <- repsum[,1]
  repsum[,1] <- NULL
  repsum
  R=nrow(repsum)
  repldf=R-1
  trtsum= aggregate(DATA[,3:p], by=list(TRT=DATA[,1]), FUN=sum)
  rownames(trtsum) <- trtsum[,1]
  trtsum[,1] <- NULL
  trtsum
  trtmean= aggregate(DATA[,3:p], by=list(DATA[,1]), FUN=mean)
  rownames(trtmean) <- trtmean[,1]
  trtmean[,1] <- NULL
  trtmean
  T =nrow(trtsum)
  trtdf =T-1
  errordf=totdf-trtdf-repldf
  totsum= colSums(DATA[,3:p])
  grandmean= totsum/N
  A=as.matrix((DATA[,3:p]))
  TS=as.matrix(t(totsum))
  B=as.matrix(trtsum)
  C=as.matrix(repsum)
  Xj= crossprod(TS,TS)/N
  Mj=crossprod(B,B)
  Nj=crossprod(C,C)
  Oj=crossprod(A,A)
  tp =(Mj/R)-Xj
  rp=(Nj/T)-Xj
  totp=Oj-Xj
  ep=totp-tp-rp
  SPmat=rbind(tp,rp,ep,totp)
  trtMsp = tp/trtdf
  rMsp = rp/repldf
  eMsp = ep /errordf
  totMsp =totp/totdf
  MSPmat=rbind(trtMsp,rMsp,eMsp,totMsp)
  trtmsp=diag(trtMsp)
  replmsp= diag(rMsp)
  errormsp=diag(eMsp)
  totmsp= diag(totMsp)
  Ftrt= trtmsp/errormsp
  Frepl= replmsp/errormsp
  F=rbind(Ftrt,Frepl)
  SEm=round((errormsp/R)^(1/2) , digits=4)
  SEd= round((2*errormsp/R)^(1/2) , digits=4)
  ttable=abs(qt(0.025,errordf))
  ftrttable=abs(qf(.95, df1=trtdf, df2=errordf) )
  sig=  ifelse( Ftrt> ftrttable, "*", "NS")
  CD= round((SEd*ttable), digits=4)
  result =  ifelse( sig== "*",CD , "NS")
  Genvar=(trtMsp-eMsp)/R
  Phenvar=Genvar+eMsp
  genvar=diag(Genvar)
  phenvar=diag(Phenvar)
  gv=as.matrix(t(genvar))
  pv=as.matrix(t(phenvar))
  sd=as.matrix(t(errormsp))
  GVj=crossprod(gv,gv)
  SDj=crossprod(sd,sd)
  PVj=crossprod(pv,pv)
  CV= round(((sqrt(errormsp)/grandmean)*100), digits=4)
  PCV= round(((sqrt(phenvar)/grandmean)*100), digits=4)
  GCV= round(((sqrt(genvar)/grandmean)*100), digits=4)
  hertability= round((genvar/phenvar), digits=4)
  Treatmean=round(trtmean,digits=2)
  final=rbind(Treatmean,"===",SEm, SEd, sig, result,"===",CV,PCV,GCV, hertability)
  lab=row.names(Treatmean)
  nams = c(lab,"=====","SEm","SEd", "Significance" ,"CD","===","CV","PCV","GCV", "Hbroad")
  output<- data.frame(final, row.names=nams)
  GenCorrel=Genvar/sqrt(GVj)
  EnvCorrel=eMsp/sqrt(SDj)
  PhenCorrel=Phenvar/sqrt(PVj)
  rg=nrow(GenCorrel)
  g1mat=GenCorrel[1:(rg-1),1:(rg-1)]
  gmat=GenCorrel[1:(rg-1),rg]
  gimat =(solve(g1mat))%*%gmat
  rsqr_g=1-colSums(gimat*gmat)
  GenCorreffect=sweep(g1mat,2,gimat,"*")
  pg=nrow(PhenCorrel)
  p1mat=PhenCorrel[1:(pg-1),1:(pg-1)]
  pmat=PhenCorrel[1:(pg-1),pg]
  pimat =(solve(p1mat))%*%pmat
  rsqr_p=1-colSums(pimat*pmat)
  PhenCorreffect=sweep(p1mat,2,pimat,"*")
    list(RBD = output,
       GenCorrel=GenCorrel,
       EnvCorrel=EnvCorrel,
       PhenCorrel=PhenCorrel,
       PatheffectGenCorr=GenCorreffect,
       ResidualeffectGenCorr=rsqr_g,
       PatheffectPhenCorr=PhenCorreffect,
       ResidualeffectPhenCorr=rsqr_p
  )
}
