#'Provides RBD analysis and Path analysis
#' @export
#' @param DATA  Experimental data from RBD
#'
Fact<- function(DATA){
  DATA=data.frame(DATA)
  N <- nrow(DATA)
  totdf =N-1
  p= ncol(DATA)
  repsum=aggregate(DATA[,4:p], by=list(DATA[,3]), FUN=sum)
  rownames(repsum) <- repsum[,1]
  repsum[,1] <- NULL
  repsum
  Mpsum= aggregate(DATA[,4:p], by=list(DATA[,1]), FUN=sum)
  rownames(Mpsum) <- Mpsum[,1]
  Mpsum[,1] <- NULL
  Mpsum
  Spsum= aggregate(DATA[,4:p], by=list(DATA[,2]), FUN=sum)
  rownames(Spsum) <- Spsum[,1]
  Spsum[,1] <- NULL
  Spsum
  intsum= aggregate(DATA[,4:p], by=c(list(DATA[,1],DATA[,2])), FUN=sum)
  names <- paste(intsum[,1], intsum[,2], sep=":")
  rownames(intsum) <- names
  intsum[,1:2] <- NULL
  intsum
  totsum=colSums(DATA[,4:p])
  M =nrow(Mpsum)
  Mdf =M-1
  S =nrow(Spsum)
  Sdf =S-1
  R=nrow(repsum)
  repldf=R-1
  intdf= Mdf*Sdf
  errordf=totdf-Mdf-Sdf-intdf-repldf
  A=as.matrix((DATA[,4:p]))
  TS=as.matrix(t(totsum))
  MB=as.matrix(Mpsum)
  SB=as.matrix(Spsum)
  C=as.matrix(repsum)
  D=as.matrix(intsum)
  Xj= crossprod(TS,TS)/N
  Mj=crossprod(MB,MB)
  Sj=crossprod(SB,SB)
  Nj=crossprod(C,C)
  Oj=crossprod(A,A)
  Dj=crossprod(D,D)
  SR=S*R
  MR=M*R
  Mp =(Mj/SR)-Xj
  Sp =(Sj/MR)-Xj
  MS=M*S
  rp=(Nj/MS)-Xj
  intp= (Dj/R)-Xj- Mp-Sp
  totp=Oj-Xj
  ep=totp-rp-Mp-Sp-intp
  rMsp = rp/repldf
  MMsp = Mp/Mdf
  SMsp = Sp/Sdf
  intMSP= intp/intdf
  eMsp = ep /errordf
  totMsp =totp/totdf
  MSPmat=rbind(rMsp,MMsp,SMsp,intMSP,eMsp,totMsp)
  replmsp= diag(rMsp)
  Mmsp=diag(MMsp)
  Smsp=diag(SMsp)
  intmsp=diag(intMSP)
  errormsp=diag(eMsp)
  totmsp= diag(totMsp)
  Frepl= replmsp/errormsp
  FM= Mmsp/errormsp
  FS= Smsp/errormsp
  FMS = intmsp/errormsp
  F=rbind(Frepl,FM,FS,FMS)
  fMtable=abs(qf(.95, df1=Mdf, df2=errordf) )
  fStable=abs(qf(.95, df1=Sdf, df2=errordf) )
  fMStable=abs(qf(.95, df1=intdf, df2=errordf) )
  ttable=abs(qt(0.025,errordf))
  sigM=  ifelse( FM> fMtable, "*", "NS")
  sigS=  ifelse( FS> fStable, "*", "NS")
  sigMS=  ifelse( FMS> fMStable, "*", "NS")
  SEdM= round(sqrt(2*errormsp/(S*R)) , digits=3)
  SEdS= round(sqrt(2*errormsp/(M*R)) , digits=3)
  SEdMS= round(sqrt(2*errormsp/R) , digits=3)
  CDM= round((SEdM*ttable), digits=3)
  CDS= round((SEdS*ttable), digits=3)
  CDMS= round((SEdMS*ttable), digits=3)
  CD_M=ifelse( sigM== "*",CDM , "NS")
  CD_S=ifelse( sigS== "*",CDS , "NS")
  CD_MS=ifelse( sigMS== "*",CDMS , "NS")
  Treatmean=round(intsum/R,digits=3)
  final=rbind(Treatmean,"====",SEdM, SEdS, SEdMS,"______" ,CD_M, CD_S,CD_MS)
  lab=row.names(Treatmean)
  nams = c(lab,"=====","SEdA", "SEdB", "SEdAB", "______","CD_A", "CD_B","CD_AB")
  output<- data.frame(final, row.names=nams)
  intmean1= aggregate(DATA[,4:p], by=c(list(DATA[,1],DATA[,2])), FUN=mean)

  df=data.frame(A=intmean1[,1],B=intmean1[,2], intmean1[,-c(1:2)])

  predictorlist<- as.list(colnames(df[,-c(1:2)]))
  for (i in predictorlist){
    intabl <- matrix(df[,i], nrow =length(unique(df$A)))
    rownames(intabl)<-unique(df$A)
    colnames(intabl)<-unique(df$B)
    intabl1=rbind(intabl,COLmean=colMeans(intabl))
    intabl2=cbind(intabl1,ROWmean=rowMeans(intabl1))
    intabl2<- round(intabl2,3)
    cat("\n Mean table for Factor A x Factor B",
        "\nResponse variable:",i[[1]],"\n")
    print(intabl2)
    cat("\n===========+++============")
  }
  cat( "\n===========Factorial RBD output ============","\n","\n")
  print( output)

}
