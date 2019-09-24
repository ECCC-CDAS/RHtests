#**** License agreement ****
#       All users of any code in this package agree to the terms and conditions described in the 
#file License_RHtests__RClimDex_SoftwarePackages.pdf, which is also included in this package.
#*****************************

# last updated on 2016-06-21
# fixed ploting bug for IBC transfered series in FindU(), FindUD() and
#   StepSize(), use insert-back-to-original-date series (odata) instead
#   of non-missing series (Y1 and meanhat)

# last updated on 2016-03-29
# in FindU/FindUD/StepSize, once final Mq==1, set QMadjusted value as missing and bypass 
# correponding plot

# last updated on 2014-02-13
# in FindU(), FindUD, StepSize, changed odata initial value for column 6+, they were NA, changed to 
# original data from input data, this way, those 0 prcp will remain 0 in output.

# last updated on 2013-12-05
# fixed Ids issue in FindUD() function, make copy of _fCs file to _mCs file in StepSize()

# last updated on 2013-07-19
# changed all as.real() to as.numeric() since R2.15 no longer support as.real()

# last updated on 2013-4-30
# changed Mq upper limit from 20 to 100 in adjDLYp() and StartGUI()

# last updated on 2013-02-22
# added StartGUI() function
# modified FindU(), FindUD() and StepSize(), output pdf file even there has no change points

# last updated on 2012-03-30
# set pthr=0.0 as default for FindU, FindUD and StepSize
# change Ti from 365/yr to 365.25/yr in ReadDLY()
# Changed in Fstat output, removed trend estimate for transfered daily PRCP data,
#   estimate trend from monthly total PRCP for both original series and QMadjusted
#   and IBC adjusted series.

# last updated on 2010-06-08
#   In FindU(), FindUD(), StepSize() and 3 corresponding .wRef() functions,
#   added a text indication "yes", "no " or "?  " in changepoint list
#   output file, also changed corresponding read-in part.

# last modified: 2010-05-06
# in QMadjGaussian(), add Nadj option, set empirical prob. segment as length
# Nadj instead of whole segment.
#
# last modified: 2010-01-14
# changed lambda choices from (0.4,0.2,0.1,0.05,0.02,0.01) to (1,0.5,0.2,0.1)
FindU.dlyPrcp<-function(InSeries,output,MissingValueCode,GUI=FALSE,pthr=0.0,
                Mq=10,p.lev=0.95,Iadj=10000,Ny4a=0){
  ErrorMSG<-NA
  assign("ErrorMSG",ErrorMSG,envir=.GlobalEnv)
  Nmin<-10
  Ncat.min<-20
  if(!p.lev%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
      ErrorMSG<<-paste("FindU.dlyPrcp: input p.lev",p.lev,"error\n",
                 get("ErrorMSG",env=.GlobalEnv),"\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
  }
  plev<-p.lev
  pkth<-match(p.lev,c(0.75,0.8,0.9,0.95,0.99,0.9999))
  assign("Nmin",Nmin,envir=.GlobalEnv)
  itmp<-ReadDLY(InSeries,MissingValueCode,pthr)
  if(itmp==(-1)){
    ErrorMSG<<-paste("FindU.dlyPrcp: Error in read data from",InSeries,"\n",
               get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }

  if(Ny4a>0&Ny4a<=5) Ny4a<-5
  Nall<-dim(ori.itable)[1]
  N<-length(Y0); Pfreq<-N/Nall; Nadj<-Ny4a*366*Pfreq
  readPFtable(N,pkth)
# smallestP<-min(Y0)
  smallestP<-min(ori.itable[ori.itable[,"data"]>0,"data"],na.rm=T)
  smallestP<-max(c(smallestP,pthr))
  sumLogPi<-sum(log(Y0))

  ofileSout<-paste(output,"_Ustat.txt",sep="")
  file.create(ofileSout)
  cat(paste("Input data filename:", InSeries,"; N=",N,"; Nday=",Nall,
      "; pthr=",pthr,"\n"), file=ofileSout)
  cat(paste("The nominal level of confidence (1-alpha)=",plev,"\n"),
      file=ofileSout,append=T)
  cat(paste("The smallest non-zero dailyP is:",smallestP,"\n"),
      file=ofileSout,append=T)

  P<-Y0
  aa<-(-1)
  bb<-1
  lambdas<-seq(-1,1,.01)
  LHs<-rep(NA,length(lambdas))
  Nss<-rep(NA,length(lambdas))
  Ipss<-matrix(NA,length(lambdas),100)
  Niter0<-0
  LH.max<-(-9999999.)
  lr.max<-(-1)
  for(dr in c(1,0.5,0.2,0.1)){
    lrs<-seq(aa,bb,by=dr)
  for(lr in 1:length(lrs)){
    lambda<-aa+(lr-1)*dr
    ind<-round((lambda+1)*100+1)
    if(abs(lambdas[ind]-lambda)>1e-8) {
      print(c(lambdas[ind],lambda,lambdas[ind]-lambda))
      stop(paste("ind=",ind,"lambda=",lambda,"error"))
    }
    if(is.na(LHs[ind])==F){  # already done before
      LH<-LHs[ind]
      Ns<-Nss[ind]
      Ips<-Ipss[ind,1:(Ns+1)]
      Y1<-BCtrans(P,lambda)
    }
    else{
      Niter0<-Niter0+1
      Y1<-BCtrans(P,lambda)

      Ip0<-N
      Pk0<-Pk.PMFT(N)

      oout<-PMFT(Y1,Ti,Pk0)
      I0<-0
      I2<-oout$KPx
      I4<-N
      I1<-PMFxKxI0I2(Y1,Ti,I0,I2)$Ic
      I3<-PMFxKxI0I2(Y1,Ti,I2,I4)$Ic
      I2<-PMFxKxI0I2(Y1,Ti,I1,I3)$Ic

      Ns<-1
      Ips<-c(I1,N)
      if(I1>0){
        otmp<-LSmultiple(Y1,Ti,Ips)
        resi<-otmp$resi
        fitted<-otmp$fitted
        otmp<-Rphi(resi,Ips,Ns)
        cor1<-otmp$cor
        corL1<-otmp$corl
        corU1<-otmp$corh
        W<-otmp$W+fitted
        otmp<-PMFxKc(Y1,Ti,I0,I4,I1)
        PFx1<-otmp$PFc
        otmp<-PMFxKc(W,Ti,I0,I4,I1)
        prob1<-otmp$prob
      }
      else{
        prob1<-0
        PFx1<-0
        cor1<-0
        corL1<-0
        corU1<-0
      }

      Ips<-c(I2,N)
      if(I2>0){
        otmp<-LSmultiple(Y1,Ti,Ips)
        resi<-otmp$resi
        fitted<-otmp$fitted
        otmp<-Rphi(resi,Ips,Ns)
        cor2<-otmp$cor
        corL2<-otmp$corl
        corU2<-otmp$corh
        W<-otmp$W+fitted
        otmp<-PMFxKc(Y1,Ti,I0,I4,I2)
        PFx2<-otmp$PFc
        otmp<-PMFxKc(W,Ti,I0,I4,I2)
        prob2<-otmp$prob
      }
      else{
        prob2<-0
        PFx2<-0
        cor2<-0
        corL2<-0
        corU2<-0
      }

      Ips<-c(I3,N)
      if(I3>0){
        otmp<-LSmultiple(Y1,Ti,Ips)
        resi<-otmp$resi
        fitted<-otmp$fitted
        otmp<-Rphi(resi,Ips,Ns)
        cor3<-otmp$cor
        corL3<-otmp$corl
        corU3<-otmp$corh
        W<-otmp$W+fitted
        otmp<-PMFxKc(Y1,Ti,I0,I4,I3)
        PFx3<-otmp$PFc
        otmp<-PMFxKc(W,Ti,I0,I4,I3)
        prob3<-otmp$prob
      }
      else{
        prob3<-0
        PFx3<-0
        cor3<-0
        corL3<-0
        corU3<-0
      }
 
      tmp<-sort(c(PFx1,PFx2,PFx3),decreasing=T,index.return=T)
      PFx.mx<-tmp$x[1]
      prob.mx<-c(prob1,prob2,prob3)[tmp$ix[1]]
      Imx<-c(I1,I2,I3)[tmp$ix[1]]
      cor.mx<-c(corL1,corL2,corL3)[tmp$ix[1]]
      PFx95L<-getPFx95(cor.mx,N)
      if(PFx.mx<PFx95L){
        Ns<-0
        Ips<-N
      }
      else{
        Ns<-1
        Ips<-c(Imx,N)

        tt<-TRUE
        Niter<-0
        while(tt){  # condition on is there any more Bps to insert in Ips?
          Niter<-Niter+1
          tt<-FALSE
          Ips0<-NULL
          for(i in 1:(Ns+1)){ # search for new break points
            I0<- if(i==1) 0 else Ips[i-1]
            I2<-Ips[i]
            otmp<-PMFxKxI0I2(Y1,Ti,I0,I2)
            if(otmp$prob>0) Ips0<-sort(c(Ips0,otmp$Ic))
          }
 # finished search for possible new changepoints, start estimate the p-value
 # of the Ips0 series and find the most significant changepoint Iseg.mx
          tt1<-TRUE
          while(tt1){
            PFx.mx<-(-9999)
            Iseg.mx<-0
            PFx95L.mx<-0
            if(length(Ips0)==0) tt1<-FALSE
            else{
              for(i in 1:length(Ips0)){
	        Ips1<-sort(c(Ips,Ips0[i]))
	        ith<-match(Ips0[i],Ips1)
	        otmp<-PMFxIseg(Y1,Ti,Ips1,ith)
	        if(otmp$PFx<otmp$PFx95L) Ips0[i]<-0
	        else if(otmp$PFx>PFx.mx){
	          PFx.mx<-otmp$PFx
	          Iseg.mx<-Ips0[i]
	          PFx95L.mx<-otmp$PFx95L
	        }
	      }
	      if(PFx.mx>=PFx95L.mx){ # Ic is significant
	        Ips<-sort(c(Ips,Iseg.mx))
	        Ns<-Ns+1
	        Ips0<-Ips0[Ips0!=Iseg.mx]
	        tt<-TRUE
	      }
	      else tt1<-FALSE
	      Ips0<-Ips0[Ips0!=0]
            }
          }
        }
# finish finding any possible new changepoints

# start to delete the least significant changepoint if it is insignificant
        tt<-TRUE
        while(tt){
          tt<-FALSE
          Iseg.mn<-0
          PFx.mn<-99999
          PFx95L.mn<-99999
          for(i in 1:Ns){
            otmp<-PMFxIseg(Y1,Ti,Ips,i)
            if(otmp$PFx<PFx.mn){
              Iseg.mn<-i
	      PFx.mn<-otmp$PFx
	      PFx95L.mn<-otmp$PFx95L
            }
          }
          if(Iseg.mn>0&PFx.mn<PFx95L.mn){
            Ips<-Ips[-Iseg.mn]
            Ns<-Ns-1
            if(Ns>0) tt<-TRUE
          }
        }
      }
      otmp<-LSmultipleRed(Y1,Ti,Ips)
      SSEf<-sum(otmp$resi**2)
      LH<--log(SSEf/(N-2-Ns))*N/2+(lambda-1)*sumLogPi
#     ind<-floor((lambda+1)*100+1)
#     if(lambdas[ind]!=lambda) stop(paste("ind=",ind,"lambda=",lambda,"error"))
#     print(paste("##1##",round(lambda,2),Ns,Niter0,Ns,Ips[1],round(LH,2)))
      LHs[ind]<-LH
      Nss[ind]<-Ns
      if(Ns>0) Ipss[ind,1:(Ns+1)]<-Ips
    }
    if(LH>LH.max){
      LH.max<-LH
      Ips.max<-Ips
      Ns.max<-Ns
      lambda.max<-lambda
    }
  }
# end of detection
  if(dr>0.1){
    aa<-max(c(lambda.max-dr,-1))
    bb<-min(c(lambda.max+dr,1))
    lrs<-seq(aa,bb,by=dr)
  }
  }

# final output
  lambda<-lambda.max
  Ns<-Ns.max
  Ips<-Ips.max

  cat(paste("Best lambda=",round(lambda,2),"; Ns=",Ns,"; likelihood=",
      round(LH.max,2),"\n\n"),file=ofileSout,append=T)
  Y1<-BCtrans(Y0,lambda)

  if(Ns>0) {
    Nsegs<-Ips-c(0,Ips[1:Ns])
    Iseg.longest<-sort(Nsegs,index=T,decreasing=T)$ix[1]
  }
  else Iseg.longest<-0

  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

  Ipd<-c(N)
  dtmp<-LSmultipleRed(Y1,Ti,Ipd)
  beta0<-dtmp$trend
# cat(paste("Ignore changepoints -> trend0 =",round(beta0,6),
#   "(p=", round(dtmp$p.tr,4), "); cor=", round(dtmp$cor,4),
#   "(",round(dtmp$corl,4),",",round(dtmp$corh,4),")\n\n"),
#   file=ofileSout,append=TRUE)

  otmp<-LSmultipleRed(Y1,Ti,Ips)
  cor<-otmp$cor
  corL<-otmp$corl
  corU<-otmp$corh
  W<-otmp$W
  WL<-otmp$WL
  WU<-otmp$WU
  Rresi<-otmp$resi
  pcor<-pt(abs(cor)*sqrt((N-2)/(1-cor^2)),N-2)
  meanhat<-otmp$meanhat
  Pmeanhat<-IVBCtrans(meanhat,lambda)
  meanhatA<-rep(NA,N)

# cat("Common trend TPR fit to the transformed Base series:\n",
#     file=ofileSout,append=T)
# cat(paste("#steps= ",Ns,"; trend=",round(otmp$trend,6),"unit/yr (p=",
#           round(otmp$p.tr,4),"); lambda=",round(lambda,2),"; cor=",
#           round(cor,4),"(",round(corL,4),",",round(corU,4),") p=",
#           round(pcor,4),"\n"),
#           file=ofileSout,append=T)
  
  if(Ns>0) {
    cat(paste("Step-sizes of the transfromed/original series:\n"),file=ofileSout,append=T)
    C<-rep(NA,Ns); sumC<-0; E<-rep(0,Ns+1)
    for(i in 1:Ns){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i+1]
      Ic<-Ips[i]
      meanhat0<-otmp$mu[i]+otmp$trend*Ti[I1:I2]
      Y<-otmp$mu[i+1]+otmp$trend*Ti[I1:I2]
      R0<-IVBCtrans(meanhat0,lambda)
      R<-IVBCtrans(Y,lambda)
      C[i]<-R[Ic-I1+1]-R0[Ic-I1+1]
      sumC<-sumC+C[i]
      E[i+1]<-sumC
      stepsize<-otmp$mu[i+1]-otmp$mu[i]
      cat(paste(Ips[i],IY0[Ips[i]],
                "transfered stepsize=",round(stepsize,4),
                "original stepsize=",round(C[i],4),
		"\n"), file=ofileSout,append=T)
    }

# calculate IBC adjustment    
    Y<-rep(NA,N)
    for(i in 1:(Ns+1)){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      Delta<-otmp$mu[Iseg.adj]-otmp$mu[i]
      Y[I1:I2]<-Y1[I1:I2]+Delta
      meanhatA[I1:I2]<-meanhat[I1:I2]+Delta
    }
#   PmeanhatA<-IVBCtrans(meanhatA,lambda)
    PmeanhatA<-rep(NA,N)
    dPmu<-PmeanhatA[1]-Pmeanhat[1]
    PdA<-rep(NA,N)
    for(i in 1:(Ns+1)){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      if(i==Iseg.adj){
        PdA[I1:I2]<-P[I1:I2]
        DeltaP<-0
	PmeanhatA[I1:I2]<-Pmeanhat[I1:I2]
      }
      else{
        DeltaP<-E[Iseg.adj]-E[i]
	DeltaP0<-DeltaP
        PdA[I1:I2]<-P[I1:I2]+DeltaP
        PmeanhatA[I1:I2]<-Pmeanhat[I1:I2]+DeltaP
	PdA[I1:I2][PdA[I1:I2]<smallestP]<-smallestP
	Diff<-mean(PdA[I1:I2]-P[I1:I2])
	Delta<-DeltaP-Diff
	tflg<-(Delta<=(-0.01))
	Niter<-1
	while(tflg){
	  Niter<-Niter+1
	  Delta0<-Delta
	  PdA[I1:I2]<-PdA[I1:I2]+Delta
	  PmeanhatA[I1:I2]<-PmeanhatA[I1:I2]+Delta
	  PdA[I1:I2][PdA[I1:I2]<smallestP]<-smallestP
	  Diff<-mean(PdA[I1:I2]-P[I1:I2])
	  Delta<-DeltaP-Diff
	  tflg<-(Delta<=(-0.01)&abs(Delta)<abs(Delta0)&Niter<5)
	}
      }
    }
#   PdA[PdA<smallestP]<-smallestP
    Ptr.mx<-max(PmeanhatA)

    Pdtr<-P+Ptr.mx-PmeanhatA
    
    QMout<-adjDLYp(Pdtr,Ips,Mq,Iseg.adj,Ptr.mx,PmeanhatA,Ncat.min,Nadj)
    PA<-QMout$PA
#   PA[PA<smallestP]<-max(c(smallestP,pthr))
    PA[PA<smallestP]<-smallestP
    cat(paste("Nseg_shortest =",QMout$Nseg.mn,"; Mq = ",QMout$Mq,"\n"),
        file=ofileSout,append=T)
    cat(paste("\n Adjust to segment", Iseg.adj,": from", 
        if(Iseg.adj==1) 1 else Ips[Iseg.adj-1]+1,
        "to",Ips[Iseg.adj],"\n"),file=ofileSout,append=T)

    if(QMout$Mq>1){
      cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
          append=T)
      write.table(round(QMout$osmean,4),file=ofileSout,append=T,
                  row.names=F,col.names=F)
      for(i in 1:(Ns+1)){
        I1<-if(i==1) 1 else Ips[i-1]+1
        I2<-Ips[i]
        cat(paste("Seg. ",i,": mean of QM-adjustments =",round(mean(QMout$W[I1:I2]),4),
            "\n",sep=""),file=ofileSout,append=T)
      }
    }

    if(QMout$Mq==1) PA<-rep(NA,length(QMout$PA))
    odata<-cbind(1:nrow(ori.itable),ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3],
#                ori.itable[,4],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
                 ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],
                 ori.itable[,4],ori.itable[,4],ori.itable[,4],0,0)
    Imd<-itable[,1]
    odata[Imd,4]<-round(Pmeanhat,4)
    odata[Imd,5]<-round(PA,4)
    odata[Imd,6]<-round(PdA,4)
    odata[Imd,7]<-round(PmeanhatA,4)
    odata[Imd,8]<-round(Y1,4)
    odata[Imd,9]<-round(meanhat,4)
    odata[Imd,10]<-round(Y,4)
    odata[Imd,11]<-round(meanhatA,4)
    odata[Imd,12]<-round(PA-itable[,5],4)
    odata[Imd,13]<-round(PdA-itable[,5],4)
#   odata[Imd,5]<-round(PA,4)
#   odata[Imd,6]<-round(PdA,4)
#   odata[Imd,7]<-round(PA-ori.itable[Imd,4],4)
#   odata[Imd,8]<-round(PdA-ori.itable[Imd,4],4)
    odataP<-cbind(itable[,1],IY0,itable[,5],Pmeanhat,PA,PdA,PmeanhatA,
                  Y1,meanhat,Y,meanhatA,PA-itable[,5],PdA-itable[,5])
    otrend.ori<-getMtrendFdly(itable[,2:5])
    otrend.IBC<-getMtrendFdly(cbind(itable[,2:4],PdA))
    cat(paste("Linear trend in the original monthly total P is trend0=",
        round(otrend.ori,5),"mm/month\n"),file=ofileSout,append=T)
    if(QMout$Mq>1){
      otrend.QM<-getMtrendFdly(cbind(itable[,2:4],PA))
      cat(paste("Linear trend in QMadjusted monthly total P is trendQM=",
          round(otrend.QM,5),"mm/month\n"),file=ofileSout,append=T)
    }
    cat(paste("Linear trend in mean adjusted(IBC) monthly total P is trend(IBC)=",
        round(otrend.IBC,5),"mm/month\n\n"),file=ofileSout,append=T)
  }
  else{
#   odata<-cbind(ori.itable[,c(1,2,3,4,4,4)],0,0)
    odataP<-cbind(itable[,1],IY0,itable[,5],Pmeanhat,itable[,5],itable[,5],
                  Pmeanhat,Y1,meanhat,Y1,meanhat,0,0)
    odata<-cbind(1:nrow(ori.itable),ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3],
                 ori.itable[,4],ori.itable[,4],0,0,0,0,0,0,0,0)
    otrend.ori<-getMtrendFdly(itable[,2:5])
    cat(paste("Linear trend in the original monthly total P is trend0=",
        round(otrend.ori,5),"mm/month\n"),file=ofileSout,append=T)
  }

  ofilePout<-paste(output,"_U.dat",sep="")
# write.table(round(odataP,4),file=ofilePout,na=MissingValueCode,col.names=F,row.names=F)
  write.table(odata,file=ofilePout,na=MissingValueCode,col.names=F,row.names=F)

  ofilePdf<-paste(output,"_U.pdf",sep="")
  pdf(file=ofilePdf,onefile=T)
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow=c(2,1))
  par(mar=c(2,4,3,1)+.1)
  par(cex=.8,cex.main=.8,cex.lab=.8,cex.axis=.8)
  if(Ns>0){
#   par(mar=c(3,4,3,2)+.1)

    p1data<-cbind(c(1:Nall),ori.itable[,4],NA)
    p1data[Imd,2]<-P
    p1data[Imd,3]<-odataP[,4]

    for(Iseg in 1:Ns){ # plot P~Pmeanhat
      I1<-max(c(1,(itable[Ips[Iseg],1]-180)))
      Ic<-itable[Ips[Iseg],1]
      I2<-min(c(Nall,(itable[Ips[Iseg],1]+180)))
      plot(I1:I2,ori.itable[I1:I2,4],type="l",col="black",ylab="prcp(mm)",
           main=paste("Original dailyP>pthr series\n","Changepoint ",Iseg,
	   " at:",IY0[Ips[Iseg]],sep=""),xaxt="n")
      ats<-c(I1,Ic,I2)
      labels=ori.itable[ats,1]*10000+ori.itable[ats,2]*100+ori.itable[ats,3]
      axis(side=1,at=c(I1,Ic,I2),labels=labels)

      tflg<-is.na(p1data[I1:I2,3])==F
      lines(c(I1:I2)[tflg],p1data[I1:I2,3][tflg],col="red",lwd=1.5)
    }
    for(Iseg in 1:Ns){ # plot Y1~meanhat
#     I1<-max(c(1,Ips[Iseg]-180))
#     Ic<-Ips[Iseg]
#     I2<-min(c(N,Ips[Iseg]+180))
      I1<-max(c(1,(itable[Ips[Iseg],1]-180)))
      Ic<-itable[Ips[Iseg],1]
      I2<-min(c(Nall,(itable[Ips[Iseg],1]+180)))
      plot(I1:I2,odata[I1:I2,8],type="l",col="black",ylab="",xaxt="n",
           main=paste("Box-Cox transformed dailyP>pthr series\n",
	        "Changepoint ",Iseg," at:",IY0[Ips[Iseg]],sep=""))
      axis(side=1,at=c(I1,Ic,I2),labels=IY0[c(I1,Ic,I2)])
      lines(c(I1:I2),odata[I1:I2,9],col="red",lwd=1.5)
    }
  }

  yrs<-unique(itable[,2])
  yrs<-as.integer(seq(yrs[1],yrs[length(yrs)],length=8))
  ymd<-yrs*10000+101
  ats<-rep(NA,length(yrs))
  IY=itable[,2]*10000+itable[,3]*100+itable[,4]
  for(i in 1:length(yrs)){
    it<-match(ymd[i],IY)
    if(!is.na(it)) ats[i]<-itable[it,1]
    else ats[i]<-itable[which.max(IY>ymd[i]),1]
  }
  pdata=rep(NA,nrow(ori.itable))
  pdata[itable[,1]]<-P
  plot(1:nrow(ori.itable),pdata,type="l",col="black",ylab="prcp(mm)",
       main="Original dailyP>pthr series",xaxt="n")
  axis(side=1,at=ats,labels=yrs)
# lines(1:nrow(ori.itable),odataP[,4],col="red")

  if(Ns>0) if(QMout$Mq>1){
    plot(1:nrow(ori.itable),odata[,5],type="l",col="black",ylab="prcp(mm)",
         main="QM adjusted dailyP>pthr series",xaxt="n")
    axis(side=1,at=ats,labels=yrs)
    plot(1:nrow(ori.itable),odata[,6],type="l",col="black",ylab="prcp(mm)",
         main="IBC adjusted dailyP>pthr series",xaxt="n")
    axis(side=1,at=ats,labels=yrs)
  }

  # test plot
  if(Ns>0) if(QMout$Mq>1){
    par(mar=c(4,5,3,2)+.1,cex=.8,mfrow=c(1,1))
    col=0
    np<-0
    osp<-QMout$osp
    osmean<-QMout$osmean
    for(i in 1:(Ns+1)){
      Fd<-.5/QMout$Mq
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      ymax<-max(osp[,2:3],na.rm=T); ymin<-min(osp[,2:3],na.rm=T)
      if(i!=Iseg.adj){
        np<-np+1
        if(col==0) { 
          col<-2
	  plot(osp[I1:I2,2],osp[I1:I2,3],xlim=c(0,1),ylim=c(ymin,ymax),
	       type="l",lwd=1,col=col,xlab="Cumulative Frequency",
	       ylab="QM Adjustment",
	       main=paste("distribution of QM adjustments with Mq=",QMout$Mq))
	  icol<-2*np
	  for(j in 1:QMout$Mq){
	    lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
	          c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
	    if(j>=1&j<QMout$Mq) lines(rep(osmean[(j+1),icol]+Fd,2),
	          c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=.5)
	  }
        }
        else{
          col<-col+1
	  lines(osp[I1:I2,2],osp[I1:I2,3],lwd=1,col=col)
	  icol<-2*np
	  for(j in 1:QMout$Mq){
	    lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
	          c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
	    if(j>=1&j<QMout$Mq) lines(rep(osmean[(j+1),icol]+Fd,2),
	          c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=.5)
	  }
        }
        text(.05,ymax-np*(ymax-ymin)/(Ns*3),paste("Seg.",i))
        lines(c(.15,.20),rep(ymax-np*(ymax-ymin)/(Ns*3),2),lwd=2,col=col)
      }
      else np<-np+1
    }
  }

  par(op)
  dev.off()

  otmp<-LSmultiple(Y1,Ti,Ips)
  resi<-otmp$resi
  otmpW<-LSmultiple(W,Ti,Ips)
  resiW<-otmpW$resi
  otmpWL<-LSmultiple(WL,Ti,Ips)
  resiWL<-otmpWL$resi
  otmpWU<-LSmultiple(WU,Ti,Ips)
  resiWU<-otmpWU$resi

  ofileIout<-paste(output,"_1Cs.txt",sep="")
  ofileMout<-paste(output,"_mCs.txt",sep="")
  file.create(ofileIout)

  if(Ns==0) {
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),file=ofileIout)
    cat("PMF finds no Type-1 changepoints in the series!\n")
  }
  else{
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),
        file=ofileIout)
    for(i in 1:Ns){
      I1<- if(i==1) 1 else Ips[i-1]+1
      Ic<-Ips[i]
      I3<-Ips[i+1]
      Nseg<-I3-I1+1
      PFx95<-getPFx95(cor,Nseg)
      PFx95l<-getPFx95(corL,Nseg)
      PFx95h<-getPFx95(corU,Nseg)
      SSEf.Iseg<-sum(resi[I1:I3]^2)
      Ips1<-Ips[-i]
      otmp1<-LSmultiple(Y1,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      Pk1<-Pk.PMFT(Nseg)
      PFx<-Fx*Pk1[Ic-I1+1]

      SSEf.Iseg<-sum(resiW[I1:I3]^2)
      otmp1<-LSmultiple(W,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(Fx>0) prob<-pf(Fx,1,Nseg-3)
      else{
        Fx<-0
        prob<-0
        PFx<-0
      }

      SSEf.Iseg<-sum(resiWL[I1:I3]^2)
      otmp1<-LSmultiple(WL,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      FxL<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(FxL>0) probL<-pf(Fx,1,Nseg-3)
      else probL<-0

      SSEf.Iseg<-sum(resiWU[I1:I3]^2)
      otmp1<-LSmultiple(WU,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      FxU<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(FxU>0) probU<-pf(Fx,1,Nseg-3)
      else probU<-0

#     else if(Id==1) { # type-1 changepoints
        if(PFx<PFx95l) Idc<-"No  "
        else if(PFx>=PFx95l&PFx<PFx95h) Idc<-"?   "
        else if(PFx>=PFx95h) Idc<-"Yes "
#     }

      cat(paste(sprintf("%1.0f",1)," ",
                sprintf("%-4.4s",Idc),
                sprintf("%10.0f", IY0[Ic])," (",
	        sprintf("%10.4f",probL),"-",
	        sprintf("%10.4f",probU),")",
		sprintf("%6.3f",plev),
	        sprintf("%10.4f",PFx)," (",
	        sprintf("%10.4f",PFx95l),"-",
                sprintf("%10.4f",PFx95h),")\n",sep=""),
		file=ofileIout,
	        append=TRUE)
      cat(paste("PMF : c=", sprintf("%4.0f",Ic), 
      		"; (Time ", sprintf("%10.0f",IY0[Ic]), 
		"); Type= 1; p=",sprintf("%10.4f",prob),"(",
		sprintf("%10.4f",probL),"-",
		sprintf("%10.4f",probU),")",
		"; PFmax=", sprintf("%10.4f",PFx), 
		"; CV95=", sprintf("%10.4f",PFx95), 
		"(", sprintf("%10.4f",PFx95l), 
		"-", sprintf("%10.4f",PFx95h),
		"); Nseg=", sprintf("%4.0f",Nseg),"\n",sep=""), 
		file=ofileSout, append=T)
    }
  }
  if(GUI)
    return(0)
  else{
    file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    cat("FindU.dlyPrcp finished successfully...\n")
  }
}

FindUD.dlyPrcp<-function(InSeries,InCs,output,MissingValueCode,GUI=FALSE,
                pthr=0.0,Mq=10,p.lev=0.95,Iadj=10000,Ny4a=0){
  ErrorMSG<-NA
  assign("ErrorMSG",ErrorMSG,envir=.GlobalEnv)
  Nmin<-10
  if(Ny4a>0&Ny4a<=5) Ny4a<-5
  Ncat.min<-20
  if(!p.lev%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
      ErrorMSG<<-paste("FindU.dlyPrcp: input p.lev",p.lev,"error\n",
                 get("ErrorMSG",env=.GlobalEnv),"\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
  }
  plev<-p.lev
  pkth<-match(p.lev,c(0.75,0.8,0.9,0.95,0.99,0.9999))
  assign("Nmin",Nmin,envir=.GlobalEnv)
  itmp<-ReadDLY(InSeries,MissingValueCode,pthr)
  if(itmp==(-1)){
    ErrorMSG<<-paste("FindUD.dlyprcp: Error in read data from",InSeries,"\n",
               get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }

  Nall<-dim(ori.itable)[1]
  N<-length(Y0); Pfreq<-N/Nall; Nadj<-Ny4a*366*Pfreq
  readPFtable(N,pkth)

  itmp<-readLines(InCs)
  if(length(itmp)<2){
    Ns.ini<-0
    Ips.ini<-N
    Ids.ini<-0
  }
  else{
    Ns.ini<-length(itmp)-1
    Ips.ini<-c(rep(0,Ns.ini),N)
    Ids.ini<-rep(0,Ns.ini)
    for(i in 1:Ns.ini){ # using YYYYMMDD as index, searching for the largest
                    # date less or equal to given YYYYMMDD
      ymdtmp<-as.numeric(substr(itmp[i+1],7,16))
      it<-match(ymdtmp,IY0)
      if(!is.na(it)) Ips.ini[i]<-it
      else Ips.ini[i]<-max(c(1:N)[IY0<=ymdtmp])
      Ids.ini[i]<-as.numeric(substr(itmp[i+1],1,1))
    }
    if(sum(is.na(Ips.ini))>0|!identical(Ips.ini,sort(Ips.ini))){
      ErrorMSG<<-paste("FindUD.dlyprcp: Ips read in from ",InCs,"error:")
      for(i in 1:Ns.ini)
        ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),Ips[i])
      ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),"\n\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
    }
  }

  smallestP<-min(ori.itable[ori.itable[,"data"]>0,"data"],na.rm=T)
  smallestP<-max(c(smallestP,pthr))
  sumLogPi<-sum(log(Y0))

  ofileSout<-paste(output,"_UDstat.txt",sep="")
  file.create(ofileSout)
  cat(paste("Input data filename:", InSeries,"; N=",N,"; Nday=",Nall,
      "; pthr=",pthr,"\n"), file=ofileSout)
  cat(paste("The nominal level of confidence (1-alpha)=",plev,"\n"),
      file=ofileSout,append=T)
  cat(paste("The smallest non-zero dailyP is:",smallestP,"\n"),
      file=ofileSout,append=T)

  P<-Y0
  aa<-(-1)
  bb<-1
  lambdas<-seq(-1,1,.01)
  LHs<-rep(NA,length(lambdas))
  Nss<-rep(NA,length(lambdas))
  Ipss<-matrix(NA,length(lambdas),500)
  Niter0<-0
  LH.max<-(-9999999.)
  lr.max<-(-1)
  for(dr in c(1,0.5,0.2,0.1)){
    lrs<-seq(aa,bb,by=dr)
  for(lr in 1:length(lrs)){
    lambda<-aa+(lr-1)*dr
    ind<-round((lambda+1)*100+1)
    if(abs(lambdas[ind]-lambda)>1e-8) {
      print(c(lambdas[ind],lambda,lambdas[ind]-lambda))
      stop(paste("ind=",ind,"lambda=",lambda,"error"))
    }
    if(is.na(LHs[ind])==F){  # already done before
      LH<-LHs[ind]
      Ns<-Nss[ind]
      Ips<-Ipss[ind,1:(Ns+1)]
      Y1<-BCtrans(P,lambda)
    }
    else{
      tt<-TRUE
      Ns<-Ns.ini
      Ips<-Ips.ini
      Niter0<-Niter0+1
      Y1<-BCtrans(P,lambda)

      if(Ns==0){
        Ip0<-N
        Pk0<-Pk.PMFT(N)

        oout<-PMFT(Y1,Ti,Pk0)
        I0<-0
        I2<-oout$KPx
        I4<-N
        I1<-PMFxKxI0I2(Y1,Ti,I0,I2)$Ic
        I3<-PMFxKxI0I2(Y1,Ti,I2,I4)$Ic
        I2<-PMFxKxI0I2(Y1,Ti,I1,I3)$Ic

        Ns<-1
        Ips<-c(I1,N)
        if(I1>0){
          otmp<-LSmultiple(Y1,Ti,Ips)
          resi<-otmp$resi
          fitted<-otmp$fitted
          otmp<-Rphi(resi,Ips,Ns)
          cor1<-otmp$cor
          corL1<-otmp$corl
          corU1<-otmp$corh
          W<-otmp$W+fitted
          otmp<-PMFxKc(Y1,Ti,I0,I4,I1)
          PFx1<-otmp$PFc
          otmp<-PMFxKc(W,Ti,I0,I4,I1)
          prob1<-otmp$prob
        }
        else{
          prob1<-0
          PFx1<-0
          cor1<-0
          corL1<-0
          corU1<-0
        }

        Ips<-c(I2,N)
        if(I2>0){
          otmp<-LSmultiple(Y1,Ti,Ips)
          resi<-otmp$resi
          fitted<-otmp$fitted
          otmp<-Rphi(resi,Ips,Ns)
          cor2<-otmp$cor
          corL2<-otmp$corl
          corU2<-otmp$corh
          W<-otmp$W+fitted
          otmp<-PMFxKc(Y1,Ti,I0,I4,I2)
          PFx2<-otmp$PFc
          otmp<-PMFxKc(W,Ti,I0,I4,I2)
          prob2<-otmp$prob
        }
        else{
          prob2<-0
          PFx2<-0
          cor2<-0
          corL2<-0
          corU2<-0
        }

        Ips<-c(I3,N)
        if(I3>0){
          otmp<-LSmultiple(Y1,Ti,Ips)
          resi<-otmp$resi
          fitted<-otmp$fitted
          otmp<-Rphi(resi,Ips,Ns)
          cor3<-otmp$cor
          corL3<-otmp$corl
          corU3<-otmp$corh
          W<-otmp$W+fitted
          otmp<-PMFxKc(Y1,Ti,I0,I4,I3)
          PFx3<-otmp$PFc
          otmp<-PMFxKc(W,Ti,I0,I4,I3)
          prob3<-otmp$prob
        }
        else{
          prob3<-0
          PFx3<-0
          cor3<-0
          corL3<-0
          corU3<-0
        }
 
        tmp<-sort(c(PFx1,PFx2,PFx3),decreasing=T,index.return=T)
        PFx.mx<-tmp$x[1]
        prob.mx<-c(prob1,prob2,prob3)[tmp$ix[1]]
        Imx<-c(I1,I2,I3)[tmp$ix[1]]
        cor.mx<-c(corL1,corL2,corL3)[tmp$ix[1]]
        PFx95L<-getPFx95(cor.mx,N)
        if(prob.mx<plev){
          Ns<-0
	  tt<-FALSE
          Ids<-0
        }
        else{
          Ns<-1
          Ips<-c(Imx,N)
	  Ids<-0
        }
      }
      else{  # there are input changepoints
	Ns<-Ns.ini
	Ips<-Ips.ini
        Ids<-Ids.ini
      }
      Ips.i<-Ips
      Niter<-0
      while(tt){  # condition on is there any more Bps to insert in Ips?
        Niter<-Niter+1
        tt<-FALSE
        Ips0<-NULL
        for(i in 1:(Ns+1)){ # search for new break points
          I0<- if(i==1) 0 else Ips[i-1]
          I2<-Ips[i]
          otmp<-PMFxKxI0I2(Y1,Ti,I0,I2)
          if(otmp$prob>0) Ips0<-sort(c(Ips0,otmp$Ic))
        }
 # finished search for possible new changepoints, start estimate the p-value
 # of the Ips0 series and find the most significant changepoint Iseg.mx
        tt1<-TRUE
        while(tt1){
          PFx.mx<-(-9999)
          Iseg.mx<-0
	  prob.mx<-(-99999)
	  probL.mx<-(-99999)
          if(length(Ips0)==0) tt1<-FALSE
          else{
            for(i in 1:length(Ips0)){
	      Ips1<-sort(c(Ips,Ips0[i]))
	      ith<-match(Ips0[i],Ips1)
	      otmp<-PMFxIseg(Y1,Ti,Ips1,ith)
	      probL<-min(c(otmp$probL,otmp$probU,otmp$prob))
	      probU<-max(c(otmp$probL,otmp$probU,otmp$prob))
#	      if(otmp$PFx<otmp$PFx95L) Ips0[i]<-0
	      if(probU<plev) Ips0[i]<-0
	      else if(otmp$PFx>PFx.mx){
	        prob.mx<-otmp$prob
		probL.mx<-probL
	        Iseg.mx<-Ips0[i]
	        PFx.mx<-otmp$PFx
	      }
	    }
	    if(probL.mx>=plev){ # Ic is significant
	      Ips<-sort(c(Ips,Iseg.mx))
	      Ns<-Ns+1
	      Ips0<-Ips0[Ips0!=Iseg.mx]
	      tt<-TRUE
	    }
	    else tt1<-FALSE
	    Ips0<-Ips0[Ips0!=0]
          }
        }
      }
# finish finding any possible new changepoints

      Ids0<-rep(NA,length(Ips))
      for(i in 1:length(Ips)){
        if(Ips[i]%in%Ips.i) Ids0[i]<-Ids[Ips.i==Ips[i]]
	else Ids0[i]<-0
      }
      Ids<-Ids0
# start to delete the least significant changepoint if it is insignificant
      while(tt){
        tt<-FALSE
        Iseg.mn<-0
        probL.mn<-99999
        for(i in 1:Ns){
	  if(Ids[i]==0){ # check those un-documented
            otmp<-PMFxIseg(Y1,Ti,Ips,i)
	    probL<-min(otmp$probL,otmp$probU)
            if(probL<probL.mn){
              Iseg.mn<-i
	      probL.mn<-probL
            }
          }
	}
        if(Iseg.mn>0&probL.mn<plev){
          Ips<-Ips[-Iseg.mn]
	  Ids<-Ids[-Iseg.mn]
          Ns<-Ns-1
          if(Ns>0) tt<-TRUE
        }
      }
    }
# all changepoints are significant, calculate stats and output
    if(length(Ips)==1) if(is.na(Ips)) Ips<-N

    otmp<-LSmultipleRed(Y1,Ti,Ips)
    SSEf<-sum(otmp$resi**2)
    LH<--log(SSEf/(N-2-Ns))*N/2+(lambda-1)*sumLogPi
    LHs[ind]<-LH
    Nss[ind]<-Ns
    if(Ns>0) Ipss[ind,1:(Ns+1)]<-Ips

    if(LH>LH.max){
      LH.max<-LH
      Ips.max<-Ips
      Ids.max<-Ids
      Ns.max<-Ns
      lambda.max<-lambda
    }
#   cat(paste("ind=",ind,"lambda=",lambda,"Ns=",Ns,"likelihood=",round(LH,1),
#             "lambda.max=",lambda.max,"\n"))
  }
# end of detection
  if(dr>0.1){
    aa<-max(c(lambda.max-dr,-1))
    bb<-min(c(lambda.max+dr,1))
    lrs<-seq(aa,bb,by=dr)
  }
  }

# final output
  lambda<-lambda.max
  Ns<-Ns.max
  Ips<-Ips.max
# Ids<-Ids.max
  if(Ns>0){
    Ids0<-rep(NA,length(Ips))
    for(i in 1:length(Ips)){
      if(Ips[i]%in%Ips.ini&(i!=Ns)) Ids0[i]<-Ids.ini[which(Ips.ini==Ips[i])]
      else Ids0[i]<-0
    }
    Ids<-Ids0
  }
  else Ids<-0


  cat(paste("Best lambda=",round(lambda,2),"; Ns=",Ns,"; likelihood=",
      round(LH.max,2),"\n\n"),file=ofileSout,append=T)
  Y1<-BCtrans(Y0,lambda)

  if(Ns>0) {
    Nsegs<-Ips-c(0,Ips[1:Ns])
    Iseg.longest<-sort(Nsegs,index=T,decreasing=T)$ix[1]
  }
  else Iseg.longest<-0

  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

  Ipd<-c(N)
  dtmp<-LSmultipleRed(Y1,Ti,Ipd)
  beta0<-dtmp$trend
  cat(paste("Ignore changepoints -> trend0 =",round(beta0,6),
    "(p=", round(dtmp$p.tr,4), "); cor=", round(dtmp$cor,4),
    "(",round(dtmp$corl,4),",",round(dtmp$corh,4),")\n\n"),
    file=ofileSout,append=TRUE)

  otmp<-LSmultipleRed(Y1,Ti,Ips)
  cor<-otmp$cor
  corL<-otmp$corl
  corU<-otmp$corh
  W<-otmp$W
  WL<-otmp$WL
  WU<-otmp$WU
  Rresi<-otmp$resi
  pcor<-pt(abs(cor)*sqrt((N-2)/(1-cor^2)),N-2)
  meanhat<-otmp$meanhat
  Pmeanhat<-IVBCtrans(meanhat,lambda)
  meanhatA<-rep(NA,N)

# cat("Common trend TPR fit to the transformed Base series:\n",
#     file=ofileSout,append=T)
# cat(paste("#steps= ",Ns,"; trend=",round(otmp$trend,6),"unit/yr (p=",
#           round(otmp$p.tr,4),"); lambda=",round(lambda,2),"; cor=",
#           round(cor,4),"(",round(corL,4),",",round(corU,4),") p=",
#           round(pcor,4),"\n"),
#           file=ofileSout,append=T)
  
  if(Ns>0) {
    cat(paste("Step-sizes of the transfromed/original series:\n"),file=ofileSout,append=T)
    C<-rep(NA,Ns); sumC<-0; E<-rep(0,Ns+1)
    for(i in 1:Ns){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i+1]
      Ic<-Ips[i]
      meanhat0<-otmp$mu[i]+otmp$trend*Ti[I1:I2]
      Y<-otmp$mu[i+1]+otmp$trend*Ti[I1:I2]
      R0<-IVBCtrans(meanhat0,lambda)
      R<-IVBCtrans(Y,lambda)
      C[i]<-R[Ic-I1+1]-R0[Ic-I1+1]
      sumC<-sumC+C[i]
      E[i+1]<-sumC
      stepsize<-otmp$mu[i+1]-otmp$mu[i]
      cat(paste(Ips[i],IY0[Ips[i]],
                "transfered stepsize=",round(stepsize,4),
                "original stepsize=",round(C[i],4),
		"\n"), file=ofileSout,append=T)
    }

# calculate IBC adjustment    
    Y<-rep(NA,N)
    for(i in 1:(Ns+1)){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      Delta<-otmp$mu[Iseg.adj]-otmp$mu[i]
      Y[I1:I2]<-Y1[I1:I2]+Delta
      meanhatA[I1:I2]<-meanhat[I1:I2]+Delta
    }
#   PmeanhatA<-IVBCtrans(meanhatA,lambda)
    PmeanhatA<-rep(NA,N)
    dPmu<-PmeanhatA[1]-Pmeanhat[1]
    PdA<-rep(NA,N)
    for(i in 1:(Ns+1)){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      if(i==Iseg.adj){
        PdA[I1:I2]<-P[I1:I2]
        DeltaP<-0
	PmeanhatA[I1:I2]<-Pmeanhat[I1:I2]
      }
      else{
        DeltaP<-E[Iseg.adj]-E[i]
	DeltaP0<-DeltaP
        PdA[I1:I2]<-P[I1:I2]+DeltaP
        PmeanhatA[I1:I2]<-Pmeanhat[I1:I2]+DeltaP
	PdA[I1:I2][PdA[I1:I2]<smallestP]<-smallestP
	Diff<-mean(PdA[I1:I2]-P[I1:I2])
	Delta<-DeltaP-Diff
	tflg<-(Delta<=(-0.01))
	Niter<-1
	while(tflg){
	  Niter<-Niter+1
	  Delta0<-Delta
	  PdA[I1:I2]<-PdA[I1:I2]+Delta
	  PmeanhatA[I1:I2]<-PmeanhatA[I1:I2]+Delta
	  PdA[I1:I2][PdA[I1:I2]<smallestP]<-smallestP
	  Diff<-mean(PdA[I1:I2]-P[I1:I2])
	  Delta<-DeltaP-Diff
	  tflg<-(Delta<=(-0.01)&abs(Delta)<abs(Delta0)&Niter<5)
	}
      }
    }
#   PdA[PdA<smallestP]<-smallestP
    Ptr.mx<-max(PmeanhatA)

    Pdtr<-P+Ptr.mx-PmeanhatA

    QMout<-adjDLYp(Pdtr,Ips,Mq,Iseg.adj,Ptr.mx,PmeanhatA,Ncat.min,Nadj)
    PA<-QMout$PA
#   PA[PA<smallestP]<-max(c(smallestP,pthr))
    PA[PA<smallestP]<-smallestP
    cat(paste("Nseg_shortest =",QMout$Nseg.mn,"; Mq = ",QMout$Mq,"\n"),
        file=ofileSout,append=T)
    cat(paste("\n Adjust to segment", Iseg.adj,": from", 
        if(Iseg.adj==1) 1 else Ips[Iseg.adj-1]+1,
        "to",Ips[Iseg.adj],"\n"),file=ofileSout,append=T)

    if(QMout$Mq>1){
      cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
          append=T)
      write.table(round(QMout$osmean,4),file=ofileSout,append=T,
                  row.names=F,col.names=F)
      for(i in 1:(Ns+1)){
        I1<-if(i==1) 1 else Ips[i-1]+1
        I2<-Ips[i]
        cat(paste("Seg. ",i,": mean of QM-adjustments =",round(mean(QMout$W[I1:I2]),4),
            "\n",sep=""),file=ofileSout,append=T)
      }
    }
    if(QMout$Mq==1) PA<-rep(NA,length(QMout$PA))

    odata<-cbind(1:nrow(ori.itable),ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3],
#                ori.itable[,4],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
                 ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],
                 ori.itable[,4],ori.itable[,4],ori.itable[,4],0,0)
    Imd<-itable[,1]
    odata[Imd,4]<-round(Pmeanhat,4)
    odata[Imd,5]<-round(PA,4)
    odata[Imd,6]<-round(PdA,4)
    odata[Imd,7]<-round(PmeanhatA,4)
    odata[Imd,8]<-round(Y1,4)
    odata[Imd,9]<-round(meanhat,4)
    odata[Imd,10]<-round(Y,4)
    odata[Imd,11]<-round(meanhatA,4)
    odata[Imd,12]<-round(PA-itable[,5],4)
    odata[Imd,13]<-round(PdA-itable[,5],4)
    odataP<-cbind(itable[,1],IY0,itable[,5],Pmeanhat,PA,PdA,PmeanhatA,
                  Y1,meanhat,Y,meanhatA,PA-itable[,5],PdA-itable[,5])
    otrend.ori<-getMtrendFdly(itable[,2:5])
    otrend.IBC<-getMtrendFdly(cbind(itable[,2:4],PdA))
    cat(paste("Linear trend in the original monthly total P is trend0=",
        round(otrend.ori,5),"mm/month\n"),file=ofileSout,append=T)
    if(QMout$Mq>1){
      otrend.QM<-getMtrendFdly(cbind(itable[,2:4],PA))
      cat(paste("Linear trend in QMadjusted monthly total P is trendQM=",
          round(otrend.QM,5),"mm/month\n"),file=ofileSout,append=T)
    }
    cat(paste("Linear trend in mean adjusted(IBC) monthly total P is trend(IBC)=",
        round(otrend.IBC,5),"mm/month\n\n"),file=ofileSout,append=T)

  }
  else{
    odataP<-cbind(itable[,1],IY0,itable[,5],Pmeanhat,itable[,5],itable[,5],
                  Pmeanhat,Y1,meanhat,Y1,meanhat,0,0)
    odata<-cbind(1:nrow(ori.itable),ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3],
                 ori.itable[,4],ori.itable[,4],0,0,0,0,0,0,0,0)
    otrend.ori<-getMtrendFdly(itable[,2:5])
    cat(paste("Linear trend in the original monthly total P is trend0=",
        round(otrend.ori,5),"mm/month\n"),file=ofileSout,append=T)
  }
# ofileDout<-paste(output,"_adj_UD.dat",sep="")
# write.table(odata,file=ofileDout,na=MissingValueCode,col.names=F,row.names=F)

  ofilePout<-paste(output,"_UD.dat",sep="")
# write.table(round(odataP,4),file=ofilePout,na=MissingValueCode,col.names=F,row.names=F)
  write.table(odata,file=ofilePout,na=MissingValueCode,col.names=F,row.names=F)

  ofilePdf<-paste(output,"_UD.pdf",sep="")
  pdf(file=ofilePdf,onefile=T)
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow=c(2,1))
# par(mar=c(3,4,3,2)+.1)
  par(mar=c(2,4,3,1)+.1)
  par(cex=.8,cex.main=.8,cex.lab=.8,cex.axis=.8)

  if(Ns>0){
    p1data<-cbind(c(1:Nall),ori.itable[,4],NA)
    p1data[Imd,2]<-P
    p1data[Imd,3]<-odataP[,4]

    for(Iseg in 1:Ns){ # plot P~Pmeanhat
      I1<-max(c(1,(itable[Ips[Iseg],1]-180)))
      Ic<-itable[Ips[Iseg],1]
      I2<-min(c(Nall,(itable[Ips[Iseg],1]+180)))
      plot(I1:I2,ori.itable[I1:I2,4],type="l",col="black",ylab="prcp(mm)",
           main=paste("Original dailyP>pthr series\n","Changepoint ",Iseg,
	   " at:",IY0[Ips[Iseg]],sep=""),xaxt="n")
      ats<-c(I1,Ic,I2)
      labels=ori.itable[ats,1]*10000+ori.itable[ats,2]*100+ori.itable[ats,3]
      axis(side=1,at=c(I1,Ic,I2),labels=labels)

      tflg<-is.na(p1data[I1:I2,3])==F
      lines(c(I1:I2)[tflg],p1data[I1:I2,3][tflg],col="red",lwd=1.5)
    }
    for(Iseg in 1:Ns){ # plot Y1~meanhat
      I1<-max(c(1,(itable[Ips[Iseg],1]-180)))
      Ic<-itable[Ips[Iseg],1]
      I2<-min(c(Nall,(itable[Ips[Iseg],1]+180)))
#     I1<-max(c(1,Ips[Iseg]-180))
#     Ic<-Ips[Iseg]
#     I2<-min(c(N,Ips[Iseg]+180))
      plot(I1:I2,odata[I1:I2,8],type="l",col="black",ylab="",xaxt="n",
           main=paste("Box-Cox transformed dailyP>pthr series\n",
	        "Changepoint ",Iseg," at:",IY0[Ips[Iseg]],sep=""))
      axis(side=1,at=c(I1,Ic,I2),labels=IY0[c(I1,Ic,I2)])
      lines(c(I1:I2),odata[I1:I2,9],col="red",lwd=1.5)
    }
  }

  yrs<-unique(itable[,2])
  yrs<-as.integer(seq(yrs[1],yrs[length(yrs)],length=8))
  ymd<-yrs*10000+101
  ats<-rep(NA,length(yrs))
  IY=itable[,2]*10000+itable[,3]*100+itable[,4]
  for(i in 1:length(yrs)){
    it<-match(ymd[i],IY)
    if(!is.na(it)) ats[i]<-itable[it,1]
    else ats[i]<-itable[which.max(IY>ymd[i]),1]
  }
  pdata<-rep(NA,nrow(ori.itable))
  pdata[itable[,1]]<-P
  plot(1:nrow(ori.itable),pdata,type="l",col="black",ylab="prcp(mm)",
       main="Original dailyP>pthr series",xaxt="n")
  axis(side=1,at=ats,labels=yrs)
  lines(1:nrow(ori.itable),odata[,4],col="red")
  if(Ns>0) if(QMout$Mq>1){
    plot(1:nrow(ori.itable),odata[,5],type="l",col="black",ylab="prcp(mm)",
         main="QM adjusted dailyP>pthr series",xaxt="n")
    axis(side=1,at=ats,labels=yrs)
    plot(1:nrow(ori.itable),odata[,6],type="l",col="black",ylab="prcp(mm)",
         main="IBC adjusted dailyP>pthr series",xaxt="n")
    axis(side=1,at=ats,labels=yrs)
  }

  # test plot
  if(Ns>0) if(QMout$Mq>1){
    par(mar=c(4,5,3,2)+.1,cex=.8,mfrow=c(1,1))
    col=0
    np<-0
    osp<-QMout$osp
    osmean<-QMout$osmean
    for(i in 1:(Ns+1)){
      Fd<-.5/QMout$Mq
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      ymax<-max(osp[,2:3],na.rm=T); ymin<-min(osp[,2:3],na.rm=T)
      if(i!=Iseg.adj){
        np<-np+1
        if(col==0) { 
          col<-2
	  plot(osp[I1:I2,2],osp[I1:I2,3],xlim=c(0,1),ylim=c(ymin,ymax),
	       type="l",lwd=1,col=col,xlab="Cumulative Frequency",
	       ylab="QM Adjustment",
	       main=paste("distribution of QM adjustments with Mq=",QMout$Mq))
	  icol<-2*np
	  for(j in 1:QMout$Mq){
	    lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
	          c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
	    if(j>=1&j<QMout$Mq) lines(rep(osmean[(j+1),icol]+Fd,2),
	          c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=.5)
	  }
        }
        else{
          col<-col+1
	  lines(osp[I1:I2,2],osp[I1:I2,3],lwd=1,col=col)
	  icol<-2*np
	  for(j in 1:QMout$Mq){
	    lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
	          c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
	    if(j>=1&j<QMout$Mq) lines(rep(osmean[(j+1),icol]+Fd,2),
	          c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=.5)
	  }
        }
        text(.05,ymax-np*(ymax-ymin)/(Ns*3),paste("Seg.",i))
        lines(c(.15,.20),rep(ymax-np*(ymax-ymin)/(Ns*3),2),lwd=2,col=col)
      }
      else np<-np+1
    }
  }
  par(op)
  dev.off()

  otmp<-LSmultiple(Y1,Ti,Ips)
  resi<-otmp$resi
  otmpW<-LSmultiple(W,Ti,Ips)
  resiW<-otmpW$resi
  otmpWL<-LSmultiple(WL,Ti,Ips)
  resiWL<-otmpWL$resi
  otmpWU<-LSmultiple(WU,Ti,Ips)
  resiWU<-otmpWU$resi

  ofileIout<-paste(output,"_pCs.txt",sep="")
  ofileMout<-paste(output,"_mCs.txt",sep="")
  file.create(ofileIout)

  if(Ns==0) {
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),file=ofileIout)
    cat("PMF finds no Type-1 changepoints in the series!\n")
  }
  else{
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),
        file=ofileIout)
    for(i in 1:Ns){
      I1<- if(i==1) 1 else Ips[i-1]+1
      Ic<-Ips[i]
      I3<-Ips[i+1]
      Nseg<-I3-I1+1
      PFx95<-getPFx95(cor,Nseg)
      PFx95l<-getPFx95(corL,Nseg)
      PFx95h<-getPFx95(corU,Nseg)
      SSEf.Iseg<-sum(resi[I1:I3]^2)
      Ips1<-Ips[-i]
      otmp1<-LSmultiple(Y1,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      Pk1<-Pk.PMFT(Nseg)
      PFx<-Fx*Pk1[Ic-I1+1]

      SSEf.Iseg<-sum(resiW[I1:I3]^2)
      otmp1<-LSmultiple(W,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(Fx>0) prob<-pf(Fx,1,Nseg-3)
      else{
        Fx<-0
        prob<-0
        PFx<-0
      }

      SSEf.Iseg<-sum(resiWL[I1:I3]^2)
      otmp1<-LSmultiple(WL,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      FxL<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(FxL>0) probL<-pf(Fx,1,Nseg-3)
      else probL<-0

      SSEf.Iseg<-sum(resiWU[I1:I3]^2)
      otmp1<-LSmultiple(WU,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      FxU<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(FxU>0) probU<-pf(Fx,1,Nseg-3)
      else probU<-0

      if(Ids[i]==0) { # type-0 changepoints
        if(probU<plev) Idc<-"No  "
	else if(probL<plev&probU>=plev) Idc<-"?   "
	else if(probL>=plev) Idc<-"YifD"
	if(PFx>=PFx95h) Idc<-"Yes "
      }
      else if(Ids[i]==1){ # type-1 changepoints
        if(PFx<PFx95l) Idc<-"No  "
        else if(PFx>=PFx95l&PFx<PFx95h) Idc<-"?   "
        else if(PFx>=PFx95h) Idc<-"Yes "
      }

      cat(paste(sprintf("%1.0f",Ids[i])," ",
                sprintf("%-4.4s",Idc),
                sprintf("%10.0f", IY0[Ic])," (",
	        sprintf("%10.4f",probL),"-",
	        sprintf("%10.4f",probU),")",
		sprintf("%6.3f",plev),
	        sprintf("%10.4f",PFx)," (",
	        sprintf("%10.4f",PFx95l),"-",
                sprintf("%10.4f",PFx95h),")\n",sep=""),
		file=ofileIout,
	        append=TRUE)
      cat(paste("PMF : c=", sprintf("%4.0f",Ic), 
      		"; (Time ", sprintf("%10.0f",IY0[Ic]), 
		"); Type= 1; p=",sprintf("%10.4f",prob),"(",
		sprintf("%10.4f",probL),"-",
		sprintf("%10.4f",probU),")",
		"; PFmax=", sprintf("%10.4f",PFx), 
		"; CV95=", sprintf("%10.4f",PFx95), 
		"(", sprintf("%10.4f",PFx95l), 
		"-", sprintf("%10.4f",PFx95h),
		"); Nseg=", sprintf("%4.0f",Nseg),"\n",sep=""), 
		file=ofileSout, append=T)
    }
  }
  if(GUI)
    return(0)
  else{
    file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    cat("FindUD.dlyprcp finished successfully...\n")
  }
}

StepSize.dlyPrcp<-function(InSeries,InCs,output,MissingValueCode,GUI=FALSE,
                  pthr=0,Mq=10,p.lev=0.95,Iadj=10000,Ny4a=0){
  ErrorMSG<-NA
  Ncat.min<-20
  if(Ny4a>0&Ny4a<=5) Ny4a<-5
  assign("ErrorMSG",ErrorMSG,envir=.GlobalEnv)
  if(!p.lev%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
      ErrorMSG<<-paste("StepSize.dlyprcp: input p.lev",p.lev,"error\n",
                 get("ErrorMSG",env=.GlobalEnv),"\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
  }
  plev<-p.lev
  pkth<-match(p.lev,c(0.75,0.8,0.9,0.95,0.99,0.9999))
  itmp<-ReadDLY(InSeries,MissingValueCode,pthr)
  N<-length(Y0)
  Nall<-dim(ori.itable)[1]; Pfreq<-N/Nall; Nadj<-Ny4a*366*Pfreq
  readPFtable(N,pkth)
  if(itmp==(-1)){
    ErrorMSG<<-paste("StepSize.dlyprcp: Error in read data from",InSeries,"\n",
               get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }
  itmp<-readLines(InCs)
  if(length(itmp)<2){
    Ns<-0
    Ips<-N
  }  else{
    Ns<-length(itmp)-1
    Ips<-c(rep(0,Ns),N)
    Ids<-rep(0,Ns)
    for(i in 1:Ns){ # using YYYYMMDD as index, searching for the largest
                    # date less or equal to given YYYYMMDD
      ymdtmp<-as.numeric(substr(itmp[i+1],7,16))
      it<-match(ymdtmp,IY0)
      if(!is.na(it)) Ips[i]<-it
      else Ips[i]<-max(c(1:N)[IY0<=ymdtmp])
      Ids[i]<-as.numeric(substr(itmp[i+1],1,1))
    }
    if(sum(is.na(Ips))>0|!identical(Ips,sort(Ips))){
      ErrorMSG<<-paste("StepSize.dlyprcp: Ips read in from ",InCs,"error:")
      for(i in 1:Ns)
        ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),Ips[i])
      ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),"\n\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
    }
  }
# smallestP<-min(Y0)
  smallestP<-min(ori.itable[ori.itable[,"data"]>0,"data"],na.rm=T)
  smallestP<-max(c(smallestP,pthr))
  sumLogPi<-sum(log(Y0))

  ofileSout<-paste(output,"_Fstat.txt",sep="")
  file.create(ofileSout)

  cat(paste("Input data filename:", InSeries,"; N=",N,"; Nday=",Nall,
      "; pthr=",pthr,"\n"), file=ofileSout)
  cat(paste("The nominal level of confidence (1-alpha)=",plev,"\n"),
      file=ofileSout,append=T)
  cat(paste("The smallest non-zero dailyP is:",smallestP,"\n"),
      file=ofileSout,append=T)

  Nsegs<-Ips-c(0,Ips[1:Ns])
  Iseg.longest<-sort(Nsegs,index=T,decreasing=T)$ix[1]
  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

  P<-Y0
  aa<-(-1)
  bb<-1
  lambdas<-seq(-1,1,.01)
  LHs<-rep(NA,length(lambdas))
  Niter0<-0
  LH.max<-(-9999999.)
  lr.max<-(-1)
  for(dr in c(1,0.5,0.2,0.1)){
    lrs<-seq(aa,bb,by=dr)
    for(lr in 1:length(lrs)){
      lambda<-aa+(lr-1)*dr
      ind<-round((lambda+1)*100+1)
      if(abs(lambdas[ind]-lambda)>1e-8) {
        print(c(lambdas[ind],lambda,lambdas[ind]-lambda))
        stop(paste("ind=",ind,"lambda=",lambda,"error"))
      }
      if(is.na(LHs[ind])==F){  # already done before
        LH<-LHs[ind]
        Y1<-BCtrans(P,lambda)
      }
      else{
        Y1<-BCtrans(P,lambda)
	otmp<-LSmultipleRed(Y1,Ti,Ips)
	SSEf<-sum(otmp$resi**2)
	LH<--log(SSEf*(N-3)/(N-2-Ns))*N/2+(lambda-1)*sumLogPi
        LHs[ind]<-LH
      }
      if(LH>LH.max){
        LH.max<-LH
	lambda.max<-lambda
      }
    }
    if(dr>0.1){
      aa<-max(c(lambda.max-dr,-1))
      bb<-min(c(lambda.max+dr,1))
      lrs<-seq(aa,bb,by=dr)
    }
  }
  lambda<-lambda.max
  Y1<-BCtrans(Y0,lambda)

  Ipd<-c(N)
  dtmp<-LSmultipleRed(Y1,Ti,Ipd)
  beta0<-dtmp$trend
  cat(paste("Ignore changepoints -> trend0 =",round(beta0,6),
    "(p=", round(dtmp$p.tr,4), "); cor=", round(dtmp$cor,4),
    "(",round(dtmp$corl,4),",",round(dtmp$corh,4),")\n\n"),
    file=ofileSout,append=TRUE)
  
  otmp<-LSmultipleRed(Y1,Ti,Ips)
  cor<-otmp$cor
  corL<-otmp$corl
  corU<-otmp$corh
  W<-otmp$W
  WL<-otmp$WL
  WU<-otmp$WU

  Rresi<-otmp$resi
  pcor<-pt(abs(cor)*sqrt((N-2)/(1-cor^2)),N-2)
  meanhat<-otmp$meanhat
  meanhatA<-rep(NA,N)
  Pmeanhat<-IVBCtrans(meanhat,lambda)

# cat("Common trend TPR fit to the transformed Base series:\n",
#     file=ofileSout,append=T)
# cat(paste("#steps= ",Ns,"; trend=",round(otmp$trend,6),"unit/yr (p=",
#           round(otmp$p.tr,4),"); lambda=",round(lambda,2),"; cor=",
#           round(cor,4),"(",round(corL,4),",",round(corU,4),") p=",
#           round(pcor,4),"\n"),file=ofileSout,append=T)
  if(Ns>0) {
    cat(paste("Step-sizes of the transfromed/original series:\n"),file=ofileSout,append=T)
    C<-rep(NA,Ns); sumC<-0; E<-rep(0,Ns+1)
    for(i in 1:Ns){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i+1]
      Ic<-Ips[i]
      meanhat0<-otmp$mu[i]+otmp$trend*Ti[I1:I2]
      Y<-otmp$mu[i+1]+otmp$trend*Ti[I1:I2]
      R0<-IVBCtrans(meanhat0,lambda)
      R<-IVBCtrans(Y,lambda)
      C[i]<-R[Ic-I1+1]-R0[Ic-I1+1]
      sumC<-sumC+C[i]
      E[i+1]<-sumC
      stepsize<-otmp$mu[i+1]-otmp$mu[i]
      cat(paste(Ips[i],IY0[Ips[i]],
                "transfered stepsize=",round(stepsize,4),
		"original stepsize=",round(C[i],4),
		"\n"), file=ofileSout,append=T)
    }
# calculate IBC adjustment
    Y<-rep(NA,N)
    for(i in 1:(Ns+1)){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      Delta<-otmp$mu[Iseg.adj]-otmp$mu[i]
      Y[I1:I2]<-Y1[I1:I2]+Delta
      meanhatA[I1:I2]<-meanhat[I1:I2]+Delta
    }
#   PmeanhatA<-IVBCtrans(meanhatA,lambda)
    PmeanhatA<-rep(NA,N)
    dPmu<-PmeanhatA[1]-Pmeanhat[1]
    PdA<-rep(NA,N)
    for(i in 1:(Ns+1)){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      if(i==Iseg.adj){
        PdA[I1:I2]<-P[I1:I2]
        DeltaP<-0
	PmeanhatA[I1:I2]<-Pmeanhat[I1:I2]
      }
      else{
        DeltaP<-E[Iseg.adj]-E[i]
	DeltaP0<-DeltaP
        PdA[I1:I2]<-P[I1:I2]+DeltaP
        PmeanhatA[I1:I2]<-Pmeanhat[I1:I2]+DeltaP
	PdA[I1:I2][PdA[I1:I2]<smallestP]<-smallestP
	Diff<-mean(PdA[I1:I2]-P[I1:I2])
	Delta<-DeltaP-Diff
	tflg<-(Delta<=(-0.01))
	Niter<-1
	while(tflg){
	  Niter<-Niter+1
	  Delta0<-Delta
	  PdA[I1:I2]<-PdA[I1:I2]+Delta
	  PmeanhatA[I1:I2]<-PmeanhatA[I1:I2]+Delta
	  PdA[I1:I2][PdA[I1:I2]<smallestP]<-smallestP
	  Diff<-mean(PdA[I1:I2]-P[I1:I2])
	  Delta<-DeltaP-Diff
	  tflg<-(Delta<=(-0.01)&abs(Delta)<abs(Delta0)&Niter<5)
	}
      }
    }
#   PdA[PdA<smallestP]<-smallestP
    Ptr.mx<-max(PmeanhatA)

    Pdtr<-P+Ptr.mx-PmeanhatA
    
    QMout<-adjDLYp(Pdtr,Ips,Mq,Iseg.adj,Ptr.mx,PmeanhatA,Ncat.min,Nadj)
#   assign("QMout",QMout,env=.GlobalEnv)
    PA<-QMout$PA
#   PA[PA<smallestP]<-max(c(smallestP,pthr))
    PA[PA<smallestP]<-smallestP
    cat(paste("Nseg_shortest =",QMout$Nseg.mn,"; Mq = ",QMout$Mq,"\n"),
        file=ofileSout,append=T)
    cat(paste("\n Adjust to segment", Iseg.adj,": from",
        if(Iseg.adj==1) 1 else Ips[Iseg.adj-1]+1,
        "to",Ips[Iseg.adj],"\n"),file=ofileSout,append=T)

    if(QMout$Mq>1){
      cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
          append=T)
      write.table(round(QMout$osmean,4),file=ofileSout,append=T,
                  row.names=F,col.names=F)
      for(i in 1:(Ns+1)){
        I1<-if(i==1) 1 else Ips[i-1]+1
        I2<-Ips[i]
        if(i!=Iseg.adj)
        cat(paste("Seg. ",i,": mean of QM-adjustments =",round(mean(QMout$W[I1:I2]),4),
            "\n",sep=""),file=ofileSout,append=T)
      }
    }
    if(QMout$Mq==1) PA<-rep(NA,length(QMout$PA))

    odata<-cbind(1:nrow(ori.itable),ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3],
#                ori.itable[,4],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
                 ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],
                 ori.itable[,4],ori.itable[,4],ori.itable[,4],0,0)
    Imd<-itable[,1]
    odata[Imd,4]<-round(Pmeanhat,4)
    odata[Imd,5]<-round(PA,4)
    odata[Imd,6]<-round(PdA,4)
    odata[Imd,7]<-round(PmeanhatA,4)
    odata[Imd,8]<-round(Y1,4)
    odata[Imd,9]<-round(meanhat,4)
    odata[Imd,10]<-round(Y,4)
    odata[Imd,11]<-round(meanhatA,4)
    odata[Imd,12]<-round(PA-itable[,5],4)
    odata[Imd,13]<-round(PdA-itable[,5],4)
    odataP<-cbind(itable[,1],IY0,itable[,5],Pmeanhat,PA,PdA,PmeanhatA,
                  Y1,meanhat,Y,meanhatA,PA-itable[,5],PdA-itable[,5])
    otrend.ori<-getMtrendFdly(itable[,2:5])
    otrend.QM<-getMtrendFdly(cbind(itable[,2:4],PA))
    otrend.IBC<-getMtrendFdly(cbind(itable[,2:4],PdA))
    cat(paste("Linear trend in the original monthly total P is trend0=",
        round(otrend.ori,5),"mm/month\n"),file=ofileSout,append=T)
    cat(paste("Linear trend in QMadjusted monthly total P is trendQM=",
        round(otrend.QM,5),"mm/month\n"),file=ofileSout,append=T)
    cat(paste("Linear trend in mean adjusted(IBC) monthly total P is trend(IBC)=",
        round(otrend.IBC,5),"mm/month\n\n"),file=ofileSout,append=T)
  }
  else{
    odataP<-cbind(itable[,1],IY0,itable[,5],Pmeanhat,itable[,5],itable[,5],
                  Pmeanhat,Y1,meanhat,Y1,meanhat,0,0)
    odata<-cbind(1:nrow(ori.itable),ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3],
                 ori.itable[,4],ori.itable[,4],0,0,0,0,0,0,0,0)
    otrend.ori<-getMtrendFdly(itable[,2:5])
    cat(paste("Linear trend in the original monthly total P is trend0=",
        round(otrend.ori,5),"mm/month\n"),file=ofileSout,append=T)
  }
# ofileDout<-paste(output,"_adj_F.dat",sep="")
# write.table(odata,file=ofileDout,na=MissingValueCode,col.names=F,row.names=F)

  ofilePout<-paste(output,"_F.dat",sep="")
# write.table(round(odataP,4),file=ofilePout,na=MissingValueCode,col.names=F,row.names=F)
  write.table(odata,file=ofilePout,na=MissingValueCode,col.names=F,row.names=F)

  ofilePdf<-paste(output,"_F.pdf",sep="")
  pdf(file=ofilePdf,onefile=T)
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow=c(2,1))
# par(mar=c(3,4,3,2)+.1)
  par(mar=c(2,4,3,1)+.1)
  par(cex=.8,cex.main=.8,cex.lab=.8,cex.axis=.8)

  if(Ns>0){
    p1data<-cbind(c(1:Nall),ori.itable[,4],NA)
    p1data[Imd,2]<-P
    p1data[Imd,3]<-odataP[,4]

    for(Iseg in 1:Ns){ # plot P~Pmeanhat
      I1<-max(c(1,(itable[Ips[Iseg],1]-180)))
      Ic<-itable[Ips[Iseg],1]
      I2<-min(c(Nall,(itable[Ips[Iseg],1]+180)))
      plot(I1:I2,ori.itable[I1:I2,4],type="l",col="black",ylab="prcp(mm)",
           main=paste("Original dailyP>pthr series\n","Changepoint ",Iseg,
	   " at:",IY0[Ips[Iseg]],sep=""),xaxt="n")
      ats<-c(I1,Ic,I2)
      labels=ori.itable[ats,1]*10000+ori.itable[ats,2]*100+ori.itable[ats,3]
      axis(side=1,at=c(I1,Ic,I2),labels=labels)

      tflg<-is.na(p1data[I1:I2,3])==F
      lines(c(I1:I2)[tflg],p1data[I1:I2,3][tflg],col="red",lwd=1.5)
    }
    for(Iseg in 1:Ns){ # plot Y1~meanhat
      I1<-max(c(1,(itable[Ips[Iseg],1]-180)))
      Ic<-itable[Ips[Iseg],1]
      I2<-min(c(Nall,(itable[Ips[Iseg],1]+180)))
#     I1<-max(c(1,Ips[Iseg]-180))
#     Ic<-Ips[Iseg]
#     I2<-min(c(N,Ips[Iseg]+180))
      plot(I1:I2,odata[I1:I2,8],type="l",col="black",ylab="",xaxt="n",
           main=paste("Box-Cox transformed dailyP>pthr series\n",
	        "Changepoint ",Iseg," at:",IY0[Ips[Iseg]],sep=""))
      axis(side=1,at=c(I1,Ic,I2),labels=IY0[c(I1,Ic,I2)])
#     tflg<-is.na(p2data[I1:I2,3])==F
#     lines(c(I1:I2)[tflg],p2data[I1:I2,3][tflg],col="red",lwd=1.5)
      lines(c(I1:I2),odata[I1:I2,9],col="red",lwd=1.5)
    }
  }

  yrs<-unique(IY0/10000)
  yrs<-as.integer(seq(yrs[1],yrs[length(yrs)],length=8))
  ymd<-yrs*10000+101
  ats<-rep(NA,length(yrs))
  IY=itable[,2]*10000+itable[,3]*100+itable[,4]
  for(i in 1:length(yrs)){
    it<-match(ymd[i],IY)
    if(!is.na(it)) ats[i]<-itable[it,1]
    else ats[i]<-itable[which.max(IY>ymd[i]),1]
  }
  pdata<-rep(NA,nrow(ori.itable))
  pdata[itable[,1]]<-P
  plot(1:nrow(ori.itable),pdata,type="l",col="black",ylab="prcp(mm)",
       main="Original dailyP>pthr series",xaxt="n")
  axis(side=1,at=ats,labels=yrs)
  lines(1:nrow(ori.itable),odata[,4],col="red")
  if(Ns>0) if(QMout$Mq>1){
    plot(1:nrow(ori.itable),odata[,5],type="l",col="black",ylab="prcp(mm)",
         main="QM adjusted dailyP>pthr series",xaxt="n")
    axis(side=1,at=ats,labels=yrs)
    plot(1:nrow(ori.itable),odata[,6],type="l",col="black",ylab="prcp(mm)",
         main="IBC adjusted dailyP>pthr series",xaxt="n")
    axis(side=1,at=ats,labels=yrs)
  }
  # test plot
  if(Ns>0) if(QMout$Mq>1){
    par(mar=c(4,5,3,2)+.1,cex=.8,mfrow=c(1,1))
    col=0
    np<-0
    osp<-QMout$osp
    osmean<-QMout$osmean
    for(i in 1:(Ns+1)){
      Fd<-.5/QMout$Mq
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      ymax<-max(osp[,2:3],na.rm=T); ymin<-min(osp[,2:3],na.rm=T)
      if(i!=Iseg.adj){
        np<-np+1
        if(col==0) { 
          col<-2
	  plot(osp[I1:I2,2],osp[I1:I2,3],xlim=c(0,1),ylim=c(ymin,ymax),
	       type="l",lwd=1,col=col,xlab="Cumulative Frequency",
	       ylab="QM Adjustment",
	       main=paste("distribution of QM adjustments with Mq=",QMout$Mq))
	  icol<-2*np
	  for(j in 1:QMout$Mq){
	    lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
	          c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
	    if(j>=1&j<QMout$Mq) lines(rep(osmean[(j+1),icol]+Fd,2),
	          c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=.5)
	  }
        }
        else{
          col<-col+1
	  lines(osp[I1:I2,2],osp[I1:I2,3],lwd=1,col=col)
	  icol<-2*np
	  for(j in 1:QMout$Mq){
	    lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
	          c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
	    if(j>=1&j<QMout$Mq) lines(rep(osmean[(j+1),icol]+Fd,2),
	          c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=.5)
	  }
        }
        text(.05,ymax-np*(ymax-ymin)/(Ns*3),paste("Seg.",i))
        lines(c(.15,.20),rep(ymax-np*(ymax-ymin)/(Ns*3),2),lwd=2,col=col)
      }
      else np<-np+1
    }
  }
  par(op)
  dev.off()

  otmp<-LSmultiple(Y1,Ti,Ips)
  resi<-otmp$resi
  otmpW<-LSmultiple(W,Ti,Ips)
  resiW<-otmpW$resi
  otmpWL<-LSmultiple(WL,Ti,Ips)
  resiWL<-otmpWL$resi
  otmpWU<-LSmultiple(WU,Ti,Ips)
  resiWU<-otmpWU$resi

  ofileIout<-paste(output,"_fCs.txt",sep="")
  ofileMout<-paste(output,"_mCs.txt",sep="")
  file.create(ofileIout)

  if(Ns==0) {
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),file=ofileIout)
    cat("PMF finds no Type-1 changepoints in the series!\n")
  }
  else{
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),
        file=ofileIout)
    for(i in 1:Ns){
      I1<- if(i==1) 1 else Ips[i-1]+1
      Ic<-Ips[i]
      I3<-Ips[i+1]
      Id<-Ids[i]
      Nseg<-I3-I1+1
      PFx95<-getPFx95(cor,Nseg)
      PFx95l<-getPFx95(corL,Nseg)
      PFx95h<-getPFx95(corU,Nseg)
      SSEf.Iseg<-sum(resi[I1:I3]^2)
      Ips1<-Ips[-i]
      otmp1<-LSmultiple(Y1,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      Pk1<-Pk.PMFT(Nseg)
      PFx<-Fx*Pk1[Ic-I1+1]

      SSEf.Iseg<-sum(resiW[I1:I3]^2)
      otmp1<-LSmultiple(W,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(Fx>0) prob<-pf(Fx,1,Nseg-3)
      else{
        Fx<-0
        prob<-0
        PFx<-0
      }

      SSEf.Iseg<-sum(resiWL[I1:I3]^2)
      otmp1<-LSmultiple(WL,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      FxL<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(FxL>0) probL<-pf(Fx,1,Nseg-3)
      else probL<-0

      SSEf.Iseg<-sum(resiWU[I1:I3]^2)
      otmp1<-LSmultiple(WU,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      FxU<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(FxU>0) probU<-pf(Fx,1,Nseg-3)
      else probU<-0
      if(Id==0) { # type-0 changepoints
        if(probU<plev) Idc<-"No  "
	else if(probL<plev&probU>=plev) Idc<-"?   "
	else if(probL>=plev) Idc<-"YifD"
	if(PFx>=PFx95h) Idc<-"Yes "
      }
      else if(Id==1){ # type-1 changepoints
        if(PFx<PFx95l) Idc<-"No  "
        else if(PFx>=PFx95l&PFx<PFx95h) Idc<-"?   "
        else if(PFx>=PFx95h) Idc<-"Yes "
      }

      cat(paste(sprintf("%1.0f",Id)," ",
                sprintf("%-4.4s",Idc),
                sprintf("%10.0f", IY0[Ic])," (",
	        sprintf("%10.4f",probL),"-",
	        sprintf("%10.4f",probU),")",
		sprintf("%6.3f",plev),
	        sprintf("%10.4f",PFx)," (",
	        sprintf("%10.4f",PFx95l),"-",
                sprintf("%10.4f",PFx95h),")\n",sep=""),
		file=ofileIout,
	        append=TRUE)
      cat(paste("PMF : c=", sprintf("%4.0f",Ic), 
      		"; (Time ", sprintf("%10.0f",IY0[Ic]), 
		"); Type= 1; p=",sprintf("%10.4f",prob),"(",
		sprintf("%10.4f",probL),"-",
		sprintf("%10.4f",probU),")",
		"; PFmax=", sprintf("%10.4f",PFx), 
		"; CV95=", sprintf("%10.4f",PFx95), 
		"(", sprintf("%10.4f",PFx95l), 
		"-", sprintf("%10.4f",PFx95h),
		"); Nseg=", sprintf("%4.0f",Nseg),"\n",sep=""), 
		file=ofileSout, append=T)
    }
  }
  if(GUI) return(0) else {
    file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    cat("StepSize.dlyprcp finished successfully...\n")
  }
}

adjDLYp<-function(P,Ips,Mq,Iseg.adj,Ptr.mx,R,Ncat.min,Nadj){
  Ns<-length(Ips)-1
  N<-length(P)
  Nseg.mn<-N
  for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Nseg<-I2-I1+1
    if(Nseg<Nseg.mn) Nseg.mn<-Nseg
  }
  if(Mq<=0) Mq<-min(floor(Nseg.mn/Ncat.min),100)
    else Mq<-min(floor(Nseg.mn/Ncat.min),Mq)
  if(Mq>100) Mq<-100
  if(Mq<=0) Mq<-1
  Fd<-.5/Mq
  Fcat<-matrix(NA,(Ns+1),(Mq+2))
  F<-matrix(NA,(Ns+1),N)
  EPb<-matrix(NA,(Ns+1),(Mq+2))
  EPa<-matrix(NA,(Ns+1),(Mq+2))
  for(i in 1:(Ns+1)) Fcat[i,]<-seq(0,by=1/Mq,length=(Mq+2))
  for(i in 1:Ns){
    I1<-if(i==1) 0 else Ips[i-1]
    I2<-Ips[i]
    if(Nadj>0) I1<-max(c(I1,I2-Nadj))
    Nseg<-I2-I1
    Y<-P[(I1+1):I2]
    iindx<-sort(Y,index=T)$ix
    irank<-sort(iindx,index=T)$ix
    F[i,1:Nseg]<-irank/Nseg
    EPb[i,]<-0
    for(k in 2:(Mq+1)){
      Mp1<-floor(Nseg*Fcat[i,(k-1)])
      Mp2<-floor(Nseg*Fcat[i,k])
      EPb[i,k]<-mean(Y[iindx[(Mp1+1):Mp2]])
    }
    EPb[i,(Mq+2)]<-EPb[i,(Mq+1)]
    I3<-Ips[i]
    I4<-Ips[i+1]
    if(Nadj>0) I4<-min(c(I4,I3+Nadj))
    if(Nadj>0|i==Ns){
      Nseg<-I4-I3
      Y<-P[(I3+1):I4]
      iindx<-sort(Y,index=T)$ix
      irank<-sort(iindx,index=T)$ix
#     F[i,1:Nseg]<-irank/Nseg
      EPa[i,]<-0
      for(k in 2:(Mq+1)){
        Mp1<-floor(Nseg*Fcat[i,(k-1)])
        Mp2<-floor(Nseg*Fcat[i,k])
	EPa[i,k]<-mean(Y[iindx[(Mp1+1):Mp2]])
      }
      EPa[i,(Mq+2)]<-EPa[i,(Mq+1)]
    }
  }
# if(Mq==1) EPb[,1]<-EPb[,2]
  if(Nadj==0) if(Ns>1) EPa[1:(Ns-1),]<-EPb[2:Ns,]
  Adj<-matrix(0,(Ns+1),(Mq+2))
  for(k in 1:(Mq+2)){
    if(Iseg.adj>1) for(i in 1:(Iseg.adj-1)) Adj[i,k]<-sum(EPa[i:(Iseg.adj-1),k])-sum(EPb[i:(Iseg.adj-1),k])
    if(Iseg.adj<=Ns) for(i in (Iseg.adj+1):(Ns+1)) Adj[i,k]<-sum(EPb[Iseg.adj:(i-1),k])-sum(EPa[Iseg.adj:(i-1),k])
  }

  osmean<-c(1:(Mq+2)) # for plot purpose
  for(Iseg in c(1:(Ns+1)))
    osmean<-cbind(osmean,Fcat[Iseg,]-Fd,Adj[Iseg,])
  # output osmean is a 2*(Ns+1)+1 by Mq+2 matrix

  PA<-P
  W<-rep(NA,N)
  osp<-NULL
  for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Nseg<-I2-I1+1
    Y<-P[I1:I2]
    iindx<-sort(Y,index=T)$ix
    irank<-sort(iindx,index=T)$ix
    F[i,1:Nseg]<-irank/Nseg
    if(i==Iseg.adj) PA[I1:I2]<-P[I1:I2]-Ptr.mx+R[I1:I2] else{
      dx<-Fcat[i,]-Fd
#     fdx<-EP[Iseg.adj,]-EP[i,]
      fdx<-Adj[i,]
      if(Mq==1) fdx[1]<-fdx[2]
      fdx2<-splineN(dx,fdx,2E30,2E30)
      for(j in I1:I2) W[j]<-splintN(dx,fdx,fdx2,F[i,(j-I1+1)])
      PA[I1:I2]<-P[I1:I2]+W[I1:I2]-Ptr.mx+R[I1:I2]
    }

    Rs<-F[i,1:Nseg]
    ors<-sort(Rs,index=T)$ix
    osp<-rbind(osp,cbind(I1:I2,Rs[ors],W[I1:I2][ors]))
  }

  oout<-list()
  oout$PA<-PA
  oout$W<-W
  oout$Nseg.mn<-Nseg.mn
  oout$Mq<-Mq
  oout$osmean<-osmean
  oout$osp<-osp
  return(oout)
}

ReadDLY<-function(idata,MissingValue,pthr){
  if(!file.exists(idata)) {
    ErrorMSG<<-paste("Input datafile",idata,"does not exist!\n")
    return(-1)
  }
  if(is.csv(idata)){
    itmp<-try(read.table(idata,sep=",",header=F,na.strings=MissingValue,
            colClasses=rep("real",4)),silent=T)
    if(inherits(itmp,"try-error")){
      ErrorMSG<<-geterrmessage()
      return(-1)
    }
    else itable<-itmp
  }
  else{
    itmp<-try(read.table(idata,sep="",header=F,na.strings=MissingValue,
            colClasses=rep("real",4)),silent=T)
    if(inherits(itmp,"try-error")){
      ErrorMSG<<-geterrmessage()
      return(-1)
    }
    else itable<-itmp
  }
  if(ncol(itable)!=4){
    ErrorMSG<<-paste(idata,"has",ncol(itable),"columns. The number of columns should be 4\n")
    return(-1)
  }
  colnames(itable)<-c("id1","id2","id3","data")

# keep input base data as ori.itable
  ori.itable<-itable
  Nall<-dim(itable)[1]
  ind<-c(1:Nall)
  otable<-cbind(ind,itable)[is.na(itable[,4])==F&itable[,4]>pthr,]
  Nday<-otable[,1]
  IY0<-otable[,2]*10000+otable[,3]*100+otable[,4]
  IY0flg<-rep(0,length(IY0))
  Y0<-otable[,5]
  Iyr<-otable[,2]
  Imd<-otable[,3]*100+otable[,4]
  Ti<-Nday/365.25
  for(i in 1:(length(IY0)-1)){
    if(Nday[i+1]-Nday[i]==1) IY0flg[i]<-1
  }
  if(sum(IY0flg)<10){  # too few available data for autocorlh 
    ErrorMSG<<-paste("Too many missing values in ", idata, "to estimate autocorrelation\n")
    return(-1)
  }
  itable<-otable
  assign("ori.itable",ori.itable,envir=.GlobalEnv)
  assign("itable",itable,envir=.GlobalEnv)
  assign("Ti",Ti,envir=.GlobalEnv) # Time index for LS fitting
  assign("Y0",Y0,envir=.GlobalEnv) # Data series for Base
  assign("IY0",IY0,envir=.GlobalEnv) # Cycle index for Base
  assign("Imd",Imd,envir=.GlobalEnv) # Cycle index for Base
  assign("IY0flg",IY0flg,envir=.GlobalEnv) # continuous flag for Base
  return(0)
}

is.csv<-function(names){
  nlen<-nchar(names)
  if(substr(names,nlen-2,nlen)=="csv"|substr(names,nlen-2,nlen)=="CSV") return(T) 
  else return(F)
} 

readPFtable<-function(Nx,pkth=4,ifile="PFmax31red_Nmin10_6CVs"){
# read in PTmax table, assign PTmax table as global variable;
# phi -- vector for cor catalog -- as global variable
# itmp<-matrix(scan(ifile,skip=2,quiet=T),ncol=6,byrow=T)
  itmp<-matrix(c(
20,1.05,1.31,2.24,3.29,6.2,18.71,-0.2, 
20,1.13,1.41,2.4,3.53,6.66,20.06,-0.15, 
20,1.21,1.51,2.57,3.79,7.15,21.49,-0.1, 
20,1.3,1.62,2.76,4.07,7.68,23.08,-0.05, 
20,1.39,1.74,2.97,4.37,8.25,24.92,0, 
20,1.44,1.81,3.08,4.53,8.55,25.8,0.025, 
20,1.5,1.87,3.19,4.7,8.87,26.75,0.05, 
20,1.55,1.94,3.31,4.87,9.19,27.68,0.075, 
20,1.61,2.02,3.43,5.05,9.52,28.68,0.1, 
20,1.67,2.09,3.56,5.23,9.87,29.63,0.125, 
20,1.73,2.17,3.69,5.43,10.23,30.68,0.15, 
20,1.8,2.25,3.83,5.63,10.6,31.83,0.175, 
20,1.86,2.33,3.97,5.84,10.99,32.96,0.2, 
20,1.93,2.42,4.12,6.05,11.38,34.06,0.225, 
20,2.01,2.51,4.27,6.28,11.8,35.18,0.25, 
20,2.08,2.61,4.44,6.51,12.23,36.47,0.275, 
20,2.16,2.71,4.6,6.75,12.67,37.74,0.3, 
20,2.24,2.81,4.77,7,13.13,39.13,0.325, 
20,2.33,2.92,4.95,7.26,13.6,40.46,0.35, 
20,2.42,3.03,5.14,7.53,14.07,41.82,0.375, 
20,2.51,3.14,5.33,7.8,14.56,43.13,0.4, 
20,2.6,3.26,5.52,8.08,15.07,44.45,0.425, 
20,2.7,3.38,5.73,8.37,15.58,45.81,0.45, 
20,2.8,3.51,5.93,8.67,16.1,47.2,0.475, 
20,2.91,3.64,6.14,8.97,16.63,48.32,0.5, 
20,3.01,3.77,6.36,9.28,17.17,49.56,0.525, 
20,3.12,3.9,6.58,9.59,17.7,50.95,0.55, 
20,3.23,4.03,6.8,9.9,18.23,52.25,0.575, 
20,3.34,4.17,7.02,10.21,18.75,53.52,0.6, 
20,3.45,4.3,7.24,10.52,19.27,54.61,0.625, 
20,3.55,4.44,7.46,10.83,19.79,55.62,0.65, 
20,3.66,4.57,7.68,11.13,20.26,56.38,0.675, 
20,3.76,4.7,7.88,11.42,20.73,57.15,0.7, 
20,3.86,4.82,8.08,11.69,21.16,57.69,0.725, 
20,3.95,4.93,8.27,11.94,21.57,58.43,0.75, 
20,4.03,5.03,8.43,12.17,21.94,58.96,0.775, 
20,4.1,5.12,8.57,12.37,22.23,59.37,0.8, 
20,4.16,5.19,8.69,12.53,22.48,59.74,0.825, 
20,4.2,5.24,8.78,12.66,22.67,60,0.85, 
20,4.22,5.27,8.83,12.74,22.79,60.04,0.875, 
20,4.22,5.28,8.85,12.77,22.84,60.2,0.9, 
20,4.21,5.26,8.84,12.77,22.85,59.86,0.925, 
20,4.18,5.24,8.81,12.73,22.83,59.86,0.95, 
20,4.16,5.21,8.78,12.69,22.78,60.09,0.975, 
25,2.93,3.32,4.55,5.85,9.18,21.75,-0.2, 
25,3.15,3.57,4.9,6.31,9.92,23.5,-0.15, 
25,3.38,3.83,5.28,6.81,10.71,25.32,-0.1, 
25,3.63,4.12,5.69,7.34,11.57,27.44,-0.05, 
25,3.9,4.43,6.13,7.92,12.51,29.63,0, 
25,4.04,4.6,6.37,8.23,13.01,30.84,0.025, 
25,4.19,4.77,6.62,8.56,13.53,32.13,0.05, 
25,4.35,4.95,6.87,8.9,14.07,33.47,0.075, 
25,4.51,5.14,7.14,9.25,14.64,34.8,0.1, 
25,4.68,5.34,7.42,9.62,15.23,36.23,0.125, 
25,4.86,5.54,7.72,10.01,15.85,37.62,0.15, 
25,5.05,5.76,8.03,10.41,16.5,39.13,0.175, 
25,5.25,5.99,8.35,10.84,17.18,40.72,0.2, 
25,5.45,6.22,8.69,11.28,17.89,42.29,0.225, 
25,5.67,6.47,9.05,11.75,18.63,43.95,0.25, 
25,5.89,6.74,9.42,12.24,19.41,45.87,0.275, 
25,6.13,7.01,9.82,12.75,20.23,47.75,0.3, 
25,6.38,7.3,10.23,13.29,21.08,49.77,0.325, 
25,6.65,7.61,10.66,13.86,21.97,51.76,0.35, 
25,6.92,7.93,11.11,14.45,22.9,53.98,0.375, 
25,7.22,8.26,11.59,15.07,23.86,56.23,0.4, 
25,7.53,8.62,12.09,15.72,24.87,58.32,0.425, 
25,7.85,8.99,12.62,16.4,25.91,60.64,0.45, 
25,8.19,9.38,13.17,17.11,26.99,62.72,0.475, 
25,8.55,9.79,13.74,17.85,28.1,65.16,0.5, 
25,8.93,10.23,14.34,18.61,29.25,67.44,0.525, 
25,9.32,10.67,14.96,19.41,30.42,69.76,0.55, 
25,9.73,11.15,15.61,20.23,31.64,72.09,0.575, 
25,10.16,11.63,16.27,21.07,32.88,74.67,0.6, 
25,10.61,12.14,16.96,21.93,34.13,76.88,0.625, 
25,11.07,12.66,17.66,22.81,35.37,78.97,0.65, 
25,11.54,13.2,18.38,23.69,36.6,80.98,0.675, 
25,12.02,13.74,19.1,24.58,37.82,82.86,0.7, 
25,12.5,14.28,19.82,25.46,39.02,84.77,0.725, 
25,12.98,14.81,20.52,26.3,40.18,86.76,0.75, 
25,13.45,15.33,21.19,27.11,41.26,88.96,0.775, 
25,13.9,15.83,21.83,27.87,42.25,90.19,0.8, 
25,14.31,16.29,22.42,28.56,43.12,91.13,0.825, 
25,14.68,16.69,22.92,29.15,43.85,92.32,0.85, 
25,14.99,17.03,23.34,29.63,44.41,92.79,0.875, 
25,15.23,17.29,23.65,29.98,44.82,93.13,0.9, 
25,15.38,17.46,23.85,30.2,45.02,93.46,0.925, 
25,15.47,17.54,23.93,30.31,45.14,93.64,0.95, 
25,15.49,17.56,23.95,30.31,45.15,93.52,0.975, 
30,3.53,3.93,5.2,6.5,9.75,21.39,-0.2, 
30,3.8,4.24,5.62,7.04,10.56,23.18,-0.15, 
30,4.1,4.58,6.08,7.62,11.45,25.18,-0.1, 
30,4.42,4.95,6.58,8.26,12.42,27.33,-0.05, 
30,4.78,5.34,7.12,8.95,13.48,29.64,0, 
30,4.96,5.56,7.41,9.32,14.06,30.89,0.025, 
30,5.16,5.78,7.72,9.71,14.66,32.17,0.05, 
30,5.36,6.01,8.03,10.12,15.28,33.62,0.075, 
30,5.58,6.25,8.37,10.55,15.94,35.05,0.1, 
30,5.8,6.51,8.72,11,16.63,36.58,0.125, 
30,6.04,6.78,9.09,11.47,17.36,38.15,0.15, 
30,6.29,7.06,9.48,11.97,18.13,39.81,0.175, 
30,6.55,7.36,9.89,12.49,18.93,41.58,0.2, 
30,6.83,7.67,10.32,13.04,19.79,43.41,0.225, 
30,7.12,8,10.78,13.63,20.68,45.29,0.25, 
30,7.43,8.35,11.26,14.24,21.62,47.4,0.275, 
30,7.75,8.72,11.77,14.9,22.62,49.63,0.3, 
30,8.1,9.12,12.32,15.59,23.67,51.91,0.325, 
30,8.46,9.53,12.89,16.32,24.79,54.32,0.35, 
30,8.85,9.97,13.49,17.09,25.97,56.86,0.375, 
30,9.26,10.44,14.14,17.92,27.22,59.59,0.4, 
30,9.7,10.94,14.82,18.79,28.53,62.18,0.425, 
30,10.17,11.47,15.55,19.71,29.92,64.97,0.45, 
30,10.67,12.04,16.33,20.69,31.38,67.9,0.475, 
30,11.2,12.64,17.15,21.73,32.92,70.92,0.5, 
30,11.76,13.29,18.02,22.83,34.55,74.28,0.525, 
30,12.37,13.97,18.95,24,36.25,77.41,0.55, 
30,13.01,14.7,19.94,25.23,38.03,80.94,0.575, 
30,13.7,15.48,20.98,26.52,39.88,84.17,0.6, 
30,14.43,16.3,22.08,27.88,41.8,87.81,0.625, 
30,15.21,17.17,23.24,29.3,43.79,91.4,0.65, 
30,16.03,18.1,24.45,30.78,45.81,95.1,0.675, 
30,16.89,19.06,25.72,32.3,47.89,99.07,0.7, 
30,17.8,20.07,27.01,33.86,49.96,102.38,0.725, 
30,18.73,21.12,28.34,35.44,52.05,105.95,0.75, 
30,19.69,22.18,29.68,37.01,54.11,109.03,0.775, 
30,20.66,23.25,31.02,38.55,56.06,111.99,0.8, 
30,21.63,24.3,32.31,40.03,57.92,114.81,0.825, 
30,22.56,25.32,33.52,41.42,59.58,116.67,0.85, 
30,23.43,26.26,34.63,42.66,61.08,118.21,0.875, 
30,24.2,27.08,35.59,43.71,62.33,120.02,0.9, 
30,24.84,27.76,36.36,44.54,63.21,121.18,0.925, 
30,25.32,28.26,36.9,45.09,63.82,122.14,0.95, 
30,25.6,28.56,37.19,45.39,64.12,122.12,0.975, 
35,3.81,4.22,5.49,6.79,9.94,20.74,-0.2, 
35,4.12,4.57,5.95,7.36,10.8,22.6,-0.15, 
35,4.46,4.94,6.45,7.99,11.74,24.59,-0.1, 
35,4.82,5.35,7,8.68,12.77,26.82,-0.05, 
35,5.22,5.8,7.6,9.43,13.9,29.21,0, 
35,5.43,6.04,7.92,9.83,14.51,30.51,0.025, 
35,5.65,6.28,8.25,10.26,15.15,31.88,0.05, 
35,5.89,6.55,8.6,10.7,15.82,33.28,0.075, 
35,6.13,6.82,8.97,11.17,16.52,34.81,0.1, 
35,6.39,7.11,9.37,11.67,17.27,36.36,0.125, 
35,6.66,7.41,9.78,12.19,18.05,38,0.15, 
35,6.94,7.73,10.21,12.74,18.89,39.8,0.175, 
35,7.24,8.07,10.67,13.32,19.77,41.73,0.2, 
35,7.56,8.43,11.15,13.93,20.7,43.66,0.225, 
35,7.9,8.81,11.67,14.59,21.69,45.74,0.25, 
35,8.25,9.21,12.21,15.28,22.73,47.99,0.275, 
35,8.63,9.64,12.79,16.01,23.83,50.36,0.3, 
35,9.03,10.09,13.41,16.8,25.01,52.9,0.325, 
35,9.46,10.57,14.07,17.63,26.27,55.55,0.35, 
35,9.91,11.09,14.77,18.52,27.6,58.31,0.375, 
35,10.4,11.64,15.52,19.47,29.03,61.33,0.4, 
35,10.92,12.23,16.32,20.48,30.55,64.47,0.425, 
35,11.47,12.86,17.18,21.57,32.17,67.66,0.45, 
35,12.07,13.54,18.1,22.73,33.9,71.09,0.475, 
35,12.72,14.27,19.09,23.98,35.74,74.8,0.5, 
35,13.41,15.06,20.15,25.31,37.7,78.45,0.525, 
35,14.16,15.9,21.3,26.75,39.79,82.7,0.55, 
35,14.97,16.82,22.53,28.28,42.01,87.19,0.575, 
35,15.84,17.8,23.86,29.93,44.38,91.91,0.6, 
35,16.78,18.86,25.28,31.69,46.88,96.59,0.625, 
35,17.8,20.01,26.8,33.56,49.52,101.33,0.65, 
35,18.9,21.24,28.43,35.55,52.3,106.26,0.675, 
35,20.08,22.57,30.17,37.66,55.22,110.93,0.7, 
35,21.34,23.99,32.01,39.87,58.24,115.86,0.725, 
35,22.69,25.49,33.96,42.19,61.34,120.44,0.75, 
35,24.12,27.08,35.98,44.59,64.48,125.51,0.775, 
35,25.61,28.73,38.06,47.03,67.6,130.65,0.8, 
35,27.14,30.43,40.16,49.48,70.67,135.25,0.825, 
35,28.7,32.13,42.25,51.88,73.62,139.28,0.85, 
35,30.25,33.81,44.29,54.16,76.42,143.24,0.875, 
35,31.73,35.41,46.15,56.23,78.86,146.53,0.9, 
35,33.07,36.86,47.8,58.02,80.85,148.94,0.925, 
35,34.19,38.05,49.13,59.41,82.35,150.77,0.95, 
35,34.98,38.87,50,60.3,83.26,151.71,0.975, 
40,3.98,4.39,5.66,6.94,10.02,20.25,-0.2, 
40,4.32,4.77,6.15,7.54,10.91,22.02,-0.15, 
40,4.68,5.17,6.68,8.2,11.88,23.97,-0.1, 
40,5.07,5.61,7.26,8.93,12.94,26.12,-0.05, 
40,5.5,6.09,7.89,9.72,14.1,28.52,0, 
40,5.73,6.34,8.23,10.14,14.73,29.86,0.025, 
40,5.97,6.61,8.59,10.59,15.39,31.17,0.05, 
40,6.22,6.89,8.96,11.06,16.09,32.65,0.075, 
40,6.49,7.19,9.36,11.55,16.82,34.19,0.1, 
40,6.76,7.5,9.77,12.07,17.59,35.79,0.125, 
40,7.06,7.83,10.21,12.62,18.41,37.48,0.15, 
40,7.36,8.17,10.68,13.2,19.29,39.26,0.175, 
40,7.69,8.54,11.17,13.82,20.2,41.15,0.2, 
40,8.04,8.93,11.68,14.47,21.18,43.17,0.225, 
40,8.4,9.34,12.24,15.17,22.23,45.33,0.25, 
40,8.79,9.78,12.83,15.91,23.33,47.6,0.275, 
40,9.2,10.24,13.45,16.69,24.51,49.96,0.3, 
40,9.64,10.74,14.12,17.53,25.77,52.57,0.325, 
40,10.11,11.27,14.83,18.43,27.1,55.38,0.35, 
40,10.61,11.83,15.6,19.4,28.54,58.3,0.375, 
40,11.15,12.44,16.42,20.43,30.07,61.38,0.4, 
40,11.73,13.09,17.3,21.54,31.72,64.89,0.425, 
40,12.35,13.8,18.25,22.74,33.49,68.44,0.45, 
40,13.02,14.55,19.27,24.03,35.4,72.32,0.475, 
40,13.75,15.37,20.38,25.41,37.45,76.39,0.5, 
40,14.54,16.26,21.58,26.92,39.66,80.86,0.525, 
40,15.39,17.22,22.88,28.54,42.03,85.57,0.55, 
40,16.32,18.27,24.28,30.29,44.59,90.56,0.575, 
40,17.33,19.41,25.81,32.2,47.34,96.17,0.6, 
40,18.43,20.65,27.48,34.27,50.3,101.77,0.625, 
40,19.64,22.01,29.29,36.51,53.47,107.61,0.65, 
40,20.96,23.49,31.26,38.92,56.88,113.66,0.675, 
40,22.4,25.11,33.38,41.52,60.49,120.04,0.7, 
40,23.97,26.87,35.68,44.31,64.33,126.68,0.725, 
40,25.68,28.78,38.17,47.29,68.38,133.9,0.75, 
40,27.53,30.85,40.82,50.45,72.61,140.82,0.775, 
40,29.52,33.06,43.62,53.75,76.95,147.38,0.8, 
40,31.64,35.39,46.55,57.18,81.32,154.16,0.825, 
40,33.86,37.84,49.56,60.66,85.65,160.72,0.85, 
40,36.15,40.33,52.57,64.1,89.88,166.1,0.875, 
40,38.44,42.82,55.54,67.39,93.84,171.16,0.9, 
40,40.65,45.2,58.29,70.37,97.34,175.74,0.925, 
40,42.65,47.33,60.67,72.88,100.16,179.42,0.95, 
40,44.17,48.93,62.39,74.68,102.02,181.6,0.975, 
50,4.19,4.6,5.86,7.11,10.08,19.51,-0.2, 
50,4.56,5.01,6.38,7.75,11,21.28,-0.15, 
50,4.95,5.45,6.95,8.45,12,23.25,-0.1, 
50,5.38,5.92,7.57,9.21,13.1,25.46,-0.05, 
50,5.85,6.44,8.25,10.05,14.31,27.82,0, 
50,6.11,6.72,8.61,10.5,14.97,29.12,0.025, 
50,6.37,7.02,8.99,10.97,15.66,30.5,0.05, 
50,6.65,7.32,9.4,11.47,16.38,31.94,0.075, 
50,6.94,7.65,9.82,12,17.15,33.46,0.1, 
50,7.24,7.99,10.27,12.55,17.96,35.05,0.125, 
50,7.57,8.35,10.74,13.14,18.82,36.8,0.15, 
50,7.91,8.73,11.24,13.76,19.73,38.63,0.175, 
50,8.27,9.13,11.78,14.42,20.7,40.68,0.2, 
50,8.65,9.56,12.34,15.12,21.72,42.74,0.225, 
50,9.06,10.01,12.94,15.86,22.82,45.01,0.25, 
50,9.49,10.49,13.58,16.65,23.99,47.43,0.275, 
50,9.95,11.01,14.26,17.5,25.24,49.95,0.3, 
50,10.44,11.55,14.99,18.41,26.57,52.65,0.325, 
50,10.96,12.14,15.77,19.38,28.01,55.6,0.35, 
50,11.53,12.77,16.61,20.43,29.55,58.83,0.375, 
50,12.13,13.45,17.51,21.56,31.21,62.24,0.4, 
50,12.79,14.19,18.49,22.78,33.01,65.95,0.425, 
50,13.49,14.98,19.54,24.1,34.95,69.91,0.45, 
50,14.26,15.84,20.69,25.54,37.06,74.16,0.475, 
50,15.09,16.77,21.94,27.1,39.35,78.55,0.5, 
50,16,17.79,23.3,28.8,41.85,83.41,0.525, 
50,16.99,18.91,24.79,30.65,44.56,88.67,0.55, 
50,18.08,20.13,26.42,32.68,47.53,94.47,0.575, 
50,19.28,21.48,28.22,34.91,50.77,100.53,0.6, 
50,20.61,22.97,30.2,37.37,54.33,107.37,0.625, 
50,22.07,24.62,32.4,40.08,58.2,114.89,0.65, 
50,23.7,26.45,34.83,43.07,62.47,122.52,0.675, 
50,25.52,28.48,37.52,46.37,67.13,130.91,0.7, 
50,27.55,30.76,40.5,50.02,72.19,140.25,0.725, 
50,29.81,33.29,43.79,54.03,77.71,149.73,0.75, 
50,32.35,36.12,47.44,58.43,83.69,159.72,0.775, 
50,35.17,39.25,51.46,63.22,90.09,170.73,0.8, 
50,38.29,42.71,55.86,68.39,96.9,181.72,0.825, 
50,41.74,46.51,60.59,73.92,104.05,192.76,0.85, 
50,45.49,50.61,65.62,79.69,111.35,203.34,0.875, 
50,49.49,54.97,70.84,85.6,118.64,213.88,0.9, 
50,53.65,59.44,76.05,91.43,125.53,223.61,0.925, 
50,57.76,63.83,81,96.77,131.75,231.78,0.95, 
50,61.26,67.51,84.97,100.97,136.44,235.97,0.975, 
60,4.33,4.73,5.97,7.21,10.09,19.17,-0.2, 
60,4.71,5.16,6.52,7.87,11.03,20.93,-0.15, 
60,5.13,5.62,7.11,8.59,12.06,22.87,-0.1, 
60,5.58,6.12,7.76,9.38,13.18,25.05,-0.05, 
60,6.08,6.67,8.47,10.24,14.41,27.49,0, 
60,6.35,6.96,8.85,10.71,15.08,28.82,0.025, 
60,6.63,7.27,9.25,11.2,15.78,30.17,0.05, 
60,6.92,7.6,9.66,11.71,16.52,31.64,0.075, 
60,7.23,7.94,10.11,12.26,17.3,33.16,0.1, 
60,7.55,8.3,10.58,12.83,18.13,34.76,0.125, 
60,7.89,8.68,11.07,13.44,19.01,36.47,0.15, 
60,8.26,9.08,11.59,14.08,19.94,38.31,0.175, 
60,8.64,9.5,12.14,14.76,20.92,40.27,0.2, 
60,9.05,9.95,12.73,15.49,21.97,42.35,0.225, 
60,9.48,10.43,13.36,16.26,23.09,44.57,0.25, 
60,9.94,10.94,14.03,17.09,24.29,46.92,0.275, 
60,10.43,11.49,14.74,17.97,25.57,49.51,0.3, 
60,10.95,12.07,15.51,18.91,26.95,52.26,0.325, 
60,11.51,12.69,16.33,19.93,28.43,55.27,0.35, 
60,12.11,13.37,17.21,21.02,30.02,58.47,0.375, 
60,12.76,14.09,18.16,22.2,31.75,61.99,0.4, 
60,13.46,14.87,19.19,23.48,33.61,65.68,0.425, 
60,14.22,15.71,20.31,24.87,35.62,69.77,0.45, 
60,15.04,16.63,21.52,26.38,37.83,74.13,0.475, 
60,15.94,17.64,22.85,28.03,40.23,78.91,0.5, 
60,16.92,18.73,24.3,29.83,42.86,83.83,0.525, 
60,18,19.94,25.9,31.81,45.75,89.59,0.55, 
60,19.18,21.27,27.66,34,48.94,95.8,0.575, 
60,20.5,22.74,29.6,36.43,52.45,102.46,0.6, 
60,21.96,24.37,31.77,39.11,56.35,110,0.625, 
60,23.59,26.2,34.18,42.11,60.67,118.21,0.65, 
60,25.42,28.25,36.89,45.45,65.45,127.25,0.675, 
60,27.48,30.55,39.93,49.2,70.77,137.53,0.7, 
60,29.8,33.16,43.37,53.41,76.69,148.77,0.725, 
60,32.45,36.11,47.24,58.13,83.31,160.56,0.75, 
60,35.47,39.48,51.62,63.44,90.65,173.35,0.775, 
60,38.91,43.31,56.56,69.39,98.78,187.46,0.8, 
60,42.83,47.68,62.12,76.03,107.7,201.79,0.825, 
60,47.3,52.61,68.35,83.38,117.4,216.93,0.85, 
60,52.34,58.14,75.2,91.33,127.6,232.29,0.875, 
60,57.96,64.28,82.62,99.83,138.3,247.95,0.9, 
60,64.09,70.88,90.43,108.62,148.98,263.07,0.925, 
60,70.48,77.69,98.24,117.15,159.11,274.42,0.95, 
60,76.39,83.92,105.02,124.42,167.29,284.22,0.975, 
70,4.42,4.82,6.06,7.27,10.1,18.69,-0.2, 
70,4.82,5.27,6.62,7.95,11.05,20.48,-0.15, 
70,5.25,5.74,7.23,8.68,12.08,22.42,-0.1, 
70,5.73,6.26,7.89,9.49,13.22,24.58,-0.05, 
70,6.24,6.83,8.62,10.37,14.47,26.92,0, 
70,6.52,7.14,9.01,10.85,15.14,28.2,0.025, 
70,6.81,7.46,9.42,11.35,15.85,29.57,0.05, 
70,7.11,7.79,9.85,11.88,16.6,30.95,0.075, 
70,7.43,8.15,10.31,12.43,17.4,32.47,0.1, 
70,7.77,8.52,10.79,13.02,18.23,34.07,0.125, 
70,8.13,8.91,11.3,13.64,19.12,35.76,0.15, 
70,8.5,9.33,11.84,14.3,20.06,37.62,0.175, 
70,8.9,9.77,12.41,15,21.06,39.56,0.2, 
70,9.33,10.24,13.01,15.74,22.12,41.66,0.225, 
70,9.78,10.74,13.66,16.53,23.26,43.84,0.25, 
70,10.25,11.27,14.34,17.37,24.47,46.2,0.275, 
70,10.77,11.83,15.08,18.28,25.77,48.74,0.3, 
70,11.31,12.44,15.87,19.24,27.17,51.46,0.325, 
70,11.9,13.08,16.71,20.28,28.67,54.43,0.35, 
70,12.53,13.78,17.62,21.41,30.29,57.63,0.375, 
70,13.2,14.53,18.6,22.62,32.05,60.99,0.4, 
70,13.93,15.35,19.66,23.93,33.95,64.76,0.425, 
70,14.73,16.23,20.82,25.36,36.02,68.82,0.45, 
70,15.59,17.19,22.08,26.91,38.28,73.3,0.475, 
70,16.53,18.24,23.45,28.61,40.73,78.2,0.5, 
70,17.56,19.39,24.96,30.48,43.45,83.46,0.525, 
70,18.69,20.65,26.62,32.55,46.45,89.46,0.55, 
70,19.95,22.05,28.47,34.82,49.76,95.84,0.575, 
70,21.34,23.6,30.51,37.35,53.44,102.83,0.6, 
70,22.88,25.33,32.79,40.18,57.51,110.68,0.625, 
70,24.62,27.27,35.35,43.35,62.07,119.65,0.65, 
70,26.58,29.46,38.23,46.92,67.22,129.39,0.675, 
70,28.8,31.94,41.51,50.96,73.02,140.56,0.7, 
70,31.34,34.77,45.24,55.55,79.55,152.85,0.725, 
70,34.25,38.02,49.51,60.77,86.92,166.49,0.75, 
70,37.61,41.77,54.4,66.75,95.35,181.97,0.775, 
70,41.51,46.12,60.05,73.61,104.86,198.67,0.8, 
70,46.06,51.18,66.57,81.43,115.53,216.63,0.825, 
70,51.37,57.05,74.04,90.34,127.46,236.33,0.85, 
70,57.54,63.87,82.57,100.36,140.56,257.25,0.875, 
70,64.68,71.66,92.13,111.41,154.65,278.41,0.9, 
70,72.75,80.4,102.55,123.22,169.37,298.89,0.925, 
70,81.57,89.82,113.45,135.3,183.83,317.9,0.95, 
70,90.2,98.9,123.47,146.07,196.19,331.84,0.975, 
80,4.49,4.9,6.12,7.32,10.09,18.39,-0.2, 
80,4.9,5.35,6.69,8.01,11.05,20.18,-0.15, 
80,5.35,5.84,7.31,8.76,12.1,22.11,-0.1, 
80,5.84,6.37,7.99,9.58,13.25,24.25,-0.05, 
80,6.37,6.96,8.73,10.48,14.51,26.57,0, 
80,6.65,7.27,9.13,10.96,15.19,27.83,0.025, 
80,6.95,7.6,9.55,11.47,15.91,29.15,0.05, 
80,7.27,7.94,9.99,12.01,16.66,30.56,0.075, 
80,7.6,8.31,10.46,12.57,17.46,32.06,0.1, 
80,7.95,8.69,10.95,13.17,18.3,33.69,0.125, 
80,8.31,9.1,11.47,13.8,19.2,35.41,0.15, 
80,8.7,9.52,12.02,14.47,20.15,37.19,0.175, 
80,9.11,9.98,12.6,15.18,21.16,39.08,0.2, 
80,9.55,10.46,13.22,15.93,22.23,41.13,0.225, 
80,10.01,10.97,13.87,16.74,23.37,43.34,0.25, 
80,10.5,11.51,14.58,17.59,24.6,45.63,0.275, 
80,11.03,12.09,15.33,18.51,25.91,48.14,0.3, 
80,11.59,12.72,16.13,19.49,27.32,50.88,0.325, 
80,12.2,13.38,16.99,20.55,28.84,53.82,0.35, 
80,12.84,14.1,17.92,21.69,30.47,56.95,0.375, 
80,13.54,14.87,18.93,22.92,32.25,60.4,0.4, 
80,14.3,15.71,20.01,24.26,34.18,64.12,0.425, 
80,15.11,16.62,21.19,25.71,36.26,68.15,0.45, 
80,16,17.61,22.48,27.3,38.55,72.63,0.475, 
80,16.98,18.69,23.89,29.03,41.06,77.45,0.5, 
80,18.04,19.88,25.44,30.94,43.82,82.63,0.525, 
80,19.22,21.18,27.15,33.05,46.87,88.65,0.55, 
80,20.52,22.63,29.04,35.39,50.25,95,0.575, 
80,21.97,24.24,31.15,38,54.02,102.62,0.6, 
80,23.58,26.03,33.52,40.92,58.23,110.71,0.625, 
80,25.4,28.06,36.18,44.21,62.98,119.98,0.65, 
80,27.45,30.35,39.19,47.94,68.35,130.33,0.675, 
80,29.79,32.96,42.63,52.18,74.44,142.06,0.7, 
80,32.47,35.96,46.57,57.05,81.4,155.2,0.725, 
80,35.58,39.42,51.12,62.64,89.35,170.41,0.75, 
80,39.2,43.46,56.41,69.15,98.52,187.58,0.775, 
80,43.45,48.21,62.59,76.7,109.1,206.58,0.8, 
80,48.49,53.82,69.85,85.48,121.25,227.32,0.825, 
80,54.5,60.47,78.39,95.74,135.18,250.97,0.85, 
80,61.65,68.38,88.39,107.58,150.96,277.65,0.875, 
80,70.15,77.71,99.96,121.06,168.45,305.14,0.9, 
80,80.11,88.54,113.05,136.01,187.4,332.73,0.925, 
80,91.43,100.68,127.19,151.85,206.59,358.78,0.95, 
80,103.09,112.99,140.9,166.69,223.68,378.64,0.975, 
90,4.55,4.95,6.17,7.36,10.09,18.26,-0.2, 
90,4.97,5.41,6.75,8.06,11.06,19.98,-0.15, 
90,5.43,5.91,7.38,8.82,12.11,21.91,-0.1, 
90,5.93,6.46,8.07,9.65,13.27,24.05,-0.05, 
90,6.47,7.06,8.83,10.56,14.54,26.37,0, 
90,6.76,7.38,9.24,11.05,15.23,27.64,0.025, 
90,7.07,7.71,9.66,11.57,15.95,29,0.05, 
90,7.39,8.07,10.11,12.11,16.7,30.42,0.075, 
90,7.73,8.44,10.58,12.68,17.51,31.93,0.1, 
90,8.08,8.83,11.08,13.28,18.36,33.51,0.125, 
90,8.46,9.24,11.61,13.92,19.25,35.17,0.15, 
90,8.86,9.68,12.17,14.6,20.2,36.99,0.175, 
90,9.28,10.14,12.76,15.32,21.21,38.88,0.2, 
90,9.72,10.63,13.39,16.08,22.29,40.91,0.225, 
90,10.2,11.16,14.06,16.9,23.44,43.12,0.25, 
90,10.7,11.71,14.77,17.76,24.68,45.5,0.275, 
90,11.24,12.31,15.53,18.69,25.99,48.07,0.3, 
90,11.82,12.94,16.35,19.69,27.41,50.77,0.325, 
90,12.43,13.62,17.23,20.76,28.93,53.69,0.35, 
90,13.1,14.36,18.17,21.91,30.58,56.82,0.375, 
90,13.81,15.15,19.19,23.16,32.37,60.2,0.4, 
90,14.59,16,20.3,24.52,34.3,63.89,0.425, 
90,15.42,16.93,21.5,25.99,36.41,67.98,0.45, 
90,16.34,17.94,22.81,27.59,38.72,72.45,0.475, 
90,17.33,19.05,24.25,29.35,41.25,77.38,0.5, 
90,18.43,20.26,25.83,31.29,44.04,82.77,0.525, 
90,19.63,21.6,27.57,33.44,47.12,88.81,0.55, 
90,20.97,23.08,29.5,35.82,50.55,95.43,0.575, 
90,22.45,24.73,31.66,38.48,54.39,102.87,0.6, 
90,24.12,26.58,34.08,41.46,58.68,111.21,0.625, 
90,25.99,28.67,36.81,44.84,63.53,120.81,0.65, 
90,28.12,31.03,39.91,48.67,69.03,131.51,0.675, 
90,30.54,33.74,43.46,53.05,75.34,143.8,0.7, 
90,33.34,36.86,47.56,58.11,82.63,157.9,0.725, 
90,36.6,40.48,52.33,63.99,91.05,173.47,0.75, 
90,40.43,44.74,57.91,70.87,100.8,191.69,0.775, 
90,44.96,49.79,64.53,78.96,112.19,211.67,0.8, 
90,50.39,55.85,72.41,88.53,125.58,235.53,0.825, 
90,56.95,63.14,81.84,99.92,141.22,262.36,0.85, 
90,64.95,72.01,93.13,113.41,159.49,294.27,0.875, 
90,74.72,82.75,106.57,129.25,180.34,328.58,0.9, 
90,86.49,95.62,122.26,147.39,203.55,362.77,0.925, 
90,100.33,110.56,139.91,167.22,227.95,398.17,0.95, 
90,115.27,126.34,157.73,186.68,250.72,426.38,0.975, 
100,4.6,5,6.21,7.39,10.09,18.14,-0.2, 
100,5.03,5.47,6.8,8.1,11.07,19.89,-0.15, 
100,5.5,5.98,7.44,8.87,12.14,21.81,-0.1, 
100,6.01,6.54,8.14,9.71,13.3,23.93,-0.05, 
100,6.56,7.15,8.91,10.63,14.58,26.26,0, 
100,6.86,7.47,9.32,11.13,15.27,27.51,0.025, 
100,7.17,7.81,9.75,11.65,15.99,28.82,0.05, 
100,7.5,8.17,10.21,12.2,16.76,30.23,0.075, 
100,7.84,8.55,10.69,12.78,17.57,31.75,0.1, 
100,8.2,8.95,11.19,13.38,18.42,33.33,0.125, 
100,8.59,9.37,11.73,14.03,19.32,35,0.15, 
100,8.99,9.81,12.29,14.71,20.28,36.81,0.175, 
100,9.42,10.28,12.89,15.44,21.31,38.73,0.2, 
100,9.88,10.78,13.53,16.21,22.39,40.77,0.225, 
100,10.36,11.32,14.21,17.03,23.55,42.91,0.25, 
100,10.88,11.88,14.93,17.91,24.78,45.26,0.275, 
100,11.42,12.49,15.7,18.85,26.11,47.77,0.3, 
100,12.01,13.13,16.53,19.85,27.53,50.39,0.325, 
100,12.64,13.83,17.42,20.93,29.06,53.28,0.35, 
100,13.32,14.57,18.37,22.1,30.71,56.44,0.375, 
100,14.05,15.38,19.41,23.36,32.5,59.76,0.4, 
100,14.84,16.25,20.53,24.73,34.45,63.46,0.425, 
100,15.69,17.19,21.74,26.22,36.57,67.51,0.45, 
100,16.63,18.22,23.07,27.85,38.89,71.9,0.475, 
100,17.64,19.35,24.53,29.63,41.45,76.88,0.5, 
100,18.76,20.59,26.13,31.59,44.27,82.36,0.525, 
100,19.99,21.95,27.9,33.76,47.39,88.43,0.55, 
100,21.36,23.46,29.86,36.17,50.86,95.17,0.575, 
100,22.88,25.14,32.06,38.86,54.73,102.52,0.6, 
100,24.58,27.03,34.52,41.9,59.11,111.1,0.625, 
100,26.5,29.17,37.3,45.33,64.04,120.51,0.65, 
100,28.68,31.59,40.48,49.24,69.69,131.42,0.675, 
100,31.18,34.37,44.11,53.73,76.16,143.79,0.7, 
100,34.07,37.59,48.32,58.94,83.66,158.05,0.725, 
100,37.45,41.34,53.25,65.02,92.37,174.93,0.75, 
100,41.43,45.78,59.08,72.2,102.59,194.24,0.775, 
100,46.19,51.08,66.03,80.71,114.7,216.79,0.8, 
100,51.92,57.49,74.4,90.95,129.14,243.22,0.825, 
100,58.98,65.33,84.57,103.31,146.32,273.17,0.85, 
100,67.71,75.02,97,118.27,166.71,307.64,0.875, 
100,78.58,87.05,112.21,136.28,190.69,347.41,0.9, 
100,92.09,101.86,130.46,157.5,218.18,390.53,0.925, 
100,108.5,119.61,151.63,181.56,248.09,432.82,0.95, 
100,126.95,139.17,173.91,205.98,276.66,469.68,0.975, 
150,4.77,5.16,6.36,7.51,10.12,17.63,-0.2, 
150,5.22,5.66,6.98,8.24,11.12,19.39,-0.15, 
150,5.72,6.2,7.65,9.04,12.21,21.27,-0.1, 
150,6.26,6.79,8.38,9.91,13.4,23.39,-0.05, 
150,6.85,7.43,9.18,10.87,14.71,25.75,0, 
150,7.17,7.78,9.61,11.38,15.41,27.02,0.025, 
150,7.5,8.14,10.06,11.92,16.15,28.3,0.05, 
150,7.85,8.52,10.54,12.49,16.93,29.7,0.075, 
150,8.22,8.92,11.04,13.09,17.76,31.16,0.1, 
150,8.6,9.34,11.57,13.72,18.63,32.7,0.125, 
150,9.01,9.79,12.12,14.39,19.54,34.37,0.15, 
150,9.44,10.26,12.71,15.09,20.52,36.18,0.175, 
150,9.9,10.76,13.34,15.84,21.55,38.03,0.2, 
150,10.38,11.29,14.01,16.64,22.65,40.02,0.225, 
150,10.9,11.85,14.71,17.49,23.83,42.17,0.25, 
150,11.44,12.45,15.47,18.39,25.08,44.47,0.275, 
150,12.03,13.09,16.27,19.36,26.43,46.94,0.3, 
150,12.65,13.77,17.13,20.4,27.87,49.59,0.325, 
150,13.32,14.5,18.06,21.51,29.42,52.43,0.35, 
150,14.05,15.29,19.06,22.72,31.1,55.56,0.375, 
150,14.82,16.14,20.13,24.02,32.91,58.95,0.4, 
150,15.66,17.07,21.3,25.43,34.89,62.6,0.425, 
150,16.57,18.07,22.57,26.97,37.04,66.58,0.45, 
150,17.57,19.16,23.95,28.65,39.39,71.04,0.475, 
150,18.65,20.35,25.47,30.49,41.97,76.01,0.5, 
150,19.84,21.66,27.14,32.52,44.84,81.5,0.525, 
150,21.16,23.1,28.99,34.76,48.02,87.4,0.55, 
150,22.62,24.71,31.05,37.27,51.56,94.04,0.575, 
150,24.25,26.5,33.35,40.08,55.55,101.63,0.6, 
150,26.07,28.52,35.94,43.24,60.06,110.28,0.625, 
150,28.14,30.8,38.88,46.84,65.18,120.25,0.65, 
150,30.5,33.4,42.25,50.97,71.11,131.64,0.675, 
150,33.21,36.41,46.15,55.74,77.96,144.68,0.7, 
150,36.37,39.91,50.69,61.33,86.03,160.21,0.725, 
150,40.09,44.04,56.07,67.95,95.55,178.89,0.75, 
150,44.54,48.99,62.54,75.93,107.08,201.19,0.775, 
150,49.95,55.01,70.44,85.67,121.1,228.09,0.8, 
150,56.67,62.51,80.26,97.78,138.55,261.01,0.825, 
150,65.21,72.03,92.75,113.17,160.54,301.78,0.85, 
150,76.35,84.46,109.01,133.12,188.69,352.55,0.875, 
150,91.33,101.12,130.66,159.39,225.25,417.19,0.9, 
150,111.95,123.99,159.92,194.4,272.6,497.74,0.925, 
150,140.66,155.47,198.89,240,331.88,591.79,0.95, 
150,178.95,196.68,247.44,294.61,398.88,681.65,0.975, 
200,4.87,5.27,6.45,7.59,10.16,17.4,-0.2, 
200,5.34,5.78,7.08,8.34,11.17,19.15,-0.15, 
200,5.86,6.34,7.77,9.16,12.28,21.05,-0.1, 
200,6.42,6.95,8.53,10.05,13.49,23.14,-0.05, 
200,7.03,7.61,9.35,11.03,14.81,25.48,0, 
200,7.36,7.97,9.79,11.55,15.53,26.73,0.025, 
200,7.7,8.34,10.26,12.1,16.28,28.07,0.05, 
200,8.07,8.74,10.75,12.68,17.07,29.48,0.075, 
200,8.45,9.15,11.26,13.3,17.9,30.98,0.1, 
200,8.85,9.59,11.8,13.94,18.78,32.54,0.125, 
200,9.27,10.05,12.38,14.62,19.71,34.2,0.15, 
200,9.72,10.54,12.98,15.35,20.7,35.96,0.175, 
200,10.19,11.05,13.62,16.11,21.75,37.82,0.2, 
200,10.69,11.6,14.3,16.93,22.86,39.82,0.225, 
200,11.23,12.18,15.03,17.79,24.05,41.96,0.25, 
200,11.8,12.8,15.8,18.72,25.32,44.22,0.275, 
200,12.4,13.46,16.63,19.7,26.68,46.61,0.3, 
200,13.05,14.17,17.52,20.76,28.13,49.24,0.325, 
200,13.75,14.93,18.47,21.9,29.7,52.02,0.35, 
200,14.5,15.74,19.49,23.13,31.39,55.05,0.375, 
200,15.3,16.63,20.6,24.45,33.24,58.39,0.4, 
200,16.18,17.58,21.8,25.89,35.23,62.06,0.425, 
200,17.12,18.61,23.1,27.46,37.4,66.1,0.45, 
200,18.15,19.74,24.52,29.17,39.78,70.41,0.475, 
200,19.28,20.98,26.07,31.05,42.4,75.17,0.5, 
200,20.52,22.33,27.79,33.12,45.29,80.49,0.525, 
200,21.89,23.83,29.68,35.41,48.5,86.46,0.55, 
200,23.4,25.49,31.79,37.96,52.09,93.24,0.575, 
200,25.09,27.35,34.15,40.82,56.1,100.89,0.6, 
200,27,29.44,36.81,44.05,60.66,109.44,0.625, 
200,29.15,31.81,39.84,47.73,65.89,119.26,0.65, 
200,31.61,34.51,43.31,51.95,71.88,130.66,0.675, 
200,34.45,37.64,47.32,56.85,78.86,143.88,0.7, 
200,37.76,41.29,52.03,62.61,87.1,159.56,0.725, 
200,41.67,45.62,57.62,69.47,96.91,178.53,0.75, 
200,46.38,50.83,64.39,77.77,108.85,201.25,0.775, 
200,52.14,57.22,72.7,88,123.57,230.16,0.8, 
200,59.36,65.24,83.18,100.93,142.19,266.26,0.825, 
200,68.66,75.6,96.76,117.68,166.33,312.16,0.85, 
200,81.09,89.47,114.94,140.1,198.45,373.17,0.875, 
200,98.41,108.79,140.27,171.2,242.58,454.88,0.9, 
200,123.76,137.01,177.01,215.86,304.75,564.6,0.925, 
200,162.61,179.96,231.61,281.06,392.18,707.1,0.95, 
200,222.22,244.78,310.1,370.95,506.46,872.55,0.975, 
250,4.95,5.34,6.52,7.65,10.21,17.41,-0.2, 
250,5.43,5.87,7.17,8.41,11.24,19.18,-0.15, 
250,5.96,6.44,7.87,9.25,12.35,21.06,-0.1, 
250,6.54,7.06,8.64,10.15,13.58,23.18,-0.05, 
250,7.16,7.75,9.48,11.15,14.92,25.51,0, 
250,7.5,8.11,9.93,11.68,15.64,26.77,0.025, 
250,7.86,8.49,10.4,12.24,16.4,28.08,0.05, 
250,8.23,8.9,10.9,12.83,17.2,29.49,0.075, 
250,8.62,9.32,11.43,13.46,18.04,30.99,0.1, 
250,9.03,9.77,11.98,14.11,18.93,32.57,0.125, 
250,9.46,10.24,12.57,14.8,19.87,34.19,0.15, 
250,9.92,10.74,13.18,15.54,20.86,35.96,0.175, 
250,10.41,11.27,13.84,16.32,21.92,37.84,0.2, 
250,10.93,11.83,14.54,17.14,23.05,39.83,0.225, 
250,11.47,12.43,15.28,18.02,24.24,41.95,0.25, 
250,12.06,13.06,16.06,18.96,25.53,44.19,0.275, 
250,12.68,13.74,16.91,19.97,26.9,46.64,0.3, 
250,13.35,14.46,17.81,21.04,28.37,49.27,0.325, 
250,14.06,15.24,18.78,22.2,29.96,52.11,0.35, 
250,14.83,16.08,19.82,23.45,31.67,55.15,0.375, 
250,15.65,16.98,20.95,24.79,33.52,58.53,0.4, 
250,16.55,17.96,22.17,26.25,35.53,62.09,0.425, 
250,17.52,19.02,23.49,27.84,37.71,66.05,0.45, 
250,18.58,20.17,24.94,29.57,40.1,70.44,0.475, 
250,19.74,21.44,26.53,31.47,42.73,75.26,0.5, 
250,21.01,22.83,28.28,33.57,45.64,80.67,0.525, 
250,22.42,24.36,30.21,35.89,48.87,86.7,0.55, 
250,23.98,26.07,32.36,38.48,52.48,93.29,0.575, 
250,25.72,27.97,34.77,41.38,56.54,100.74,0.6, 
250,27.68,30.11,37.48,44.66,61.13,109.37,0.625, 
250,29.89,32.55,40.57,48.4,66.37,119.1,0.65, 
250,32.43,35.33,44.11,52.69,72.42,130.49,0.675, 
250,35.36,38.55,48.21,57.66,79.47,143.78,0.7, 
250,38.78,42.31,53.03,63.52,87.78,159.77,0.725, 
250,42.82,46.77,58.76,70.51,97.76,178.83,0.75, 
250,47.7,52.16,65.69,78.99,109.88,202.17,0.775, 
250,53.69,58.79,74.26,89.5,124.98,231.59,0.8, 
250,61.23,67.15,85.12,102.86,144.19,268.96,0.825, 
250,71.02,78.01,99.3,120.35,169.41,317.29,0.85, 
250,84.24,92.72,118.58,144.14,203.7,383.97,0.875, 
250,103.01,113.66,146.03,177.96,252.29,474.62,0.9, 
250,131.49,145.38,187.54,228.94,324.47,607.02,0.925, 
250,178.05,197.1,254.35,309.58,435.27,800.81,0.95, 
250,258.09,284.88,362.85,436,599.26,1055.74,0.975, 
300,4.99,5.39,6.56,7.69,10.23,17.36,-0.2, 
300,5.49,5.92,7.22,8.46,11.27,19.13,-0.15, 
300,6.02,6.5,7.93,9.31,12.4,21.08,-0.1, 
300,6.61,7.14,8.71,10.22,13.63,23.21,-0.05, 
300,7.25,7.83,9.57,11.23,14.99,25.55,0, 
300,7.59,8.2,10.02,11.77,15.72,26.8,0.025, 
300,7.95,8.59,10.5,12.34,16.48,28.14,0.05, 
300,8.33,9,11.01,12.94,17.29,29.54,0.075, 
300,8.73,9.44,11.54,13.57,18.14,31.01,0.1, 
300,9.15,9.89,12.1,14.23,19.04,32.59,0.125, 
300,9.59,10.37,12.7,14.93,19.99,34.26,0.15, 
300,10.06,10.88,13.32,15.68,21,35.97,0.175, 
300,10.56,11.42,13.99,16.46,22.07,37.88,0.2, 
300,11.08,11.99,14.7,17.3,23.2,39.89,0.225, 
300,11.64,12.6,15.45,18.19,24.41,42.01,0.25, 
300,12.24,13.24,16.25,19.14,25.71,44.31,0.275, 
300,12.87,13.93,17.1,20.16,27.09,46.74,0.3, 
300,13.55,14.67,18.02,21.25,28.58,49.36,0.325, 
300,14.28,15.47,19.01,22.42,30.17,52.15,0.35, 
300,15.06,16.32,20.06,23.68,31.9,55.25,0.375, 
300,15.91,17.24,21.21,25.04,33.76,58.58,0.4, 
300,16.82,18.23,22.45,26.52,35.78,62.19,0.425, 
300,17.82,19.32,23.8,28.13,37.99,66.19,0.45, 
300,18.9,20.49,25.27,29.89,40.4,70.55,0.475, 
300,20.08,21.78,26.88,31.81,43.06,75.38,0.5, 
300,21.38,23.2,28.65,33.93,45.99,80.64,0.525, 
300,22.82,24.77,30.62,36.29,49.24,86.52,0.55, 
300,24.41,26.51,32.81,38.91,52.86,93.16,0.575, 
300,26.19,28.46,35.25,41.85,56.95,100.59,0.6, 
300,28.2,30.65,38.02,45.18,61.57,109.12,0.625, 
300,30.47,33.13,41.16,48.96,66.85,118.95,0.65, 
300,33.06,35.98,44.76,53.3,72.93,130.13,0.675, 
300,36.06,39.27,48.93,58.35,80.04,143.35,0.7, 
300,39.57,43.12,53.84,64.3,88.42,159.08,0.725, 
300,43.73,47.69,59.68,71.4,98.47,177.96,0.75, 
300,48.74,53.21,66.76,80.02,110.7,201.31,0.775, 
300,54.9,60.01,75.53,90.72,125.96,230.09,0.8, 
300,62.69,68.62,86.66,104.37,145.52,267.45,0.825, 
300,72.83,79.86,101.27,122.31,171.31,316.6,0.85, 
300,86.61,95.16,121.23,146.94,206.79,384.91,0.875, 
300,106.39,117.19,150.08,182.53,258.01,482.77,0.9, 
300,136.99,151.3,194.88,237.72,336.75,630.68,0.925, 
300,189.25,209.52,270.73,330.04,465.81,862.69,0.95, 
300,288,318.19,407.04,490.97,678.2,1204.17,0.975, 
400,5.01,5.41,6.58,7.71,10.25,17.31,-0.2, 
400,5.52,5.95,7.25,8.5,11.3,19.09,-0.15, 
400,6.07,6.55,7.98,9.35,12.44,21.03,-0.1, 
400,6.66,7.19,8.77,10.28,13.69,23.17,-0.05, 
400,7.31,7.9,9.64,11.3,15.06,25.55,0, 
400,7.66,8.28,10.1,11.85,15.8,26.8,0.025, 
400,8.03,8.67,10.59,12.43,16.58,28.13,0.05, 
400,8.42,9.09,11.1,13.04,17.39,29.52,0.075, 
400,8.82,9.53,11.65,13.68,18.25,30.99,0.1, 
400,9.25,10,12.22,14.35,19.16,32.54,0.125, 
400,9.7,10.49,12.82,15.06,20.12,34.16,0.15, 
400,10.18,11.01,13.46,15.82,21.14,35.95,0.175, 
400,10.69,11.56,14.14,16.62,22.22,37.83,0.2, 
400,11.23,12.14,14.86,17.47,23.37,39.87,0.225, 
400,11.8,12.76,15.62,18.38,24.6,42.02,0.25, 
400,12.41,13.42,16.44,19.34,25.9,44.29,0.275, 
400,13.06,14.13,17.31,20.38,27.3,46.73,0.3, 
400,13.76,14.88,18.24,21.48,28.8,49.4,0.325, 
400,14.5,15.69,19.25,22.67,30.41,52.25,0.35, 
400,15.31,16.56,20.33,23.96,32.15,55.36,0.375, 
400,16.17,17.51,21.5,25.34,34.04,58.71,0.4, 
400,17.11,18.53,22.76,26.84,36.08,62.34,0.425, 
400,18.13,19.64,24.14,28.48,38.32,66.3,0.45, 
400,19.24,20.84,25.64,30.27,40.76,70.69,0.475, 
400,20.46,22.17,27.29,32.23,43.44,75.47,0.5, 
400,21.8,23.63,29.1,34.39,46.4,80.74,0.525, 
400,23.28,25.24,31.11,36.79,49.69,86.61,0.55, 
400,24.92,27.03,33.34,39.46,53.35,93.17,0.575, 
400,26.76,29.03,35.85,42.45,57.46,100.58,0.6, 
400,28.83,31.29,38.68,45.84,62.13,108.94,0.625, 
400,31.17,33.85,41.89,49.7,67.46,118.67,0.65, 
400,33.86,36.79,45.58,54.13,73.6,129.93,0.675, 
400,36.97,40.19,49.86,59.28,80.76,143.04,0.7, 
400,40.6,44.17,54.89,65.34,89.2,158.7,0.725, 
400,44.93,48.9,60.89,72.57,99.32,177.3,0.75, 
400,50.14,54.63,68.16,81.37,111.73,200.55,0.775, 
400,56.58,61.7,77.18,92.3,127.17,229.33,0.8, 
400,64.71,70.65,88.65,106.26,147,267.21,0.825, 
400,75.35,82.39,103.75,124.73,173.27,317.38,0.85, 
400,89.87,98.45,124.54,150.21,209.88,387.73,0.875, 
400,110.92,121.82,154.97,187.66,263.72,492.14,0.9, 
400,144.21,158.86,203.5,247.5,349.97,656.3,0.925, 
400,204.06,225.55,290.96,355.08,502.85,939.75,0.95, 
400,332.87,368.32,473.91,574.59,801.32,1447.65,0.975, 
600,5.06,5.45,6.62,7.75,10.27,17.21,-0.2, 
600,5.57,6.01,7.3,8.55,11.33,19,-0.15, 
600,6.13,6.61,8.04,9.41,12.49,20.95,-0.1, 
600,6.74,7.27,8.85,10.36,13.76,23.08,-0.05, 
600,7.41,8,9.73,11.4,15.15,25.42,0, 
600,7.77,8.38,10.21,11.96,15.89,26.7,0.025, 
600,8.15,8.79,10.71,12.55,16.68,28.02,0.05, 
600,8.55,9.22,11.23,13.17,17.51,29.42,0.075, 
600,8.96,9.67,11.79,13.82,18.38,30.92,0.1, 
600,9.4,10.15,12.37,14.51,19.3,32.5,0.125, 
600,9.87,10.65,12.99,15.23,20.27,34.17,0.15, 
600,10.36,11.18,13.64,16,21.3,35.93,0.175, 
600,10.88,11.75,14.33,16.82,22.4,37.78,0.2, 
600,11.43,12.35,15.07,17.69,23.57,39.78,0.225, 
600,12.02,12.98,15.85,18.61,24.81,41.9,0.25, 
600,12.65,13.66,16.68,19.59,26.13,44.18,0.275, 
600,13.32,14.39,17.58,20.65,27.55,46.64,0.3, 
600,14.04,15.16,18.53,21.78,29.07,49.25,0.325, 
600,14.81,16,19.56,22.99,30.7,52.09,0.35, 
600,15.63,16.9,20.66,24.3,32.47,55.12,0.375, 
600,16.53,17.87,21.86,25.71,34.38,58.42,0.4, 
600,17.5,18.92,23.16,27.25,36.45,62.03,0.425, 
600,18.55,20.06,24.57,28.92,38.71,65.91,0.45, 
600,19.7,21.3,26.11,30.74,41.18,70.23,0.475, 
600,20.95,22.67,27.79,32.74,43.89,74.95,0.5, 
600,22.34,24.17,29.65,34.95,46.9,80.18,0.525, 
600,23.87,25.83,31.71,37.4,50.22,85.92,0.55, 
600,25.57,27.68,34.01,40.12,53.92,92.4,0.575, 
600,27.48,29.75,36.58,43.18,58.08,99.76,0.6, 
600,29.62,32.09,39.48,46.64,62.8,108.17,0.625, 
600,32.06,34.74,42.78,50.57,68.19,117.78,0.65, 
600,34.85,37.78,46.57,55.1,74.39,128.98,0.675, 
600,38.08,41.3,50.97,60.36,81.63,141.9,0.7, 
600,41.88,45.44,56.15,66.55,90.16,157.38,0.725, 
600,46.38,50.37,62.32,73.95,100.37,175.98,0.75, 
600,51.84,56.32,69.82,82.95,112.84,198.63,0.775, 
600,58.58,63.69,79.12,94.15,128.4,227.21,0.8, 
600,67.12,73.06,90.96,108.46,148.37,264.03,0.825, 
600,78.32,85.36,106.6,127.4,174.95,314.08,0.85, 
600,93.7,102.27,128.19,153.65,212.1,384.85,0.875, 
600,116.16,127.04,160.05,192.54,267.6,490.45,0.9, 
600,152.18,166.93,211.8,256.08,358.77,666.97,0.925, 
600,219.48,241.83,309.87,376.93,532.98,1002.58,0.95, 
600,385.61,426.81,551.15,672.05,948.55,1756.31,0.975, 
1200,5.12,5.51,6.68,7.79,10.29,17.14,-0.2, 
1200,5.65,6.08,7.37,8.61,11.37,18.96,-0.15, 
1200,6.23,6.71,8.13,9.5,12.55,20.93,-0.1, 
1200,6.86,7.39,8.96,10.47,13.84,23.09,-0.05, 
1200,7.55,8.13,9.87,11.53,15.25,25.5,0, 
1200,7.92,8.53,10.36,12.11,16.01,26.77,0.025, 
1200,8.31,8.95,10.87,12.71,16.81,28.12,0.05, 
1200,8.72,9.4,11.41,13.34,17.66,29.55,0.075, 
1200,9.15,9.86,11.98,14.01,18.54,31.06,0.1, 
1200,9.61,10.36,12.58,14.72,19.48,32.65,0.125, 
1200,10.09,10.88,13.22,15.46,20.48,34.34,0.15, 
1200,10.6,11.43,13.89,16.25,21.53,36.13,0.175, 
1200,11.14,12.01,14.61,17.09,22.64,38.01,0.2, 
1200,11.72,12.64,15.36,17.98,23.83,40.04,0.225, 
1200,12.33,13.3,16.17,18.93,25.09,42.19,0.25, 
1200,12.98,14,17.03,19.94,26.44,44.48,0.275, 
1200,13.68,14.75,17.95,21.02,27.89,46.93,0.3, 
1200,14.43,15.56,18.94,22.18,29.44,49.56,0.325, 
1200,15.23,16.43,20,23.43,31.1,52.41,0.35, 
1200,16.09,17.36,21.15,24.77,32.9,55.46,0.375, 
1200,17.03,18.37,22.38,26.23,34.85,58.78,0.4, 
1200,18.04,19.47,23.72,27.81,36.96,62.4,0.425, 
1200,19.14,20.66,25.18,29.52,39.26,66.32,0.45, 
1200,20.34,21.96,26.78,31.4,41.77,70.62,0.475, 
1200,21.66,23.38,28.53,33.46,44.53,75.37,0.5, 
1200,23.11,24.95,30.45,35.73,47.59,80.59,0.525, 
1200,24.71,26.69,32.59,38.25,50.97,86.35,0.55, 
1200,26.5,28.63,34.97,41.06,54.74,92.91,0.575, 
1200,28.5,30.8,37.63,44.21,58.98,100.22,0.6, 
1200,30.76,33.24,40.65,47.77,63.78,108.59,0.625, 
1200,33.33,36.03,44.08,51.82,69.25,118.04,0.65, 
1200,36.28,39.22,48.02,56.48,75.57,129.03,0.675, 
1200,39.69,42.93,52.6,61.91,82.92,141.78,0.7, 
1200,43.7,47.29,57.98,68.29,91.59,156.97,0.725, 
1200,48.48,52.48,64.41,75.92,101.97,175.09,0.75, 
1200,54.28,58.77,72.22,85.2,114.62,197.39,0.775, 
1200,61.45,66.58,81.93,96.74,130.41,225.09,0.8, 
1200,70.57,76.51,94.31,111.5,150.68,261.41,0.825, 
1200,82.58,89.61,110.68,131.07,177.65,309.74,0.85, 
1200,99.13,107.68,133.35,158.28,215.23,378.02,0.875, 
1200,123.47,134.29,166.92,198.7,271.46,481.35,0.9, 
1200,162.91,177.6,221.92,265.3,365.07,656.31,0.925, 
1200,238.5,260.98,329.07,396.02,551.25,1011.1,0.95, 
1200,445.46,490.91,629.23,765.9,1083.53,2038.59,0.975
  ),ncol=8,byrow=T)
# itmp<-itmp[,c(1,2,2,2,3,4,5,6)] # temp, will remove later
  Nmax<-max(itmp[,1])
  Nmin<-min(itmp[,1])
  ncol<-dim(itmp)[2]
  nlevs<-unique(itmp[,1])
  PTmax<-matrix(0,max(Nmax,Nx),45)
  phi<-c(-0.975,itmp[1:44,ncol])
  for(i in 0:(length(nlevs)-1)){
    for(j in 1:44){
      if(i==0){
        PTmax1<-0
	PTmax2<-itmp[j,(pkth+1)]
	ind1<-1
	ind2<-nlevs[1]
      }
      else{
        PTmax1<-itmp[(i*44-44+j),(pkth+1)]
        PTmax2<-itmp[(i*44+j),(pkth+1)]
	ind1<-nlevs[i]
	ind2<-nlevs[i+1]
      }
      PTmax[ind1:ind2,(j+1)]<-seq(PTmax1,PTmax2,length=(ind2-ind1+1))
    }
    PTmax[ind1:ind2,1]<-PTmax[ind1:ind2,2]-(phi[2]-phi[1])*
      (PTmax[ind1:ind2,3]-PTmax[ind1:ind2,2])/(phi[3]-phi[2])
  }
  if(Nx>Nmax)
    for(j in 1:45){
      delta<-(PTmax[nlevs[length(nlevs)],j]-PTmax[nlevs[length(nlevs)-1],j])/
      	     (nlevs[length(nlevs)]-nlevs[length(nlevs)-1])
      PTmax[Nmax:Nx,j]<-seq(from=PTmax[nlevs[length(nlevs)],j],
           length=(Nx-Nmax+1),by=delta)
      }
# PTmax<-PTmax[-1,]
  assign("phi",phi,envir=.GlobalEnv)
  assign("PFmax",PTmax,envir=.GlobalEnv)
}

Pk.PMFT<-function(N){
# if(floor(N)!=N) stop("input data error in Pk")
  Nlt40<- if(N<40) TRUE else FALSE
  Nle500<- if(N<=500) TRUE else FALSE
  
  K<-seq(1,(N-1))
  Kmin<- if(floor((N-1)/2)==(N-1)/2) c(1:floor((N-1)/2),floor((N-1)/2):1)
         else c(1:(floor((N-1)/2)+1),floor((N-1)/2):1)
  W<- floor(if(Nle500) 11*N/50 else 21*N/100)
  A<-abs(1-2*K/N)
  B<-log(N)
  C<-log(B)
  D<-log(log(N+150))
  Q<-c(abs(1-50*K[1:floor(N/2)]/(N*11)),
       abs(1-(50*K[(floor(N/2)+1):(N-1)]-28*N)/(N*11)))
  tmp1<-11/B^3
  tmp2<-abs(2*C^4/(900-9*Kmin))
  flg<-tmp1<tmp2
  S<-tmp1*flg+tmp2*(!flg)
  S[Kmin==100]<-tmp1
  tmp1<-(1-Q^(B*(3-C)/6))
  tmp2<-(1+A^(B*(3-C)/6))
  F<-c(tmp1[1:W],tmp2[(W+1):(N-W-1)],tmp1[(N-W):(N-1)])
  tmp1<-(11.7-S)*B^.01-11.8
  tmp2<-((C+353)*N^.01-340)/200
  v<-c(tmp1[1:W],rep(tmp2,(N-W*2-1)),tmp1[(N-W):(N-1)])

  tmp1<-((64*N^(1/40)+35-C)*F^v)/100
  tmp2<-min(c(1,(98000+6*N)/100000))*F^v
  P0<-c(tmp1[1:W],tmp2[(W+1):(N-W-1)],tmp1[(N-W):(N-1)])

  P<-P0
  if(N>=40){
    L<-floor(1+(317*N^.75-2867)/1000)
    if(N<=100)
      delta<-rep(D^(1/3)*(P0[L+1]-P0[L])+C^3/(N*5),N-1)
    else{
      delta<-rep(NA,N-1)
      delta[1:L]<-(P0[L]-P0[1])*A[1:L]^(C^3)/(L+B-2*C-1)
      delta[(N-L):(N-1)]<-(P0[N-L]-P0[N-1])*A[(N-L):(N-1)]^(C^3)/(L+B-2*C-1)
    }
    P[1:L]<-P0[L]-(L-(1:L))*delta[1:L]
    P[(N-L):(N-1)]<-P0[N-L]-((N-L):(N-1)-N+L)*delta[(N-L):(N-1)]
  }
  return(P)
}

LSmatrix<-function(Y,T,Ic){
  Nx<-length(Y)
  D<-rep(1,Nx)
  X<-t(t(Y))
  D<-cbind(D,T)
  if(!is.na(Ic)) D<-cbind(D,c(rep(0,Ic),rep(1,Nx-Ic)))
  sig<-solve(t(D)%*%D)%*%t(D)%*%X
  fitted<-D%*%sig
  resi<-X-fitted
  SSE<-sum(resi^2)
  oout<-list()
  oout$sig<-as.vector(sig)
  oout$fitted<-as.vector(fitted)
  oout$resi<-as.vector(resi)
  oout$SSE<-SSE
  return(oout)
}

PMFT<-function(Y,T,Pk0){
  N<-length(Y)
  PFx<-(-99999.)
  Fx<-(-99999.)
  KPx<-0
  oout1<-LSmatrix(Y,T,NA)
  for(i in Nmin:(N-Nmin)){
    oout2<-LSmatrix(Y,T,i)
    Fc<-(oout1$SSE-oout2$SSE)*(N-3)/oout2$SSE
    PFc<-Fc*Pk0[i]
    if(PFc>PFx){
      PFx<-PFc
      KPx<-i
      Fx<-Fc
    }
  }
  oout<-list()
  oout$PFx<-PFx
  oout$KPx<-KPx
  oout$Fx<-Fx
  return(oout)
}

PMFxKc<-function(Y,T,I0,I2,Ic){
  Nseg<-(I2-I0)
  Ic1<-Ic-I0
  oout<-list()
  if(Ic>I0){
    Y0<-Y[(I0+1):I2]
    T0<-T[(I0+1):I2]
    oout1<-LSmatrix(Y0,T0,NA)
    oout2<-LSmatrix(Y0,T0,Ic1)
    Fc<-(oout1$SSE-oout2$SSE)*(Nseg-3)/oout2$SSE
    prob<-pf(Fc,1,(Nseg-3))
    Pk0<-Pk.PMFT(Nseg)
    PFc<-Fc*Pk0[Ic1]
    oout$Fc<-Fc
    oout$PFc<-PFc
    oout$prob<-prob
  }
  else{
    oout$Fc<-0
    oout$PFc<-0
    oout$prob<-0
  }
  return(oout)
}

PMFxKxI0I2<-function(Y,T,I0,I2){
  Nmin2<-Nmin*2
  prob<-(-1)
  Ic<-I0
  Nseg<-(I2-I0)
  if(Nseg>=Nmin2){
    Y0<-Y[(I0+1):I2]
    T0<-T[(I0+1):I2]
    Pk0<-Pk.PMFT(Nseg)
    oout<-PMFT(Y0,T0,Pk0)
    Ic<-I0+oout$KPx
    prob<-pf(oout$Fx,1,(Nseg-3))
  }
  oout<-list()
  oout$Ic<-Ic
  oout$prob<-prob
  return(oout)
}

rmCycle<-function(idata){
  tdata<-cbind(idata[,2]*100+idata[,3],idata[,4])
  inds<-sort(unique(tdata[,1]))
  nx<-length(inds)
  mu<-rep(0,nx)
  for(i in 1:nx){
    mu[i]<-mean(tdata[tdata[,1]==inds[i],2],na.rm=T)
    tdata[tdata[,1]==inds[i],2]<-tdata[tdata[,1]==inds[i],2]-mu[i]
  }
  oout<-list()
  oout$EB<-mu
  oout$Base<-tdata[,2]
  return(oout)
}

LSmultiple<-function(Y,T,Ips){
  Nx<-length(Y)
  Ns<-length(Ips)-1
  X<-t(t(Y))
  D<-rep(1,Nx)
  D<-cbind(D,T)
  if(Ns>=1){
    for(i in 1:Ns){
      tmp<-rep(0,Nx)
      tmp[(Ips[i]+1):Ips[i+1]]<-1
      D<-cbind(D,tmp)
    }
  }
  sig<-solve(t(D)%*%D)%*%t(D)%*%X
  fitted<-D%*%sig
  resi<-X-fitted
  SSE<-sum(resi^2)
  oout<-list()
  oout$SSE<-SSE
  oout$fitted<-as.vector(fitted)
  oout$resi<-as.vector(resi)
  oout$sig<-as.vector(sig)
  return(oout)
}

Rphi<-function(Y0,Ips,Ns){
# calculate auto-correlation of given data vector Y0 and breakpoints Ips
# output: cor -- autocorrelation; W -- prewhitenning vector of Y0
#  corl -- lower bound of cor; corh -- upper bound of cor
# if(Ns!=length(Ips)-1) stop("input data length error in Rphi")
  Y<-Y0
  N<-length(Y0)
  mu<-rep(0,Ns+1)
  for(i in 0:Ns){
    I1<- if(i==0) 1 else Ips[i]+1
    I2<-Ips[i+1]
    mu[i+1]<-mean(Y0[I1:I2])
    Y[I1:I2]<-Y0[I1:I2]-mu[i+1]
  }
  cor<-autocorlh(Y,IY0flg)
  W1<-Y
  W2<-Y
  W3<-Y
  W1[2:N]<-Y[2:N]-Y[1:(N-1)]*cor$cor
  W2[2:N]<-Y[2:N]-Y[1:(N-1)]*cor$corl
  W3[2:N]<-Y[2:N]-Y[1:(N-1)]*cor$corh
  W<-IY0flg*W1+(!IY0flg)*Y # if IY0flg==1 (continuous), W1; else Y
  WL<-IY0flg*W2+(!IY0flg)*Y # if IY0flg==1 (continuous), W1; else Y
  WU<-IY0flg*W3+(!IY0flg)*Y # if IY0flg==1 (continuous), W1; else Y
  for(i in 0:Ns){
    I1<- if(i==0) 1 else Ips[i]+1
    I2<-Ips[i+1]
    W[I1:I2]<-W[I1:I2]+mu[i+1]
    WL[I1:I2]<-WL[I1:I2]+mu[i+1]
    WU[I1:I2]<-WU[I1:I2]+mu[i+1]
  }
  t.stat<-abs(cor$cor)*sqrt((N-2)/(1-cor$cor^2))
  pcor<-pt(t.stat,N-1)
  oout<-list()
  oout$cor<-cor$cor
  oout$corl<-cor$corl
  oout$corh<-cor$corh
  oout$pcor<-pcor
  oout$W<-W
  oout$WL<-WL
  oout$WU<-WU
  return(oout)
}

autocorlh<-function(Y,IY){
# calculate autocorrelation of given data vector, using given time vector to
# judge continuouse
  N<-length(Y)
  cnt<-sum(IY)
  m0<-mean(Y,na.rm=T)
  xsd0<-0
  xsd1<-0
  S1<-sum(((Y-m0)^2*IY)[1:(N-1)])
  S2<-sum(((Y-m0)*(c(Y[2:N],0)-m0)*IY)[1:(N-1)])
  cor<-S2/S1
# else stop("too few available data in autocor") 
  z975<-1.96
  z<-.5*log((1+cor)/(1-cor))
  df<-sum(IY[1:(N-1)])
  zl<-z-z975/sqrt(df-3)
  zh<-z+z975/sqrt(df-3)
  cl<-tanh(zl)
  ch<-tanh(zh)
  corl<-min(c(cl,ch))
  corh<-max(c(cl,ch))
  oout<-list()
  oout$cor<-cor
  oout$corl<-corl
  oout$corh<-corh
  return(oout)
}

getPFx95<-function(cor,N){
# if(cor<phi[1]|cor>phi[length(phi)]) stop("input series autocorrelation outbound!")
  if(cor<=phi[1])
    PTx95<-PFmax[N,1]
  else if(cor>=phi[length(phi)])
    PTx95<-PFmax[N,length(phi)]
  else{
    for(i in 1:(length(phi)-1))
      if(cor>phi[i]&cor<phi[i+1]) {
        Kr1<-i
        Kr2<-i+1
        cor1<-phi[i]
        cor2<-phi[i+1]
      }
    tmp1<-PFmax[N,Kr1]
    tmp2<-PFmax[N,Kr2]
    PTx95<-tmp1+(tmp2-tmp1)*(cor-cor1)/(cor2-cor1)
  }
  return(PTx95)
}

PMFxIseg<-function(Y0,Ti,Ips,Iseg){
  Ns<-length(Ips)-1
  N<-length(Y0)
  I0<- if(Iseg==1) 0 else Ips[Iseg-1]
  I3<-Ips[Iseg+1]
  Ic<-Ips[Iseg]
  Nseg<-I3-I0
  Ip0<-Ips[-Iseg]

  otmp<-LSmultipleRed(Y0,Ti,Ips)
  cor<-otmp$cor
  corl<-otmp$corl
  corh<-otmp$corh
  W<-otmp$W
  WL<-otmp$WL
  WU<-otmp$WU
  PFx95<-getPFx95(cor,Nseg)
  PFx95L<-getPFx95(corl,Nseg)
  PFx95U<-getPFx95(corh,Nseg)

  otmp<-LSmultiple(W,Ti,Ip0)
  SSE0.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  otmp<-LSmultiple(W,Ti,Ips)
  SSEf.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
  if(Fx<0){
    Fx<-0
    PFx<-0
    prob<-0
  }
  else prob<-pf(Fx,1,(Nseg-3))

  otmp<-LSmultiple(WL,Ti,Ip0)
  SSE0.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  otmp<-LSmultiple(WL,Ti,Ips)
  SSEf.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  Fx1<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
  probL1<- if(Fx1<0) 0 else pf(Fx1,1,(Nseg-3))

  otmp<-LSmultiple(WU,Ti,Ip0)
  SSE0.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  otmp<-LSmultiple(WU,Ti,Ips)
  SSEf.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  Fx1<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
  probU1<- if(Fx1<0) 0 else pf(Fx1,1,(Nseg-3))

  probL<-min(c(probL1,probU1)); probU<-max(c(probL1,probU1))

  otmp<-LSmultiple(Y0,Ti,Ips)
# fitted<-otmp$fitted
  resi<-otmp$resi
# otmp<-Rphi(resi,Ips,Ns)
# W<-otmp$W+fitted
# cor<-otmp$cor
  SSEf.Iseg<-sum(resi[(I0+1):I3]^2)

  otmp<-LSmultiple(Y0,Ti,Ip0)
  SSE0.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
  Pk0<-Pk.PMFT(Nseg)
  PFx<-Fx*Pk0[Ic-I0]

  oout<-list()
  oout$Fx<-Fx
  oout$PFx<-PFx
  oout$prob<-prob
  oout$probL<-probL
  oout$probU<-probU
  oout$PFx95<-PFx95
  oout$PFx95L<-PFx95L
  oout$PFx95U<-PFx95U
  return(oout)
}

LSmultipleRed<-function(Y0,Ti,Ips){
  Ns<-length(Ips)-1
  N<-length(Y0)
  otmp<-LSmultiple(Y0,Ti,Ips)
  sig<-otmp$sig
  beta<-otmp$sig[2]
  resi<-otmp$resi
  otmp<-autocorlh(resi,IY0flg)
  cor<-otmp$cor
  corl<-otmp$corl
  corh<-otmp$corh
  resi<-resi+beta*Ti
  W1<-resi/(1-cor)
  W2<-c(W1[1],(resi[2:N]-cor*resi[1:(N-1)])/(1-cor))
  W<-c(1,IY0flg[1:(N-1)])*W2+(!c(1,IY0flg[1:(N-1)]))*W1
  otmp<-LSmatrix(W,Ti,NA)
  beta<-otmp$sig[2]
  St0<-sum((Ti-mean(Ti))^2)
  df<-N-2
  sigmaE2<-otmp$SSE/df
  t.stat<-abs(beta)/sqrt(sigmaE2/St0)
  p.tr<-pt(t.stat,df)
  itmp<-Y0-beta*Ti
  mu<-rep(0,Ns+1)
  meanhat<-mu
  for(i in 0:Ns){
    I0<- if(i==0) 1 else Ips[i]+1
    I2<-Ips[i+1]
    mu[i+1]<-mean(itmp[I0:I2])
    meanhat[I0:I2]<-mu[i+1]+beta*Ti[I0:I2]
    resi[I0:I2]<-Y0[I0:I2]-meanhat[I0:I2]
  }
  W1<-resi
  W2<-c(resi[1],resi[2:N]-cor*resi[1:(N-1)])
  W3<-c(resi[1],resi[2:N]-corl*resi[1:(N-1)])
  W4<-c(resi[1],resi[2:N]-corh*resi[1:(N-1)])
  W<-c(1,IY0flg[1:(N-1)])*W2+(!c(1,IY0flg[1:(N-1)]))*W1
  WL<-c(1,IY0flg[1:(N-1)])*W3+(!c(1,IY0flg[1:(N-1)]))*W1
  WU<-c(1,IY0flg[1:(N-1)])*W4+(!c(1,IY0flg[1:(N-1)]))*W1
  for(i in 0:Ns){
    I0<- if(i==0) 1 else Ips[i]+1
    I2<-Ips[i+1]
    W[I0:I2]<-W[I0:I2]+mean(itmp[I0:I2])+beta*Ti[I0:I2]
    WL[I0:I2]<-WL[I0:I2]+mean(itmp[I0:I2])+beta*Ti[I0:I2]
    WU[I0:I2]<-WU[I0:I2]+mean(itmp[I0:I2])+beta*Ti[I0:I2]
  }
  oout<-list()
  oout$W<-W
  oout$WL<-WL
  oout$WU<-WU
  oout$sig<-sig
  oout$cor<-cor
  oout$corl<-corl
  oout$corh<-corh
  oout$resi<-resi
  oout$mu<-mu
  oout$meanhat<-meanhat
  oout$trend<-beta
  oout$p.tr<-p.tr
  return(oout)
}

BCtrans<-function(P,lambda){
  if(sum(P<=0)>0) stop("non-positive value encountered, y<=0!")
  y<-if(lambda==0) log(P) else (P**lambda-1)/lambda
  return(y)
}

IVBCtrans<-function(P,lambda){
  y<-if(lambda==0) exp(P) else (P*lambda+1)**(1/lambda)
  return(y)
}

splineN<-function(x,y,yp1,yp2){
  n<-length(x)
  if(length(y)!=n) stop("input vector length differ")
  y2<-rep(NA,0)
  u<-rep(NA,0)
  if(yp1>1E30){
    y2[1]<-0
    u[1]<-0 }
  else{
    y2[1]<-(-0.5)
    u[1]<-(3/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
  }
  if(n>2) for(i in c(2:(n-1))){
    sig<-(x[i]-x[i-1])/(x[i+1]-x[i-1])
    p<-sig*y2[i-1]+2
    y2[i]<-(sig-1)/p
    u[i]<-(6*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/
          (x[i+1]-x[i-1])-sig*u[i-1])/p
  }
  if(yp2>1E30){
    qn<-0
    un<-0 }
  else{
    qn<-0.5
    un<-(3/(x[n]-x[n-1]))*(yp2-(y[n]-y[n-1])/(x[n]-x[n-1]))
  }
  y2[n]<-(un-qn*u[n-1])/(qn*y2[n-1]+1)
  for(i in c((n-1):1)) y2[i]<-y2[i]*y2[i+1]+u[i]
  return(y2)
}

splintN<-function(xa,ya,y2a,x){
  n<-length(xa)
  if(length(ya)!=n|length(y2a)!=n) stop("input vector length differ")
  klo<-1
  khi<-n
  while((khi-klo)>1){
    k<-ceiling((khi+klo)/2)
    if(xa[k]>x) khi<-k else klo<-k
  }
  h<-xa[khi]-xa[klo]
  if(h==0) stop("bad xa input in splintN")
  a<-(xa[khi]-x)/h
  b<-(x-xa[klo])/h
  y<-a*ya[klo]+b*ya[khi]+((a**3-a)*y2a[klo]+(b**3-b)*y2a[khi])*(h**2)/6
  return(y)
}

getMtrendFdly<-function(idata){
  cmon<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Agu","Sep","Oct","Nov","Dec")
  if(dim(idata)[2]!=4) stop("input data does not contain 4 columns in getMtrendFdly())")
  colnames(idata)<-c("year","month","day","data")
  yrs<-unique(idata[,"year"])
  mmean<-rep(NA,12)
  tdata<-NULL
  for(yth in 1:length(yrs)) {
    itmp<-ori.itable[ori.itable[,1]==yrs[yth],]
    for(mon in 1:12){
      mp<-NA
      mtmp<-itmp[itmp[,2]==mon,]
      if(length(mtmp)>0)
        if(dim(mtmp)[1]>20&sum(is.na(mtmp[,4]))<=3) mp<-sum(idata[idata[,1]==yrs[yth]&idata[,2]==mon,4])
      tdata<-rbind(tdata,c(yrs[yth],mon,mp))
    }
  }
  t1data<-tdata
  for(i in 1:12){
    if(sum(is.na(tdata[tdata[,2]==i,3])==F)<5) print(paste("monthly data too few at",cmon[i]))
    mmean[i]<-mean(tdata[tdata[,2]==i,3],na.rm=T)
    t1data[t1data[,2]==i,3]<-t1data[t1data[,2]==i,3]-mmean[i]
  }
  return(summary(lm(t1data[,3]~c(1:dim(tdata)[1])))$coef[2,1])
}

Read.file<-function(){
  ifname<-tclvalue(tkgetOpenFile())
  if(!nchar(ifname)){
    tkinsert(txt,"end","No file selected in Data Transform!\n")
    return()
  }
  outdirtmp<-strsplit(ifname,"/")[[1]]
  if(length(outdirtmp)<=2){
    curdir<-paste(outdirtmp[1])
    outdir<-paste(outdirtmp[1],"log",sep=":/")
  }
  else{
    curdir<-outdirtmp[1]
    for(i in 2:(length(outdirtmp)-1))
      curdir<-paste(curdir,outdirtmp[i],sep="/")
  }
  # if(!file.exists(outdir)) dir.create(outdir)
  ofname<-outdirtmp[length(outdirtmp)]
  if(is.csv(ofname)) csv<-1
  else csv<-0
  if(csv){
    itmp<-try(read.table(ifname,header=F,
                         col.names=c("year","month","day","prcp","tmax","tmin"),
                         sep=",",na.strings="-99.9",colClasses=rep("real",6)),
              silent=T)
    if(inherits(itmp,"try-error")){
      ErrorMSG<<-geterrmessage()
      tkinsert(txt,"end",paste(ErrorMSG,"\n"))
      return()
    }
    else iidata<-itmp
  }
  else{
    itmp<-try(read.table(ifname,header=F,
                         col.names=c("year","month","day","prcp","tmax","tmin"),
                         na.strings="-99.9",colClasses=rep("real",6)),
              silent=T)
    if(inherits(itmp,"try-error")){
      ErrorMSG<<-geterrmessage()
      tkinsert(txt,"end",paste(ErrorMSG,"\n"))
      return()
    }
    else iidata<-itmp
  }
  if(ncol(iidata)!=6){
    ErrorMSG<-paste(ifname,"has",ncol(iidata),"columns. The number of columns should be 6\n")
    tkinsert(txt,"end",paste(ErrorMSG,"\n"))
    return()
  }
  nlen<-dim(iidata)[1]
  syear<-iidata[1,1]
  eyear<-iidata[nlen,1]
  smonth<-iidata[1,2]
  emonth<-iidata[nlen,2]
  if(eyear<(syear+1)) {
    ErrorMSG<-paste("Time series",ifname,"too short for Transform Data\n")
    return(-1)
  }
  nyrs<-eyear-syear+1
  vars<-c("tmax","tmin","prcp")
  ivars<-length(vars)
  
  tkinsert(txt,"end",paste("Data from:",syear,"to:",eyear,sep=" "))
  tkinsert(txt,"end","\n")
  
  for(i in 1:ivars){
    mdata<-NULL
    if(vars[i]=="prcp") mdata1mm<-NULL
    tmpdata<-iidata[iidata[,"year"]==syear,]
    for(k in smonth:12){
      if(sum(is.na(tmpdata[tmpdata[,"month"]==k,vars[i]]))<=3)
        if(vars[i]=="prcp"){
          itmp<-sum(tmpdata[tmpdata[,"month"]==k,vars[i]],na.rm=T)
          itmp1mm<-sum(tmpdata[tmpdata[,"month"]==k&
                                 tmpdata[,vars[i]]>=1,vars[i]],na.rm=T)
        }
      else
        itmp<-mean(tmpdata[tmpdata[,"month"]==k,vars[i]],na.rm=T)
      else{
        itmp<-NA
        itmp1mm<-NA
      }
      mdata<-rbind(mdata,c(syear,k,0,itmp))
      if(vars[i]=="prcp") mdata1mm<-rbind(mdata1mm,c(syear,k,0,itmp1mm))
    }
    for(j in (syear+1):(eyear-1)){
      year<-j
      tmpdata<-iidata[iidata[,"year"]==year,]
      for(k in 1:12){
        if(sum(is.na(tmpdata[tmpdata[,"month"]==k,vars[i]]))<=3)
          if(vars[i]=="prcp"){
            itmp<-sum(tmpdata[tmpdata[,"month"]==k,vars[i]],na.rm=T)
            itmp1mm<-sum(tmpdata[tmpdata[,"month"]==k&
                                   tmpdata[,vars[i]]>=1,vars[i]],na.rm=T)
          }
        else
          itmp<-mean(tmpdata[tmpdata[,"month"]==k,vars[i]],na.rm=T)
        else{
          itmp<-NA
          itmp1mm<-NA
        }
        mdata<-rbind(mdata,c(year,k,0,itmp))
        if(vars[i]=="prcp") mdata1mm<-rbind(mdata1mm,c(year,k,0,itmp1mm))
      }
    }
    tmpdata<-iidata[iidata[,"year"]==eyear,]
    for(k in 1:emonth){
      if(sum(is.na(tmpdata[tmpdata[,"month"]==k,vars[i]]))<=3)
        if(vars[i]=="prcp"){
          itmp<-sum(tmpdata[tmpdata[,"month"]==k,vars[i]],na.rm=T)
          itmp1mm<-sum(tmpdata[tmpdata[,"month"]==k&
                                 tmpdata[,vars[i]]>=1,vars[i]],na.rm=T)
        }
      else
        itmp<-mean(tmpdata[tmpdata[,"month"]==k,vars[i]],na.rm=T)
      else{
        itmp<-NA
        itmp1mm<-NA
      }
      mdata<-rbind(mdata,c(eyear,k,0,itmp))
      if(vars[i]=="prcp") mdata1mm<-rbind(mdata1mm,c(eyear,k,0,itmp1mm))
    }
    if(vars[i]=="prcp") {prcp<-mdata; prcp1mm<-mdata1mm}
    else if(vars[i]=="tmax") tmax<-mdata
    else tmin<-mdata
  }
  logprcp<-prcp
  if(min(prcp[,4],na.rm=T)>0) logprcp[,4]<-log(prcp[,4])
  else logprcp[,4]<-log(prcp[,4]+1)
  logprcp1mm<-prcp1mm
  if(min(prcp1mm[,4],na.rm=T)>0) logprcp1mm[,4]<-log(prcp1mm[,4])
  else logprcp1mm[,4]<-log(prcp1mm[,4]+1)
  otmp<-strsplit(ofname,"\\.")[[1]]
  if(length(otmp)>1){
    ind<-length(otmp)-1
    ofmax2<-paste(otmp[ind],"_tmaxMLY",sep="")
    ofmin2<-paste(otmp[ind],"_tminMLY",sep="")
    ofprcp2<-paste(otmp[ind],"_prcpMLY",sep="")
    ofprcpL2<-paste(otmp[ind],"_LogprcpMLY",sep="")
    ofprcp1mm2<-paste(otmp[ind],"_prcpMLY1mm",sep="")
    ofprcp1mmL2<-paste(otmp[ind],"_LogprcpMLY1mm",sep="")
    ofmaxD2<-paste(otmp[ind],"_tmaxDLY",sep="")
    ofminD2<-paste(otmp[ind],"_tminDLY",sep="")
    ofprcpD2<-paste(otmp[ind],"_prcpDLY",sep="")
    if(ind==1){
      ofmax<-ofmax2
      ofmin<-ofmin2
      ofprcp<-ofprcp2
      ofprcpL<-ofprcpL2
      ofprcp1mm<-ofprcp1mm2
      ofprcp1mmL<-ofprcp1mmL2
      ofmaxD<-ofmaxD2
      ofminD<-ofminD2
      ofprcpD<-ofprcpD2
    }
    else{
      ofmax<-otmp[1]
      ofmin<-otmp[1]
      ofprcp<-otmp[1]
      ofprcpL<-otmp[1]
      ofprcp1mm<-otmp[1]
      ofprcp1mmL<-otmp[1]
      ofmaxD<-otmp[1]
      ofminD<-otmp[1]
      ofprcpD<-otmp[1]
    }
    for(i in 2:length(otmp)){
      if(i==ind){
        ofmax<-paste(ofmax,ofmax2,sep=".")
        ofmin<-paste(ofmin,ofmin2,sep=".")
        ofprcp<-paste(ofprcp,ofprcp2,sep=".")
        ofprcpL<-paste(ofprcpL,ofprcpL2,sep=".")
        ofprcp1mm<-paste(ofprcp1mm,ofprcp1mm2,sep=".")
        ofprcp1mmL<-paste(ofprcp1mmL,ofprcp1mmL2,sep=".")
        ofmaxD<-paste(ofmaxD,ofmaxD2,sep=".")
        ofminD<-paste(ofminD,ofminD2,sep=".")
        ofprcpD<-paste(ofprcpD,ofprcpD2,sep=".")
      }
      else{
        ofmax<-paste(ofmax,otmp[i],sep=".")
        ofmin<-paste(ofmin,otmp[i],sep=".")
        ofprcp<-paste(ofprcp,otmp[i],sep=".")
        ofprcpL<-paste(ofprcpL,otmp[i],sep=".")
        ofprcp1mm<-paste(ofprcp1mm,otmp[i],sep=".")
        ofprcp1mmL<-paste(ofprcp1mmL,otmp[i],sep=".")
        ofmaxD<-paste(ofmaxD,otmp[i],sep=".")
        ofminD<-paste(ofminD,otmp[i],sep=".")
        ofprcpD<-paste(ofprcpD,otmp[i],sep=".")
      }
    }
  }
  else{
    ofmax<-paste(otmp,"_tmaxMLY",sep="")
    ofmin<-paste(otmp,"_tminMLY",sep="")
    ofprcp<-paste(otmp,"_prcpMLY",sep="")
    ofprcpL<-paste(otmp,"_LogprcpMLY",sep="")
    ofprcp1mm<-paste(otmp,"_prcpMLY1mm",sep="")
    ofprcp1mmL<-paste(otmp,"_LogprcpMLY1mm",sep="")
    ofmaxD<-paste(otmp,"_tmaxDLY",sep="")
    ofminD<-paste(otmp,"_tminDLY",sep="")
    ofprcpD<-paste(otmp,"_prcpDLY",sep="")
  }
  ofmax<-paste(curdir,ofmax,sep="/")
  ofmin<-paste(curdir,ofmin,sep="/")
  ofprcp<-paste(curdir,ofprcp,sep="/")
  ofprcpL<-paste(curdir,ofprcpL,sep="/")
  ofprcp1mm<-paste(curdir,ofprcp1mm,sep="/")
  ofprcp1mmL<-paste(curdir,ofprcp1mmL,sep="/")
  ofmaxD<-paste(curdir,ofmaxD,sep="/")
  ofminD<-paste(curdir,ofminD,sep="/")
  ofprcpD<-paste(curdir,ofprcpD,sep="/")
  rownames(tmax)<-NULL
  rownames(tmin)<-NULL
  rownames(prcp)<-NULL
  write.table(tmax,file=ofmax,sep=" ",na=MissingStr,col.names=F,row.names=F)
  write.table(tmin,file=ofmin,sep=" ",na=MissingStr,col.names=F,row.names=F)
  write.table(prcp,file=ofprcp,sep=" ",na=MissingStr,col.names=F,row.names=F)
  write.table(logprcp,file=ofprcpL,sep=" ",na=MissingStr,col.names=F,row.names=F)
  write.table(prcp1mm,file=ofprcp1mm,sep=" ",na=MissingStr,col.names=F,row.names=F)
  write.table(logprcp1mm,file=ofprcp1mmL,sep=" ",na=MissingStr,col.names=F,row.names=F)
  write.table(iidata[,c("year","month","day","prcp")],file=ofprcpD,
              na=MissingStr,col.names=F,row.names=F)
  write.table(iidata[,c("year","month","day","tmax")],file=ofmaxD,
              na=MissingStr,col.names=F,row.names=F)
  write.table(iidata[,c("year","month","day","tmin")],file=ofminD,
              na=MissingStr,col.names=F,row.names=F)
  tkinsert(txt,"end","Data transform finished, monthly series output:\n")
  tkinsert(txt,"end",paste(ofmax,"\n"))
  tkinsert(txt,"end",paste(ofmin,"\n"))
  tkinsert(txt,"end",paste(ofprcp,"\n"))
  tkinsert(txt,"end",paste(ofprcpL,"\n"))
  tkinsert(txt,"end",paste(ofprcp1mm,"\n"))
  tkinsert(txt,"end",paste(ofprcp1mmL,"\n"))
  tkinsert(txt,"end","Daily series output:\n")
  tkinsert(txt,"end",paste(ofmaxD,"\n"))
  tkinsert(txt,"end",paste(ofminD,"\n"))
  tkinsert(txt,"end",paste(ofprcpD,"\n\n"))
  return(0)
}

str40<-function(vari.name,vari.value){
  olen<-40
  b20<-"                    "
  if(!exists(vari.name,env=sys.parent())) 
    otmp<-paste(b20,"                  NA",sep="")
  else{
    if(nchar(vari.value)>olen) otmp<-paste("...",substr(vari.value,
                                                        nchar(vari.value)-olen+4,nchar(vari.value)),sep="")
    else{
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(vari.value),olen)<-vari.value
    }
  }
  return(otmp)
}

OnFindU.dlyPrcp<-function(){
  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No file selected in FindU.dlyPrcp!\n\n")
      return()
    }
    otmp<-str40("ifname",ifname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=2,sticky="e")
    outdirtmp<-strsplit(ifname,"/")[[1]]
    if(length(outdirtmp)<=2){
      curdir<-paste(outdirtmp[1])
      outdir<-paste(outdirtmp[1],"output",sep=":/")
    }
    else{
      curdir<-outdirtmp[1]
      for(i in 2:(length(outdirtmp)-1))
        curdir<-paste(curdir,outdirtmp[i],sep="/")
      outdir<-paste(curdir,"output",sep="/")
    }
    setwd(curdir)
    if(!file.exists(outdir)) dir.create(outdir)
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    if(length(itmp)<2) ofbody<-itmp[1]
    else{
      ofbody<-itmp[1]
      if(length(itmp)>2) for(i in 2:(length(itmp)-1))
        ofbody<-paste(ofbody,itmp[i],sep="/")
    }
    ofname<-paste(outdir,ofbody,sep="/")
    assign("curdir",curdir,envir=.GlobalEnv)
    assign("outdir",outdir,envir=.GlobalEnv)
    assign("ifname",ifname,envir=.GlobalEnv)
    assign("ofbody",ofbody,envir=.GlobalEnv)
    assign("ofname",ofname,envir=.GlobalEnv)
  }
  tt<-tktoplevel()
  button.chg1<-tkbutton(tt,text="Change",command=getfile1)
  tkwm.title(tt,"FindU.dlyPrcp")
  oifname<-str40("ifname",ifname)
  fontLable<-tkfont.create(family="times",size=20,weight="bold")
  #tkgrid(tklabel(tt,text="!!Do not",font=fontLable),sticky="e",column=1,row=1)
  #tkgrid(tklabel(tt,text="choose daily precipitation data!!",
  #                font=fontLable),row=1,column=2)
  tkgrid(tklabel(tt,text="Input Data filename:"),column=1,row=2,sticky="w")
  tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=2,sticky="e")
  tkgrid(button.chg1,column=3,row=2)
  
  OnOk2<-function(){
    oflg<-1
    if(!file.exists(ifname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Data file ",ifname," does not exist!\n",sep="")
    }
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      outdirtmp<-strsplit(ifname,"/")[[1]]
      if(length(outdirtmp)<=2){
        curdir<-paste(outdirtmp[1])
        outdir<-paste(outdirtmp[1],"output",sep=":/")
      }
      else{
        curdir<-outdirtmp[1]
        for(i in 2:(length(outdirtmp)-1))
          curdir<-paste(curdir,outdirtmp[i],sep="/")
        outdir<-paste(curdir,"output",sep="/")
      }
      setwd(curdir)
      if(!file.exists(outdir)) dir.create(outdir)
      ofname<-outdirtmp[length(outdirtmp)]
      itmp<-strsplit(ofname,"\\.")[[1]]
      if(length(itmp)<2) ofbody<-itmp[1]
      else{
        ofbody<-itmp[1]
        if(length(itmp)>2) for(i in 2:(length(itmp)-1))
          ofbody<-paste(ofbody,itmp[i],sep="/")
      }
      ofname<-paste(outdir,ofbody,sep="/")
      assign("curdir",curdir,envir=.GlobalEnv)
      assign("outdir",outdir,envir=.GlobalEnv)
      assign("ofbody",ofbody,envir=.GlobalEnv)
      assign("ofname",ofname,envir=.GlobalEnv)
      if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
        tkinsert(txt,"end","P_lev setting error, reset P_lev...")
        return()
      }
      itmp<-FindU.dlyPrcp(InSeries=ifname,output=ofname,
                  MissingValueCode=MissingStr,p.lev=as.numeric(PlevStr),
                  Iadj=as.numeric(AdjStr),Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),
                  pthr=as.numeric(pthrstr),GUI=TRUE)
      if(itmp<0){
        tkinsert(txt,"end",ErrorMSG)
        return()
      }
      else{
        UIpsName1<-paste(ofname,"_1Cs.txt",sep="")
        UIpsName<-paste(ofname,"_mCs.txt",sep="")
        file.copy(from=UIpsName1,to=UIpsName,overwrite=TRUE)
        assign("iDIpsName",UIpsName,envir=.GlobalEnv)
        oact<-str40("ofbody",ofbody)
        ocurdir<-str40("curdir",curdir)
        ooutdir<-str40("outdir",outdir)
        if(exists("ofref")) rm("ofref",envir=.GlobalEnv)
        b20<-"                    "
        oref<-paste(b20,"                  NA",sep="")
        tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=7,sticky="e")
        #tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=8,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9,sticky="e")
        tkinsert(txt,"end","FindU.dlyPrcp finished successfully...\n")
        tkinsert(txt,"end",paste("Modify",UIpsName,"for further calculation...\n\n"))
      }
      tkdestroy(tt)
      tkfocus(main)
      return()
    }
    
  }
  Ok2.but<-tkbutton(tt,text="   OK   ",command=OnOk2)
  tkgrid(Ok2.but,column=3,row=3)
  tkfocus(tt)
}

OnFindUD.dlyPrcp<-function(){
  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No file selected in FindUD.dlyPrcp!\n\n")
      return()
    }
    otmp<-str40("ifname",ifname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=2,sticky="e")
    outdirtmp<-strsplit(ifname,"/")[[1]]
    if(length(outdirtmp)<=2){
      curdir<-paste(outdirtmp[1])
      outdir<-paste(outdirtmp[1],"output",sep=":/")
    }
    else{
      curdir<-outdirtmp[1]
      for(i in 2:(length(outdirtmp)-1))
        curdir<-paste(curdir,outdirtmp[i],sep="/")
      outdir<-paste(curdir,"output",sep="/")
    }
    setwd(curdir)
    if(!file.exists(outdir)) dir.create(outdir)
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    if(length(itmp)<2) ofbody<-itmp[1]
    else{
      ofbody<-itmp[1]
      if(length(itmp)>2) for(i in 2:(length(itmp)-1))
        ofbody<-paste(ofbody,itmp[i],sep="/")
    }
    ofname<-paste(outdir,ofbody,sep="/")
    assign("curdir",curdir,envir=.GlobalEnv)
    assign("outdir",outdir,envir=.GlobalEnv)
    assign("ifname",ifname,envir=.GlobalEnv)
    assign("ofbody",ofbody,envir=.GlobalEnv)
    assign("ofname",ofname,envir=.GlobalEnv)
  }
  getfile2<-function(){
    if(!exists("iDIpsName")) iDIpsName<-tclvalue(tkgetOpenFile())
    else iDIpsName<-tclvalue(tkgetOpenFile(initialfile=iDIpsName))
    otmp<-str40("iDIpsName",iDIpsName)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=3,sticky="e")
    assign("iDIpsName",iDIpsName,env=.GlobalEnv)
  }
  tt<-tktoplevel()
  button.chg1<-tkbutton(tt,text="Change",command=getfile1)
  button.chg2<-tkbutton(tt,text="Change",command=getfile2)
  tkwm.title(tt,"FindUD.dlyPrcp")
  oifname<-str40("ifname",ifname)
  oIpsName<-str40("iDIpsName",iDIpsName)
  fontLable<-tkfont.create(family="times",size=20,weight="bold")
  #tkgrid(tklabel(tt,text="!!Do not",font=fontLable),sticky="e",column=1,row=1)
  #tkgrid(tklabel(tt,text="choose daily precipitation data!!",
  #               font=fontLable),row=1,column=2)
  tkgrid(tklabel(tt,text="Input Data filename:"),column=1,row=2,sticky="w")
  tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=2,sticky="e")
  tkgrid(button.chg1,column=3,row=2)
  
  tkgrid(tklabel(tt,text="Input changepoints filename:"),column=1,row=3,sticky="w")
  tkgrid(tklabel(tt,text=oIpsName,width=40),column=2,row=3,sticky="e")
  tkgrid(button.chg2,column=3,row=3)
  
  OnOk2<-function(){
    oflg<-1
    if(!file.exists(ifname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Data file ",ifname," does not exist!\n",sep="")
    }
    if(!file.exists(iDIpsName)) {
      oflg<-0
      GuiErrorMSG<-paste(GuiErrorMSG,"Input changepoint file ",iDIpsName,
                         " does not exist!\n",sep="")
    }
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
        tkinsert(txt,"end","P_lev setting error, reset P_lev...")
        return()
      }
      outdirtmp<-strsplit(ifname,"/")[[1]]
      if(length(outdirtmp)<=2){
        curdir<-paste(outdirtmp[1])
        outdir<-paste(outdirtmp[1],"output",sep=":/")
      }
      else{
        curdir<-outdirtmp[1]
        for(i in 2:(length(outdirtmp)-1))
          curdir<-paste(curdir,outdirtmp[i],sep="/")
        outdir<-paste(curdir,"output",sep="/")
      }
      setwd(curdir)
      if(!file.exists(outdir)) dir.create(outdir)
      ofname<-outdirtmp[length(outdirtmp)]
      itmp<-strsplit(ofname,"\\.")[[1]]
      if(length(itmp)<2) ofbody<-itmp[1]
      else{
        ofbody<-itmp[1]
        if(length(itmp)>2) for(i in 2:(length(itmp)-1))
          ofbody<-paste(ofbody,itmp[i],sep="/")
      }
      ofname<-paste(outdir,ofbody,sep="/")
      assign("curdir",curdir,envir=.GlobalEnv)
      assign("outdir",outdir,envir=.GlobalEnv)
      assign("ofbody",ofbody,envir=.GlobalEnv)
      assign("ofname",ofname,envir=.GlobalEnv)
      itmp<-FindUD.dlyPrcp(InSeries=ifname,output=ofname,InCs=iDIpsName,
                   MissingValueCode=MissingStr,p.lev=as.numeric(PlevStr),
                   Iadj=as.numeric(AdjStr),Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),
                   pthr=as.numeric(pthrstr),GUI=TRUE)
      if(itmp<0){
        tkinsert(txt,"end",ErrorMSG)
        return()
      }
      else{
        UDIpsName1<-paste(ofname,"_pCs.txt",sep="")
        UDIpsName<-paste(ofname,"_mCs.txt",sep="")
        file.copy(from=UDIpsName1,to=UDIpsName,overwrite=TRUE)
        assign("iDIpsName",UDIpsName,envir=.GlobalEnv)
        oact<-str40("ofbody",ofbody)
        ocurdir<-str40("curdir",curdir)
        ooutdir<-str40("outdir",outdir)
        if(exists("ofref")) rm("ofref",envir=.GlobalEnv)
        b20<-"                    "
        oref<-paste(b20,"                  NA",sep="")
        tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=7,sticky="e")
        #tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=7,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9,sticky="e")
        tkinsert(txt,"end","FindUD.dlyPrcp finished successfully...\n")
        tkinsert(txt,"end",paste("Modify",UDIpsName,"for further calculation...\n\n"))
      }
      tkdestroy(tt)
      tkfocus(main)
      return()
    }
  }
  Ok2.but<-tkbutton(tt,text="   OK   ",command=OnOk2)
  tkgrid(Ok2.but,column=3,row=4)
  tkfocus(tt)
}

OnStepSize.dlyprcp<-function(){
  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No file selected in StepSize.dlyprcp!\n\n")
      return()
    }
    otmp<-str40("ifname",ifname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=2,sticky="e")
    outdirtmp<-strsplit(ifname,"/")[[1]]
    if(length(outdirtmp)<=2){
      curdir<-paste(outdirtmp[1])
      outdir<-paste(outdirtmp[1],"output",sep=":/")
    }
    else{
      curdir<-outdirtmp[1]
      for(i in 2:(length(outdirtmp)-1))
        curdir<-paste(curdir,outdirtmp[i],sep="/")
      outdir<-paste(curdir,"output",sep="/")
    }
    setwd(curdir)
    if(!file.exists(outdir)) dir.create(outdir)
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    if(length(itmp)<2) ofbody<-itmp[1]
    else{
      ofbody<-itmp[1]
      if(length(itmp)>2) for(i in 2:(length(itmp)-1))
        ofbody<-paste(ofbody,itmp[i],sep="/")
    }
    ofname<-paste(outdir,ofbody,sep="/")
    assign("curdir",curdir,envir=.GlobalEnv)
    assign("outdir",outdir,envir=.GlobalEnv)
    assign("ifname",ifname,envir=.GlobalEnv)
    assign("ofbody",ofbody,envir=.GlobalEnv)
    assign("ofname",ofname,envir=.GlobalEnv)
  }
  getfile2<-function(){
    if(!exists("iDIpsName")) iDIpsName<-tclvalue(tkgetOpenFile())
    else iDIpsName<-tclvalue(tkgetOpenFile(initialfile=iDIpsName))
    otmp<-str40("iDIpsName",iDIpsName)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=3,sticky="e")
    assign("iDIpsName",iDIpsName,env=.GlobalEnv)
  }
  tt<-tktoplevel()
  button.chg1<-tkbutton(tt,text="Change",command=getfile1)
  button.chg2<-tkbutton(tt,text="Change",command=getfile2)
  tkwm.title(tt,"StepSize.dlyprcp")
  oifname<-str40("ifname",ifname)
  oIpsName<-str40("iDIpsName",iDIpsName)
  fontLable<-tkfont.create(family="times",size=20,weight="bold")
  tkgrid(tklabel(tt,text="Input Data filename:"),column=1,row=2,sticky="w")
  tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=2,sticky="e")
  tkgrid(button.chg1,column=3,row=2)
  
  tkgrid(tklabel(tt,text="Input changepoints filename:"),column=1,row=3,sticky="w")
  tkgrid(tklabel(tt,text=oIpsName,width=40),column=2,row=3,sticky="e")
  tkgrid(button.chg2,column=3,row=3)
  
  OnOk2<-function(){
    oflg<-1
    if(!file.exists(ifname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Data file ",ifname," does not exist!\n",sep="")
    }
    if(!file.exists(iDIpsName)) {
      oflg<-0
      GuiErrorMSG<-paste(GuiErrorMSG,"Input changepoint file ",iDIpsName,
                         " does not exist!\n",sep="")
    }
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
        tkinsert(txt,"end","P_lev setting error, reset P_lev...")
        return()
      }
      outdirtmp<-strsplit(ifname,"/")[[1]]
      if(length(outdirtmp)<=2){
        curdir<-paste(outdirtmp[1])
        outdir<-paste(outdirtmp[1],"output",sep=":/")
      }
      else{
        curdir<-outdirtmp[1]
        for(i in 2:(length(outdirtmp)-1))
          curdir<-paste(curdir,outdirtmp[i],sep="/")
        outdir<-paste(curdir,"output",sep="/")
      }
      setwd(curdir)
      if(!file.exists(outdir)) dir.create(outdir)
      ofname<-outdirtmp[length(outdirtmp)]
      itmp<-strsplit(ofname,"\\.")[[1]]
      if(length(itmp)<2) ofbody<-itmp[1]
      else{
        ofbody<-itmp[1]
        if(length(itmp)>2) for(i in 2:(length(itmp)-1))
          ofbody<-paste(ofbody,itmp[i],sep="/")
      }
      ofname<-paste(outdir,ofbody,sep="/")
      assign("curdir",curdir,envir=.GlobalEnv)
      assign("outdir",outdir,envir=.GlobalEnv)
      assign("ofbody",ofbody,envir=.GlobalEnv)
      assign("ofname",ofname,envir=.GlobalEnv)
      itmp<-StepSize.dlyPrcp(InSeries=ifname,output=ofname,InCs=iDIpsName,
                     MissingValueCode=MissingStr,p.lev=as.numeric(PlevStr),
                     Iadj=as.numeric(AdjStr),Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),
                     pthr=as.numeric(pthrstr),GUI=TRUE)
      if(itmp<0){
        tkinsert(txt,"end",ErrorMSG)
        return()
      }
      else{
        UIpsName1<-paste(ofname,"_fCs.txt",sep="")
        UIpsName<-paste(ofname,"_mCs.txt",sep="")
        file.copy(from=UIpsName1,to=UIpsName,overwrite=TRUE)
        oact<-str40("ofbody",ofbody)
        ocurdir<-str40("curdir",curdir)
        ooutdir<-str40("outdir",outdir)
        if(exists("ofref")) rm("ofref",envir=.GlobalEnv)
        b20<-"                    "
        oref<-paste(b20,"                  NA",sep="")
        tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=6,sticky="e")
        tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=7,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9,sticky="e")
        tkinsert(txt,"end","StepSize.dlyprcp finished successfully...\n")
        tkinsert(txt,"end",paste("Final output at ",outdir,"/",ofbody,"_*\n\n",sep=""))
      }
      tkdestroy(tt)
      tkfocus(main)
      return()
    }
  }
  Ok2.but<-tkbutton(tt,text="   OK   ",command=OnOk2)
  tkgrid(Ok2.but,column=3,row=4)
  tkfocus(tt)
}

OnadjDLYp<-function(){
  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No file selected in adjDLYp!\n\n")
      return()
    }
    otmp<-str40("ifname",ifname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=2,sticky="e")
    outdirtmp<-strsplit(ifname,"/")[[1]]
    if(length(outdirtmp)<=2){
      curdir<-paste(outdirtmp[1])
      outdir<-paste(outdirtmp[1],"output",sep=":/")
    }
    else{
      curdir<-outdirtmp[1]
      for(i in 2:(length(outdirtmp)-1))
        curdir<-paste(curdir,outdirtmp[i],sep="/")
      outdir<-paste(curdir,"output",sep="/")
    }
    setwd(curdir)
    if(!file.exists(outdir)) dir.create(outdir)
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    if(length(itmp)<2) ofbody<-itmp[1]
    else{
      ofbody<-itmp[1]
      if(length(itmp)>2) for(i in 2:(length(itmp)-1))
        ofbody<-paste(ofbody,itmp[i],sep="/")
    }
    ofname<-paste(outdir,ofbody,sep="/")
    assign("curdir",curdir,envir=.GlobalEnv)
    assign("outdir",outdir,envir=.GlobalEnv)
    assign("ifname",ifname,envir=.GlobalEnv)
    assign("ofbody",ofbody,envir=.GlobalEnv)
    assign("ofname",ofname,envir=.GlobalEnv)
  }
  getfile2<-function(){
    if(!exists("iDIpsName")) iDIpsName<-tclvalue(tkgetOpenFile())
    else iDIpsName<-tclvalue(tkgetOpenFile(initialfile=iDIpsName))
    otmp<-str40("iDIpsName",iDIpsName)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=3,sticky="e")
    assign("iDIpsName",iDIpsName,env=.GlobalEnv)
  }
  tt<-tktoplevel()
  button.chg1<-tkbutton(tt,text="Change",command=getfile1)
  button.chg2<-tkbutton(tt,text="Change",command=getfile2)
  tkwm.title(tt,"QMadjDLYp")
  oifname<-str40("ifname",ifname)
  oIpsName<-str40("iDIpsName",iDIpsName)
  fontLable<-tkfont.create(family="times",size=20,weight="bold")
  #tkgrid(tklabel(tt,text="!!Do not",font=fontLable),sticky="e",column=1,row=1)
  #tkgrid(tklabel(tt,text="choose daily precipitation data!!",
  #               font=fontLable),row=1,column=2)
  tkgrid(tklabel(tt,text="Input Data filename:"),column=1,row=2,sticky="w")
  tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=2,sticky="e")
  tkgrid(button.chg1,column=3,row=2)
  
  tkgrid(tklabel(tt,text="Input changepoints filename:"),column=1,row=3,sticky="w")
  tkgrid(tklabel(tt,text=oIpsName,width=40),column=2,row=3,sticky="e")
  tkgrid(button.chg2,column=3,row=3)
  
  rb1 <- tkradiobutton(tt)
  rb2 <- tkradiobutton(tt)
  rb3 <- tkradiobutton(tt)
  rbValue<-tclVar("1")
  tkconfigure(rb1,variable=rbValue,value="1")
  tkconfigure(rb2,variable=rbValue,value="4")
  tkconfigure(rb3,variable=rbValue,value="12")
  tkgrid(tklabel(tt,text="Choose one of the following    "),column=1,row=4,sticky="w")
  tkgrid(tklabel(tt,text="No seasonality in trend or distribution"),column=1,row=5,sticky="w")
  tkgrid(tklabel(tt,text="Seasonal trends and distributions (DJF,MAM...)"),column=1,row=6,sticky="w")
  tkgrid(tklabel(tt,text="Monthly trends and distributions (Jan,Feb...) "),column=1,row=7,sticky="w")
  tkgrid(rb1,column=2,row=5)
  tkgrid(rb2,column=2,row=6)
  tkgrid(rb3,column=2,row=7)
  
  OnOk2<-function(){
    oflg<-1
    if(!file.exists(ifname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Data file ",ifname," does not exist!\n",sep="")
    }
    if(!file.exists(iDIpsName)) {
      oflg<-0
      GuiErrorMSG<-paste(GuiErrorMSG,"Input changepoint file ",iDIpsName,
                         " does not exist!\n",sep="")
    }
    rbVal<-as.character(tclvalue(rbValue))
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      if(rbVal=="1"){
        itmp<-adjDLYp(InSeries=ifname,output=ofname,InCs=iDIpsName,
                                MissingValueCode=MissingStr,Iadj=as.numeric(AdjStr),
                                Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),GUI=TRUE)
        if(itmp<0){
          tkinsert(txt,"end",ErrorMSG)
          tkdestroy(tt)
          tkfocus(main)
          return()
        }
      }
      else{
        if(rbVal=="4"){
          Snames<-c("DJF","MAM","JJA","SON")
          Sstrts<-c(1201,301,601,901)
          Sends<-c(229,531,831,1130)
          Nseas<-4
        }
        else if(rbVal=="12"){
          Snames<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
          Sstrts<-c((1:12)*100+1)
          Sends<-c(131,229,331,430,531,630,731,831,930,1031,1130,1231)
          Nseas<-12
        }
        if(is.csv(ifname)) idata<-read.csv(ifname)
        else idata<-read.table(ifname)
        mmdd<-idata[,2]*100+idata[,3]
        for(ith in 1:Nseas){
          if(Sstrts[ith]>Sends[ith])
            odata<-idata[mmdd>=Sstrts[ith]|mmdd<=Sends[ith],]
          else
            odata<-idata[mmdd>=Sstrts[ith]&mmdd<=Sends[ith],]
          iftmp<-paste(ofname,"_",Snames[ith],".dat",sep="")
          oftmp<-paste(ofname,"_",Snames[ith],sep="")
          write.table(odata,file=iftmp,col.names=F,row.names=F)
          itmp<-adjDLYp(InSeries=iftmp,output=oftmp,InCs=iDIpsName,
                                  MissingValueCode=MissingStr,Iadj=as.numeric(AdjStr),
                                  Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),GUI=TRUE)
          if(itmp<0){
            tkinsert(txt,"end",ErrorMSG)
            tkdestroy(tt)
            tkfocus(main)
            return()
          }
        }
        ofileSout<-paste(ofname,"_adjDLYpstat.txt",sep="")
        ofileAout<-paste(ofname,"_adjDLYp.dat",sep="")
        if(file.exists(ofileSout)) file.remove(ofileSout)
        odata<-NULL
        for(ith in 1:Nseas){
          ifileSout<-paste(ofname,"_",Snames[ith],"_adjDLYpstat.txt",sep="")
          ifileAout<-paste(ofname,"_",Snames[ith],"_adjDLYp.dat",sep="")
          iftmp<-paste(ofname,"_",Snames[ith],".dat",sep="")
          if(ith>1) cat("\n\n",file=ofileSout,append=T)
          cat(paste("#  ",ifname,"season:",Snames[ith],"\n"),file=ofileSout,append=T)
          file.append(ofileSout,ifileSout)
          file.remove(ifileSout)
          i1data<-read.table(ifileAout)
          odata<-rbind(odata,i1data)
          file.remove(ifileAout)
          file.remove(iftmp)
        }
        ymd<-odata[,1]*10000+odata[,2]*100+odata[,3]
        o1data<-odata[order(ymd),]
        write.table(o1data,file=ofileAout,col.names=F,row.names=F)
      }
      
      oact<-str40("ofbody",ofbody)
      ocurdir<-str40("curdir",curdir)
      ooutdir<-str40("outdir",outdir)
      if(exists("ofref")) rm("ofref",envir=.GlobalEnv)
      b20<-"                    "
      oref<-paste(b20,"                  NA",sep="")
      tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=5,sticky="e")
      tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=6,sticky="e")
      tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=7,sticky="e")
      tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=8,sticky="e")
      tkinsert(txt,"end","adjDLYp finished successfully...\n")
      tkinsert(txt,"end",paste("Final output at ",outdir,"/",ofbody,"_*\n\n",sep=""))
      tkdestroy(tt)
      tkfocus(main)
      return()
    }
  }
  Ok2.but<-tkbutton(tt,text="   OK   ",command=OnOk2)
  tkgrid(Ok2.but,column=3,row=7)
  tkfocus(tt)
}

OnQuit<-function(){
  if(exists("curdir")) rm("curdir",envir=.GlobalEnv)
  if(exists("outdir")) rm("outdir",envir=.GlobalEnv)
  if(exists("ofbody")) rm("ofbody",envir=.GlobalEnv)
  if(exists("ofref")) rm("ofref",envir=.GlobalEnv)
  if(exists("ofname")) rm("ofname",envir=.GlobalEnv)
  if(exists("ifname")) rm("ifname",envir=.GlobalEnv)
  if(exists("ifrname")) rm("ifrname",envir=.GlobalEnv)
  if(exists("iIpsName")) rm("iIpsName",envir=.GlobalEnv)
  if(exists("iDIpsName")) rm("iDIpsName",envir=.GlobalEnv)
  if(exists("MissingStr")) rm("MissingStr",envir=.GlobalEnv)
  if(exists("PlevStr")) rm("PlevStr",envir=.GlobalEnv)
  if(exists("AdjStr")) rm("AdjStr",envir=.GlobalEnv)
  if(exists("Mq0Str")) rm("Mq0Str",envir=.GlobalEnv)
  if(exists("Ny4a")) rm("Ny4a",envir=.GlobalEnv)
  if(exists("pthrstr")) rm("pthrstr",envir=.GlobalEnv)
  if(exists("xscr")) rm("xscr",envir=.GlobalEnv)
  if(exists("yscr")) rm("yscr",envir=.GlobalEnv)
  if(exists("txt")) rm("txt",envir=.GlobalEnv)
  if(exists("ErrorMSG")) rm("ErrorMSG",envir=.GlobalEnv)
  if(exists("textMissing")) rm("textMissing",envir=.GlobalEnv)
  tkdestroy(main)
  if(exists("main")) rm("main",envir=.GlobalEnv)
}

Chg.Para<-function(){
  tt<-tktoplevel()
  tkwm.title(tt,"Change Parameters")
  
  textMissing<<-tclVar(paste(MissingStr))
  Entry.Missing<-tkentry(tt,width="10",textvariable=textMissing)
  tkgrid(tklabel(tt,text="Please enter the Missing Value Code."),sticky="w",
         column=1,row=1)
  tkgrid(Entry.Missing,column=2,row=1)
  
  textPlev<<-tclVar(paste(PlevStr))
  Entry.Plev<-tkentry(tt,width="10",textvariable=textPlev)
  tkgrid(tklabel(tt,text="Please enter the nominal conf. level p.lev value."),
         sticky="w",column=1,row=2)
  tkgrid(Entry.Plev,column=2,row=2)
  
  textAdj<<-tclVar(paste(AdjStr))
  Entry.Adj<-tkentry(tt,width="10",textvariable=textAdj)
  tkgrid(tklabel(tt,text="Please enter integer Iadj (0 to 10000 inclusive)"),
         sticky="w",column=1,row=3)
  tkgrid(Entry.Adj,column=2,row=3)
  
  textMq0<<-tclVar(paste(Mq0Str))
  Entry.Mq0<-tkentry(tt,width="10",textvariable=textMq0)
  tkgrid(tklabel(tt,text="Please enter integer Mq (# of points for evaluating PDF)"),
         sticky="w",column=1,row=4)
  tkgrid(Entry.Mq0,column=2,row=4)
  
  textNy4a<<-tclVar(paste(Ny4aStr))
  Entry.Ny4a<-tkentry(tt,width="10",textvariable=textNy4a)
  tkgrid(tklabel(tt,text="Please enter integer Ny4a (>=5, or 0 for choosing the whole segment)"),
         sticky="w",column=1,row=5)
  tkgrid(Entry.Ny4a,column=2,row=5)
  
  textpthr<<-tclVar(paste(pthrstr))
  Entry.pthr<-tkentry(tt,width="10",textvariable=textpthr)
  tkgrid(tklabel(tt,text="Please enter the lower precipitation threshold pthr (>=0)"),
         sticky="w",column=1,row=6) 
  tkgrid(Entry.pthr,column=2,row=6)
  
  OnOk1<-function(){
    oflg<-1
    GuiErrorMSG<-NULL
    MissingStr<-tclvalue(textMissing)
    olen<-40
    assign("MissingStr",MissingStr,envir=.GlobalEnv)
    if(nchar(MissingStr)>olen) {
      GuiErrorMSG<-"MissingCode length error!\n"
      oflg<-0
    }
    
    PlevStr<-tclvalue(textPlev)
    if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
      GuiErrorMSG<-paste(GuiErrorMSG,"p.lev must be one of these: 0.75,0.80,0.90,0.95,0.99,0.9999. Please re-enter\n")
      oflg<-0
    }
    
    AdjStr<-tclvalue(textAdj)
    if(!as.numeric(AdjStr)%in%c(0:10000)){
      GuiErrorMSG<-paste(GuiErrorMSG,"Integer Iadj must be between 0 and 10000 inclusive, please re-enter\n")
      oflg<-0
    }
    
    Mq0Str<-tclvalue(textMq0)
    if(!as.numeric(Mq0Str)%in%c(0:100)){
      GuiErrorMSG<-paste(GuiErrorMSG,"Mq setting must be an integer between 0 and 100, please re-enter\n")
      oflg<-0
    }
    
    Ny4aStr<-tclvalue(textNy4a)
    if(as.numeric(Ny4aStr)!=as.integer(Ny4aStr)){
      GuiErrorMSG<-paste(GuiErrorMSG,"Ny4a must be an integer, please re-enter\n")
      oflg<-0
    }
    
    pthrstr<-tclvalue(textpthr)
    if(as.numeric(pthrstr)!=as.numeric(pthrstr)){
      GuiErrorMSG<-paste(GuiErrorMSG,"pthr must be an real number, please re-enter\n")
      oflg<-0
    }
    
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      assign("MissingStr",MissingStr,envir=.GlobalEnv)
      assign("PlevStr",PlevStr,envir=.GlobalEnv)
      assign("AdjStr",AdjStr,envir=.GlobalEnv)
      assign("Mq0Str",Mq0Str,envir=.GlobalEnv)
      assign("Ny4aStr",Ny4aStr,envir=.GlobalEnv)
      assign("pthrstr",pthrstr,envir=.GlobalEnv)
      
      tkinsert(txt,"end",paste("MissingValueCode set to:",MissingStr,"..\n",sep=" "))
      tkinsert(txt,"end",paste("The nominal level p.lev = ",PlevStr,"..\n",sep=" "))
      tkinsert(txt,"end",paste("Current Iadj value is",AdjStr,"..\n",sep=" "))
      tkinsert(txt,"end",paste("Current Mq value is",Mq0Str,"..\n",sep=" "))
      tkinsert(txt,"end",paste("Current Ny4a value is",Ny4aStr,"..\n",sep=" "))
      tkinsert(txt,"end",paste("Current pthr value is",pthrstr,"..\n",sep=" "))
      
      b20<-"                    "
      olen<-40
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(MissingStr),olen)<-MissingStr
      omiss<-otmp
      tkgrid(tklabel(frameMiddle,text=omiss,width=40),column=2,row=1,sticky="e")
      
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(PlevStr),olen)<-PlevStr
      oplev<-otmp
      tkgrid(tklabel(frameMiddle,text=oplev,width=40),column=2,row=2,sticky="e")
      
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(AdjStr),olen)<-AdjStr
      oadj<-otmp
      tkgrid(tklabel(frameMiddle,text=oadj,width=40),column=2,row=3,sticky="e")
      
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(Mq0Str),olen)<-Mq0Str
      omq0<-otmp
      tkgrid(tklabel(frameMiddle,text=omq0,width=40),column=2,row=4,sticky="e")
      
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(Ny4aStr),olen)<-Ny4aStr
      ony0<-otmp
      tkgrid(tklabel(frameMiddle,text=ony0,width=40),column=2,row=5,sticky="e")
      
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(pthrstr),olen)<-pthrstr
      onp0<-otmp
      tkgrid(tklabel(frameMiddle,text=onp0,width=40),column=2,row=6,sticky="e")
      
      tkdestroy(tt)
      tkfocus(main)
      return()
    }
  }
  
  tkbind(Entry.Ny4a,"<Return>",OnOk1)
  
  Ok1.but<-tkbutton(tt,text="   OK   ",command=OnOk1)
  tkbind(Entry.Missing,"<Return>",OnOk1)
  tkgrid(Ok1.but,column=2,sticky="w",row=7)
  #tkgrid(Ok1.but,column=1,sticky="e",row=7)
  tkfocus(tt)
}

StartGUI<-function(){
  require(tcltk)

GUI<-function(){
  if(!exists("MissingStr")) MissingStr<-"-99.9"
  if(!exists("PlevStr")) PlevStr<-"0.95"
  if(!exists("AdjStr")) AdjStr<-"10000"
  if(!exists("Mq0Str")) Mq0Str<-"100"
  if(!exists("Ny4aStr")) Ny4aStr<-"10"
  if(!exists("pthrstr")) pthrstr<-"0.0"
  assign("MissingStr",MissingStr,envir=.GlobalEnv)
  assign("PlevStr",PlevStr,envir=.GlobalEnv)
  assign("AdjStr",AdjStr,envir=.GlobalEnv)
  assign("Mq0Str",Mq0Str,envir=.GlobalEnv)
  assign("Ny4aStr",Ny4aStr,envir=.GlobalEnv)
  assign("pthrstr",pthrstr,envir=.GlobalEnv)
  
  main<-tktoplevel()
  assign("main",main,envir=.GlobalEnv)
  tkwm.title(main,"RHtests_dlyPrcp")
  fontHeading<-tkfont.create(family="times",size=20,weight="bold",
                             slant="italic")
  fontLable<-tkfont.create(family="times",size=15,weight="bold")
  
  frameOverall<-tkframe(main)
  frameUpper<-tkframe(frameOverall,relief="groove",borderwidth=2)
  assign("frameUpper",frameUpper,env=.GlobalEnv)
  frameMiddle<-tkframe(frameOverall,relief="groove",borderwidth=2)
  assign("frameMiddle",frameMiddle,env=.GlobalEnv)
  tkgrid(tklabel(frameUpper,text="RHtests for daily",font=fontHeading),row=1,column=1)
  tkgrid(tklabel(frameUpper,text="precipitation data",font=fontHeading),row=1,column=2)
  
  ChgPara.but<-tkbutton(frameUpper,text="Change Pars",width=14,command=Chg.Para)
  #Rfile.but<-tkbutton(frameUpper,text="Transform Data",width=14,command=Read.file)
  FindU.but<-tkbutton(frameUpper,text= "FindU",width=14,command=OnFindU.dlyPrcp)
  FindUD.but<-tkbutton(frameUpper,text="FindUD",width=14,command=OnFindUD.dlyPrcp)
  StepSize.but<-tkbutton(frameUpper,text="StepSize",width=14,command=OnStepSize.dlyprcp)
  Cancel.but<-tkbutton(frameUpper,text="Quit",width=14,command=OnQuit)
#   QMadjGaussian.wRef.but<-tkbutton(frameUpper,text= "QMadjGaussian.wRef",width=19,
#                            command=OnQMadjGaussian.wRef)
  #FindUD.wRef.but<-tkbutton(frameUpper,text="FindUD.wRef",width=14,
  #                          command=OnFindUD.wRef)
  #StepSize.wRef.but<-tkbutton(frameUpper,text="StepSize.wRef",width=14,
  #                            command=OnStepSize.wRef)
  
  adjDLYp.but<-tkbutton(frameUpper,text="QMadjDLYp",width=14,
                         command=OnadjDLYp)
  #adjDLYp.wRef.but<-tkbutton(frameUpper,text="adjDLYp.wRef",width=14,
  #                      command=OnadjDLYp.wRef)
  
  olen<-40
  b20<-"                    "
  if(nchar(MissingStr)>olen) stop("MissingValueCode length error!")
  otmp<-paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(MissingStr),olen)<-MissingStr
  omiss<-otmp
  if(nchar(PlevStr)>olen) stop("P_lev length error!")
  otmp<-paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(PlevStr),olen)<-PlevStr
  oplev<-otmp
  otmp<-paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(AdjStr),olen)<-AdjStr
  oadj<-otmp
  otmp<-paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(Mq0Str),olen)<-Mq0Str
  omq0<-otmp
  otmp<-paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(Ny4aStr),olen)<-Ny4aStr
  oNy4a<-otmp
  otmp<-paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(pthrstr),olen)<-pthrstr
  opthr<-otmp
    
  
  if(!exists("ofbody")) otmp<-paste(b20,"                  NA",sep="")
  else{
    if(nchar(ofbody)>olen) otmp<-paste("...",substr(ofbody,nchar(ofbody)-olen+4,nchar(ofbody)),sep="")
    else{
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(ofbody),olen)<-ofbody
    }
  }
  oact<-otmp
#  if(!exists("ofref")) otmp<-paste(b20,"                  NA",sep="")
#  else{
#    if(nchar(ofref)>olen) otmp<-paste("...",substr(ofref,nchar(ofref)-olen+4,nchar(ofref)),sep="")
#    else{
#      otmp<-paste(b20,b20,sep="")
#      substr(otmp,olen+1-nchar(ofref),olen)<-ofref
#    }
#  }
#  oref<-otmp
  if(!exists("curdir")) otmp<-paste(b20,"                  NA",sep="")
  else{
    if(nchar(curdir)>olen) otmp<-paste("...",substr(curdir,nchar(curdir)-olen+4,nchar(curdir)),sep="")
    else{
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(curdir),olen)<-curdir
    }
  }
  ocurdir<-otmp
  if(!exists("outdir")) otmp<-paste(b20,"                  NA",sep="")
  else{
    if(nchar(outdir)>olen) otmp<-paste("...",substr(outdir,nchar(outdir)-olen+4,nchar(outdir)),sep="")
    else{
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(outdir),olen)<-outdir
    }
  }
  ooutdir<-otmp
  # arrange menu buttons and labels
  tkgrid(tklabel(frameUpper,text="          "),ChgPara.but,tklabel(frameUpper,text="          "),Cancel.but)
  tkgrid(tklabel(frameUpper,text="PMF and F tests:",font=fontLable),sticky="w")
  tkgrid(FindU.but,column=1,row=3)
  tkgrid(FindUD.but,column=2,row=3)
  tkgrid(StepSize.but,column=3,row=3)
  #tkgrid(tklabel(frameUpper,text="PMT and t tests:",font=fontLable),sticky="w")
  #tkgrid(Cancel.but,column=3,row=4)
  #tkgrid(FindUD.wRef.but,column=2,row=4)
  #tkgrid(StepSize.wRef.but,column=3,row=4)
  #tkgrid(tklabel(frameUpper,text="  To adjust daily",font=fontLable),sticky="w")
  #tkgrid(tklabel(frameUpper,text="Gaussian data:",font=fontLable),sticky="w",column=1,row=4)
  
  #tkgrid(adjDLYp.but,column=1,row=4)####jan21
  #tkgrid(adjDLYp.wRef.but,column=3,row=4)
  
  tkgrid(tklabel(frameMiddle,text="Current Missing Value Code:"),
         column=1,row=1,sticky="w") 
  tkgrid(tklabel(frameMiddle,text=omiss,width=40),column=2,row=1,sticky="w")
  
  tkgrid(tklabel(frameMiddle,text="Current nominal level of confidence (p.lev):"),
         column=1,row=2,sticky="w") 
  tkgrid(tklabel(frameMiddle,text=oplev,width=40),column=2,row=2,sticky="e") 
  
  tkgrid(tklabel(frameMiddle,text="Segment to which to adjust the series (Iadj):"),
         column=1,row=3,sticky="w")
  tkgrid(tklabel(frameMiddle,text=oadj,width=40),column=2,row=3,sticky="e")
  
  tkgrid(tklabel(frameMiddle,text="Current Mq (# of points for evaluating PDF):"),
         column=1,row=4,sticky="w")
  tkgrid(tklabel(frameMiddle,text=omq0,width=40),column=2,row=4,sticky="e")
  
  tkgrid(tklabel(frameMiddle,text="Current Ny4a (max # of years of data for estimating PDF):"),
         column=1,row=5,sticky="w")
  tkgrid(tklabel(frameMiddle,text=oNy4a,width=40),column=2,row=5,sticky="e")
  
  tkgrid(tklabel(frameMiddle,text="Current pthr (Lower threshold of precipitation):"),
         column=1,row=6,sticky="w")
  tkgrid(tklabel(frameMiddle,text=opthr,width=40),column=2,row=6,sticky="e")
  
  tkgrid(tklabel(frameMiddle,text="Current input Base series filename:"),
         column=1,row=7,sticky="w")
  tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=7,sticky="w")
  
  tkgrid(tklabel(frameMiddle,text="Current data directory:    "),
         column=1,row=8,sticky="w")
  tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
  
  tkgrid(tklabel(frameMiddle,text="Current output directory:    "),
         column=1,row=9,sticky="w")
  tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9)
  
  # frameMiddle<-tkframe(frameOverall,relief="groove",borderwidth=2)
  # assign("frameMiddle",frameMiddle,env=.GlobalEnv)
  frameLower<-tkframe(frameOverall,relief="groove",borderwidth=2)
  assign("xscr",tkscrollbar(frameLower,repeatinterval=5,orient="horizontal",
                            command=function(...)tkxview(txt,...)),envir=.GlobalEnv)
  assign("yscr",tkscrollbar(frameLower,repeatinterval=5,
                            command=function(...)tkyview(txt,...)),envir=.GlobalEnv)
  assign("txt",tktext(frameLower,bg="white",font="courier",
                      xscrollcommand=function(...)tkset(xscr,...),
                      yscrollcommand=function(...)tkset(yscr,...) 
  ),envir=.GlobalEnv)
  tkgrid(frameUpper)
  tkgrid(frameMiddle)
  tkgrid(txt,yscr)
  tkgrid(xscr)
  tkgrid.configure(yscr,sticky="ns")
  tkgrid.configure(xscr,sticky="ew")
  tkfocus(txt)
  tkgrid(frameLower)
  tkgrid(frameOverall)
}

UA1<-'RClimDex and RHtests software packages (all versions included), herein after called "The Library" \n'
UA2<-'©  Copyright, Environment Canada, 2012, \n'
UA3<-"The Library was created by the Climate Research Division of Environment Canada and as such all intellectual property rights (including copyright) that may exist in The Library are owned by the Government of Canada.  The Library software code is provided free under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, version 3.0 of the License. It is distributed under the terms of this license 'as-is' and has not been designed or prepared to meet any Licensee's particular requirements. Environment Canada makes no warranty, either express or implied, including but not limited to, warranties of merchantability or fitness for a particular purpose. In no event will Environment Canada be liable for any indirect, special, consequential or other damages attributed to the Licensee's use of The Library. In downloading The Library you understand and agree to these terms and those of the associated LGP License. See the GNU Lesser General Public License for more details.\n"
UA4<-"You should have received a copy of the GNU Lesser General Public License along with This Library; if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA."
UserAgrement <- paste(UA1,UA2,UA3,UA4,sep='\n')
yesAgree<-function(){
tkdestroy(topFrame)
GUI()
}
noAgree<-function(){
tkdestroy(topFrame)
}

topFrame <- tktoplevel(pady=2, bg="gray94")
outFrame <- tkframe(topFrame, padx=2, bg="gray94")
buttonFrame <- tkframe(topFrame, padx=2, bg="gray94")
tkwm.title(topFrame, "RHtestsV4 User Agreement")
tkwm.resizable(topFrame, 0, 0)
tkgrid(outFrame)
tkgrid(buttonFrame)
tkgrid.configure(outFrame, sticky="nsew")
scrollBar <- tkscrollbar(outFrame, repeatinterval=5, command=function(...)tkyview(textFrame,...))
textFrame <- tktext(outFrame,bg="gray94",font="courier",yscrollcommand=function(...)tkset(scrollBar,...))
tkgrid(textFrame,scrollBar)
tkgrid.configure(scrollBar,sticky="ns")
tkinsert(textFrame,"end",paste(UserAgrement,sep=""))
tkconfigure(textFrame, state="disabled")
yesButton <- tkbutton(buttonFrame, text="I Agree", command=function()yesAgree())
noButton <- tkbutton(buttonFrame, text="I Do Not Agree", command=function()noAgree())
tkgrid(yesButton,noButton)
tkbind(topFrame,"<Destroy>",function()noAgree())
tkfocus(textFrame)
}


#StartGUI()
