RClimDex2RHtest<-function(ifname,MissingStr="-99.9"){
# this function transfer RClimDex data format (daily) to RHtest format (monthly)
# input data file name -- ifname, and MissingCodeValue -- MissingStr
# output 4 data files: *_prcp *_tmax *_tmin *_dtr, those files will be at 
# the same directory as input data file, and using the same prex as input file
# name.
  if(!nchar(ifname)){
    cat("No file selected in Data Transform!\n")
    return()
  }
  outdirtmp<-strsplit(ifname,"/")[[1]]
  if(length(outdirtmp)==1) curdir<-""
  else if(length(outdirtmp)<=2) curdir<-paste(outdirtmp[1])
  else{
    curdir<-outdirtmp[1]
    for(i in 2:(length(outdirtmp)-1))
      curdir<-paste(curdir,outdirtmp[i],sep="/")
  }
  ofname<-outdirtmp[length(outdirtmp)]
  if(is.csv(ofname)) csv<-1
  else csv<-0
  if(csv){
    itmp<-try(read.table(ifname,header=F,
            col.names=c("year","month","day","prcp","tmax","tmin"),
            sep=",",na.strings=MissingStr,colClasses=rep("real",6)),
	    silent=T)
    if(inherits(itmp,"try-error")){
      ErrorMSG<<-geterrmessage()
      cat(paste(ErrorMSG,"\n"))
      return()
    }
    else iidata<-itmp
  }
  else{
    itmp<-try(read.table(ifname,header=F,
            col.names=c("year","month","day","prcp","tmax","tmin"),
            na.strings=MissingStr,colClasses=rep("real",6)),
	    silent=T)
    if(inherits(itmp,"try-error")){
      ErrorMSG<<-geterrmessage()
      cat(paste(ErrorMSG,"\n"))
      return()
    }
    else iidata<-itmp
  }
  if(ncol(iidata)!=6){
    ErrorMSG<-paste(ifname,"has",ncol(iidata),"columns. The number of columns should be 6")
    cat(paste(ErrorMSG,"\n"))
    return()
  }
  nlen<-dim(iidata)[1]
  syear<-iidata[1,1]
  eyear<-iidata[nlen,1]
  smonth<-iidata[1,2]
  emonth<-iidata[nlen,2]
  if(eyear<(syear+1)) {
    ErrorMSG<-paste("Time series",ifname,"too short for Transform Data")
    cat(paste(ErrorMSG,"\n"))
    return()
  }
  nyrs<-eyear-syear+1
  vars<-c("tmax","tmin","prcp")
  ivars<-length(vars)
  
  cat(paste("Data from:",syear,"to:",eyear,"\n",sep=" "))

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
    if(vars[i]=="tmax") tmax<-mdata
    else if(vars[i]=="tmin") tmin<-mdata
    else { prcp<-mdata; prcp1mm<-mdata1mm}
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
  if(curdir!=""){
    ofmax<-paste(curdir,ofmax,sep="/")
    ofmin<-paste(curdir,ofmin,sep="/")
    ofprcp<-paste(curdir,ofprcp,sep="/")
    ofprcpL<-paste(curdir,ofprcpL,sep="/")
    ofprcp1mm<-paste(curdir,ofprcp1mm,sep="/")
    ofprcp1mmL<-paste(curdir,ofprcp1mmL,sep="/")
    ofmaxD<-paste(curdir,ofmaxD,sep="/")
    ofminD<-paste(curdir,ofminD,sep="/")
    ofprcpD<-paste(curdir,ofprcpD,sep="/")
  }
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
  cat(paste("Data transform finished, monthly output:\n",ofmax,"\n",ofmin,"\n",
            ofprcp,"\n",ofprcpL,"\n",ofprcp1mm,"\n",ofprcp1mmL,"\n",
	    "daily output:\n",ofmaxD,"\n",ofminD,"\n",ofprcpD,"\n"))
}

RHtest2RClimDex<-function(ifprcp,icolprcp,iftmax,icoltmax,iftmin,icoltmin,
  ofname,MissingStr="-99.9"){
# this function will combine 3 columns from 3 different input file into a single
# output file, the input files should be output file from RHtest. 3 input file
# name should be full file name with directory, togather with the column number
# for this file. Input prcp file -- ifprcp, column number for prcp -- icolprcp;
# input tmax file -- iftmax, column number for tmax -- icoltmax; input tmin 
# file -- iftmin, column number for tmin. Column numbers should greater than 3
# since first 3 columns will be YYYY MM DD
  if(!nchar(ifprcp)){
    cat("No input prcp filename in Data Transform!\n")
    return()
  }
  if(!nchar(iftmax)){
    cat("No input tmax filename in Data Transform!\n")
    return()
  }
  if(!nchar(iftmin)){
    cat("No input tmin filename in Data Transform!\n")
    return()
  }
  prcp<-try(read.table(ifprcp,header=F,na.strings=MissingStr,
            colClasses="real"),silent=T)
  if(inherits(prcp,"try-error")){
      ErrorMSG<<-geterrmessage()
      cat(paste(ErrorMSG,"\n"))
      return()
    }
  if(dim(prcp)[2]<icolprcp){
    cat("Input prcp file has less column than selected icol for prcp")
    return()
  }
  tmax<-try(read.table(iftmax,header=F,na.strings=MissingStr,
            colClasses="real"),silent=T)
  if(inherits(tmax,"try-error")){
      ErrorMSG<<-geterrmessage()
      cat(paste(ErrorMSG,"\n"))
      return()
    }
  if(dim(tmax)[2]<icoltmax){
    cat("Input tmax file has less column than selected icol for tmax")
    return()
  }
  tmin<-try(read.table(iftmin,header=F,na.strings=MissingStr,
            colClasses="real"),silent=T)
  if(inherits(tmin,"try-error")){
      ErrorMSG<<-geterrmessage()
      cat(paste(ErrorMSG,"\n"))
      return()
    }
  if(dim(tmin)[2]<icoltmin){
    cat("Input tmin file has less column than selected icol for tmin")
    return()
  }
  otmp<-merge(prcp[,c(2,icolprcp)],tmax[,c(2,icoltmax)],by.x=1,
              by.y=1,all=T,sort=F)
  otmp2<-merge(otmp,tmin[,c(2,icoltmin)],by.x=1,by.y=1,all=T,sort=F)
  ymd<-otmp2[,1]
  yyyy<-as.integer(ymd/10000)
  mm<-as.integer((ymd-yyyy*10000)/100)
  dd<-ymd-as.integer(ymd/100)*100
  write.table(cbind(yyyy,mm,dd,otmp2[,2:4]),file=ofname,col.names=F,row.names=F,na=MissingStr)
}

is.csv<-function(names){
  nlen<-nchar(names)
  if(substr(names,nlen-2,nlen)=="csv"|substr(names,nlen-2,nlen)=="CSV") return(T)
  else return(F)
}
