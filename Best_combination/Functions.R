
kdmTopwm<-function(motifs,nmat=1,nflank=5,base=2,alpha=1e-6,cutoff=2,bkg=c(0.25,0.25,0.25,0.25),pseudo=1,um=TRUE,fraction=0.9,estimate=TRUE,threshold=0){
        map<-list()
        map[["c"]]<-c(0,1,0,0)
        map[["g"]]<-c(0,0,1,0)
        map[["t"]]<-c(0,0,0,1)
        map[["n"]]<-bkg	
        map[["A"]]<-c(1,0,0,0)
        map[["C"]]<-c(0,1,0,0)
        map[["G"]]<-c(0,0,1,0)
        map[["T"]]<-c(0,0,0,1)
        map[["N"]]<-bkg		

        getPWM<-function(seq){
            matrix(unlist(map[strsplit(seq,split="")[[1]]]),nrow=4)
        }

        update_ss<-function(kmers,PWM,cutoff,pseudo,bkg,base,klength,alpha){
            ###ALIGN
            res<-data.frame(seq=kmers[,1],weight=kmers$weight,pos=NA,score=NA)
            PWMM<-matrix(0,4,ncol(PWM))
            scol<-which(apply(PWM,2,sum)>0)
            
            
            PWMM<-apply(PWM*1000+pseudo,2,function(a){log((a/sum(a)/bkg),base)})
            maxpos<-ncol(PWMM)-klength+1
            for(i in 1:nrow(kmers)){
                pwm<-getPWM(kmers[i,1])
                scores<-rep(0,maxpos)
                for(pos in 1:maxpos){
                    scores[pos]<-sum(apply(pwm*PWMM[,pos:(pos+klength-1)],2,sum))
                }
                
                res$pos[i]<-which.max(scores)
                res$score[i]<-scores[res$pos[i]]
            }
            
            ###UPDATE
            sel<-which(res$score>cutoff)
            newPWM<-matrix(0,4,ncol(PWMM))
            for(i in 1:length(sel)){
                k<-sel[i]
                start<-res$pos[k]
                newPWM[,start:(start+klength-1)]<-newPWM[,start:(start+klength-1)] + getPWM(res$seq[k])*res$weight[k] + alpha
            }
            newPWM<-newPWM[,which(apply(newPWM,2,sum)>0)]
            newPWM<-apply(newPWM,2,function(a){a/sum(a)})
            rownames(newPWM)<-c("A","C","G","T")
            list(PWM=newPWM,selected=sel,res=res)
        }
        
        cbkg<-(bkg*1000+pseudo)/sum(bkg*1000+pseudo)
        PWMS<-list();
        if(motifs@k>0 & estimate){
            motifs<-kdmEstimateMotifs(motifs)
        }
        
        if(estimate){
            kk1<-kdmGkmersString(kmerlength=motifs@l,ngaps=0,double.strand=FALSE)
        }else{
            kk1<-kdmGkmersString(kmerlength=motifs@l,ngaps=motifs@k,double.strand=FALSE)
        }
        if(motifs@doubleStrand==FALSE){
            #if(estimate){
                #kk2<-kdmGkmersString(kmerlength=motifs@l,ngaps=0,double.strand=TRUE)
            #}else{
                #kk2<-kdmGkmersString(kmerlength=motifs@l,ngaps=motifs@k,double.strand=TRUE)
            #}
            for(k in 1:ncol(motifs)){
                PWMS[[k]]<-list()
                kk1$score<-as.numeric(motifs[,k]^2)
    # 			kk2$weight<-alpha*as.numeric(motifs[,k]^2)
                kk1$weight<-as.numeric(motifs[,k]^2)
                
                oo<-order(kk1$score,decreasing=TRUE)
                kmers<-kk1[oo,]
                #kmers<-kmers[1:floor(nrow(kk2)*fraction),]
    # 			print(cumsum(kk2$score[oo]))
                last<-min(which(cumsum(kk1$score[oo])>fraction))
                #cat("PPPPPP---------------------->",last,"\n")
                kmers<-kmers[1:last,]
    # 			browser()
                
                
                for(m in 1:nmat){
                    #PWM0<-cbind(matrix(0,4,nflank),getPWM(kmers$fwd[1]),matrix(0,4,nflank))
                    PWM0<-cbind(matrix(0,4,nflank),getPWM(kmers[,1]),matrix(0,4,nflank))
                    r<-1000000
    # 				while(r>0.0001){
                    while(r>threshold){
                        PWM<-update_ss(kmers,PWM0,cutoff=cutoff,pseudo=pseudo,bkg=cbkg,base=base,motifs@l,alpha=alpha)
                        if(ncol(PWM$PWM)==ncol(PWM0)){
                            r<-sqrt(sum( (PWM0-PWM$PWM)^2))
                        }
                        PWM0<-PWM$PWM
                        #print(PWM0)
                        #cat("matrix=",m,"\tchange=",r,"\tnkmers=",length(PWM$selected),"\n",sep="")
                    }
                    PWMS[[k]][[m]]<-PWM$PWM
    #  				print(PWM$PWM)
    # 				print(PWM$res)
                    
                    kmers<-kmers[-PWM$selected,]
                }
                if(um){
                    #aa<-list()
                    for(i in 1:length(PWMS[[k]])){
                        aa<-create_motif(PWMS[[k]][[i]],alphabet="DNA",type="PPM",name=paste("GKM.",k,sep=""),bkg=bkg)
                        aa@alphabet<-"DNA"
                    }
                    PWMS[[k]]<-aa
                }
            }
        }else{
            stop("Not yet implemented for single strand motifs")
        }
        PWMS
    }




from_exons_to_bed<-function(input){
    new_rows<-list()
    
    chr<-input["chr"]
    strand<-input["strand"]
    start<-as.numeric(input["exonStart_0base"])
    end<-as.numeric(input["exonEnd"])
    score<-as.numeric(input["IncLevelDifference"])
    name<-paste0(chr,":",start,"-",end)
    
    if(strand=="+"){
    new_rows<-append(new_rows,list(c(chr,"start"=start-(L+5),"end"=start-5,"name"=name,"score"=score,strand,"region"="upstream")))
    new_rows<-append(new_rows,list(c(chr,"start"=end+5,"end"=end+(L+5),"name"=name,"score"=score,strand,"region"="downstream")))
    }
    
    if(strand=="-"){
    new_rows<-append(new_rows,list(c(chr,"start"=start-(L+5),"end"=start-5,"name"=name,"score"=score,strand,"region"="downstream")))
    new_rows<-append(new_rows,list(c(chr,"start"=end+5,"end"=end+(L+5),"name"=name,"score"=score,strand,"region"="upstream")))
    
    }
    
    return(new_rows)
}


from_bed_to_counts<-function(input){
    bed<-new("kbed",input%>%mutate(start=as.numeric(start),end=as.numeric(end),score=as.numeric(score)),type="bed6")
    colnames(bed)<-c("chrom","start","end","name","score","strand")
    seq<-kdmGetSequence(bed,paste0(ref,"/",genome,"/",genome,".2bit"))
    return(kdmDistr(seq,length,gaps,double.strand=FALSE,strict=TRUE))
    }

calculate_coeff<-function(upstream,downstream,W){
    bkg_counts_upstream<-from_bed_to_counts(bkg_upstream_bed[1:ncol(upstream),])
    bkg_counts_downstream<-from_bed_to_counts(bkg_downstream_bed[1:ncol(downstream),])

    pos_features<-cbind(kdmFeatures(W,upstream),kdmFeatures(W,downstream))
    neg_features<-cbind(kdmFeatures(W,bkg_counts_upstream),kdmFeatures(W,bkg_counts_downstream))

    N<-nrow(pos_features)
    mot<-ncol(W)
    ids<-split(sample(1:N),ceiling(seq_along(seq(1:N))/(N/10)))
    coefficients_all<-data.frame()
    for(k in 1:10){
        id<-seq(1:N)[-ids[[k]]]
        train_features<-rbind(pos_features[id,],neg_features[id,])
        train_label<-c(rep(1,length(id)),rep(0,length(id)))
        mm<-cv.glmnet(y=train_label,x=train_features, family="binomial", type.measure="auc",alpha=alpha)
        coefficients_all<-rbind(coefficients_all,data.frame(model=c("intercept",rep(seq(1,mot),2)),window=c("intercept",rep("upstream",mot),rep("downstream",mot)),cv=k,coeff=coef(mm,s="lambda.min")[,1]))
    }
    coeff<-coefficients_all%>%filter (model!="intercept") %>% group_by(model,window) %>% summarise(Coeff=median(coeff))
    return(coeff)
}

annotate_motif<-function(pos.bed,neg.bed){
    
    
    #pos.bed<-U_upstream_bed
    #neg.bed<-bkg_upstream_bed

    
    pos.bed<-pos.bed%>%mutate(start=as.numeric(start),end=as.numeric(end),score=as.numeric(score))%>%mutate(start=start-(hwin+1),end=end+(hwin+1))
    pos.bed<-new("kbed",pos.bed,type="bed6")
    colnames(pos.bed)<-c("chrom","start","end","name","score","strand")
    centers<-floor((pos.bed$end+pos.bed$start)/2)-pos.bed$start
    temp<-pos.bed
    temp$start<-pos.bed$start+centers-flanks
    temp$end<-pos.bed$start+centers+flanks+1 
    pseq<-kdmGetSequence(temp,paste0(ref,"/",genome,"/",genome,".2bit"))
        
    
    #pseq<-kdmGetSequence(pos.bed,paste0(ref,"/",genome,"/",genome,".2bit"))
    
    neg.bed<-neg.bed[1:length(pseq),]
    neg.bed<-neg.bed%>%mutate(start=as.numeric(start),end=as.numeric(end),score=as.numeric(score))%>%mutate(start=start-(hwin+1),end=end+(hwin+1))
    neg.bed<-new("kbed",neg.bed,type="bed6")
    colnames(neg.bed)<-c("chrom","start","end","name","score","strand")
    centers<-floor((neg.bed$end+neg.bed$start)/2)-neg.bed$start
    temp<-neg.bed
    temp$start<-neg.bed$start+centers-flanks
    temp$end<-neg.bed$start+centers+flanks+1 
    nseq<-kdmGetSequence(temp,paste0(ref,"/",genome,"/",genome,".2bit"))
    
    #nseq<-kdmGetSequence(neg.bed,paste0(ref,"/",genome,"/",genome,".2bit"))
    
    pinfo<-kdmSeqProfileInfo(pseq,nseq,W,half.interval=flanks,half.win=hwin)
    enr<-kdmEnrichment(pinfo)
    
    
    return(list(profile=pinfo,enrichment=enr))
}

calculate_score<-function(input){
pos_features<-lapply(1:70,function(i){do.call(cbind,lapply(input$profile$PP$pfeat,function(x) x[,i]))})
neg_features<-lapply(1:70,function(i){do.call(cbind,lapply(input$profile$PP$nfeat,function(x) x[,i]))})
profile<-rbind()
for(m in 1:70){
#profile<-rbind(profile,data.frame(Position=seq(1:93),Score=rowMeans(pos_features[[m]]),Type="Positive",Motif=m))
profile<-rbind(profile,data.frame(Position=seq(1:93),Score=apply(pos_features[[m]],1,median),Type="Positive",Motif=m))
#profile<-rbind(profile,data.frame(Position=seq(1:93),Score=rowMeans(neg_features[[m]]),Type="Negative",Motif=m))}
profile<-rbind(profile,data.frame(Position=seq(1:93),Score=apply(neg_features[[m]],1,median),Type="Negative",Motif=m))
#colnames(profile)<-c("Position","Score","Type")
}
return(profile)
}



plot_summary<-function(exons,th){
if(th=="FDR"){
g0<-ggplot(exons,aes(x=IncLevelDifference,y=-log10(FDR+1e-10),col=Type))+theme_bw()+
geom_point(show.legend=FALSE)+scale_color_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
theme(plot.title=element_text(hjust=0.5),axis.text=element_text(size=10))+ggtitle("FDR vs InclusionLevel")+
geom_hline(yintercept=-log10(0.1+1e-10),lty="dashed")+geom_vline(xintercept=c(0.1,-0.1),lty="dashed")
}
if(th=="pValue"){
g0<-ggplot(exons,aes(x=IncLevelDifference,y=-log10(PValue+1e-10),col=Type))+theme_bw()+
geom_point(show.legend=FALSE)+scale_color_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
theme(plot.title=element_text(hjust=0.5),axis.text=element_text(size=10))+ggtitle("pValue vs InclusionLevel")+
geom_hline(yintercept=-log10(0.05+1e-10),lty="dashed")+geom_vline(xintercept=c(0.1,-0.1),lty="dashed")
}

n<-exons%>%group_by(Type)%>%summarise(n=n())
n2<-exons%>%filter(upstreamI>(L+10),downstreamI>(L+10))%>%group_by(Type)%>%summarise(n=n())%>%mutate(vjust=c(0,1,-1))
g1<-ggplot(exons,aes(x=exonLength,y=Type,col=Type,fill=Type))+
geom_density_ridges(alpha=0.7,show.legend=FALSE,scale = 0.8)+scale_x_log10()+theme_bw()+
scale_fill_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
scale_color_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
theme(axis.title=element_blank(),plot.title=element_text(hjust=0.5),axis.text=element_text(size=10))+ggtitle("Exon Length")+
geom_text(data=n,aes(label=n),x=1,col="black",vjust=-1)


g2<-ggplot(exons,aes(x=upstreamI,y=Type,col=Type,fill=Type))+
geom_density_ridges(alpha=0.7,show.legend=FALSE,scale = 0.8)+scale_x_log10()+theme_bw()+
scale_fill_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
scale_color_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
theme(axis.title=element_blank(),plot.title=element_text(hjust=0.5),axis.text=element_text(size=10))+ggtitle("Upstream Introns Length")+
geom_vline(xintercept=(L+10),lty="dashed")


g3<-ggplot(exons,aes(x=downstreamI,y=Type,col=Type,fill=Type))+
geom_density_ridges(alpha=0.7,show.legend=FALSE,scale = 0.8)+scale_x_log10()+theme_bw()+
scale_fill_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
scale_color_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
theme(axis.title=element_blank(),plot.title=element_text(hjust=0.5),axis.text=element_text(size=10))+ggtitle("Downstream Introns Length")+
geom_vline(xintercept=(L+10),lty="dashed")

g4<-ggplot(exons%>%filter(Type!="Control"),aes(x=upstreamI,y=downstreamI,fill=Type))+
geom_point(shape=21,col="black",show.legend=FALSE)+scale_fill_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
theme_bw()+geom_vline(xintercept=(L+10),lty="dashed")+geom_hline(yintercept=(L+10),lty="dashed")+
scale_x_log10()+scale_y_log10()+xlab("Upstream Intron length")+ylab("Downstream Intron length")+
geom_text(data=n2%>%filter(Type!="Control"),x=2.5,y=5,aes(label=n,col=Type,hjust=vjust),show.legend=FALSE)+
scale_color_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))


g6<-ggplot(exons%>%filter(Type!="Control")%>%mutate(chr=factor(chr,levels= c(paste0("chr",seq(1:22)),"chrX","chrY"))),
  aes(y=chr,fill=Type,col=Type))+geom_bar(alpha=0.6,width=0.6,show.legend=FALSE)+
facet_grid(~Type)+theme_bw()+scale_y_discrete(limits=rev)+
scale_color_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
scale_fill_manual(values=c("Control"="grey","Enhanced"="red","Silenced"="blue"))+
theme(axis.title=element_blank(),plot.title=element_text(hjust=0.5),axis.text=element_text(size=10))+
ggtitle("Chromosome frequency")

pdf(paste0(wd,"/Figure/Summary_",rbp,"_",cell_line,"_",genome,".pdf"),height = 6,width = 10)
print(plot_grid(g0,g1,g6,g2,g4,g3,nrow=2))
dev.off()
}



plot_motif<-function(total_data){

scores<-total_data$profile%>%
mutate(Annotation=ifelse(Type=="Negative","Negative",ifelse(Alt=="UP","UP-regulated","DOWN-regulated")))%>%
mutate(Val=ifelse(Alt=="UP",Score,-Score))%>%filter(Position%in%seq(2,(L+2)))%>%
mutate(Pos=ifelse(window=="upstream",Position-(L+7),Position+3))
RES<-total_data$test%>%mutate(Alt = recode(Alt, "UP" = "ENHANCED", "DOWN" = "SILENCED"),Alt=factor(Alt,levels=c("SILENCED","ENHANCED")))
total_coeff<-total_data$coeff%>%mutate(Alt = recode(Alt, "UP" = "ENHANCED", "DOWN" = "SILENCED"),Alt=factor(Alt,levels=c("SILENCED","ENHANCED")))
pwms_converted<-total_data$pwms
kdmtom<-total_data$kdmtom
sign<-total_data$sel

lim_coeff<-range(total_coeff$Coeff)
lim_over<-range(RES$Value)
lim_scores<-range(scores$Val)

#order<-unique(total_coeff%>%arrange(desc(Coeff))%>%pull(model))
#delta<-scores%>%group_by(Motif,Alt,window,Pos)%>%summarise(delta=Score[1]-Score[2])
delta<-scores%>%group_by(Motif,Alt,window,Pos)%>%summarise(delta=log10(Score[1]/Score[2]))
max.d<-delta%>%group_by(Motif)%>%summarise(max=max(delta))
order<-max.d%>%arrange(desc(max))%>%pull(Motif)
max.d<-max.d%>%mutate(r=round(max/max(max),3))

plots<-list()

for(i in 1:length(order)){
if(!order[i]%in%sign$model){next}
#motif=as.numeric(sign$model[i])
motif=as.numeric(order[i])
print(motif)

p0<-ggplot(max.d%>%filter(Motif==motif),aes(x="",y="",fill=max))+
geom_point(shape=21,col="black",size=20,show.legend=FALSE)+theme_void()+geom_text(aes(label=r),size=5)+
scale_fill_gradient2(low="white",high="darkorange1",midpoint=median(max.d$max),mid="white",limits=range(max.d$max))+
coord_fixed()


p1<-ggplot(total_coeff%>%filter(model==motif))+geom_tile(aes(y=Alt,x=window,fill=Coeff),height=0.6,width=0.6,col="black",show.legend=FALSE)+
theme_minimal()+scale_fill_gradient2(low="grey30",high="gold",midpoint=0,mid="white",limits=lim_coeff)+
theme(panel.grid=element_blank(),axis.title=element_blank(),plot.title=element_text(hjust=0.5))+ggtitle("GLM coeffcient")+
scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+
coord_fixed()


p2<-ggplot(RES%>%filter(model==motif))+geom_tile(aes(y=Alt,x=window,fill=Value),height=0.6,width=0.6,col="black",show.legend=FALSE)+
theme_minimal()+scale_fill_gradient2(low="purple",high="forestgreen",midpoint=0,mid="white",limits=lim_over)+
theme(panel.grid=element_blank(),axis.title=element_blank(),plot.title=element_text(hjust=0.5))+ggtitle("Overrepresentation")+
scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+
coord_fixed()

# p3<-ggplot(data=NULL,aes(x=Pos,y=Val,col=Annotation,fill=Annotation))+
#  geom_area(data=scores%>%filter(Motif==motif,Alt=="DOWN",Type=="Negative"),alpha=0.5,show.legend=FALSE)+
#  geom_area(data=scores%>%filter(Motif==motif,Alt=="UP",Type=="Negative"),alpha=0.5,show.legend=FALSE)+
#  geom_area(data=scores%>%filter(Motif==motif,Type=="Positive"),alpha=0.5,show.legend=FALSE)+
#  scale_color_manual(values = c("Negative"="grey30","UP-regulated"="red","DOWN-regulated"="blue"))+
#  scale_fill_manual(values = c("Negative"="grey30","UP-regulated"="red","DOWN-regulated"="blue"))+
#  facet_grid(~window,scales="free_x")+
#  theme_bw()+ggtitle(paste0("Motif: ",motif))+theme(axis.title=element_blank(),panel.grid=element_blank())+
#  scale_x_continuous(expand=c(0,0),breaks=seq(-95,95,10))+ylim(lim_scores)

p3<-ggplot(delta%>%filter(Motif==motif),aes(x=Pos,y=delta,col=Alt))+
geom_line(show.legend=FALSE)+scale_color_manual(values=c("UP"="red","DOWN"="blue"))+theme_bw()+
facet_grid(~window,scale="free")+scale_x_continuous(expand=c(0,0))+
theme(axis.title=element_blank(),panel.grid=element_blank())+geom_hline(yintercept=0,lty="dashed")+
ggtitle(paste0("Motif: ",motif))

# p3<-ggplot(data=NULL,aes(x=Position,y=Score,col=Annotation,fill=Annotation))+
  # geom_stream(data=scores%>%filter(Type=="Negative",Motif==motif,Alt=="UP"),type = "ridge",
              # geom="polygon",size = 0.5,
              # bw=0.4,alpha=0.4,show.legend = FALSE
  # )+
  # geom_stream(data=scores%>%filter(Type=="Negative",Motif==motif,Alt=="DOWN"),type = "ridge",
              # geom="polygon",size = 0.5,
              # bw=0.4,alpha=0.4,show.legend = FALSE
  # )+
  # geom_stream(data=scores%>%filter(Type=="Positive",Motif==motif,Alt=="UP"),type = "ridge",
              # geom="polygon",size = 0.5,
             # bw=0.4,alpha=0.4,show.legend = FALSE
  # )+
  # geom_stream(data=scores%>%filter(Type=="Positive",Motif==motif,Alt=="DOWN"),type = "ridge",
              # geom="polygon",size = 0.5,
              # bw=0.4,alpha=0.7,show.legend = FALSE
  # )+
  # scale_color_manual(values = c("Negative"="black","UP-regulated"="red","DOWN-regulated"="blue"))+
  # scale_fill_manual(values = c("Negative"="black","UP-regulated"="red","DOWN-regulated"="blue"))+
  #scale_color_manual(values = c("red"))+
  #scale_fill_manual(values = c("red"))+
  #facet_wrap(Alt~window,strip.position = "left",scales="free_x")+
  # facet_grid(Alt~window,scales="free_x")+
  # theme_bw()+ggtitle(paste0("Motif: ",motif))

 
 p4<-ggseqlogo(pwms_converted[[motif]])
 
 if(!motif%in%kdmtom$query){p5<-NULL}
 if(motif%in%kdmtom$query){
p5<-ggplot(kdmtom%>%filter(query==motif),aes(x="",y=""))+theme_void()+geom_tile(aes(fill=class),alpha=0.5,col="black",show.legend = FALSE)+
scale_fill_manual(values=c("1"="white","2"="forestgreen","3"="goldenrod","4"="magenta4"))+geom_text(aes(label=label))
}

plots[[(length(plots)+1)]]<-plot_grid(p0,p3,p1,p2,p4,p5,align="h",axis="tb",rel_widths=c(1,8,2,2,4,1),nrow=1)

}


pdf(paste0(wd,"/Figure/",rbp,"_",cell_line,"_",genome,"_",tt,".pdf"),height = 70,width = 30)
#ggarrange(p1,p2,p3,p4,p5,p6,nrow = 1,align = "h",widths = c(4,1,1,1,2,4))
print(plot_grid(plotlist=plots,ncol = 1,align = "v",axis ="tb"))
dev.off()
}

