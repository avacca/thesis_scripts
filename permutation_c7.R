#/usr/bin/Rscript

###
difference_matrix<-function(w,mino){ #compute the difference tp1 -tp2 between gene 1 and gene2 pairwise across all the genes in the list w , in the matrix mino
  	mx<-matrix(0,nrow=length(w),ncol=length(w))
  	colnames(mx)<-w
  	rownames(mx)<-w
   	for(i in 1:length(mino$gene_name)){
     		for(j in 1:length(mino$gene_name)){
       			if(j>=i){
         			indi<-mino[i,]$shuffled
         			indj<-mino[j,]$shuffled
         			mx[as.vector(mino[i,]$gene_name),as.vector(mino[j,]$gene_name)]=indi-indj;
 				mx[as.vector(mino[j,]$gene_name),as.vector(mino[i,]$gene_name)]=indj-indi}
       }
   }
  return(mx)	}

create_sig_mx <- function(a,aa){ #count in how many datasets the difference between gene 1 and gene 2 is negative or positive, assign the sign + if the majority of the datasets has positive difference and - if they have negative difference
sigmatrix<-matrix(NA,nrow=length(a),ncol=length(a))
colnames(sigmatrix)<-a
rownames(sigmatrix)<-a
   for (i in 1: length(a)){
 	for (j in 1: length(a)){
 		if(j>=i){
			signeg<-0
			sigpos<-0
 			for (k in 1:8){
 				if((a[i] %in% colnames(aa[[k]]))&(a[j] %in% colnames(aa[[k]]))){
					if((aa[[k]][a[i],a[j]]>0)|(aa[[k]][a[j],a[i]]<0)){ ###if gene in row peak after gene in column   
						signeg<-signeg-1}
					if((aa[[k]][a[i],a[j]]<0)|(aa[[k]][a[j],a[i]]>0)){###if gene in row peak before gene in column 
						sigpos<-sigpos+1}
			}}
		if ((sigpos==0) &(signeg==0)){
  			sigmatrix[a[i],a[j]]<-0
  			sigmatrix[a[j],a[i]]<-0} 		
 		if(sigpos>-signeg){
 			sigmatrix[a[i],a[j]]<-signeg
 			sigmatrix[a[j],a[i]]<- sigpos}
               if(-signeg>sigpos){
 			sigmatrix[a[i],a[j]]<-sigpos
 			sigmatrix[a[j],a[i]]<- signeg}
		if(-signeg==sigpos){
 			sigmatrix[a[i],a[j]]<-signeg
 			sigmatrix[a[j],a[i]]<-sigpos}
 	}}}
return(sigmatrix)
 }



create_connectivity_mx<-function(ax,sigmatrix,num){
sx<-matrix(0,nrow=length(ax),ncol=length(ax))
colnames(sx)<-ax
rownames(sx)<-ax
   for (i in 1: length(ax)){
 	for (j in 1:length(ax)){
		if(j>=i){
 			if(sigmatrix[i,j]<= - as.numeric(num)){
				sx[i,j]<-1}
			else if(sigmatrix[i,j]>=as.numeric(num)){
				sx[j,i]<-1}
			else if(sigmatrix[j,i]>=as.numeric(num)){
				sx[i,j]<-1}
			else if (sigmatrix[j,i]<= - as.numeric(num)){
				sx[j,i]<-1}}}}
return(sx)
}
###
args<-commandArgs(trailingOnly = TRUE)
library(plyr)

r<-read.table("all_RNAs_tp_log2FC_dataframe_allpeakinggenes_onlyfirstpeakingtime_alldatasets_nomitochondrialgenes.txt",sep="\t",h=T)
  a<-read.table("all_PCs_tp_log2FC_dataframe_allpeakinggenes_onlyfirstpeakingtime_alldatasets_nomitochondrialgenes.txt",sep="\t",h=T)
  ar<- rbind(a,r)
  minoar1<-ar[ar$dataset=="PMDM_LPS",c(1,2)][order(ar[ar$dataset=="PMDM_LPS",2]),]
  minoar2<-ar[ar$dataset=="PAC_FGF2",c(1,2)][order(ar[ar$dataset=="PAC_FGF2",2]),]
   minoar3<-ar[ar$dataset=="PAC_IL1B",c(1,2)][order(ar[ar$dataset=="PAC_IL1B",2]),]
   minoar4<-ar[ar$dataset=="MCF7_HRG",c(1,2)][order(ar[ar$dataset=="MCF7_HRG",2]),]
   minoar5<-ar[ar$dataset=="MCF7_EGF1",c(1,2)][order(ar[ar$dataset=="MCF7_EGF1",2]),]
  minoar6<-ar[ar$dataset=="PEC_VEGF",c(1,2)][order(ar[ar$dataset=="PEC_VEGF",2]),]
   minoar7<-ar[ar$dataset=="MPSC_MIX",c(1,2)][order(ar[ar$dataset=="MPSC_MIX",2]),]
   minoar8<-ar[ar$dataset=="SAOS2_OST",c(1,2)][order(ar[ar$dataset=="SAOS2_OST",2]),]
   minoar<-list(minoar1,minoar2,minoar3,minoar4,minoar5,minoar6,minoar7,minoar8)
 war<-count(unique(ar[,c(1,2)]), "gene_name")
 w<-war[war$freq>=7,]
 
w<-w[w$gene_name!="U1",] #eliminate multigene families
w<-w[w$gene_name!="U2",]
w<-w[w$gene_name!="U3",]
w<-w[w$gene_name!="U4",]
w<-w[w$gene_name!="7SK",]
w<-as.vector(w$gene_name)


t <- as.numeric(Sys.time())
seed <- 1e8 * (t - floor(t))
set.seed(seed); print(seed)

sigmatrix<-NULL
cgc<-NULL
bb<-NULL
mx<-NULL
gc7<-NULL
 ptm <- proc.time()
 for (m in 1:args[1]){
 	for (j in 1:8){
 		minoar[[j]]$shuffled<-sample(minoar[[j]]$tp)
 		minoar[[j]]<-minoar[[j]][minoar[[j]]$gene_name %in% as.vector(w),]
 		bb[[j]]<-difference_matrix(as.vector(w),minoar[[j]])
 	 	bb[[j]][is.na(bb[[j]])]<-0}
 	sgx<-create_sig_mx(as.vector(w),bb)
 	c_matrix7<-create_connectivity_mx(as.vector(w),sgx,7)
 	gc7<-c(gc7,sum(c_matrix7))
 }

 proc.time() - ptm
 
 write.table(gc7,file=paste(seed,"seed_conservedorderconnections.txt",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
