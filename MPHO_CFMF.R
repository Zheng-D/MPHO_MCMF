
fitness1 <- function(r_list){
  k <- length(r_list)
  q_m <- matrix(0, 1, k+4)
  q_m[1,1:k] <- r_list
  r_list1 <- unique(r_list)
  if(length(r_list1)<length(r_list)){
    return(q_m)
    break
  }
  r <- which(r_list==0)
  if(length(r)==0){
    vect <-matrix(0,nrow=1,ncol = nrow(adj_matrix))
    vect[1,r_list] <- 1
    count1 <- 1
    q_w <- matrix(0,nrow = 1,ncol = k)
    A_adj <- adj_matrix[r_list,r_list]
    M_adj <- M[r_list,r_list]
    for (h in 1:length(r_list)) {
      q_w[1,h] <- sum(M_adj[h,])
    }
    adj_mean <- sum(A_adj)/(k*(k-1)) - 1
    q_m[1,k+2] <- (sum(q_w[1,])/(k*(k-1)))+adj_mean
    y <- mut_mtrix %*% t(vect)
    y_i <- matrix(0,nrow = nrow(y),ncol = 1)
    true_list <- which(y[,1]!=0)
    y_i[true_list,1] <- 1
    v1 <- sum(y_i[,1])/nrow(mut_mtrix)
    q_m[1,k+1] <- v1
    v2_up <- apply(y,2,j2)
    q_m[1,k+3] <-  sum(v2_up)/sum(y_i[,1])
    #q_m[1,k+3] <- sum(v2_up)/nrow(mut_mtrix)
    q_m[1,k+4] <- q_m[1,k+3]*(q_m[1,k+1] + q_m[1,k+2])
    #q_m[1,k+4] <- q_m[1,k+3]*q_m[1,k+1]  
    return(q_m)
  }else{
    return(q_m)
  }
}
j2 <- function(vect1){
  for (i in 1:length(vect1)) {
    j <- vect1[i]
    if(j>0){
      vect1[i] <- 1/exp(j-1)
    }
  }
  return(vect1)
}
 
getscore_m <- function(MS,pg1){
  for (g in 1:nrow(pg1)) {
    gscore1 <- pg1[g,K+4]
    #gscore2 <- pg1[g,K+2] + 1
    for (g1 in 1:K) {
      g2 <- pg1[g,g1]
      if(MS[4,g2]==0){
        if(MS[4,g2]<abs(gscore1)){
          MS[4,g2] <- gscore1
        }
      }else{
        if(MS[4,g2]<gscore1){
          MS[4,g2] <- gscore1
        }
      }
      
      
    }
    
  }
  return(MS)
}

NGA1 <- function(Pb1,Pb2){  
  #pnext <- matrix(0,2*popsize,K+5)
  p_all <- rbind(Pb1,Pb2)
  for (i in 1:nrow(p_all)) {
    p_all[i,1:K] <- mfo(p_all[i,])
  }
  #pnext[,1:K] <- p_all[,1:K]
  #pnext <- t(parApply(cl,pnext[,1:K],1,mfo))
  pnext <- t(parApply(cl,p_all[,1:K],1,fitness1))
  
  return(pnext)
}

mfo <- function(list1){
  CONS <- 2
  mfr1 <- 1
  #mfr1 <- list1[1]
  x1 <- list1[1:K]
  #fit_value <- list1[(K+5)]
  y1 <- target_pop[mfr1,1:K]
  v1 <- abs(y1 - x1)
  apha <- matrix(0, 1, K)
  mask1 <- matrix(0, 1, K)
  for (b1 in 1:K) {
    t <- 2*runif(1) - 1
    mask1[1,b1] <- exp(CONS*t)
    #apha[1,b1] <-  1-((runif(1))^(1/8))
    apha[1,b1] <-  cos(2*pi*t)
    n1 <- 0.8
    #n2 <- alp/2
    if(runif(1)>n1 & v1[b1]==0){
      v1[b1] <- rnorm(1,0,sqrt(n))
    }
  }
  x12 <- mask1*apha*v1
  x2 <- round(y1 + x12)
  flg <- runif(1)
  if(flg<alp){
    beixuan <- pool_list
  }else{
    beixuan <- 1:n
  }
  beixuan <- beixuan[!beixuan%in%x1] 
  
  for (b1 in 1:K) {
    if(x12[b1]>0){
      index1 <- which(beixuan>=x2[b1])
      if(length(index1)==0){
        x2[b1] <- sample(beixuan,1)
      }else{
        x2[b1] <- beixuan[index1[1]]
      }
      
    } 
    if(x12[b1]<0){
      index1 <- which(beixuan<=x2[b1])
      if(length(index1)==0){
        x2[b1] <- sample(beixuan,1)
      }else{
        x2[b1] <- beixuan[index1[length(index1)]]
      }
    }
    beixuan <- setdiff(beixuan,x2[b1])
  }
  return(x2)
  
}



crossover2 <- function(parent1,parent2,n){ 
  temp22=parent1
  temp <- matrix(0,1,n) 
  temp[parent1] <- 1
  parent1 <- temp
  temp <- matrix(0,1,n)
  temp[parent2] <- 1
  parent2 <- temp
  newpop <- matrix(0,1,n)
  index <- which((parent1+parent2)==2)  
  newpop[index] <- 1
  parent1[index] <- 0
  parent2[index] <- 0
  temp <- which(parent1+parent2==1) 
  index <- sample(1:sum(parent1+parent2==1),sum(parent1+parent2))
  newpop[temp[index[1:(K-sum(newpop))]]]=1
  newpop <- which(newpop==1)
  return(newpop)
}

GA3 <- function(tag1,temp,alp1,g_range){
  rd <- alp1 #There are two functions that generate random numbers, and they are runif(),rnorm()
  Kgene <- tag1[1:K]
  #if(runif(1)>rd){
  #beixuan <- 1:n
  #beixuan <- MS[1,which(MS[4,]!=0)]
  #}else{
  #beixuan <- pool_list
  #}
  beixuan2 <- MS[1,which(MS[4,]!=0)]
  sam_len <- g_range*length(beixuan2)
  beixuan <- sample(beixuan2,sam_len)
  if(runif(1)<rd){
    beixuan <- pool_list
  }
  
  beixuan <- beixuan[!beixuan%in%Kgene] 
  #sam_len <- floor((0.6+ jit)*length(beixuan))
  #sam_len <- min(sam_len,length(beixuan))
  #beixuan <- sample(beixuan,sam_len)
  #beixuan <- beixuan[!beixuan%in%Kgene] 
  #popTemp2=matrix(-1000,2,K+4)
  #popTemp2[1,] <- tag1[index1,]
  #temp <- sample(1:K,1)
  # ff=floor((n-k)/2)
  ff <- length(beixuan)
  # ff=n-k
  Kgene2 <- matrix(rep(Kgene,each=ff),ff,K)
  Kgene2[,temp] <- beixuan
  pop2=t(parApply(cl,Kgene2,1,fitness1))
  pop2 <- pop2[order(pop2[,K+4],decreasing = T),]
  #print(pop2[1,K+4])
  if(pop2[1,K+4]>tag1[K+4]){
    tag1 <- pop2[1,]
  }
  return(tag1)
}


crossover <- function(pc){
  pc2 <- matrix(0,popsize,(K+4))
  pknumber <- popsize*K
  beixuan1 <- MS[1,which(MS[4,]==0)]
  beixuan2 <- MS[1,which(MS[4,]!=0)]
  beixuan3 <- MS[1,which(MS[4,]>0)]
  if(length(beixuan1)<popsize){
    #beixuan <- sample(beixuan2,pknumber)
    beixuan <-  beixuan2
    if(length(beixuan3)>popsize){
      fit2 <- MS[4,beixuan3]
      index <- vector(length = K)
      popsize2 <- popsize/2
      while (length(index)<=popsize2) {
        index <- select_order_fitness2(fit2) 
      }
      beixuan <- beixuan3[index]
    }
    if(length(beixuan1)!=0){
      beixuan <- union(beixuan1,beixuan)
    }
    #beixuan4 <- sample(beixuan2,pknumber)
    #beixuan <- union(beixuan,beixuan4)
  }
  if(length(beixuan1)>=popsize){
    beixuan <- beixuan1
  }
  Kgene <- unique(as.vector(pc[,1:K]))
  
  beixuan <- beixuan[!beixuan%in%Kgene] 
  if(length(beixuan)<popsize){
    beixuan <- union(beixuan,beixuan2)
  }
  pc1 <- pc[!duplicated(pc[,K+4]),]
  if(length(unique(pc[,K+4]))==1){
    npS <- 1
  }else{
    npS <- nrow(pc1)
    
  }
  l0 <- popsize - npS
  if(runif(1)<alp){
    beixuan <- pool_list
  }
  if(l0>0){
    pc2[1:npS,] <- pc1[1:npS,]
    for (g1 in (npS+1):popsize) {
      pc2[g1,1:K] <- sample(beixuan,K)
      #pc[g1,] <- select_generate(beixuan,kenel)
    }
  }else{
    pc2 <- pc1
  }
  pc2 <- t(parApply(cl,pc2[,1:K],1,fitness1))
  pc2 <- pc2[order(pc2[,K+4],decreasing = 'T'),]
  pn1 <- matrix(0,popsize,K+4)
  for (i in 1:popsize) {
    index <- sample(1:nrow(pc2),K) 
    #index1 <- index[1]
    index1 <- min(index)
    s1 <- setdiff(1:nrow(pc2),index1)
    index2 <- min(sample(s1,K))
    #index2 <- index[2]
    pn1[i,1:K] <- crossover2(pc2[index1,(1:K)],pc2[index2,1:K],n)
    
  }
  pn1 <- t(parApply(cl,pn1[,1:K],1,fitness1))
  #pn2 <- rbind(pc,pn1)
  #pn2 <- pn2[order(pn2[,K+4],decreasing = 'T'),]
  #for (i in 1:pp2) {
  #drow <- which(pc[,K+4]<pn1[i,K+4])
  #if(length(drow)==0)next
  #dindex <- drow[1]
  #pc[dindex,] <- pn1[i,]
  #}
  return(pn1)
}



mutation_other1 <- function(tag1,tn,temp){
  #pop2 <- tag1[tn,]
  pknumber <- popsize*K
  b0 <- 0
  flg <- runif(1)
  Kgene=tag1[tn,1:K]
  rd <-  alp 
  #rd <- min(rd1,alp)
  if(flg<rd){
    beixuan <- pool_list
    index <- beixuan
  }
  if(flg>=rd){
    beixuan1 <- MS[1,which(MS[4,]>0)]
    beixuan2 <- MS[1,which(MS[4,]!=0)]
    if(length(beixuan1)<=popsize){
      beixuan <- beixuan2
    }
    if(length(beixuan1)>popsize){
      fit  <- MS[4,beixuan1]
      popsize2 <- popsize/2
      index <- numeric(popsize2)
      while(length(index)<=popsize2){
        index <- select_order_fitness2(fit) 
      }
      beixuan <- beixuan1[index]
    }
    
  } 
  if(length(beixuan)>pknumber){
    beixuan <- sample(beixuan,pknumber)
  }
  
  beixuan <- beixuan[!beixuan%in%Kgene]
  ff <- length(beixuan)
  #index1 <- index[1]
  #index1 <- sample(index,1)
  #pop2[temp] <- index1
  Kgene2 <- matrix(rep(Kgene,each=ff),ff,K)
  Kgene2[,temp] <- beixuan
  #res <- parApply(cl,Kgene2,1,fitness1)
  pop2=t(parApply(cl,Kgene2,1,fitness1))
  #I <- order(pop2[,K+4],decreasing = T)
  pop2 <- pop2[order(pop2[,K+4],decreasing = T),]
  #pop2 <- fitness1(pop2[1:K])
  if(pop2[1,K+4]>tag1[tn,K+4]){
    tag1[tn,] <- pop2[1,]
  }
  return(tag1[tn,])
}
GA1 <- function(pop,kenel){  
   
  pop_next <-  pop
  for (i in 1:nrow(pop_next)) {
    pop_next[i,] <- parell_mutation(pop_next[i,])
  }
  pop_next <- t(parApply(cl,pop_next[,1:K],1,fitness1))
  pop_next=pop_next[order(pop_next[,K+4],decreasing = 'T'),]
  return(pop_next)
}

parell_mutation <- function(list1){
  Kgene <- list1[1:K]
  list2 <- list1
  beixuan1 <- MS[1,which(MS[4,]==0)]
  beixuan2 <- MS[1,which(MS[4,]!=0)]
  if(length(beixuan1)<popsize){
    beixuan <- sample(beixuan2,popsize)
    if(length(beixuan1)!=0){
      beixuan <- union(beixuan1,beixuan)
    }
  }
  if(length(beixuan1)>=popsize){
    beixuan <- beixuan1
  }
  if(runif(1)<alp){
    beixuan <- pool_list
  } 
  beixuan <- beixuan[!beixuan%in%Kgene]
  m_number <- sample(1:K,1)
  m_number1 <- sample(1:K,m_number)
  beixuan <- sample(beixuan,m_number)
  for (pois in 1:m_number) {
    pos1 <- m_number1[pois]
    list2[pos1] <- beixuan[pois]
  }
  return(list2)
}



select_order_fitness2 <- function(fit_vector){
  n1 <- length(fit_vector) 
  p <- matrix(0,2,n1)
  p[1,] <- 1:n1
  p[2,] <- fit_vector
  
  p=p[,order(p[2,],decreasing = 'F')]
  s=sum(fit_vector)
  p[2,] <- p[2,]/s
  p[2,] <- cumsum(p[2,])
  #p_cumsum <- cumsum(p)#p=1,2,3,4;   p_cumsum=1,3,6,10
  random_data <- runif(1)
  temp <- which(p[2,]>=random_data)
  #index1 <- temp[1] 
  return(p[1,temp])
}



#####################PRAD
M2 <- read.csv("E:/dz/pathway/data/prad/PRAD_M2.csv")
rownames(M2) <- M2[,1]
M2 <- M2[,-1]
net_gene <- colnames(M2)

M <- M2
M <- as.matrix(M)

adj_matrix <- matrix(0,nrow = length(net_gene),ncol = length(net_gene))
rownames(adj_matrix) <- net_gene
colnames(adj_matrix) <- net_gene
for (i in 1:nrow(M)) {
  gene1_neib <- which(M[i,]!=0)
  adj_matrix[i,gene1_neib] <- 1
  
}
sum(adj_matrix)
G_mut <- read.csv("E:/dz/pathway/data/prad/prad_mutation.csv")
rownames(G_mut) <- G_mut[,1]
G_mut<-G_mut[,-1]
total_gene <- net_gene



###############GBM
M2 <- read.csv("E:/dz/pathway/data/gbm/GBM_M2.csv")
rownames(M2) <- M2[,1]
M2 <- M2[,-1]
net_gene <- colnames(M2)
M <- M2
M <- as.matrix(M)
adj_matrix <- matrix(0,nrow = length(net_gene),ncol = length(net_gene))
rownames(adj_matrix) <- net_gene
colnames(adj_matrix) <- net_gene
for (i in 1:nrow(M)) {
  gene1_neib <- which(M[i,]!=0)
  adj_matrix[i,gene1_neib] <- 1
  
}
sum(adj_matrix)
G_mut <- read.csv('E:/dz/pathway/data/gbm/SNVdata_440.csv')
rownames(G_mut) <- G_mut[,1]
G_mut<-G_mut[,-1]
total_gene <- net_gene

###########################OV
M2 <- read.csv("E:/dz/pathway/data/ov/OV_M2.csv")
rownames(M2) <- M2[,1]
M2 <- M2[,-1]
net_gene <- colnames(M2)
M <- M2
M <- as.matrix(M)
 
adj_matrix <- matrix(0,nrow = length(net_gene),ncol = length(net_gene))
rownames(adj_matrix) <- net_gene
colnames(adj_matrix) <- net_gene
for (i in 1:nrow(M)) {
  gene1_neib <- which(M[i,]!=0)
  adj_matrix[i,gene1_neib] <- 1
  
}
sum(adj_matrix)
G_mut <- read.csv('E:/dz/pathway/data/ov/OVCA_SNVdata_2547.csv')
rownames(G_mut) <- G_mut[,1]
G_mut<-G_mut[,-1]
total_gene <- net_gene







geneName <- total_gene

mut_mtrix <- matrix(0,nrow = nrow(G_mut),ncol = length(total_gene))
colnames(mut_mtrix) <- total_gene
rownames(mut_mtrix) <- rownames(G_mut)
for (i in 1:ncol(adj_matrix)) {
  gene1_G_coldex <- which(colnames(G_mut)==total_gene[i])
  mut_mtrix[,i] <- G_mut[,gene1_G_coldex]
}
#mut_mtrix <- as.matrix(mut_mtrix)
md_matrix <- matrix(0,3,ncol(mut_mtrix))
md_matrix[2,] <- colSums(mut_mtrix)
ol <- order(md_matrix[2,],decreasing = T)
md_matrix[1,] <- 1:ncol(mut_mtrix)
md_matrix <- md_matrix[,ol]

initial_gene <- matrix(0,nrow(adj_matrix),3)
initial_gene[,1] <- total_gene
for (i in 1:nrow(initial_gene)) {
  initial_gene[i,2] <- length(which(adj_matrix[i,]!=0))
}

write.csv(initial_gene,"initial_gene.csv")

result_gene <- read.csv("initial_gene.csv")
result_gene <- result_gene[order(result_gene[,3],decreasing = 'T'),]
initial_gene <- result_gene[1:200,2]
G_mut <- as.matrix(G_mut)








kenel <- 12
library(doParallel)
n <- ncol(adj_matrix)
var1 <- 0.2
cl <- makeCluster(getOption("cl.cores", kenel))      
registerDoParallel(cl)       #Process registration
clusterExport(cl, "adj_matrix") 
clusterExport(cl, "mut_mtrix") 
clusterExport(cl, "M")
clusterExport(cl, "total_gene")
clusterExport(cl, "n")
clusterExport(cl, varlist = c("j2"))  
clusterExport(cl, varlist = c("select_order_fitness2"))  
clusterExport(cl, "var1") 
clusterExport(cl, "G_mut") 

#cl_cores <- detectCores()
GIVEN_POOL=0
iteration <- 200
itR <- 10


#filename <- paste("brca_score_",ite,".csv") 



if(GIVEN_POOL==TRUE){
  alp=0.8
}else{
  alp=0
}
clusterExport(cl, "alp") 


for (K in 2:7) {
  #K <- 7
  clusterExport(cl, "K") 
  #3.Initial population
  popsize <- floor(log2(n^K))
  clusterExport(cl, "popsize")
  t1 <- Sys.time()   #start counting
  RS <- matrix(0,10,K+5)
  for(v in 1:10){
    #str=paste("Execute the ", v,"th time", sep = "")
    #print(str)
    j=1
    R=1
    t=0
    flag1 <- 0
    
    popsize2 <- 2*popsize
    if(GIVEN_POOL==TRUE){
      pool_list_mut <- as.numeric(md_matrix[1,1:floor(popsize/2)])
      pool_list_net <- as.numeric(result_gene[1:floor(popsize/2),1])
      pool_list <- union(pool_list_net,pool_list_mut)
      pool_list <- pool_list[order(pool_list)]
    } 
    
    #clusterExport(cl, "pool_list",envir = .GlobalEnv) 
    MS <- matrix(0,4,n)
    MS[1,] <- 1:n
    t1 <- Sys.time()  
    P1 <- matrix(0,popsize,K+4)
    if(GIVEN_POOL==TRUE){
      for(i in 1:popsize){ 
        P1[i,1:K] <- sample(pool_list,K)
      }
    }else{
      for(i in 1:popsize){ 
        P1[i,1:K] <- sample(1:n,K)
      }
    }
    
    P1 <- t(parApply(cl,P1[,1:K],1,fitness1))
    I1=P1[order(P1[,K+4],decreasing = 'T'),]
    MS <- getscore_m(MS,I1)
    P2 <- matrix(0,popsize,K+4)
    if(GIVEN_POOL==TRUE){
      for(i in 1:popsize){ 
        P2[i,1:K] <- sample(pool_list,K)
      }
    }else{
      for(i in 1:popsize){ 
        P2[i,1:K] <- sample(1:n,K)
      }
    }
    
    P2 <- t(parApply(cl,P2[,1:K],1,fitness1))
    I2=P2[order(P2[,K+4],decreasing = 'T'),]
    MS <- getscore_m(MS,I2)
    pp2 <- floor(popsize/2) + 1
    pp3 <- popsize-pp2
    target_pop <- matrix(-1000,popsize,K+4)
    target_pop[1:pp2,]=I1[1:pp2,]
    target_pop[(pp2+1):popsize,]=I2[1:pp3,]
    target_pop=target_pop[order(target_pop[,K+4],decreasing = 'T'),]
    pr1 <- rep(1,K)
    best_ind=target_pop[1,]
    while(j<=iteration & R<=itR) { 
      flag2 <- 0
      jit <- R/itR
      I1[nrow(I1),] <- best_ind
      I2[nrow(I2),] <- best_ind
      best_ind <- target_pop[which.max(target_pop[,K+4]),]
      if(GIVEN_POOL==TRUE){
        clusterExport(cl, "pool_list",envir = .GlobalEnv) 
      } 
      clusterExport(cl, "MS",envir = .GlobalEnv)
      I1 <- crossover(I1)
      MS <- getscore_m(MS,I1)
      I1=GA1(I1,kenel)
      MS <- getscore_m(MS,I1)
      I2 <- crossover(I2)
      MS <- getscore_m(MS,I2)
      I2=GA1(I2,kenel)
      MS <- getscore_m(MS,I2)
      target_pop[2:pp2,]=I1[1:(pp2-1),]
      target_pop[(pp2+1):popsize,]=I2[1:pp3,]
      
      if(jit>=0.3){
        clusterExport(cl, "target_pop",envir = .GlobalEnv) 
        I12 <- rbind(I1,I2)
        PT <- NGA1(I1,I2)
        MS <- getscore_m(MS,PT)
        for (i in 1:nrow(PT)) {
          if(PT[i,(K+4)]>I12[i,(K+4)]){
            I12[i,] <- PT[i,]
          }
        }
        I1 <- I12[1:popsize,]
        I2 <- I12[(popsize+1):nrow(I12),]
        I1 <- I1[order(I1[,K+4],decreasing = 'T'),]
        #MS <- getscore_m(MS,I3)
        I2 <- I2[order(I2[,K+4],decreasing = 'T'),]
        #MS <- getscore_m(MS,I4)
        target_pop[2:pp2,] <- I1[1:(pp2-1),]
        target_pop[(pp2+1):popsize,] <- I2[1:pp3,]
        target_pop <- target_pop[order(target_pop[,K+4],decreasing = 'T'),]
        #I1 <- I3
        #I2 <- I4
      }
      condition=(j%%3==0) || jit>=0.6
      if(condition){
        target_pop1 <- crossover(target_pop)
        MS <- getscore_m(MS,target_pop1)
        if(K>2){
          target_pop1 <- target_pop1[order(target_pop1[,K+4],decreasing = 'T'),]
          p0 <- which(pr1==0)
          if(length(p0)==K){
            pr1 <- rep(1,K)
          }
          psite <- sample(1:K,1,prob = pr1)
          pr1[psite] <- 0
          target_pop1[1,] <- GA3(target_pop1[1,],psite,alp,0.8)
          target_pop1 <- target_pop1[!duplicated(target_pop1[,K+4]),]
          if(length(unique(target_pop1[,K+4]))==1){
            np1 <- 1
          }else{
            np1 <- nrow(target_pop1)
          }
          MS <- getscore_m(MS,target_pop1)
          clusterExport(cl, "MS",envir = .GlobalEnv)
          M_number <-  round((1- j/iteration) * nrow(target_pop1))
          M_number <- max(sample(M_number,1),K)
          for (i4 in 2:M_number) {
            #var2 <- max(var1,tar_value1)
            var2 <-  var1
            #var3 <- max(var1,tar_value2)
            var3 <- var1
            if(runif(1)<var2){
              psite1 <- sample(1:K,1)
              target_pop1[i4,] <- mutation_other1(target_pop1,i4,psite1)
              flag2 <- 1
            }
            MS <- getscore_m(MS,target_pop1)
            
          }
          
          MS <- getscore_m(MS,target_pop1)
        }
        target_pop2 <- rbind(target_pop,target_pop1)
        target_pop2 <- target_pop2[order(target_pop2[,K+4],decreasing = 'T'),]
        target_pop <- target_pop2[1:popsize,]
        
        #I1 <- I3
        #I2 <- I4
      }
      pjit  <- K
      for (pj in 1:pjit) {
        rdindex <- sample(1:popsize,1)
        I1[rdindex,] <- target_pop[rdindex,]
      }
      for (pj in 1:pjit) {
        rdindex <- sample(1:popsize,1)
        I2[rdindex,] <- target_pop[rdindex,]
      }
      
      target_pop <- target_pop[order(target_pop[,K+4],decreasing = 'T'),]
      # MS <- getscore_m(MS,target_pop)
      if(GIVEN_POOL==TRUE){
        MS_max <- MS[4,]
        pool_1 <- MS[,order(MS_max,decreasing = 'T')]
        pool_1 <- pool_1[,which(pool_1[2,]>0)]
        pool_list  <- union(pool_list,pool_1[1,])
      } 
      
      #ll=length(which(best_ind[K+4]==best_ind1[K+4])==TRUE) #Determine whether it is equal
      if(best_ind[K+4]==target_pop[1,(K+4)]){
        R <- R+1
      }else{
        R <- 0
        pr1 <- rep(1,K)
      }
      #print(R)
      #print(geneName[target_pop[1,1:K]])
      j=j+1
    }
    t3 <- Sys.time() 
    print(t3-t1)
    
    maxpop <- target_pop
    maxpop[,1:K] <- apply(maxpop[,1:K,drop=T],2,function(x) {geneName[x]})
    
    #output the optimal solution.
    print(maxpop[1,1:K])
    print(maxpop[1,(K+4)])
    RS[v,1:(K+4)] <- maxpop[1,1:(K+4)]
    RS[v,(K+5)] <- t3-t1
    v=v+1
  }
  filename=paste("E:/dz/pathway/MPHO_MCMF_nopool12CORE_",K,".csv")
  write.csv(RS,file = filename)
}


stopCluster(cl)


