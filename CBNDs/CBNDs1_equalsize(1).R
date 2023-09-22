
################################################################################
# CQRNDs2_equalsize: Circular Quasi Rees neighbor designs for block of equal 
# size(K)
################################################################################

# Algorithm from paper:
# Akbar Firdos,Mahmood Ul Hassan,Farrukh Jamal,Hurria Ali,Khadija Noreen and Rashid Ahmed Algorithms to Construct Minimal Circular Strongly. 
# Coded by Fardos Akbar et al., 2023
################################################################################

################################################################
# Division of adjusted A in i-1 groups of size k and one
# of size K-2 to get the set(s) of shifts
################################################################
grouping1<-function(A,k,v,i){
  bs<-c()
  z=0;f=1
  A1=A
  while(f<=(i-1)){
    
    for(y in 1:5000){
      comp<-sample(1:length(A1),k)
      com<-A1[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs<-rbind(bs,com)
        A1<-A1[-comp]
        z<-z+1
        f=f+1
      }
      if(z==i) break
    }
    if(z<i) {bs<-c();z=0;f=1;A1=A}  
    
  }
  
  
  bs1<-t(apply(bs,1,sort))
  bs1<-cbind(bs1,rowSums(bs),rowSums(bs)/v)
  rownames(bs1)<-paste("G",1:i, sep="")
  colnames(bs1)<-c(paste(1:k, sep=""),"sum" ,"sum/v")
  
  bs2<-t(apply(bs,1,sort))
  bs2<-delmin(bs2)
  list(B1=list(bs2),B2=list(bs1),B3=A1)
}


#######################################################################
# Obtaing set(s) of shifts by deleting smallest value of each group
#######################################################################

delmin<-function(z){
  fs<-c()
  n<-nrow(z)
  c<-ncol(z)-1
  for(i in 1:n){
    z1<-z[i,]
    z2<-z1[z1!=min(z1)]
    fs<-rbind(fs,z2)
  }
  rownames(fs)<-paste("S",1:n, sep="")
  colnames(fs)<-rep("",c)
  return(fs)
}


################################################################################
# Selection of adjusted A and the set(s) of shifs to obtain Circular  
# balance neighbor design for block of equal size. 
################################################################################

# D=1: Circular Balanced Neighbor Designs
# D=2: Circular Strongly Balanced Neighbor Designs
#   K: Block sizes
#   i: Number of set of shifts for K


CBND_equalsize<-function(k,i,D=1){
  
  if(k<=2) stop("k= Block size: Block size must be greater than 2")
  if(i<=0) stop("i= Must be a positive integer")
  
  setClass( "stat_test", representation("list"))
  
  setMethod("show", "stat_test", function(object) {
    
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    if(D==2){
      cat("Following are required sets of shifts to obtain the 
minimal CSBND for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")}
    
    if(D==1){
      cat("Following are required sets of shifts to obtain the 
minimal CBND for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")}
    
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
  })
  
  if(D==2){  
    
    v=2*i*k-1; m=(v-1)/2
    
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(j-1),(j+1):m,(v-j))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
  }
  
  if(D==1){
    v=2*i*k+1; m=(v-1)/2
    
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(j-1),(j+1):m,(v-j))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    
    
    
  }
  new("stat_test", x)
  
}

##################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CWBND_equalsize
# The output is called using the design_CWBND to generate design

design_CBND<-function(H){
  
  setClass( "CBND_design", representation("list"))
  setMethod("show", "CBND_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal CBND for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    for(i in 1:length(ss)){
      W<-ss[[i]]
      nr<-dim(W)[1]
      for(j in 1:nr){
        print(object$Design[[i]][[j]])
        cat("\n\n")
      }}
  })  
  
  v<-H$R[1]
  k<-H$R[2]
  ss<-H$S  
  treat<-(1:v)-1
  fn<-(1:v)
  G<-list()
  
  
  for(j in 1:length(ss)){ 
    W<-ss[[j]]
    nr<-dim(W)[1]
    nc<-dim(W)[2]
    D<-list()
    
    for(i in 1:nr){
      dd<-c()
      d1<-matrix(treat,(nc+1),v,byrow = T)
      ss1<-cumsum(c(0,W[i,]))
      dd2<-d1+ss1
      dd<-rbind(dd,dd2)
      rr<-dd[which(dd>=v)]%%v
      dd[which(dd>=v)]<-rr
      colnames(dd)<-paste("B",fn, sep="")
      rownames(dd)<-rep("",(nc+1))
      fn<-fn+v
      D[[i]]<-dd
    }
    G[[j]]<-D
    
  }
  
  x<-list(Design=G,R=H$R)
  new("CBND_design", x)
}



################################################################################
# Examples: Using CBND_equal size function to obtain the set(s) of shifts
# for construction of circular balance neighbor design for equal block  
# sizes (k)
################################################################################
# D=1: Circular Balanced Neighbor Designs
# D=2: Circular Strongly Balanced Neighbor Designs
#   K: Block sizes
#   i: Number of set of shifts for K
#MCBNDs and MCSBNDs, say here D1 & D2 respectively, for v odd

# example#1
(H<-CBND_equalsize(k=5,i=2,D=1))
(D<-design_CBND(H))
# example#2
(H<-CBND_equalsize(k=4,i=4,D=2))
(D<-design_CBND(H))

