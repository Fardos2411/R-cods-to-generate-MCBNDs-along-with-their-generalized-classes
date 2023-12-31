#################################################################################
# CSBGND_3diffsize: Circular Strongly Balance Generalized Neighbor Design for 
# block of three different sizes (K1,k2 and k3)
# Algorithm from paper:
# Akbar Firdos,Mahmood Ul Hassan,Farrukh Jamal,Hurria Ali,Khadija Noreen and Rashid Ahmed Algorithms to Construct Minimal Circular Strongly. 
# Coded by Fardos Akbar et al., 2023
#################################################################################




#################################################################################
# Selection of i group of size K1 from adjusted A. The set of remaining 
# (Unselected) elements are saved in object B2. 
#################################################################################

grouping1<-function(A,k,v,i){
  bs<-c()
  z=0;f=1
  A1=A
  while(f<=i){
    
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
  list(B1=bs,B2=A1)
}

#################################################################################
# Selection of i group of size K1 from adjusted A and selection of required 
# number of groups of size K2 from B2. The set of remaining (Unselected) 
# elements are saved in B3.
#################################################################################

grouping2<-function(A,k,v,i,sk2){
  bs1<-c()
  j=i+sk2
  z=0;f=1
  A1=A
  while(f<=j){
    s<-grouping1(A1,k[1],v,i)
    A2<-s$B2
    z=i;f=f+i
    for(y in 1:2000){
      comp<-sample(1:length(A2),k[2])
      com<-A2[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs1<-rbind(bs1,com)
        A2<-A2[-com]
        z<-z+1
        f=f+1
      }
      if(z==j) break
    }
    
    
    if(z<j) {bs1<-c();z=0;f=1;A1=A}  
    
  }
  
  list(B1=s$B1,B2=bs1,B3=A2)
}

#################################################################################
# Selection of i group of size K1 from adjusted A, selection of required number 
# of groups of size K2 from B2 and division of required number of groups of size 
# K3 from B3.
#################################################################################
grouping3<-function(A,k,v,i,sk2,sk3){
  bs1<-c()
  j=i+sk2+sk3
  z=0;f=1
  A1=A
  while(f<=j){
    s<-grouping2(A1,k,v,i,sk2)
    A3<-s$B3
    z=i+sk2;f=f+i+sk2
    for(y in 1:1000){
      comp<-sample(1:length(A3),k[3])
      com<-A3[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs1<-rbind(bs1,com)
        A3<-A3[-com]
        z<-z+1
        f=f+1
      }
      if(z==j) break
    }
    if(z<j) {bs1<-c();z=0;f=1;A1=A}  
  }
  
  gs1<-t(apply(s$B1,1,sort))
  gs1<-cbind(gs1,rowSums(gs1),rowSums(gs1)/v)
  rownames(gs1)<-paste("G",1:i, sep="")
  colnames(gs1)<-c(paste(1:k[1], sep=""),"sum" ,"sum/v")
  
  gs2<-t(apply(s$B2,1,sort))
  gs2<-cbind(gs2,rowSums(gs2),rowSums(gs2)/v)
  rownames(gs2)<-paste("G",(i+1):(i+sk2), sep="")
  colnames(gs2)<-c(paste(1:k[2], sep=""),"sum" ,"sum/v")
  
  
  gs3<-t(apply(bs1,1,sort))
  gs3<-cbind(gs3,rowSums(gs3),rowSums(gs3)/v)
  rownames(gs3)<-paste("G",(i+sk2+1):(i+sk2+sk3), sep="")
  colnames(gs3)<-c(paste(1:k[3], sep=""),"sum" ,"sum/v")
  
  
  fs1<-t(apply(s$B1,1,sort))
  fs1<-delmin(fs1)
  rownames(fs1)<-paste("S",1:i, sep="")
  colnames(fs1)<-rep("",(k[1])-1)
  
  fs2<-t(apply(s$B2,1,sort))
  fs2<-delmin(fs2)
  rownames(fs2)<-paste("S",(i+1):(i+sk2), sep="")
  colnames(fs2)<-rep("",(k[2])-1)
  
  
  fs3<-t(apply(bs1,1,sort))
  fs3<-delmin(fs3)
  rownames(fs3)<-paste("S",(i+sk2+1):(i+sk2+sk3), sep="")
  colnames(fs3)<-rep("",(k[3]-1))
  
  
  list(B1=list(fs1,fs2,fs3),B4=list(gs1,gs2,gs3),B5=A3)
}

#################################################################################
# Selection of adjusted A and the set(s) of shifts to obtain Circular Strongly 
# Balance Generalized Neighbor Design for three different block size.
##################################################################################

# D=9: Minimal CSBGNDs in which v/2 of the ordered pairs appear twice as neighbors
# D=10: Minimal CSBGNDs in which 3v/2 of the unordered pairs appear twice as neighbors 
#   K: Vector of three different block sizes 
#   i: Number of sets of shifts for K1
# Sk2: Number of sets of shifts for K2
# Sk3: Number of sets of shifts for K3


CSBGND_3diffsize<-function(k,i,D,sk2,sk3){
  if(length(k)>3 | length(k)<3){stop("length(k)=3")}
  if(any(k<=3)!=0) stop("k=Block size: Each block size must be greater than 3")
  if(i<=0) stop("i=must be a positive integer")
  if(k[1]<k[2] | k[2]<k[3] |  k[1]<k[3]  ) stop("k1>K2>K3")

  setClass( "stat_test", representation("list"))
  
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
cat("Following are required sets of shifts to obtain the 
minimal CPBND for", "v=" ,object$R[1], ",","k1=",object$R[2],
",","k2=",object$R[3],"and","k3=",object$R[4],"\n")
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
    cat("\n")
    print(object$S[[2]])
    cat("\n")
    print(object$S[[3]])
  })
  

if(D==9 & sk2==1 & sk3==1){  
    v=2*i*k[1]+2*k[2]+2*k[3]-2 ; m<-(v-2)/2
    
if(m%%4==2){
      A=0:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
if(m%%4==3){
      A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+2k2+2k3 and and k1,k2 k3 block sizes")}
}  
   
  
  if(D==10 & sk2==1 & sk3==1){
    v= 2*i*k[1]+2*k[2]+2*k[3]-4; m<-(v-2)/2
    
    if(m%%4==0){
      A=c(0:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)

    }
    
    if(m%%4==1){
      A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%4==2){
      A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%4==3){
      A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
}   
    
 
if(D==9 & sk2==2 & sk3==2){  
  v=2*i*k[1]+4*k[2]+4*k[3]-2 ; m=(v-2)/2
  
  if(m%%4==2){
    A=0:(m+1)
    A1<-grouping3(A,k,v,i,sk2,sk3)
    A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
    x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
  }
  
  if(m%%4==3){
    A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
    A1<-grouping3(A,k,v,i,sk2,sk3)
    A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
    x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
  }
  
  if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+4k2+4k3 and and k1,k2 k3 block sizes")}
}  


if(D==10 & sk2==2 & sk3==2){
  v= 2*i*k[1]+4*k[2]+4*k[3]-4  ;m<-(v-2)/2
  
  if(m%%4==0){
    A=c(0:(m-1),(m-1),(m+1),(m+2))
    A1<-grouping3(A,k,v,i,sk2,sk3)
    A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
    x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    
  }
  
  if(m%%4==1){
    A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
    A1<-grouping3(A,k,v,i,sk2,sk3)
    A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
    x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
  }
  
  if(m%%4==2){
    A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
    A1<-grouping3(A,k,v,i,sk2,sk3)
    A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
    x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
  }
  
  if(m%%4==3){
    A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
    A1<-grouping3(A,k,v,i,sk2,sk3)
    A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
    x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
  }
} 

  
  if(D==9 & sk2==1 & sk3==2){  
    v=2*i*k[1]+2*k[2]+4*k[3]-2 ; m=(v-2)/2
    
    if(m%%4==2){
      A=0:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%4==3){
      A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+2k2+4k3 and and k1,k2 k3 block sizes")}
  }  
  
  
  if(D==10 & sk2==1 & sk3==2){
    v= 2*i*k[1]+2*k[2]+4*k[3]-4  ;m<-(v-2)/2
    
    if(m%%4==0){
      A=c(0:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
      
    }
    
    if(m%%4==1){
      A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%4==2){
      A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%4==3){
      A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  } 
  
  if(D==9 & sk2==2 & sk3==1){  
    v=2*i*k[1]+4*k[2]+2*k[3]-2 ; m=(v-2)/2
    
    if(m%%4==2){
      A=0:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%4==3){
      A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+4k2+2k3 and and k1,k2 k3 block sizes")}
  }  
  
  
  if(D==10 & sk2==2 & sk3==1){
    v= 2*i*k[1]+4*k[2]+2*k[3]-4  ;m<-(v-2)/2
    
    if(m%%4==0){
      A=c(0:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
      
    }
    
    if(m%%4==1){
      A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%4==2){
      A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%4==3){
      A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
  } 
new("stat_test", x)  
}


###################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CSBGND_equalsize
# The output is called using the design_CSBGND to generate design
design_CSBGND<-function(H){
  
  setClass( "CSBGND_design", representation("list"))
  setMethod("show", "CSBGND_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal CSBGND for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
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
  new("CSBGND_design", x)
}


###############################################################################
# Examples: Using CSBGND_3diffsize function to obtain the set(s) of shifts
# for construction of Circular Strongly Balance Generalized Neighbor Design 
# for block of three different sizes (k1, K2 and k3)
###############################################################################
# D=9: Minimal CSBGNDs in which v/2 of the ordered pairs appear twice as neighbors
# D=10: Minimal CSBGNDs in which 3v/2 of the unordered pairs appear twice as neighbors 
#   K: Vector of three different block sizes 
#   i: Number of sets of shifts for K1
# Sk2: Number of sets of shifts for K2
# Sk3: Number of sets of shifts for K3
#	MCSBGNDs-I and MCSBGNDs-II, say D9 & D10 for v even


# Example#1
(H<-CSBGND_3diffsize(k=c(12,7,4),i=5,D=9,sk2=1,sk3=1))

#Example#2
(H<-CSBGND_3diffsize(k=c(8,5,4),i=6,D=10,sk2=1,sk3=1))
(design_CSBGND(H))










