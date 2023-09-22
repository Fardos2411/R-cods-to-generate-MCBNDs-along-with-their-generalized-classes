#################################################################################
# CGND_2diffsize:Circular Generalized neighbor design for block of two different
# sizes (K1 and k2)
# Algorithm from paper:
# Akbar Firdos,Mahmood Ul Hassan,Farrukh Jamal,Hurria Ali,Khadija Noreen and Rashid Ahmed Algorithms to Construct Minimal Circular Strongly. 
# Coded by Fardos Akbar et al., 2023
#####################################################################################
# Selection of i groups of size K1 from adjusted A. The set of remaining (Unselected)
# elements are saved in the object named as B2. 
#####################################################################################
grouping1<-function(A,k,v,i){
  bs<-c()
  z=0;f=1
  A1=A
  while(f<=i){
    for(y in 1:5000){
      com<-sample(A1,k)
      cs<-sum(com)
      if(cs%%v==0){
        bs<-rbind(bs,com)
        A1<-A1[!A1 %in% com]
        z<-z+1
        f=f+1
      }
      if(z==i) break
    }
    if(z<i) {bs<-c();z=0;f=1;A1=A}  
   }
list(B1=bs,B2=A1)
}
#######################################################################################
# Selection of i group of size K1 from adjusted A and division of required number of 
# groups of size K2 from B2. 
#######################################################################################
grouping2<-function(A,k,v,i,sk2){
  bs1<-c()
  j=i+sk2
  z=0;f=1
  A1=A
  while(f<=j){
    s<-grouping1(A1,k[1],v,i)
    A2<-s$B2
    z=i;f=f+i
    for(y in 1:1000){
      com<-sample(A2,k[2])
      cs<-sum(com)
      if(cs%%v==0){
        bs1<-rbind(bs1,com)
        A2<-A2[!A2 %in% com]
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
  gs2<-t(apply(bs1,1,sort))
  gs2<-cbind(gs2,rowSums(gs2),rowSums(gs2)/v)
  rownames(gs2)<-paste("G",(nrow(gs1)+1):(nrow(gs1)+sk2), sep="")
  colnames(gs2)<-c(paste(1:k[2], sep=""),"sum" ,"sum/v")
  fs1<-t(apply(s$B1,1,sort))
  fs1<-delmin(fs1)
  rownames(fs1)<-paste("S",1:i, sep="")
  colnames(fs1)<-rep("",(k[1])-1)
  fs2<-t(apply(bs1,1,sort))
  fs2<-delmin(fs2)
  rownames(fs2)<-paste("S",(i+1):(i+sk2), sep="")
  colnames(fs2)<-rep("",(k[2]-1))
  list(B1=fs1,B2=fs2,B3=list(gs1,gs2),B4=A2)
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
  return(fs)
}
#################################################################################
# Selection of adjusted A and the set(s) of shifts to obtain Minimal 
# Circular Generalized neighbor design for two different block size.
#################################################################################
# D=3: minimal CGNDs in which v/2 unordered pairs do not appear
# D=4: minimal CGNDs in which 3v/2 unordered pairs do not appear 
#   K: Vector of two different block sizes
#   i: Number of sets of shifts for K1
# Sk2: Number of sets of shifts for K2
CPND_2diffsize<-function(k,i,D=3,sk2=1){
 if(length(k)>2 | length(k)<2){stop("length(k)=2 ")}
if(any(k<=3)!=0) stop("k=Block size: Each block size must be greater than 3")
if(i<=0) stop("i= Must be a positive integer")
if(k[1]<k[2]) stop("k1>K2")
setClass( "stat_test", representation("list"))
setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
cat("Following are required sets of shifts to obtain the 
minimal CPND for", "v=" ,object$R[1], ",","k1=",object$R[2],"and","k2=",object$R[3],"\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object[[1]])
    cat("\n")
    print(object[[2]])
  })
  if(D==3 & sk2==1){  
    v=2*i*k[1]+2*k[2]+2 ; m=(v-2)/2
    if(m%%4==0){
      A<-1:m
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+2 and k1= block size and k2= block size")}
  }  
  if(D==3 & sk2==2){
    v=2*i*k[1]+4*k[2]+2 ; m=(v-2)/2
    if(m%%4==0){
      A<-1:m
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
    if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+2 and k1= block size and k2= block size")}
  }  
  if(D==4 & sk2==1){
    v= 2*i*k[1]+2*k[2]+4 ; m=(v-2)/2 
    if(m%%4==0){
      A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(2:(m-2),m,(2*m+1))
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
  }
  if(D==4 & sk2==2){
    v= 2*i*k[1]+4*k[2]+4  ; m=(v-2)/2 
    if(m%%4==0){
      A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
        if(m%%4==2){
      A=c(2:(m-2),m,(2*m+1))
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
      if(m%%4==3){
      A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping2(A,k,v,i,sk2)
      A2<-c(v,k);names(A2)<-c("V","K1","K2")
      x<-list(S1=A1$B1,S2=A1$B2,G=A1$B3,R=A2,A=A)
    }
  }
  new("stat_test", x) 
} 
###############################################################################
# Examples: Using CGND_2diffsize function to obtain the set(s) of shifts
# for construction of Minimal Circular Generalized neighbor design for block of 
# two different sizes (k1 and k2)
###############################################################################
# D=3: minimal CGNDs in which v/2 unordered pairs do not appear
# D=4: minimal CGNDs in which 3v/2 unordered pairs do not appear 
#   K: Vector of two different block sizes
#   i: Number of sets of shifts for K1
# Sk2: Number of sets of shifts for K2
#MCPBNDs-I and MCPBNDs-II, say D3 & D4 for v even
# Example#1
(H<-CPND_2diffsize(k=c(6,3),i=5,D=3,sk2=2))
# Example#2
(H<-CPND_2diffsize(k=c(10,6),i=3,D=4,sk2=2))

