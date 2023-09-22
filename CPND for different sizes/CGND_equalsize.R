##################################################################################
# CGND_equalsize: Circular Generalized Neighbor Design for block of equal size(K)

# Algorithm from paper:
# Akbar Firdos,Mahmood Ul Hassan,Farrukh Jamal,Hurria Ali,Khadija Noreen and Rashid Ahmed Algorithms to Construct Minimal Circular Strongly. 
# Coded by Fardos Akbar et al., 2023
#################################################################################




################################################################
# Division of adjusted A in i groups to get the set(s) of shifts
################################################################
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
  bs1<-t(apply(bs,1,sort))
  bs1<-cbind(bs1,rowSums(bs),rowSums(bs)/v)
  rownames(bs1)<-paste("G",1:i, sep="")
  colnames(bs1)<-c(paste(1:k, sep=""),"sum" ,"sum/v")
  bs2<-t(apply(bs,1,sort))
  bs2<-delmin(bs2)
  list(B1=bs2,B2=bs1,B3=A1)
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
####################################################################################
# Selection of adjusted A and the set(s) of shifts to obtain Minimal Circular Generalized 
# neighbor design for block of equal size. 
##################################################################################
# D=3: minimal CGNDs in which v/2 unordered pairs do not appear
# D=4: minimal CGNDs in which 3v/2 unordered pairs do not appear 
#   K: Block sizes
#   i: Number of sets of shifts for K
CGND_equalsize<-function(k,i,D=3){
if(k<=3) stop("k= Block size: Block size must be greater than 3")
if(i<=0) stop("i= Must be a positive integer")
setClass( "stat_test", representation("list"))
setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
cat("Following are required sets of shifts to obtain the 
minimal CGND for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object[[1]])
  })

if(D==3){  

v=2*i*k+2; m=(v-2)/2

if(m%%4==0){
   A=1:m
   A1<-grouping1(A,k,v,i)
   A2<-c(v,k);names(A2)<-c("V","K")
   x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%4==3){
  A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
  A1<-grouping1(A,k,v,i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
 }

if(m%%4==1 |  m%%4==2){return("The minimal CGNDs in which v/2 unordered pairs cannot be constructed for v=2ik+2 and k=block size")}
}

if(D==4){
v=2*i*k+4; m=(v-2)/2


if(m%%4==0){
  A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
  A1<-grouping1(A,k,v,i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%4==1){
  A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
  A1<-grouping1(A,k,v,i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%4==2){
  A=c(2:(m-2),m,(2*m+1))
  A1<-grouping1(A,k,v,i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%4==3){
  A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
  A1<-grouping1(A,k,v,i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

}
new("stat_test", x)

}


##################################################################################
# Examples: Using CGND_equalsize function to obtain the set(s) of shifts
# for construction of Minimal circular Generalized neighbor design for equal 
# block sizes (k)
##################################################################################

#D=3: minimal CGNDs in which v/2 unordered pairs do not appear
#D=4: minimal CGNDs in which 3v/2 unordered pairs do not appear 
#   K: Block sizes
#   i: Number of sets of shifts for K
#MCPBNDs-I and MCPBNDs-II, say D3 & D4 for v even

# example#1
(H<-CGND_equalsize(k=12,i=6,D=3))
# example #2
(H<-CGND_equalsize(k=11,i=5,D=4))




