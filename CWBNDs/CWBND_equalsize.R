##########################################################################################
# CWBND_equalsize: Circular Weakly balance neighbour design for block of equal size(K)
#################################################################################
# CWBND_equalsize: Circular Weakly balance neighbour design for block of equal 
# size(K)

# Algorithm from paper:
# Akbar Firdos,Mahmood Ul Hassan,Farrukh Jamal,Hurria Ali,Khadija Noreen and Rashid Ahmed Algorithms to Construct Minimal Circular Strongly. 
# Coded by Fardos Akbar et al., 2023
#################################################################################




################################################################
# Division of adjusted A in i groups to get the set(s) of shifs
################################################################
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


####################################################################################
# Selection of adjusted A and the set(s) of shifs to obtain Circular weakly 
# balance neighbour design for block of equal size. 
##################################################################################

# D=5: Minimal CWBNDs in which v/2 of the ordered pairs appear twice as neighbors
# D=6: Minimal CWBNDs in which 3v/2 of the unordered pairs appear twice as neighbors
#   K: Block sizes
#   i: Number of set of shifts for K


CWBND_equalsize<-function(k,i,D=5){
  
if(k<3) stop("k= Block size: Block size must be greater than 2")
if(i<=0) stop("i= Must be a positive integer")

setClass( "stat_test", representation("list"))
  
setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
cat("Following are required sets of shifts to obtain the 
minimal CWBND for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
  })

if(D==5){  

v=2*i*k; m=(v-2)/2

if(m%%4==2){
   A=1:(m+1)
   A1<-grouping1(A,k,v,i)
   A2<-c(v,k);names(A2)<-c("V","K")
   x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%4==3){
  A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
  A1<-grouping1(A,k,v,i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
 }

if(m%%4==1 |  m%%4==0){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik and k=block size")}
}

if(D==6){
v=2*i*k-2; m=(v-2)/2


if(m%%4==0){
  A=c(1:(m-1),(m-1),(m+1),(m+2))
  A1<-grouping1(A,k,v,i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%4==1){
  A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
  A1<-grouping1(A,k,v,i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%4==2){
  A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
  A1<-grouping1(A,k,v,i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%4==3){
  A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
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

design_CWBND<-function(H){
  
  setClass( "CWBND_design", representation("list"))
  setMethod("show", "CWBND_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal CWBND for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
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
  new("CWBND_design", x)
}



##################################################################################
# Examples: Using CWBND_equalsize function to obtain the set(s) of shifs
# for construction of circular weakly balance neighbour design for equal block  
# sizes (k)
##################################################################################
# D=5: Minimal CWBNDs in which v/2 of the ordered pairs appear twice as neighbors
# D=6: Minimal CWBNDs in which 3v/2 of the unordered pairs appear twice as neighbors
#	MCWBNDs-I and MCWBNDs-II, say D5 & D6 for v even

# example#1
(H<-CWBND_equalsize(k=4,i=2,D=6))
(D<-design_CWBND(H))


# example #2
(H<-CWBND_equalsize(k=10,i=1,D=2))
design_CWBND(H)



# example #3
(H<-CWBND_equalsize(k=4,i=2,D=2))
design_CWBND(H)


