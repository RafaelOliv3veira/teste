# Tom Jones - It's Not Unusual
# i want you back



A <-diag(1)
B <-diag(c(1,1))

ncolA <- ncol(A)
nrowA <- nrow(A)
ncolB <- ncol(B)
nrowB <- nrow(B)

C <-matrix(0,ncol = ncolA*ncolB, nrow = nrowA*nrowB)

for(i in 1:nrowA){
  for(j in 1:ncolA){
    row_start = (i - 1) * nrowB + 1;
    row_end = (i - 1) * nrowB + nrowB;
    col_start = (j - 1) * ncolB  + 1;
    col_end = (j - 1) * ncolB + ncolB;
    
    C[row_start:row_end, col_start:col_end] <- A[i,j]*B
    
  }
}

C





#Produto de kronecker vetor matriz


A <-c(1,1)
B <-diag(c(1,1))

ncolA <- 1
nrowA <- 2
ncolB <- 2
nrowB <- 2

C <-matrix(0,ncol = ncolA*ncolB, nrow = nrowA*nrowB)

for(i in 1:nrowA){
    row_start = (i - 1) * nrowB + 1;
    row_end = (i - 1) * nrowB + nrowB;
    col_start =  1;
    col_end =  ncolB;
    
    C[row_start:row_end, col_start:col_end] <- A[i]*B
    
  
}

C

