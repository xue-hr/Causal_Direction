# Description -------------------------------------------------------------

# Calculate the limiting distribution of sample correlation matrix,
# then get the V_{Xg} and V_{Yg}
# Use Theorem 2 in "The Asymptotic Variance Matrix of the Sample Correlation Matrix"
# by H. Neudecker and A. M. Wesselman

# calculate asymptotic covariance matrix ----------------------------------
calculate_asymptotic_variance <- function(SNP,rho)
{
  # This is the main function
  # SNP is samplesize by n
  # rho is (n+1) by (n+1), the correlation matrix of n SNPs with X
  
  n = ncol(SNP)
  SNP = scale(SNP)
  SNP = SNP * sqrt(nrow(SNP)) / sqrt(nrow(SNP)-1)
  M_s = generate_M_s(n + 1)
  M_d = generate_M_d(n + 1)
  
  M_1 = diag((n+1)^2) - M_s %*% ( kronecker( diag(n+1) , rho) ) %*% M_d
  
  V = generate_V(SNP,rho)
  
  asymp_cov = M_1 %*% V %*% t(M_1)
  target_ind1 = n * (n+1) + 1
  target_ind2 = (n+1)^2-1
  
  return( asymp_cov[(target_ind1:target_ind2) , (target_ind1:target_ind2)] )
}

# Matrix M_s --------------------------------------------------------------

generate_M_s <- function(n)
{
  # For n, generate matrix M_s in formula 2.9
  K = matrix(0,n^2,n^2)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      E_ij = matrix(0,n,n)
      E_ij[i,j] = 1
      K = K + kronecker(E_ij,t(E_ij))
    }
  }
  M_s = 1/2 * (diag(n^2) + K)
  return(M_s)
}


# Matrix M_d --------------------------------------------------------------

generate_M_d <- function(n)
{
  # For n, generate matrix M_d in formula 2.13
  K = matrix(0,n^2,n^2)
  for(i in 1:n)
  {
    E_ii = matrix(0,n,n)
    E_ii[i,i] = 1
    K = K + kronecker(E_ii,E_ii)
  }
  
  return(K)
}



# Matrix V ----------------------------------------------------------------

generate_V <- function(SNP,rho)
{
  # generate matrix V from formula 3.6
  # rho is (n+1) by (n+1) correlation matrix: n SNPs and X
  # SNP is standardized, of size sample_size by n
  
  n = ncol(SNP)
  Sigma = rho[1:n,1:n]
  inv_Sigma = solve(Sigma,tol = 0)
  rho_X = as.matrix(rho[1:n,(n+1)])
  alpha =  inv_Sigma %*% rho_X
  e2 = 1 - t(rho_X) %*% inv_Sigma %*% rho_X
  
  SNP_alpha = SNP %*% diag(as.numeric(alpha))
  S = rowSums(SNP_alpha)
  
  V = matrix(0,(n+1)^2,(n+1)^2)
  for(i in 1:(n+1)^2)
  {
    for(j in i:(n+1)^2)
    {
      ind_A1 = floor( (i-0.1) / (n+1) ) + 1
      ind_A2 = floor( (j-0.1) / (n+1) ) + 1
      ind_B1 = i - (ind_A1-1) * (n+1)
      ind_B2 = j - (ind_A2-1) * (n+1)
      ind_v = sort(c(ind_A1,ind_A2,ind_B1,ind_B2))
      if( sum(ind_v == (n+1)) == 4)
      {
        V[i,j] = calculate_xxxx(ind_v, SNP, S,e2)
      }
      if( sum(ind_v == (n+1)) == 3)
      {
        V[i,j] = calculate_xxx(ind_v, SNP, S,e2)
      }
      if( sum(ind_v == (n+1)) == 2)
      {
        V[i,j] = calculate_xx(ind_v, SNP, S,e2)
      }
      if( sum(ind_v == (n+1)) == 1)
      {
        V[i,j] = calculate_x(ind_v, SNP, S,e2)
      }
      if( sum(ind_v == (n+1)) == 0)
      {
        V[i,j] = calculate_nox(ind_v, SNP, S,e2)
      }
      V[j,i] = V[i,j]
    }
  }
  
  vSigma = as.matrix(c(rho))
  V = V - vSigma%*%t(vSigma)
  
  return(V)
}



# calculate_xxxx ----------------------------------------------------------
calculate_xxxx <- function(ind_v, SNP, S,e2)
{
  return( mean(S^4) + 6 * mean(S^2) * e2 + 3 * e2^2 )
}


# calculate_xxx -----------------------------------------------------------
calculate_xxx <- function(ind_v, SNP, S,e2)
{
  ind1 = ind_v[1]
  return( mean(S^3 * SNP[,ind1]) + 3 * mean(S * SNP[,ind1]) * e2 )
}


# calculate_xx ------------------------------------------------------------
calculate_xx <- function(ind_v, SNP, S,e2)
{
  ind1 = ind_v[1]
  ind2 = ind_v[2]
  return( mean(S^2 * SNP[,ind1] * SNP[,ind2]) + mean(SNP[,ind1] * SNP[,ind2]) * e2 )
}


# calculate_x -------------------------------------------------------------
calculate_x <- function(ind_v, SNP, S,e2)
{
  ind1 = ind_v[1]
  ind2 = ind_v[2]
  ind3 = ind_v[3]
  return( mean(S * SNP[,ind1] * SNP[,ind2] * SNP[,ind3]) )
}


# calculate_nox -----------------------------------------------------------
calculate_nox <- function(ind_v, SNP, S,e2)
{
  ind1 = ind_v[1]
  ind2 = ind_v[2]
  ind3 = ind_v[3]
  ind4 = ind_v[4]
  return( mean( SNP[,ind1] * SNP[,ind2] * SNP[,ind3] * SNP[,ind4] ) )
}


