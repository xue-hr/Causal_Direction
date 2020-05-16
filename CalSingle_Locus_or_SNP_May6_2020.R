library(pracma)

# CD_3_methods -----------------------------------------------------------

CD_3_methods<-function(pruned,num_iteration = 20)
{
  SNP = pruned$loci_bed
  SNP = scale(SNP)
  sig_part = pruned$sig_part
  p = ncol(SNP)
  #
  N_T1 = sig_part[,8]
  N_T2 = sig_part[,11]
  T1_T = sig_part[,6] / sig_part[,7]
  T1_r = T1_T / sqrt(N_T1 - 2 + T1_T^2)
  T2_T = sig_part[,9] / sig_part[,10]
  T2_r = T2_T / sqrt(N_T2 - 2 + T2_T^2)
  #
  rho_T1 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T1[1:p,1:p] = cor(SNP)
  rho_T1[1:p,p+1] = T1_r
  rho_T1[p+1,1:p] = T1_r
  rho_T1[p+1,p+1] = 1
  #
  rho_T2 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T2[1:p,1:p] = cor(SNP)
  rho_T2[1:p,p+1] = T2_r
  rho_T2[p+1,1:p] = T2_r
  rho_T2[p+1,p+1] = 1
  #
  V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
  V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
  #
  CD_GLS_result = CD_GLS(pruned,V_T1 = V_T1,V_T2 = V_T2)
  
  CD_Egger_result = 
    CD_Egger(pruned,num_iteration,V_T1 = V_T1,V_T2 = V_T2)
  
  CD_Ratio_result = 
    CD_Ratio(pruned,num_iteration,V_T1 = V_T1,V_T2 = V_T2)
  
  return(list(CD_GLS_result = CD_GLS_result,
              CD_Egger_result = CD_Egger_result,
              CD_Ratio_result = CD_Ratio_result
              ))
}


# CD_GLS ------------------------------------------------------------------

CD_GLS <- function(pruned,V_T1 = NULL,V_T2 = NULL)
{
  ### "pruned" is a list, has "loci_bed" and "sig_part"
  ### loci_bed is the reference panel, n by p
  ### sig_part is the summary statistics, p by 12,
  ### columns are "chr,pos,rsid,A1,A2,beta_T1,se_T1,N_T1,beta_T2,se_T2,N_T2,loci"
  
  SNP = pruned$loci_bed
  sig_part = pruned$sig_part
  p = ncol(SNP)
  #
  N_T1 = sig_part[,8]
  N_T2 = sig_part[,11]
  T1_T = sig_part[,6] / sig_part[,7]
  T1_r = T1_T / sqrt(N_T1 - 2 + T1_T^2)
  T2_T = sig_part[,9] / sig_part[,10]
  T2_r = T2_T / sqrt(N_T2 - 2 + T2_T^2)
  #
  rho_T1 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T1[1:p,1:p] = cor(SNP)
  rho_T1[1:p,p+1] = T1_r
  rho_T1[p+1,1:p] = T1_r
  rho_T1[p+1,p+1] = 1
  #
  rho_T2 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T2[1:p,1:p] = cor(SNP)
  rho_T2[1:p,p+1] = T2_r
  rho_T2[p+1,1:p] = T2_r
  rho_T2[p+1,p+1] = 1
  #
  if(is.null(V_T1))
  {
    V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
  }
  if(is.null(V_T2))
  {
    V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
  }
  #
  
  #T1 to T2
  jacobian = cbind( diag(1/T1_r) , 
                    -diag(T2_r/T1_r^2) )
  combined_V = rbind( cbind(V_T2,matrix(0,ncol = p, nrow = p)) / mean(N_T2),
                      cbind(matrix(0,ncol = p, nrow = p),V_T1) / mean(N_T1)
  )
  V = jacobian %*% combined_V %*% t(jacobian)
  inv_V = solve(V, tol = 0)
  est_vec = T2_r/T1_r
  gls_est = sum(inv_V %*% est_vec) / sum(inv_V)
  gls_var = 1 / sum(inv_V)
  
  T1toT2 = c(gls_est,sqrt(gls_var))
  Q_T1toT2 = (est_vec - gls_est)%*%inv_V%*%(est_vec - gls_est)
  
  #T2 to T1
  jacobian = cbind( diag(1/T2_r) , 
                    -diag(T1_r/T2_r^2) )
  combined_V = rbind( cbind(V_T1,matrix(0,ncol = p, nrow = p)) / mean(N_T1),
                      cbind(matrix(0,ncol = p, nrow = p),V_T2) / mean(N_T2)
  )
  V = jacobian %*% combined_V %*% t(jacobian)
  inv_V = solve(V, tol = 0)
  est_vec = T1_r/T2_r
  gls_est = sum(inv_V %*% est_vec) / sum(inv_V)
  gls_var = 1 / sum(inv_V)
  
  T2toT1 = c(gls_est,sqrt(gls_var))
  Q_T2toT1 = (est_vec - gls_est)%*%inv_V%*%(est_vec - gls_est)
  
  return(list(T1toT2 = T1toT2,
              T2toT1 = T2toT1,
              Q_T1toT2 = Q_T1toT2,
              Q_T2toT1 = Q_T2toT1))
}






# CD_Egger ----------------------------------------------------------------

CD_Egger <- function(pruned,num_iteration = 1,V_T1 = NULL,V_T2 = NULL)
{
  SNP = pruned$loci_bed
  SNP = scale(SNP)
  sig_part = pruned$sig_part
  p = ncol(SNP)
  SIGMA_SNP = cov(SNP)
  SIGMA_INV = solve(SIGMA_SNP)
  v_coef = as.numeric(SIGMA_SNP%*%rep(1,p))
  SIGMA_SQUARE = SIGMA_SNP%*%SIGMA_SNP
  
  #
  N_T1 = sig_part[,8]
  N_T2 = sig_part[,11]
  T1_T = sig_part[,6] / sig_part[,7]
  T1_r = T1_T / sqrt(N_T1 - 2 + T1_T^2)
  T2_T = sig_part[,9] / sig_part[,10]
  T2_r = T2_T / sqrt(N_T2 - 2 + T2_T^2)
  #
  rho_T1 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T1[1:p,1:p] = cor(SNP)
  rho_T1[1:p,p+1] = T1_r
  rho_T1[p+1,1:p] = T1_r
  rho_T1[p+1,p+1] = 1
  #
  rho_T2 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T2[1:p,1:p] = cor(SNP)
  rho_T2[1:p,p+1] = T2_r
  rho_T2[p+1,1:p] = T2_r
  rho_T2[p+1,p+1] = 1
  #
  if(is.null(V_T1))
  {
    V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
  }
  if(is.null(V_T2))
  {
    V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
  }
  
  ### T1 to T2
  var_pleio = 0
  for(i in 1:num_iteration)
  {
    cov_mat = V_T2/mean(N_T2) + var_pleio*SIGMA_SQUARE
    M_tem = rbind(v_coef,T1_r)
    
    I_InfoMat = M_tem%*%solve(cov_mat,tol = 0)%*%t(M_tem)
    
    b0_K = solve(I_InfoMat,tol = 0)%*%M_tem%*%solve(cov_mat,tol = 0)%*%T2_r
    
    vv = 
      as.numeric(SIGMA_INV%*%(T2_r - b0_K[1]*v_coef - b0_K[2]*T1_r))
    var_pleio = estimate_sigma2(v = vv,
                                Vg = SIGMA_INV%*%
                                  (V_T2/mean(N_T2))%*%SIGMA_INV)
  }
  vec_tem = T2_r - b0_K[1]*v_coef - b0_K[2]*T1_r
  Q_T1toT2 = vec_tem%*%solve(cov_mat,tol = 0)%*%vec_tem
  T1toT2 = c(b0_K,sqrt(diag(solve(I_InfoMat,tol = 0))),var_pleio)
  
  ### T2 to T1
  var_pleio = 0
  for(i in 1:num_iteration)
  {
    cov_mat = V_T1/mean(N_T1) + var_pleio*SIGMA_SQUARE
    M_tem = rbind(v_coef,T2_r)
    
    I_InfoMat = M_tem%*%solve(cov_mat,tol = 0)%*%t(M_tem)
    
    b0_K = solve(I_InfoMat,tol = 0)%*%M_tem%*%solve(cov_mat,tol = 0)%*%T1_r
    
    vv = 
      as.numeric(SIGMA_INV%*%(T1_r - b0_K[1]*v_coef - b0_K[2]*T2_r))
    var_pleio = estimate_sigma2(v = vv,
                                Vg = SIGMA_INV%*%
                                  (V_T1/mean(N_T1))%*%SIGMA_INV)
  }
  vec_tem = T1_r - b0_K[1]*v_coef - b0_K[2]*T2_r
  Q_T2toT1 = vec_tem%*%solve(cov_mat,tol = 0)%*%vec_tem
  T2toT1 = c(b0_K,sqrt(diag(solve(I_InfoMat,tol = 0))),var_pleio)
  
  
  return(list(T1toT2 = T1toT2,
              T2toT1 = T2toT1,
              Q_T1toT2 = Q_T1toT2,
              Q_T2toT1 = Q_T2toT1))
}





# CD_Ratio ----------------------------------------------------------------

CD_Ratio <- function(pruned,num_iteration = 1,V_T1 = NULL,V_T2 = NULL)
{
  SNP = pruned$loci_bed
  SNP = scale(SNP)
  sig_part = pruned$sig_part
  p = ncol(SNP)
  SIGMA_SNP = cov(SNP)
  SIGMA_INV = solve(SIGMA_SNP)
  v_coef = as.numeric(SIGMA_SNP%*%rep(1,p))
  SIGMA_SQUARE = SIGMA_SNP%*%SIGMA_SNP
  #
  N_T1 = sig_part[,8]
  N_T2 = sig_part[,11]
  T1_T = sig_part[,6] / sig_part[,7]
  T1_r = T1_T / sqrt(N_T1 - 2 + T1_T^2)
  T2_T = sig_part[,9] / sig_part[,10]
  T2_r = T2_T / sqrt(N_T2 - 2 + T2_T^2)
  #
  rho_T1 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T1[1:p,1:p] = cor(SNP)
  rho_T1[1:p,p+1] = T1_r
  rho_T1[p+1,1:p] = T1_r
  rho_T1[p+1,p+1] = 1
  #
  rho_T2 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T2[1:p,1:p] = cor(SNP)
  rho_T2[1:p,p+1] = T2_r
  rho_T2[p+1,1:p] = T2_r
  rho_T2[p+1,p+1] = 1
  #
  if(is.null(V_T1))
  {
    V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
  }
  if(is.null(V_T2))
  {
    V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
  }
  #
  
  ###T1 to T2
  combined_V = 
    rbind(cbind(V_T2,matrix(0,ncol = p, nrow = p)) / mean(N_T2),
          cbind(matrix(0,ncol = p, nrow = p),V_T1) / mean(N_T1)
    )  
  b0 = 0
  var_pleio = 0
  for(i in 1:num_iteration)
  {
    jacobian = cbind( diag(1/T1_r) , 
                      -diag((T2_r-b0*v_coef)/T1_r^2) )
    Vg = jacobian%*%combined_V%*%t(jacobian)
    cov_mat = Vg + var_pleio*(diag(1/T1_r)%*%
                                SIGMA_SQUARE%*%
                                diag(1/T1_r))
    
    M_tepm = rbind(v_coef/T1_r,rep(1,p))
    
    I_InfoMat = M_tepm%*%solve(cov_mat,tol = 0)%*%t(M_tepm)
    
    b0_K = 
      c(solve(I_InfoMat,tol = 0)%*%M_tepm%*%
          solve(cov_mat,tol = 0)%*%(T2_r/T1_r))
    
    b0 = b0_K[1]
    vv = 
      as.numeric(SIGMA_INV%*%
                   diag(T1_r)%*%
                   ((T2_r-b0*v_coef)/T1_r - b0_K[2])
      )
    
    
    var_pleio = estimate_sigma2(v = vv,
                                Vg = 
                                  SIGMA_INV%*%
                                  diag(T1_r)%*%Vg%*%
                                  diag(T1_r)%*%SIGMA_INV)
    

  }
  T1toT2 = c(b0_K,sqrt(diag(solve(I_InfoMat,tol = 0))),var_pleio)
  vec_tem = (T2_r - b0_K[1]*v_coef) / T1_r - b0_K[2]
  Q_T1toT2 = vec_tem%*%solve(cov_mat,tol = 0)%*%vec_tem
  
  ###T2 to T1
  combined_V = 
    rbind(cbind(V_T1,matrix(0,ncol = p, nrow = p)) / mean(N_T1),
          cbind(matrix(0,ncol = p, nrow = p),V_T2) / mean(N_T2)
    )  
  b0 = 0
  var_pleio = 0
  for(i in 1:num_iteration)
  {
    jacobian = cbind( diag(1/T2_r) , 
                      -diag((T1_r-b0*v_coef)/T2_r^2) )
    Vg = jacobian%*%combined_V%*%t(jacobian)
    cov_mat = Vg + var_pleio*(diag(1/T2_r)%*%
                                SIGMA_SQUARE%*%
                                diag(1/T2_r))
    
    M_tepm = rbind(v_coef/T2_r,rep(1,p))
    
    I_InfoMat = M_tepm%*%solve(cov_mat,tol = 0)%*%t(M_tepm)
    
    b0_K = 
      c(solve(I_InfoMat,tol = 0)%*%M_tepm%*%
          solve(cov_mat,tol = 0)%*%(T1_r/T2_r))
    
    b0 = b0_K[1]
    vv = 
      as.numeric(SIGMA_INV%*%
                   diag(T2_r)%*%
                   ((T1_r-b0*v_coef)/T2_r - b0_K[2])
      )
    
    
    var_pleio = estimate_sigma2(v = vv,
                                Vg = 
                                  SIGMA_INV%*%
                                  diag(T2_r)%*%Vg%*%
                                  diag(T2_r)%*%SIGMA_INV)
    
    
  }
  T2toT1 = c(b0_K,sqrt(diag(solve(I_InfoMat,tol = 0))),var_pleio)
  vec_tem = (T1_r - b0_K[1]*v_coef) / T2_r - b0_K[2]
  Q_T2toT1 = vec_tem%*%solve(cov_mat,tol = 0)%*%vec_tem
  
  return(list(T1toT2 = T1toT2,
              T2toT1 = T2toT1,
              Q_T1toT2 = Q_T1toT2,
              Q_T2toT1 = Q_T2toT1))
  
}






# estimate_sigma2 ---------------------------------------------------------

estimate_sigma2 <- function(v,Vg)
{
  
  eigen_Vg = eigen(Vg)
  Q = t(eigen_Vg$vectors)
  eigen_values = eigen_Vg$values
  ###Now: Vg = t(Q)%*%diag(eigen_values)%*%Q
  
  x = Q%*%v
  left_end = max(min(c(x^2 - eigen_values)),max(-eigen_values)+1e-14)
  right_end = max(c(x^2 - eigen_values))
  f = function(r) {
    sum(-(x^2 / (eigen_values+r)^2)) + sum ( 1 / (eigen_values+r))
  }
  
  if(f(left_end)>0 | is.na(f(left_end)))
  {
    root = left_end
  } else {
    while(f(right_end)<0)
    {
      right_end = right_end+1
    }
    root_find = bisect(f, left_end, right_end,maxiter = 500)
    root = root_find$root
  }
  
  return(max(0,root))
  
  
}
