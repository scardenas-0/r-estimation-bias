lx_y <- function(x, y, ks, Rs) {
  exp(
    lgamma(y + ks*x) - 
      lgamma(y + 1) - 
      lgamma(ks*x) + 
      ks * x * log(ks/(Rs + ks)) + 
      y * log(Rs/(Rs + ks))
  )
}
l_i_n_c <- function(i, n, ks, Rs) {
  (i / n) * lx_y(n, n-i, ks, Rs)
}
p_i <- function(i, R_p) {
  ((R_p - 1)/R_p) ^ (i - 1) / R_p
}

Cs_p_fn <- function(Rp, Rs, pobs, ks, jmax = 200, nmax = 200) {
  n_vals <- 2:max(nmax, jmax)
  i_vals <- 1:(jmax - 1)

  l_mat <- outer(i_vals, n_vals, Vectorize(function(i, n) l_i_n_c(i, n, ks, Rs)))
  
  p_vec <- p_i(i_vals, Rp)
  mat <- sweep(l_mat, 1, p_vec, `*`)
  
  cs_n <- apply(mat, 1, function(row) rev(cumsum(rev(row))))
  cs_n <- t(cs_n)

  rj_vals <- sapply(2:jmax, function(j) {
    i_range <- 1:(j - 1)
    n_index <- j - min(n_vals) + 1
    sum(cs_n[i_range, n_index]) * (1 - Rs) / (Rp * Rs)
  })
  
  j_vals <- 2:jmax
  sum(rj_vals * pobs * (1 - pobs) ^ (j_vals - 1))
}
Cs_p_fn(1.5, 0.9, 0.8, 1.0, 400, 400)