`ICMg.links.iteration.Weights` <-
function(L, Niter, N, M, Lindices,
                       C, z, q, n, alpha, beta, conv) {
  ## Main iteration loop
  for (s in 1:Niter) {
    cat(".")
    ## Sample new component for each link
    for (li in 1:N) {

      l = Lindices[li]
      i = L[l,1]
      j = L[l,2]
	 weight=L[l,3]
      
      ## Subtract the contribution of the link from the counts
      q[z[l],i] <- q[z[l],i]-weight
      q[z[l],j] <- q[z[l],j]-weight
      n[z[l]] <- n[z[l]]-1

      ## Loop for computing probabilities for the components to be sampled
      uz <- vector("numeric", C)
      for (p in 1:C) {
        A <- (n[p] +alpha)
        B <- (q[p,i] +beta)*(q[p,j] +beta) / ((2*n[p] +1 +M*beta)*(2*n[p] +M*beta))
        uz[p] <- A *B
      }
            
      ## Draw a new component for the links and update the counts */
      newz <- ICMg.multinom.single(uz)
      conv[l] <- uz[newz]/sum(uz)
      n[newz] <- n[newz]+1
      q[newz,i] <- q[newz,i]+weight
      q[newz,j] <- q[newz,j]+weight
      z[l] <- newz
    }
  }
  cat("\n")
  return(list(z=z, q=q, n=n, conv=conv))
}
