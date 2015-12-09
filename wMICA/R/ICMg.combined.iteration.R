`ICMg.combined.iteration` <-
function(L, X, Niter, N, M, Lindices, Nindices,
                        D, C, z, w, n, m, q, E,
                        alpha, beta, pm0, V0, V,
                        convl, convn) {
  V0i <- V0^(-1)
  Vi <- V^(-1)
  
  ## Main iteration loop
  for (s in 1:Niter) {
    cat(".")
    
    ## Link loop
    for (li in 1:N) {
      
      l = Lindices[li]
      i = L[l,1]
      j = L[l,2]
      
      ## Subtract the contribution of the link from the counts
      q[z[l],i] <- q[z[l],i]-1
      q[z[l],j] <- q[z[l],j]-1
      n[z[l]] <- n[z[l]]-1

      ## Loop for computing probabilities for the components to be sampled
      uz <- vector("numeric", C)
      for (p in 1:C) {
        A <- (n[p] +m[p] +alpha)
        B <- (q[p,i] +beta)*(q[p,j] +beta) / ((2*n[p] +m[p] +1 +M*beta)*(2*n[p] +m[p] +M*beta))
        uz[p] <- A *B
      }
            
      ## Draw a new component for the links and update the counts */
      newz <- ICMg.multinom.single(uz)
      convl[l] <- uz[newz]/sum(uz)
      n[newz] <- n[newz]+1
      q[newz,i] <- q[newz,i]+1
      q[newz,j] <- q[newz,j]+1
      z[l] <- newz
    }

    ## Node loop
    for (ki in 1:M) {
      
      k = Nindices[ki]
      
      ## Subtract the contribution of the node from the counts
      q[w[k],k] <- q[w[k],k]-1
      m[w[k]] <- m[w[k]]-1
      E[w[k], ] <- E[w[k], ] - X[k, ]

      ## Loop for computing probabilities for the components to be sampled
      uw <- vector("numeric", C)
      AB <- vector("numeric", C)
      H1 <- vector("numeric", C)
      He <- vector("numeric", C)
      
      for (p in 1:C) {
        
        AB[p] <- (n[p] +m[p] +alpha) * (q[p,k] +beta) / (2*n[p] +m[p] +(M-1)*beta)

        S <- (V0i + (m[p]+1)*Vi)^(-1)
        Sd <-(V0i + m[p]*Vi)^(-1)     
        
        A <- as.vector(S * (V0i * pm0 + Vi * E[p, ] + Vi * X[k, ]))
        Ad <- as.vector(Sd * (V0i * pm0 + Vi * E[p, ]))
      
        H1[p] <- ( S / Sd )^(D/2)  
        
        H31e <- crossprod(A) / S
        H32e <- crossprod(Ad) / Sd
        He[p] <- 1/2*(H31e -H32e)
      }

      for (p in 1:C) 
        uw[p] <- AB[p] *H1[p] * exp(He[p] - max(He)) 
      
      ## Draw a new component for the links and update the counts */
      neww <- ICMg.multinom.single(uw)
      convn[k] <- uw[neww]/sum(uw)
      m[neww] <- m[neww]+1
      q[neww,k] <- q[neww,k]+1
      E[neww, ] <- E[neww, ] + X[k, ]
      w[k] <- neww
    }
  }
  cat("\n")
  
  return(list(z=z, w=w, n=n, m=m, q=q, E=E, convl=convl, convn=convn))
}

