`ICMg.links.wrapper.Weights` <-
function(L, Niter, N, M, Lindices,
                        C, z, q, n,
                        alpha, beta, conv, C.boost) {

  if (C.boost) {
   # W=L[,3]
    #print(W[1:10])
   # L=L[,-3]
   # print(L[1:10,])
    out <- .C("ICMgLinksIterationWeights",PACKAGE="ICMg.iterations", L=as.double(L),Niter=as.integer(Niter),#W=as.double(W),
              N=as.integer(N),M=as.integer(M),
              Lindices=as.integer(Lindices),C=as.integer(C),
              z=as.integer(z),q=as.double(q),
              n=as.integer(n),alpha=as.double(alpha),
              beta=as.double(beta), conv=as.double(conv))
  } else {
    out <- ICMg.links.iteration.Weights(L, Niter, N, M, Lindices,
                      C, z, q, n, alpha, beta, conv)
  }
  
  return(out)
}