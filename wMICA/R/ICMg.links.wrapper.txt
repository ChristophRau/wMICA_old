`ICMg.links.wrapper` <-
function(L, Niter, N, M, Lindices,
                        C, z, q, n,
                        alpha, beta, conv, C.boost) {

  if (C.boost) {
    out <- .C("ICMgLinksIteration",PACKAGE="ICMg.iterations", L=as.integer(L),Niter=as.integer(Niter),
              N=as.integer(N),M=as.integer(M),
              Lindices=as.integer(Lindices),C=as.integer(C),
              z=as.integer(z),q=as.integer(q),
              n=as.integer(n),alpha=as.double(alpha),
              beta=as.double(beta), conv=as.double(conv))
  } else {
    out <- ICMg.links.iteration(L, Niter, N, M, Lindices,
                      C, z, q, n, alpha, beta, conv)
  }
  return(out)
}

