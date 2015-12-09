`ICMg.combined.wrapper` <-
function(L, X, Niter, N, M, Lindices, Nindices,
                          D, C, z, w, n, m, q, E,
                          alpha, beta, pm0, V0, V,
                          convl, convn, C.boost) {

  if (C.boost) {
    out <- .C("ICMgCombinedIteration",PACKAGE="ICMg.iterations",L=as.integer(L),X=as.double(X),
              Niter=as.integer(Niter),N=as.integer(N),
              M=as.integer(M),Lindices=as.integer(Lindices),
              Nindices=as.integer(Nindices),D=as.integer(D),
              C=as.integer(C),z=as.integer(z),w=as.integer(w),
              n=as.integer(n), m=as.integer(m),q=as.integer(q),
              E=as.double(E),alpha=as.double(alpha),
              beta=as.double(beta), pm0=as.double(pm0),
              V0=as.double(V0), V=as.double(V),
              convl=as.double(convl),convn=as.double(convn))
  } else {
    out <- ICMg.combined.iteration(L, X, Niter, N, M, Lindices, Nindices,
                        D, C, z, w, n, m, q, E,
                        alpha, beta, pm0, V0, V,
                        convl, convn)
  }
  return(out)
}

