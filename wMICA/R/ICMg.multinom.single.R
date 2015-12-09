`ICMg.multinom.single` <-
function(prob) {
  cs <- cumsum(prob)
  which.max(runif(1) <= cs/cs[length(cs)])
}

