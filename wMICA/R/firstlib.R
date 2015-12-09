.First.lib <- function(lib,pkg)
{
   library.dynam("ICMg.iterations",pkg,lib)
   cat("ICMg loaded\n")
}
