## set and get states for testing; something like:
##
## ## keep old versions for testing; use something like
## ##     PoARX:::.test_state$push("a", "PoARXFilter")
## ## to set a state before calling the function, maybe in the global env.
## ##      PoARX:::.test_state$get("PoARXFilter")
##
##     if(PoARX:::.test_state$get("PoARXFilter") == "0.3.2"){
##         old code
##     }else{
##         new code
##     }
##     (could have multiple else if's
.test_state <- local({
  states <- list()
  push <- function(state, funname = "global"){
    if(is.null(states[[funname]]))
      states[[funname]] <<- ""
    states[[funname]] <<- c(state, states[[funname]])
  }

  get <- function(funname = "global"){
    if(is.null(states[[funname]]))
      states[[funname]] <<- ""
    states[[funname]][1]
  }

  pop <- function(state, funname = "global"){
    val <- get(funname)
    states[[funname]] <<- states[[funname]][-1]
    val
  }
  environment()
})
