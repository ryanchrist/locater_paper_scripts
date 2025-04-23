
hogs <- function(n=NULL, envir = NULL){
  if(is.null(envir)){ envir <- globalenv() }
  # n is the number of large memory objects that we wish to see
  if(is.null(n)){
    sort(sapply(mget(ls(envir = envir),envir = envir),pryr::object_size),decreasing=T)
  }else{
    head(sort(sapply(mget(ls(envir = envir),envir = envir),pryr::object_size),decreasing=T),n)
  }
}