write_log <- function(stage, success = T, f = "scAPA.script.log"){
  if(success){
    write(x = paste0(Sys.time(), "\t Finished ", stage),
          file = f, append = T)
  } else {
    write(x = paste0(Sys.time(), "\t", stage, " not end, shell script ",
                     "scAPA.script.R did NOT end"),
          file = f, append = T)
  }
}

write_log_start <- function(stage, command, f= "scAPA.script.log") {
  if(!is.na(command)){
  write(x = paste0(Sys.time(), "\t Starting ", stage, "\ncommand: ", command),
        file = f, append = T)
  } else {write(x = paste0(Sys.time(), "\t Starting ", stage),
             file = f, append = T)}
}

read_args <- function(arg = script.args, arg.string, defult = NULL){
  nchr.arg <- nchar(arg.string)
  arg.string <- paste0(arg.string, " .*")
  pos <- regexpr(pattern = arg.string, text  = arg)
  if(pos > 0){
    pos <- pos + nchr.arg + 1
    out <- substring(arg, first = pos)
    out <- substring(out, first = 1, last = (regexpr(pattern = " ",
                                                     text  = out) -1))
} else {
  out <- defult
}
  out
  }

write.bed <- function(.x,f){
  write.table(x = .x, file = f, sep ="\t",
              col.names = F, row.names = F,
              quote = F)
}
