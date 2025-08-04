getDiscordanceFlag <- function(df) {
   wcols <- match(df$ABX, names(df))
   wrows_list <- tapply(X = seq_along(wcols),
                        INDEX = wcols,
                        FUN = list)
   
   makeFlag <- function(vec) {
      stopifnot(length(vec) > 0L)
      vec_copy <- as.character(vec)
      if (any(is.na(vec)))
         vec_copy[is.na(vec)] <- 'Unknown'
      if (any(vec_copy == '0'))
         vec_copy[vec_copy == '0'] <- 'Concordant'
      if (any(vec_copy == '1'))
         vec_copy[vec_copy == '1'] <- 'Discordant'
      return(vec_copy)
   }
   
   df$FLAG <- NA_character_
   for (i in seq_along(wrows_list)) {
      abx <- names(df)[as.integer(names(wrows_list)[i])]
      w <- wrows_list[[i]]
      df$FLAG[w] <- makeFlag(df[[abx]][w])
   }
   
   w <- which(is.na(df$FLAG))
   if (length(w) > 0L)
      df$FLAG[w] <- 'Unknown'
   
   df <- df %>% relocate(FLAG, .after=ABX)
   return(df)
}
