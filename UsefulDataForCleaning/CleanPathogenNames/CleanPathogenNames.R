cleanPathogenNames <- function(df) {
   bugs <- read.csv(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/CleanPathogenNames/CleanPathogenNames.csv')
   bugs$type <- apply(bugs, 1, function(x) paste(which(x != ''), collapse=', '))
   
   #### Group X Streptococci
   w <- which(apply(bugs, 1, function(x) any(grepl('\\[', x))))
   for (i in w) {
      w_i <- intersect(intersect(grep(pattern = bugs$Pattern1[i],
                                      x = df$PATH_NAME, 
                                      ignore.case = TRUE),
                                 grep(pattern = '.*group ([ACDFG])( .*|$)',
                                      x = df$PATH_NAME, 
                                      ignore.case = TRUE)),
                       grep(pattern = 'not group|agalactiae',
                            x = df$PATH_NAME, 
                            ignore.case = TRUE,
                            invert = TRUE))
      
      df$BUG[w_i] <- paste0('Group ', 
                            toupper(gsub(pattern = '.*Group ([ACDFG])( .*|$)',  '\\1',
                                         x = df$PATH_NAME[w_i],
                                         ignore.case = TRUE)), 
                            ' Streptococci')
   }
   df$BUG[which(df$BUG == 'Group A Streptococci')] <- 'Streptococcus pyogenes'
   
   # for COAGULASE NEGATIVE STAPH - weird one
   w <- which(bugs$type == '1, 2, 3, 4, 5')
   for (i in w) {
      w_i <- intersect(intersect(which(stringi::stri_detect_regex(pattern = bugs$Pattern1[i],
                                                                  str = df$PATH_NAME,
                                                                  case_insensitive = TRUE)),
                                 which(stringi::stri_detect_regex(pattern = bugs$Pattern2[i],
                                                                  str = df$PATH_NAME, 
                                                                  case_insensitive = TRUE))),
                       intersect(which(stringi::stri_detect_regex(pattern = bugs$Pattern3[i],
                                                                  str = df$PATH_NAME, 
                                                                  case_insensitive = TRUE)),
                                 which(stringi::stri_detect_regex(pattern = bugs$DontMatch[i],
                                                                  str = df$PATH_NAME, 
                                                                  case_insensitive = TRUE,
                                                                  negate = TRUE))))
      
      df$BUG[w_i] <- bugs$Spelling[i]
   }
   
   #### searching for simple cases of genus_species
   w <- which(bugs$type == '1, 2' & !apply(bugs, 1, function(x) any(grepl('\\(|\\[', x))))
   pats <- paste0(bugs[w,1], ' ', bugs[w,2])
   for (i in seq_along(pats)) {
      print(i)
      w_i <- which(stringi::stri_detect_regex(pattern = pats[i], 
                                              str = df$PATH_NAME,
                                              case_insensitive = TRUE))
      
      df$BUG[w_i] <- pats[i]
   }
   print('Handled all simple cases.')
   
   #### for searching pieces separately
   w <- which(bugs$type == '1, 2, 5')
   for (i in w) {
      print(i)
      w_i <- intersect(which(stringi::stri_detect_regex(pattern = bugs$Pattern1[i],
                                                        str = df$PATH_NAME,
                                                        case_insensitive = TRUE)),
                       which(stringi::stri_detect_regex(pattern = bugs$Pattern2[i],
                                                        str = df$PATH_NAME,
                                                        case_insensitive = TRUE)))
      
      df$BUG[w_i] <- bugs$Spelling[i]
   }
   print('Handled all simple cases #2.')
   
   #### for searching single
   w <- which(bugs$type == '1, 5')
   for (i in w) {
      print(i)
      w_i <- which(stringi::stri_detect_regex(pattern = bugs$Pattern1[i],
                                              str = df$PATH_NAME, 
                                              case_insensitive = TRUE))
      
      df$BUG[w_i] <- bugs$Spelling[i]
   }
   print('Handled all simple cases #3.')
   
   # LEFTOVERS - fix name and do last
   w <- which(bugs$type == '1, 5, 6')
   for (i in w) {
      print(i)
      w_i <- which(stringi::stri_detect_regex(pattern = bugs$Pattern1[i],
                                              str = df$PATH_NAME,
                                              case_insensitive = TRUE))
      
      df$BUG[w_i] <- bugs$Spelling[i]
   }
   
   # LEFTOVERS - keep and do last
   w <- which(bugs$type == '1, 6')
   for (i in w) {
      print(i)
      w_i <- which(stringi::stri_detect_regex(pattern = bugs$Pattern1[i],
                                              str = df$PATH_NAME,
                                              case_insensitive = TRUE))
      
      df$BUG[w_i] <- bugs$Pattern1[i]
   }
   
   #### No species identified
   w <- grep('^[A-z]+ species$', bugs$Pattern1)
   genus <- gsub(' species$', '', bugs$Pattern1[w])
   w_df <- which(is.na(df$BUG))
   for (i in seq_along(w)) {
      print(i)
      w_i <- w_df[which(stringi::stri_detect_regex(pattern = genus[i],
                                                   str = df$PATH_NAME[w_df],
                                                   case_insensitive = TRUE))]
      
      df$BUG[w_i] <- paste0(genus[i], ' species')
   }
   print('Handled all sans species cases.')
   
   df <- df %>% mutate(BUG = gsub(' ?', '', BUG, fixed=TRUE))
   df <- df %>% mutate(BUG = gsub('?', '', BUG, fixed=TRUE))
   return(df)
}
