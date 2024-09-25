getConcordanceFlag <- function(df, abx) {
   if (length(abx) == 0) return('No empiric therapy given')
   
   check <- abx %in% col_names  # indicates which abx given have corresponding rows in empDF 
   vals <- unlist(df[, abx[check]])     # values associated with abx given
   
   # if any AST value is 0 (susceptible), the therapy is concordant
   if (any(vals == 0, na.rm = T)) return("CONCORDANT")
   
   # all are in table
   if (all(check)) {
      if (all(is.na(vals)))                   return('NOT TESTED')
      if (all(!is.na(vals)) & all(vals == 1)) return('DISCORDANT')  # all tested and all resistant
      if (any(vals == 1) & any(is.na(vals)))  return('DISCORDANT')      # >=1 not tested, >=1 resistant
      print('if Im here, something went wrong')
   }
   
   # all are not in table
   if (all(!check)) return('NOT TESTED')
   
   # NT + never
   r <- any(vals == 1, na.rm=TRUE)
   nt <- any(is.na(vals))
   nev <- any(!check)
   
   if (r & nt & !nev) return('DISCORDANT')
   if (!r & nt & nev) return('NOT TESTED')
   if (r & !nt & nev) return('DISCORDANT')
   if (r & nt & nev)  return('DISCORDANT')
}
