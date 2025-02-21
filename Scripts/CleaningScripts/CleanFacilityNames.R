getFacilityNames <- function(facility) {
   site_names <- c("CHP" = 'Childrens', "UPMCALT" = 'Altoona',  "UPMCBED" = 'Bedford', 
                   "UPMCCHA" = 'Chatauqua', "UPMCEAS" = 'East', "UPMCHAM" = 'Hamot', 
                   "UPMCHZN" = 'Horizon', "UPMCJAM" = 'Jameson', "UPMCMCK" = 'McKeesport', 
                   "UPMCMER" = 'Mercy', "UPMCMUN" = '', "UPMCMWH" = 'Magee-W', 
                   "UPMCNOR" = 'Northwest', "UPMCPAS" = 'Passavant', "UPMCPUH" = 'Presbyterian', 
                   "UPMCSHY" = 'Shadyside', "UPMCSMH" = 'St. Margaret', 
                   "UPMCSOL" = 'SOL', "UPMCSUN" = 'SUN',  'UPMCLOC' = 'LOC', 'UPMCMUN' = 'MUN',
                   "UPMCWIL" = 'Williamsport')
   # site_groups <- list(
   #    academic = c('Mercy', 'Presbyterian', 'Shadyside', 'Magee-W'),
   #    regional = c('Hamot', 'Williamsport', 'Jameson', 'Altoona'),
   #    community = c('Passavant', 'East', 'McKeesport', 'St. Margaret'),
   #    rural = c('Bedford', 'Northwest', 'Horizon', 'Chatauqua')
   # )
   sites <- unname(site_names[facility])
   x <- tibble(facility,
               site = sites) %>%
      filter(is.na(site), !is.na(facility))
   if (nrow(x) > 0L)
      print(x)
   
   return(sites)
}
