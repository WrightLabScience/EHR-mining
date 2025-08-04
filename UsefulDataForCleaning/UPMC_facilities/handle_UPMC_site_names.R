site_info <- list(
   c('CHP', 'Childrens', 'Specialty'),
   c('UPMCALT', 'Altoona', 'Regional'),
   c('UPMCBED', 'Bedford', 'Rural'),
   c('UPMCCHA', 'Chatauqua', 'Rural'),
   c('UPMCEAS', 'East', 'Community'),
   c('UPMCHAM', 'Hamot', 'Regional'),
   c('UPMCHZN', 'Horizon', 'Rural'),
   c('UPMCJAM', 'Jameson', 'Regional'),
   c('UPMCMCK', 'McKeesport', 'Community'),
   c('UPMCMER', 'Mercy', 'Academic'),
   c('UPMCMUN', 'Muncy', 'Regional'),
   c('UPMCMWH', 'MageeW', 'Specialty'),
   c('UPMCNOR', 'Northwest', 'Rural'),
   c('UPMCPAS', 'Passavant', 'Community'),
   c('UPMCPUH', 'Presbyterian', 'Academic'),
   c('UPMCSHY', 'Shadyside', 'Academic'),
   c('UPMCSMH', 'StMargaret', 'Community'),
   c('UPMCSOL', 'SOL', 'Community'),
   c('UPMCSUN', 'Sunbury', 'Rural'),
   c('UPMCLOC', 'LOC', 'Community'),
   c('UPMCWIL', 'Williamsport', 'Regional')
)

site_info <- setNames(data.frame(do.call(rbind, site_info)), 
                      c('code', 'facility', 'category'))

code2fac <- setNames(site_info$facility,
                     site_info$code)

fac2cat <- setNames(site_info$category,
                    site_info$facility)

code2cat <- setNames(site_info$category,
                     site_info$code)

code2mixed <- c(code2cat[code2cat != 'Academic'],
                setNames(site_info$facility[site_info$category == 'Academic'],
                         site_info$code[site_info$category == 'Academic']))

fac2mixed <- c(fac2cat[fac2cat != 'Academic'],
               setNames(site_info$facility[site_info$category == 'Academic'],
                        site_info$facility[site_info$category == 'Academic']))

mixed2cat <- c(fac2cat[fac2cat == 'Academic'],
               setNames(unique(site_info$category[site_info$category != 'Academic']),
                        unique(site_info$category[site_info$category != 'Academic'])))

code2fac_fxn <- function(code) {
   return(unname(code2fac[code]))
}

fac2cat_fxn <- function(fac) {
   return(unname(fac2cat[fac]))
}

code2cat_fxn <- function(code) {
   return(unname(code2cat[code]))
}

code2mixed_fxn <- function(code) {
   return(unname(code2mixed[code]))
}

fac2mixed_fxn <- function(fac) {
   return(unname(fac2mixed[fac]))
}

mixed2cat_fxn <- function(mixed) {
   return(unname(mixed2cat[mixed]))
}





