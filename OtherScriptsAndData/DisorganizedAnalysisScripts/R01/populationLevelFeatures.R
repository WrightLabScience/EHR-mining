# get UPMC facility of all AST
# get UPMC facility of all ABX admin
# with these, we can calculate for each infection
#     total number Rx1/C at the facility they are at and the facilities they have been
#     total number Res1/C at the facility they are at and the facilities they have been

source(file = '~/Desktop/EHR/EHR work/config_file.R')

load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_outcomes_variables.Rdata')

encs <- tibble(dbGetQuery(conn, "SELECT * FROM AMB_ETL.SENS_ENCOUNTER_VW WHERE PERSON_ID = 1000000122"))

