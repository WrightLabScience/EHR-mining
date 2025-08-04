## Example config file setup 
options(java.parameters = "-Xmx4g")
library(RJDBC)
library(dplyr)
library(dbplyr)
library(tidyr)

# Enter the location of your oracle JDBC driver (ojdbc7.jar or ojdbc8.jar)
drv <- JDBC("oracle.jdbc.OracleDriver", '/path_to_oracle_driver_jar/ojdbc8.jar')

# Enter your connection string and username here. You will be prompted for your password upon running the scripts.
connection_string <- '' # something like jdbc:oracle:thin:@dbmi......
username <- '' # mine was SEB203
conn <- dbConnect(drv, connection_string, username, password = getPass::getPass())

# execute simple query to get list of view names from schema
sql <- "SELECT view_name FROM all_views WHERE owner = 'AMB_ETL'" # enter schema name here
query <- DBI::dbSendQuery(conn, sql)
result <- DBI::dbFetch(query)
DBI::dbClearResult(query)