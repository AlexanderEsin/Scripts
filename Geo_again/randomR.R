x_tbl		<- dbSendQuery(conn, 'SELECT DISTINCT genome_l FROM t1 WHERE taxid = :taxid_list LIMIT 1')
dbBind(x_tbl, param = list(taxid_list = taxid_l))
x_df		<- dbFetch(x_tbl)
dbClearResult(x_tbl)