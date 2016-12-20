loadGeneLists = function( location, type = "GMT", sep = "\t", header = FALSE, col = 1 ) {
    if( toupper(type) == "GMT" ) {
        ## Load in the gene lists from the GMT file
        raw_data = strsplit( scan( location, what = "character", sep = "\n", quiet = TRUE ), sep );
        gene_lists = lapply( raw_data, function(x) { keep = x!=""; keep[1:2] = FALSE; return( unique(x[keep]) ); } );
        names( gene_lists ) = unlist( lapply( raw_data, function(x) return(x[1]) ) );
        return( gene_lists );
    } else if( toupper(type) == "DIR" ) {
        ## Identify the gene list files within the directory
        gene_list_files = list.files( location );
        gene_list_files_full = paste( sep="", location, "/", gene_list_files );

        ## Load in the gene lists from the files in the directory
        gene_lists = list();
        for( i in 1:length(gene_list_files) ) {
            raw_data = read.table( gene_list_files_full[i], header = FALSE, sep = sep );
            if( header == TRUE ) {
                gene_lists[[i]] = as.character( raw_data[-1,col] );
            } else {
                gene_lists[[i]] = as.character( raw_data[,col] );
            }
        }
        names( gene_lists ) = unlist( lapply( gene_list_files, function(x) { if( substr(x,nchar(x)-3,nchar(x)) == ".txt" ) { return( substr(x,1,nchar(x)-4) ); } else { return(x); } } ) );
        return( gene_lists );
    } else {
        stop( "Invalid type" );
    }
}
