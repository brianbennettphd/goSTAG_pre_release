performGOEnrichment = function( gene_lists, go_terms, filter_method = "pval", significance_threshold = 0.05, p.adjust_method = "BH" ) {
    gene_list_sizes = unlist( lapply( gene_lists, length ) );
    universe_size = length( unique(unlist(go_terms)) );

    pvals = matrix( nrow = length(go_terms), ncol = length(gene_lists) );
    for( i in 1:length(go_terms) ) {
        m = length( go_terms[[i]] );
        for( j in 1:length(gene_lists) ) {
            x = length( intersect( toupper(go_terms[[i]]), toupper(gene_lists[[j]]) ) );
            k = gene_list_sizes[j];
            N = universe_size;
            pvals[i,j] = phyper( x-1, m, N-m, k, lower.tail=FALSE ); # hypergeometric test
        }
    }
    rownames(pvals) = names(go_terms);
    colnames(pvals) = names(gene_lists);

    adjusted_pvals = matrix( nrow = length(go_terms), ncol = length(gene_lists) );
    for( i in 1:length(gene_lists) ) {
        adjusted_pvals[,i] = p.adjust( pvals[,i], method = p.adjust_method );
    }

    if( toupper( filter_method ) == "PVAL" ) {
        final_matrix = pvals[ apply( pvals<significance_threshold, 1, sum ) > 0, ];
        final_matrix[ final_matrix>=significance_threshold ] = 1;
        final_matrix = -log10( final_matrix );
    } else if( toupper( filter_method ) == "P.ADJUST" ) {
        final_matrix = pvals;
        final_matrix[ adjusted_pvals>=significance_threshold ] = 1;
        final_matrix = final_matrix[ apply( adjusted_pvals<significance_threshold, 1, sum ) > 0, ];
        final_matrix = -log10( final_matrix );
    } else {
        stop( "Invalid filter method" );
    }

    if( dim(final_matrix)[1] == 0 ) {
        stop( "Zero significant GO terms (try less stringent signficance threshold)" );
    } else if( dim(final_matrix)[1] == 1 ) {
        stop( "Only one significant GO term (try less stringent signficance threshold)" );
    }

    final_matrix[is.infinite(final_matrix)] = max(final_matrix[!is.infinite(final_matrix)]);

    return( final_matrix );
}
