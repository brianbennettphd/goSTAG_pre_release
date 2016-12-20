annotateClusters = function( clusters ) {
    ## Load the GO parents from the GO.db package
    if( Ontology( GOTERM[ clusters[[1]][1] ] ) == "BP" ) {
        parents = as.list(GOBPPARENTS);
    } else if( Ontology( GOTERM[ clusters[[1]][1] ] ) == "CC" ) {
        parents = as.list(GOCCPARENTS);
    } else if( Ontology( GOTERM[ clusters[[1]][1] ] ) == "MF" ) {
        parents = as.list(GOMFPARENTS);
    } else {
        stop( "Invalid domain" );
    }

    ## Cycle through all the clusters
    cluster_labels = c();
    for( i in 1:length(clusters) ) {
        ## Define the recursive function to get the number of paths to the root
        get_num_paths = function( go_id, parents ) {
            this_parents = unlist(as.list(parents[go_id]));
            this_parents = this_parents[ unlist( lapply( strsplit(names(this_parents),"[.]"), function(x) return( x[2] == "is_a" ) ) ) ];
            if( this_parents[1] == "all" ) {
                return(1);
            } else {
                return( sum( unlist( lapply( this_parents, function(x) return( get_num_paths( x, parents ) ) ) ) ) );
            }
        }

        ## Get the number of paths to the root for each GO term in the cluster
        num_paths = unlist( lapply( clusters[[i]], function(x) return(get_num_paths(x,parents)) ) );

        ## Select the GO term with the largest number of paths to the root
        cluster_labels[i] = Term( GOTERM[ clusters[[i]][ which.max(num_paths) ] ] );
    }

    return( cluster_labels );
}
