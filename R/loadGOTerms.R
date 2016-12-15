loadGOTerms = function( species = "human", domain = "BP", min_num_genes = 5, use_archived = TRUE ) {
    if( ! toupper(domain) %in% c("BP","CC","MF") ) {
        stop( "Invalid domain" );
    }

    domains = Ontology(GOTERM);

    if( use_archived == TRUE ) {
        if( toupper(species) == "HUMAN" ) {
            go_genes_all = goSTAG::goSTAG_go_genes_human;
        } else if( toupper(species) == "MOUSE" ) {
            go_genes_all = goSTAG::goSTAG_go_genes_mouse;
        } else if( toupper(species) == "RAT" ) {
            go_genes_all = goSTAG::goSTAG_go_genes_rat;
        } else {
            stop( "Invalid species" );
        }

        go_genes = go_genes_all[ domains[names(go_genes_all)]==toupper(domain) ];
    } else {
        if( toupper(species) == "HUMAN" ) {
            ensembl_dataset = "hsapiens_gene_ensembl";
            gene_symbol_attribute = "hgnc_symbol";
        } else if( toupper(species) == "MOUSE" ) {
            ensembl_dataset = "mmusculus_gene_ensembl";
            gene_symbol_attribute = "mgi_symbol";
        } else if( toupper(species) == "RAT" ) {
            ensembl_dataset = "rnorvegicus_gene_ensembl";
            gene_symbol_attribute = "rgd_symbol";
        } else {
            stop( "Invalid species" );
        }

        if( ! toupper(domain) %in% c("BP","CC","MF") ) {
            stop( "Invalid domain" );
        }

        ensembl_mart = useMart( biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset=ensembl_dataset );
        query = getBM( attributes = c(gene_symbol_attribute,"go_id"), mart = ensembl_mart );

        go_ids = unique(query[,2]);
        go_ids = go_ids[go_ids!=""];
        go_ids_filtered = go_ids[ domains[go_ids]==toupper(domain) ];
        query = query[ query[,1]!="", ];
        query_filtered = query[ query[,2] %in% go_ids_filtered, ];

        go_genes = lapply( go_ids_filtered, function(x) return( unique(query_filtered[query_filtered[,2]==x,1]) ) );
        names(go_genes) = go_ids_filtered;
    }

    go_terms = go_genes[ lapply(go_genes,length) >= min_num_genes ];

    return( go_terms );
}
