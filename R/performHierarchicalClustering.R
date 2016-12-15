performHierarchicalClustering = function( enrichment_matrix, feature = "row", distance_method = "correlation", clustering_method = "average" ) {
  if( toupper( feature ) == "ROW" ) {
    final_matrix = enrichment_matrix;
  } else if( toupper( feature ) == "COL" ) {
    final_matrix = t(enrichment_matrix);
  } else {
    stop( "Invalid feature" );
  }

  if( toupper( distance_method ) == "CORRELATION" ) {
    return( hclust( as.dist( 1-abs(cor(t(final_matrix))) ), method = clustering_method ) );
  } else {
    return( hclust( dist( final_matrix, method = distance_method ), method = clustering_method ) );
  }
}
