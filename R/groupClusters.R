groupClusters = function( hclust_results, distance_threshold = 0.2 ) {
  get_root = function( hclust_results, step ) {
    if( hclust_results$merge[ step, 1 ] < 0 ) {
      return( hclust_results$merge[ step, 1 ] );
    } else {
      return( get_root( hclust_results, hclust_results$merge[ step, 1 ] ) );
    }
  }

  hclust_clusters = 1:length(hclust_results$order);
  if( sum(hclust_results$height<=distance_threshold) >= 1 ) {
    for( i in 1:sum(hclust_results$height<=distance_threshold) ) {
      if( hclust_results$merge[i,1] < 0 && hclust_results$merge[i,2] < 0 ) {
        hclust_clusters[ -hclust_results$merge[i,1] ] = hclust_clusters[ -hclust_results$merge[i,2] ];
      } else if( hclust_results$merge[i,1] < 0 && hclust_results$merge[i,2] > 0 ) {
        hclust_clusters[ -hclust_results$merge[i,1] ] = hclust_clusters[ -get_root( hclust_results, hclust_results$merge[i,2] ) ];
      } else if( hclust_results$merge[i,1] > 0 && hclust_results$merge[i,2] > 0 ) {
        hclust_clusters[ hclust_clusters == hclust_clusters[ -get_root( hclust_results, hclust_results$merge[i,1] ) ] ] = hclust_clusters[ -get_root( hclust_results, hclust_results$merge[i,2] ) ];
      }
    }
  }

  hclust_final_clusters = rep( 1, length(hclust_results$order) );
  for( i in 2:length(hclust_results$order) ) {
    if( hclust_clusters[ hclust_results$order[i] ] == hclust_clusters[ hclust_results$order[i-1] ] ) {
      hclust_final_clusters[i] = hclust_final_clusters[i-1];
    } else {
      hclust_final_clusters[i] = hclust_final_clusters[i-1] + 1;
    }
  }

  clusters = list();
  for( i in 1:max(hclust_final_clusters) ) {
    clusters[[i]] = hclust_results$labels[ hclust_results$order[ hclust_final_clusters==i ] ];
  }
  names(clusters) = paste( sep="", "Cluster", seq( 1, max(hclust_final_clusters) ) );

  return( clusters );
}
