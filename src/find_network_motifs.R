find_network_motifs <- function(networks, motifs, plot_individual_motifs = T,
                                test_significance = T) 
{
  #'
  #' ___________________________________________________________________________
  
  
  
  # Plotting the individual motifs as reference
  
  if (plot_individual_motifs)
  {
    plot_dir <- paste(getwd(), "Results", "Plots", "Individual_motifs", 
                      "Networks", "Motifs", sep = "/")
    
    dir.create(plot_dir)
    
    apply(names(motifs), function(x)
      {
      pdf(file = paste0(plot_dir, "/", x, "_Motif.pdf"), height = 4, width = 4)
      
      plot(motifs[[x]])
      
      dev.off()
    })
  }
  
  return(results)
}