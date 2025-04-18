box(title = "2. Parameters for molar ellipticity / extinction", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(8, p(HTML("<b>Shared experiment parameters</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-2"),
                  checkboxInput("sharedExperimentParameters", NULL, value = FALSE),
                  tippy::tippy_this(
                    elementId = "info_uu1-2",
                    tooltip = "Select this option if all the CD curves share the same parameters:
                    Molecular weight, Concentration, Number of chromophore units, Path length and Input units.
                    The Molecular weight, Concentration and Path length are required to 
                    calculate the molar ellipticity/extinction. 
                    The 'number of chromophore units' is useful for polypeptides, 
                    because it allows normalisation by the number of peptide bonds.
                     ",placement = "right"))),
                       
     column(12,DT::dataTableOutput('cdFilesInfo')),
                       tags$script(HTML("Shiny.addCustomMessageHandler('unbind-DT', function(id) {
          Shiny.unbindAll($('#'+id).find('table').DataTable().table().node());
        })"))
                       
    ))