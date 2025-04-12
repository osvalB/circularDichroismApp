box(title = "Samples dataset - Plot settings", width = 12, solidHeader = T, status = "primary", 
    
    fluidRow(
      
      column(3, p(HTML("<b>PC - Horizontal axis</b>"),
                  numericInput("PC_ha_samples", NULL,1,min = 1,max = 8))),
      
      column(3, p(HTML("<b>PC - Vertical axis</b>"),
                  numericInput("PC_va_samples", NULL,2,min = 2,max = 8))),            
      
      column(3, p(HTML("<b>Number of clusters</b>"),
                  selectInput("GQnClust_samples", NULL,
                              paste0(3:10)))),
      #c('Auto',paste0(3:10))))),
      
      column(3, p(HTML("<b>Show graph options</b>"),
                  selectInput("showAdvanced_GQ_samples", NULL,
                              choices = c('None','PCA','Cluster'))))),
    
    conditionalPanel(
      "input.showAdvanced_GQ_samples == 'PCA'",
      
      fluidRow(
        
        column(3, p(HTML("<b>Graph title</b>"),
                    textInput("pca_graph_title_samples", NULL,""))),
        
        column(3, p(HTML("<b>Graph label size</b>"),
                    numericInput("pca_graph_label_samples", NULL,6,min = 0,max = 15))),
        
        column(3, p(HTML("<b>Axis label size</b>"),
                    numericInput("pca_axis_label_samples", NULL,16,min = 0,max = 30))), 
        
        column(3, p(HTML("<b>Point size</b>"),
                    numericInput("pca_point_size_samples", NULL,3,min = 0,max = 10))),
        
        column(3, p(HTML("<b>Arrow width</b>"),
                    numericInput("pca_arrow_width_samples", NULL,1,min = 0,max = 10))),
        
        column(3, p(HTML("<b>Ellipse confidence level</b>"),
                    numericInput("pca_ellipse_confidence_samples", NULL,85,min = 0,max = 100))),
        
        column(3, p(HTML("<b>Show ellipse</b>"),
                    checkboxInput("show_ellipse_pcaGQ_samples", NULL,TRUE))),
        
        column(3, p(HTML("<b>Show legend</b>"),
                    checkboxInput("show_legend_pcaGQ_samples", NULL,FALSE)))),
      
      fluidRow(
        
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                    downloadButton("download_samples_png", "PCA samples plot (.png)" ))),
        
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                    downloadButton("download_samples_pdf", "PCA samples plot (.pdf)" ))),
        
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                    downloadButton("download_combined_png", "PCA combined plot (.png)" ))),
        
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                    downloadButton("download_combined_pdf", "PCA combined plot (.pdf)" )))
      )
      
    ),
    
    conditionalPanel(
      "input.showAdvanced_GQ_samples == 'Cluster'",
      
      fluidRow(
        
        column(3, p(HTML("<b>Graph title</b>"),
                    textInput("cluster_graph_title_samples", NULL,""))),
        
        column(3, p(HTML("<b>Text size</b>"),
                    numericInput("cluster_text_size_samples", NULL,9,min = 0,max = 20)))),
        
      fluidRow(
      
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                    downloadButton("download_cluster_png_samples", "Samples clustering (.png)"))),
        
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                    downloadButton("download_cluster_pdf_samples", "Samples clustering (.pdf)"))),
        
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                    downloadButton("download_cluster_png_combined", "Combined clustering (.png)"))),
        
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                    downloadButton("download_cluster_pdf_combined", "Combined clustering (.pdf)")))
      )
      
    )
)
