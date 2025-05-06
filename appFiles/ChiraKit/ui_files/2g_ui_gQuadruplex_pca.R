box(title = "2. PCA Clustering", width = 12, solidHeader = T, status = "primary",
    fluidRow(

        # Select input to select the sample of interest
        column(3, p(HTML("<b>Sample</b>"),
            span(shiny::icon("info-circle"), id = "info_uu-sampleSelectGquad_pca"),
            selectInput("sampleSelectGquad_pca", NULL, choices = c("None"),
            selected = "Sample 1", multiple = FALSE,
            selectize = FALSE, width = "100%"),
            tippy::tippy_this(
                elementId = "info_uu-sampleSelectGquad_pca",
                tooltip = "If no spectrum is available, import data in the 'Import data' section.
                The selected spectrum will be used together with the reference spectra to
                perform PCA clustering.",
                placement = "right")))
    ),

    fluidRow(

        # Button to proceed to sample PCA (green button with arrow)
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            actionButton(
            inputId = "go_to_references_ui",label = "Back to references",
            icon("arrow-left"),
            style="color: #fff; background-color: #337ab7;
            border-color: #2e6da4"))),

        # Button to proceed to sample PCA (green button with arrow)
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            actionButton(
            inputId = "go_to_sample_svd_ui",label = "2b. Proceed to sec. str. estimation",
            icon("arrow-right"),
            style="color: #fff; background-color: #337ab7;
            border-color: #2e6da4")))

    )

)