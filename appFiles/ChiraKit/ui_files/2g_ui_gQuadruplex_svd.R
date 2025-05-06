box(title = "3. Structure estimation", width = 12, solidHeader = T, status = "primary",
    fluidRow(

        column(3, p(HTML("<b>Sample</b>"),
            span(shiny::icon("info-circle"), id = "info_uu-sampleSelectGquad_svd"),
            selectInput("sampleSelectGquad_svd", NULL, choices = c("None"),
            selected = "Sample 1", multiple = FALSE,
            selectize = FALSE, width = "100%"),
            tippy::tippy_this(
                elementId = "info_uu-sampleSelectGquad_svd",
                tooltip = "If no spectrum is available, import data in the 'Import data' section.
                The reference spectra will be used to perform the SVD-based structure estimation on the selected 'Sample' spectrum.",
                placement = "right")))
    ),

    fluidRow(

        # Button to proceed to sample PCA (green button with arrow)
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            actionButton(
            inputId = "go_to_samples_pca",label = "Back to PCA clustering",
            icon("arrow-left"),
            style="color: #fff; background-color: #337ab7;
            border-color: #2e6da4")))

    )

)