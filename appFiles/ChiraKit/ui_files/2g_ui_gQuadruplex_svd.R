box(title = "3. Secondary structure estimation", width = 12, solidHeader = T, status = "primary",
    fluidRow(

        #column(4, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
        #    actionButton(
        #    inputId = "launchSamplesPCAGquad",label = "2a. Run PCA on samples",
        #    icon("chart-line"),
        #    style="color: #fff; background-color: #337ab7;
        #    border-color: #2e6da4")))#,

        # Select input to select the sample of interest, add tooltip

        column(3, p(HTML("<b>Sample</b>"),
            span(shiny::icon("info-circle"), id = "info_uu-sampleSelectGquad_svd"),
            selectInput("sampleSelectGquad_svd", NULL, choices = c("None"),
            selected = "Sample 1", multiple = FALSE,
            selectize = FALSE, width = "100%"),
            tippy::tippy_this(
                elementId = "info_uu-sampleSelectGquad_svd",
                tooltip = "Select custom if you want to process the spectra manually.
                Select 'Automatic baseline subtraction' if you want to load the CD scans of the sample,
                the CD scans of the baseline  and process them automatically.
                This includes 1)
                averaging the scans, 2) subtracting the baseline, 3) zeroing the spectrum, and 4)
                normalising the final spectrum (e.g., by converting it to molar ellipticity), if desired.
                ",placement = "right")))#,

        #column(4, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
        #    actionButton(
        #    inputId = "launchSVD_GQ",label = "2b. Run SVD structure estimation",
        #    icon("tornado"),
        #    style="color: #fff; background-color: #337ab7;
        #    border-color: #2e6da4")))
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