# Protein molar extinction coefficient calculator
observeEvent(input$epsilon_calculator_click, {

    showModal(modalDialog(

        tags$h3('Please enter the protein sequence using 1 digit code
            (e.g., MALSNVIKYIQQQQKSKRLLLAAQK):'),
            textInput("proteinSequence",label = "",value = ""
        ),

        footer=tagList(
            actionButton('submitCalcEpsilon', 'Submit'),
            modalButton('Cancel')
        )

    ))

})

observeEvent(input$submitCalcEpsilon,{

    removeModal()
    sequence <- input$proteinSequence

    # Convert to uppercase
    sequence <- toupper(sequence)

    # Remove empty spaces
    sequence <- gsub(" ", "", sequence)

    # Check that we have a valid sequence of aminoacids
    valid_aminoacids <- c(
        'A', 'R', 'N', 'D', 'C', 'Q',
        'E', 'G', 'H', 'I', 'L', 'K',
        'M', 'F', 'P', 'S', 'T', 'W',
        'Y', 'V'
    )

    if (!all(unlist(strsplit(sequence, "")) %in% valid_aminoacids)) {

        popUpWarning("Invalid aminoacid sequence.
        Please enter a valid sequence.")
        return(NULL)

    }

    epsilon_dfs_lst <- calculate_epsilon(sequence)

    epsilon_table   <- epsilon_dfs_lst[[1]]

    # Round all columns 0 decimal digits
    for (col in 1:ncol(epsilon_table)) {
        epsilon_table[,col] <- as.integer(epsilon_table[,col])
    }

    # Nice names for the columns
    colnames(epsilon_table)[1:3] <- c(
        "&#949;<sub>280nm</sub> [M<sup>-1</sup>cm<sup>-1</sup>]",
        "&#949;<sub>215nm</sub> [M<sup>-1</sup>cm<sup>-1</sup>]",
        "&#949;<sub>204nm</sub> [M<sup>-1</sup>cm<sup>-1</sup>]"
        )

    output$epsilon_table <- renderTable({epsilon_table}, sanitize.text.function = function(x) x)
    # Show a Table with the results

    absorbance_table   <- epsilon_dfs_lst[[2]]

    # Nice names for the columns
    colnames(absorbance_table)[1:3] <- c(
        "Abs<sub>280nm</sub> 0.1% (=1 g/l)",
        "Abs<sub>214nm</sub> 0.1% (=1 g/l)",
        "Abs<sub>205nm</sub> 0.1% (=1 g/l)"
        )

    output$absorbance_table <- renderTable({absorbance_table}, sanitize.text.function = function(x) x)

    molecular_weight       <- aa_sequence_to_weight(sequence)
    molecular_weight_table <- data.frame(molecular_weight)
    colnames(molecular_weight_table) <- c("Molecular weight (Da)")

    output$molecular_weight_table <- renderTable({molecular_weight_table}, sanitize.text.function = function(x) x)

    number_of_residues       <- nchar(sequence)
    number_of_residues_table <- data.frame(number_of_residues)
    colnames(number_of_residues_table) <- c("Number of residues")

    output$number_of_residues_table <- renderTable({number_of_residues_table}, sanitize.text.function = function(x) x)

    showModal(modalDialog(
        tags$h3('Results:'),
        tableOutput("epsilon_table"),
        tags$h3(''),
        tableOutput("absorbance_table"),
        tags$h3(''),
        tableOutput("molecular_weight_table"),
        tags$h3(''),
        tableOutput("number_of_residues_table"),
        footer=tagList(
            modalButton('Close')
        ),
        size = "l"
    ))

})