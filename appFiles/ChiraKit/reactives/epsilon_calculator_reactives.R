# Protein molar extinction coefficient calculator
observeEvent(input$epsilon_calculator_click, {

  showModal(modalDialog(

    tags$h3('Please enter the protein sequence using 1 digit code
    (e.g., ALAS):'),
    textInput("proteinSequence",label = "",value = ""),

    footer=tagList(
      actionButton('submitCalcEpsilon', 'Submit'),
      modalButton('Cancel')
    )
  ))

})

observeEvent(input$submitCalcEpsilon,{

  removeModal()
  sequence <- input$proteinSequence
  print(sequence)

})