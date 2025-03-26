# Generate and render
# 1) The Input parameters Table,
# 2) The legends Table
# 3) The wavelength range slider
renderInputData <- function() {
  
    output$cdFilesInfo    <- NULL

    wls  <- cdAnalyzer$get_experiment_properties('wavelength')

    if (length(wls) == 0) return(NULL)

    Sys.sleep(length(wls)*0.015)

    minV <- min(unlist(wls))
    maxV <- max(unlist(wls))

    updateSliderInput(session,inputId="wavelengthRange",label = NULL,min=minV,max=maxV,value = c(minV,maxV))
    updateSelectInput(session,"experiment2delete",choices     = c("ALL",cdAnalyzer$experimentNames))
    updateSelectInput(session,"selected_cd_exp",choices       = cdAnalyzer$experimentNames)

    updateCDFilesInfoTable()

    cdAnalyzer$sharedParameters <- FALSE
    updateCheckboxInput(session,'sharedExperimentParameters',label = NULL,FALSE)

    legendDf              <- getPlottingDF(cdAnalyzer)

    output$legendInfo     <- helperRenderRHandsontable(legendDf)
    updateSESCA_ref(legendDf)

    Sys.sleep(length(wls)*0.015)
  
}

observeEvent(input$cdFiles,{
  
    req(input$cdFiles)

    reactives$data_loaded <- NULL
    output$legendInfo     <- NULL

    cd_data_files   <- input$cdFiles$datapath
    names           <- input$cdFiles[[1]]

    sorted_indices  <- order(names)

    cd_data_files   <- cd_data_files[sorted_indices]
    names           <- names[sorted_indices]

    i <- 0

    # iterate over the files
    for (name in names) {

        i              <- i + 1
        cd_data_file   <- cd_data_files[i]

        fileExtension <- getFileNameExtension(name)

        if (fileExtension == 'zip') {

            # Create a temporary directory to unzip files
            unzipDir <- tempfile()
            dir.create(unzipDir)

            unzip(cd_data_file,exdir = unzipDir)

            # List all possible types of files
            file_list <- list.files(unzipDir,
                              pattern = "\\.xlsx$|\\.gen$|\\.dat$|\\.txt$|\\.csv$|d\\d+$",
                              ignore.case = TRUE,
                              recursive = TRUE,
                              full.names = TRUE)

            for (file in file_list) load_one_experiment(file,basename(file))

        } else {

            load_one_experiment(cd_data_file,name)
        }
    }

    Sys.sleep(0.5)

    # Check that we have at least one experiment
    if (length(cdAnalyzer$experimentNames) > 0) {

        reactives$urea_example_data_loaded <- FALSE

        # Copy experiments to allow modifying the wavelength range
        cdAnalyzer$initialize_experiment_modif()
        updateMaxVoltageValue()
        renderInputData()

        reactives$show_example_data <- FALSE
        reactives$data_loaded       <- TRUE
    }

},priority = 10)

# Modal dialog to ask for the input and desired units
observeEvent(input$automaticProcess,{
  
    req(input$cdFilesSample)
    req(input$cdFilesBaseline)

    reactives$urea_example_data_loaded <- FALSE

    showModal(modalDialog(

        tags$h3('Please enter the units of the CD data files,
            the desired units of the final processed spectrum,
            and if the final spectrum should be zeroed:'),
        selectInput("inputUnitsAutomatic"  , 'CD input files unit',
            global_cd_units_choice_short),
        selectInput("workingUnitsAutomatic", 'Desired units of the final spectrum',
            global_cd_units_choice),
        checkboxInput('zeroFinalSpectrum',
            'Subtract the average signal of the highest
            10 nm wavelength interval from the final CD spectrum?',
            FALSE),

        footer=tagList(
            actionButton('submitUnits', 'Submit'),
            modalButton('Cancel')
        )

    ))
  
})

# Modal dialog to ask for the molecular weight, concentration, path length, and
# number of chromophore units, if required.
observeEvent(input$submitUnits, {
  
  removeModal()

  c1 <- grepl('molar'  , input$workingUnitsAutomatic, ignore.case = TRUE)
  c2 <- grepl('unit'   , input$workingUnitsAutomatic, ignore.case = TRUE)

  if (c1 & c2) {

    showModal(modalDialog(

      tags$h3('Please enter the following experimental details:'),
      numericInput('pathLength'   , 'Path length (mm)',     1,min = 0,max = 100,step = 1),
      numericInput('molWeight'    , 'Molecular weight (Da)',1,min = 0,max = 1e6,step = 1),
      numericInput('concentration', 'Concentration (mg/ml)',1,min = 0,max = 100,step = 0.2),
      numericInput('numberOfC'   ,  'Number of chromophore units'   ,1,min = 0,max = 1e6,step = 1),
      footer=tagList(
        actionButton('submitExperimentalDetails', 'Submit'),
        modalButton('Cancel')
      )
    ))

  } else if (c1 & !c2) {

    showModal(modalDialog(

      tags$h3('Please enter the following experimental details:'),
      numericInput('pathLength'   , 'Path length (mm)',     1,min = 0,max = 100,step = 1),
      numericInput('molWeight'    , 'Molecular weight (Da)',1,min = 0,max = 1e6,step = 1),
      numericInput('concentration', 'Concentration (mg/ml)',1,min = 0,max = 100,step = 0.2),
      footer=tagList(
        actionButton('submitExperimentalDetails', 'Submit'),
        modalButton('Cancel')
      )
    ))

  } else {

    selectedIDx <- which(global_cd_units_choice == input$workingUnitsAutomatic)

    showModal(modalDialog(

      tags$h3(paste('Selected units for the final spectrum are: ' ,names(global_cd_units_choice)[selectedIDx])),
      footer=tagList(
        actionButton('submitExperimentalDetails', 'Okay'),
        modalButton('Cancel')
      )
    ))

  }
  
})

# Generate programatically one CD experiment,
# Add this experiment to the python class cdAnalyzer
append_cd_experiment <- function(expName,spectraNames,wlNew,signalNew,inputUnits,
                                 metadata_list,temperature=NULL,
                                 signalHT = NULL,isFakeExperiment=FALSE,
                                 fakeExperimentSignal=NULL) {
  
  signalNew <- as.matrix(signalNew)

  cdAnalyzer$experimentsOri[[expName]]                   <- CdExperimentGeneral()
  cdAnalyzer$experimentsOri[[expName]]$wavelength        <- np_array(wlNew)
  cdAnalyzer$experimentsOri[[expName]]$spectraNames      <- np_array(c(spectraNames))
  cdAnalyzer$experimentsOri[[expName]]$internalID        <- np_array(c(spectraNames))
  cdAnalyzer$experimentsOri[[expName]]$signalInput       <- np_array(signalNew)

  if (!is.null(temperature)) {
    cdAnalyzer$experimentsOri[[expName]]$temperature     <- np_array(temperature)
  }
  
  na_matrix <- matrix(NA, nrow = nrow(signalNew), ncol = ncol(signalNew))
  if (is.null(signalHT)) {
    cdAnalyzer$experimentsOri[[expName]]$signalHT          <- na_matrix
  } else {
    cdAnalyzer$experimentsOri[[expName]]$signalHT          <- np_array(signalHT)
  }
  
  # Assign the desired signal
  # For fake experiments, it can't be changed once it's created!!
  if (isFakeExperiment) {
    cdAnalyzer$experimentsOri[[expName]]$signalDesiredUnit    <- np_array(signalNew)
    cdAnalyzer$experimentsOri[[expName]]$fakeExperimentSignal <- fakeExperimentSignal
    cdAnalyzer$experimentsOri[[expName]]$signalAbs            <- na_matrix
    # Trick to avoid further processing with other spectra using
    # other operation units different from 'fakeExperimentSignal'
    cdAnalyzer$experimentsOri[[expName]]$concentration        <- -1
  } else {
    cdAnalyzer$experimentsOri[[expName]]$experiment_to_absorbance_units(inputUnits)
  }
  
  cdAnalyzer$experimentsOri[[expName]]$isGenerated       <- TRUE
  
  newExperimentNames <- c(expName)
  
  cdAnalyzer$experimentsOri[[expName]]$isFakeExperiment  <- isFakeExperiment
  
  # Start the lines vector, to fill the output file
  
  lines <- c(paste0('New generated experiment: ',expName))
  lines <- c(lines,
             paste0('Number of generated CD curves:  ',length(spectraNames)))
  
  if (!is.null(metadata_list)) {
    
    metadataFeature   <- names(metadata_list)
    metadataFeature   <- paste0('#',metadataFeature,' :')
    metadataFeature   <- identate(metadataFeature,48)
    metadataValue     <- unlist(metadata_list)
    metadata_df       <- data.frame(metadataFeature,metadataValue)
    
    for (i in 1:nrow(metadata_df)) {
      lines <- c(lines,paste0(metadata_df[i,],collapse = ''))
    }
  }
  
  append_record_to_logbook(lines)

  # end of add new experiment to the logbook
  
  metadata_list[['Timestamp']]       <- format(Sys.time(),usetz = TRUE)
  metadata_list[['CD curve origin']] <- 'ChiraKit (https://spc.embl-hamburg.de/)'
  
  cdAnalyzer$experimentsOri[[expName]]$metadata          <- metadata_list
  
  cdAnalyzer$experimentNames <- c(cdAnalyzer$experimentNames,expName)
  
  return(NULL)
}

observeEvent(input$submitExperimentalDetails,{
  
  req(input$submitExperimentalDetails)
  
  removeModal()
  
  sample_data_files     <- input$cdFilesSample$datapath
  scans_sample_names    <- sapply(input$cdFilesSample[[1]],remove_file_extension)
  
  baseline_data_files   <- input$cdFilesBaseline$datapath
  scans_bl_names        <- sapply(input$cdFilesBaseline[[1]],remove_file_extension)
  
  # Check that all the files are different
  intersectFiles <- intersect(scans_sample_names,scans_bl_names)
  
  if (length(intersectFiles) > 0) {
    
    shinyalert(text = paste("<b>Please use different files for the sample and baseline</b>"),
               type = "warning",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    
    return(NULL)
  }
  
  reactives$data_loaded      <- NULL

  names <- c(scans_sample_names,scans_bl_names)
  files <- c(sample_data_files,baseline_data_files)
  
  i <- 0
  # Iterate over the files
  for (name in names) {
    
    i <- i + 1
    cd_data_file   <- files[i]

    load_one_experiment(cd_data_file,name,input$inputUnitsAutomatic)

  }

  updateMaxVoltageValue()
  ## Average sample files

  # List of signal dataframes 
  signalSel <- lapply(scans_sample_names, function(exp) {
    cdAnalyzer$experimentsOri[[exp]]$signalInput})
  
  htSel <- lapply(scans_sample_names, function(exp) {
    cdAnalyzer$experimentsOri[[exp]]$signalHT})
  
  # List of temperatures
  temperatureSel <- lapply(scans_sample_names, function(exp) {
    cdAnalyzer$experimentsOri[[exp]]$temperature
  })
  
  # Assign the mean of the sample scans temperature to the rest of the curves
  temps          <- unlist(temperatureSel)
  
  if (all(is.null(temps))) {
    newTemp <- NA
  } else {
    newTemp        <- mean(temps,na.rm = T)
  }
  
  signalSelVectorized <- list()
  htSelVectorized     <- list()
  wlSelVectorized     <- list()

  spectra_counter    <- 0
  experiment_counter <- 0
  
  # Transform list of signal dataframes into a list of signal vectors
  for (signalDF in signalSel) {
    
    experiment_counter <- experiment_counter + 1
    expName            <- names[experiment_counter]
    htDF               <- htSel[[experiment_counter]]
      
    for (columnID in 1:ncol(signalDF)) {
      spectra_counter                        <- spectra_counter + 1
      signalSelVectorized[[spectra_counter]] <- signalDF[,columnID]
      htSelVectorized[[spectra_counter]]     <- htDF[,columnID]
      wlSelVectorized[[spectra_counter]]     <- cdAnalyzer$experimentsOri[[expName]]$wavelength
    }
    
  }
  
  interpolatedAverage <- averageListOfVectors(wlSelVectorized,signalSelVectorized)
  signalNewSample     <- matrix(interpolatedAverage$average,ncol = 1)
  
  wlNewSample         <- interpolatedAverage$x
  htNewSample         <- get_ht_average(wlSelVectorized,htSelVectorized)
  
  # Find if we have already used the spectra names we want to use 
  currentSpectraNames <- getSpectraNames(cdAnalyzer)
    
  time_str <- get_hour_minute_sec()
  
  sampleAvgName <- 'Sample average'
  if ('Sample average' %in% currentSpectraNames) {
    sampleAvgName <- paste0('Sample average ',time_str)
  } 
  
  processingUnitsString <- clean_html_sup_tag(
    workingUnits2ProperLabel(input$inputUnitsAutomatic)
  )
  
  metadata_list <- list(
    'Processing'           = 'Generated by averaging many CD curves',
    'Original CD files'    = paste0(scans_sample_names,collapse = ' '),
    'Units for processing' = processingUnitsString
  )
  
  append_cd_experiment(sampleAvgName,sampleAvgName,wlNewSample,
                       signalNewSample,input$inputUnitsAutomatic,
                       metadata_list,newTemp,htNewSample)
  
  ## End of average sample files
  
  ## Average baseline files
  
  # List of signal dataframes 
  signalSel <- lapply(scans_bl_names, function(exp) {
    cdAnalyzer$experimentsOri[[exp]]$signalInput})
  
  htSel    <- lapply(scans_bl_names, function(exp) {
    cdAnalyzer$experimentsOri[[exp]]$signalHT})
  
  signalSelVectorized <- list()
  htSelVectorized     <- list()
  wlSelVectorized     <- list()
  
  spectra_counter     <- 0

  n_loaded_experiments <- experiment_counter
  
  # Transform list of signal dataframes into a list of signal vectors
  for (signalDF in signalSel) {
    
    experiment_counter <- experiment_counter + 1
    expName            <- names[experiment_counter]
    htDF               <- htSel[[experiment_counter - n_loaded_experiments]]
    
    for (columnID in 1:ncol(signalDF)) {
      spectra_counter                        <- spectra_counter + 1
      signalSelVectorized[[spectra_counter]] <- signalDF[,columnID]
      htSelVectorized[[spectra_counter]]     <- htDF[,columnID]
      wlSelVectorized[[spectra_counter]]     <- cdAnalyzer$experimentsOri[[expName]]$wavelength
    }
    
  }
  
  interpolatedAverage <- averageListOfVectors(wlSelVectorized,signalSelVectorized)
  signalNewBaseline   <- matrix(interpolatedAverage$average,ncol = 1)
  wlNewBaseline       <- interpolatedAverage$x
  
  htNewBaseline       <- get_ht_average(wlSelVectorized,htSelVectorized)
  
  baselineAvgName <- 'Baseline average'
  if ('Baseline average' %in% currentSpectraNames) {
    baselineAvgName <- paste0('Baseline average ',time_str)
  }
  
  metadata_list <- list(
    'Processing'           = 'Generated by averaging many CD curves',
    'Original CD files'    = paste0(scans_sample_names,collapse = ' '),
    'Units for processing' = processingUnitsString
  )
  
  append_cd_experiment(baselineAvgName,baselineAvgName,wlNewBaseline,
                       signalNewBaseline,input$inputUnitsAutomatic,
                       metadata_list,newTemp,htNewBaseline)
  
  ## End of average baseline files
  
  ## Subtract the baseline
  
  interpolatedSignals <- interpolateVectors(wlNewSample,signalNewSample,wlNewBaseline,signalNewBaseline)
  
  signalNew <- interpolatedSignals$y1 - interpolatedSignals$y2
  wlNew     <- interpolatedSignals$x1
  
  metadata_list <- list(
    'Processing'           = 'Generated by subtracting the sample average from the baseline average',
    'Original CD files_sample'    = paste0(scans_sample_names,collapse = ' '),
    'Original CD files_baseline'  = paste0(scans_bl_names,collapse = ' '),
    'Units for processing'        = processingUnitsString
  )
  
  baselineSubName <- 'Baseline subtracted'
  if ('Baseline subtracted' %in% currentSpectraNames) {
    baselineSubName <- paste0('Baseline subtracted ',time_str)
  }
  
  append_cd_experiment(baselineSubName,baselineSubName,wlNew,
                       signalNew,input$inputUnitsAutomatic,
                       metadata_list,newTemp,htNewSample)
  
  ## End of subtract the baseline
  
  ## Zeroing to obtain the final spectra
  ## TO DO - add modal dialog!!
  zeroFinalSpectrum <- input$zeroFinalSpectrum
  if (zeroFinalSpectrum) {
    
    maxV <- max(unlist(wlSelVectorized))
    
    signalNew <- signalNew - mean(signalNew[wlNew > (maxV - 10)])
    
    wlRange <- paste0(maxV - 10,' - ',maxV, ' nm')
    
    metadata_list <- list(
      'Processing' = 'Generated by 1) subtracting the sample average from the baseline average and 2) subtracting the mean signal from the last 10 nm',
      'Original CD files of the sample'    = paste0(scans_sample_names,collapse = ' '),
      'Original CD files of the baseline'  = paste0(scans_bl_names,collapse = ' '),
      'Units for processing'               = processingUnitsString,
      'Wavelength range to calculate offset' = wlRange
    )
    
    baselineZeroedName <- 'Baseline subtracted with offset'
    if ('Baseline subtracted with offset' %in% currentSpectraNames) {
      baselineZeroedName <- paste0('Baseline subtracted with offset ',time_str)
    }
    
    append_cd_experiment(baselineZeroedName,baselineZeroedName,wlNew,
                         signalNew,input$inputUnitsAutomatic,
                         metadata_list,newTemp,htNewSample)
    
  }
  
  ## End of zeroing to obtain the final spectra
  
  newExperimentNames <- c(baselineAvgName,sampleAvgName,baselineSubName)
  
  if (zeroFinalSpectrum) {
    newExperimentNames <- c(newExperimentNames,baselineZeroedName)
  }

  # We need to use lapply instead of a for-loop to avoid R lazy evaluation!
  lapply(newExperimentNames, function (exp) {
    
    saneExpName <- gsub(':','',exp)
    
    tabP <- tabPanel(saneExpName,value=saneExpName,fluidRow(column(12,tableOutput(paste0('metadata_',saneExpName)))))
    
    appendTab("metadata",tabP,select=TRUE)
    
    metadata_info     <- cdAnalyzer$experimentsOri[[exp]]$metadata
    metadataFeature   <- (names(metadata_info))
    metadataValue     <- unlist(metadata_info)
    
    metadata_df <- data.frame(metadataFeature,metadataValue)
    colnames(metadata_df) <- c('Metadata feature','Value') 
    
    output[[paste0('metadata_',saneExpName)]] <- renderTable({metadata_df})
    return(NULL)
    
  })
  
  # Add the information about the molecular weight, concentration, path length
  # and number of chromophore units, if required.
  
  c1 <- grepl('molar'  , input$workingUnitsAutomatic, ignore.case = TRUE)
  c2 <- grepl('unit'   , input$workingUnitsAutomatic, ignore.case = TRUE)
  
  cdAnalyzer$initialize_experiment_modif()
  
  # Spectra to convert into molar ellipticity / extinction
  experimentsToModifyParams    <- c(baselineSubName)
  
  if (zeroFinalSpectrum) {
    experimentsToModifyParams    <- c(experimentsToModifyParams,baselineZeroedName)
  }
 
  if (c1) {
    
    for (name in experimentsToModifyParams) {
      
      cdAnalyzer$set_experiment_properties(name,'concentration',   input$concentration)
      cdAnalyzer$set_experiment_properties(name,'molecularWeight', input$molWeight)
      cdAnalyzer$set_experiment_properties(name,'pathLength',      input$pathLength / 10) # Divide by 10 (from mm to cm)
      
      append_record_to_logbook(paste0('Setting concentration (mg/ml) of '    ,name,' to ',input$concentration))
      append_record_to_logbook(paste0('Setting molecular weight (Dalton) of ',name,' to ',input$molWeight))
      append_record_to_logbook(paste0('Setting path length (cm) of '         ,name,' to ',input$pathLength / 10))
      
    }
    
  }
  
  if (c2) {
    
    for (name in experimentsToModifyParams) {
      
      cdAnalyzer$set_experiment_properties(name,'numberOfCroms', input$numberOfC)
      append_record_to_logbook(paste0('Setting chromophore units of '    ,name,' to ',input$numberOfC))
    }
  }
  
  # Convert to selected units
  
  updateSelectInput(session,'workingUnits',NULL,
                    choices = getChoices(input$workingUnitsAutomatic))
  
  cdAnalyzer$all_experiments_absorbance_units_to_other_units(input$workingUnitsAutomatic)
  
  Sys.sleep(0.5)
  renderInputData()
  
  Sys.sleep(0.5)
  reactives$data_loaded     <- TRUE
  
})


observeEvent(input$triggerDeletion,{
  req(reactives$data_loaded)

  selected <- input$experiment2delete
  
  # Do not delete anything
  if (selected == "None") return(NULL)
  # Delete all experiments
  if (selected == "ALL") selected <- cdAnalyzer$experimentNames
  
  # Delete the legend table
  output$legendInfo <- NULL
  cdAnalyzer$delete_experiment(selected)
  
  # Remove the metadata Tabs
  for (sel in c(selected)) {
    
    target <-  gsub(':','',sel) # Remove the ':' character which can't be used in the Shiny Tabs names
    removeTab(inputId = "metadata", target = target)
    
    append_record_to_logbook(paste0('Deleting experiment ',sel))
    
  }
  
  exps <- cdAnalyzer$experimentNames

  if (length(exps) > 0) {
    
    updateSelectInput(session,"experiment2delete",selected = "ALL",
                      choices = c("ALL",cdAnalyzer$experimentNames))
    
    updateSelectInput(session,"selected_cd_exp",choices     = cdAnalyzer$experimentNames)
    
    # Create the legend table
    legendDf              <- getPlottingDF(cdAnalyzer)
    output$legendInfo     <- helperRenderRHandsontable(legendDf)
    updateSESCA_ref(legendDf)

    # If the parameters are not shared, then
    # update the Table with the parameters MW, #AA, etc.
    if (!cdAnalyzer$sharedParameters) {
      
      updateCDFilesInfoTable()
      
    }
      
    # Update the Table with the possible spectra to be processed
    updateProcessingTable()
    
  } else {
    output$cdFilesInfo    <- NULL
    output$proccesingInfo <- NULL

    updateSliderInput(session,'wavelengthRange',label=NULL,min=180,max=300,value = c(180,300))

    updateSelectInput(session,"experiment2delete",choices = c('None'))
    reactives$data_loaded <- FALSE

    reactives$show_example_data <- TRUE

  }
  
})

convertExperimentToWorkingUnits <- function() {
  
  # Convert everything (except 'fake' experiments) to absorbance units (from input units)
  namesAll  <- cdAnalyzer$experimentNames
  
  trueExperiments   <- !unlist(cdAnalyzer$get_experiment_properties('isFakeExperiment'))
  names             <- namesAll[trueExperiments]
  namesN            <- length(names)
  
  # Check if the input units are shared and convert to absorbance
  if (cdAnalyzer$sharedParameters) {
    
    cdAnalyzer$all_experiments_to_absorbance_units(input[["inputUnitsAll"]])
    append_record_to_logbook(paste0('Setting input units to ',input[["inputUnitsAll"]]))
    
  } else {
    
    rList <- list()
    for(i in 1:namesN) {
      selectedInputUnits <- input[[(paste0("inputUnits", i))]]
      req(selectedInputUnits)
      rList[[i]]         <- selectedInputUnits
    }
    if (length(rList) != namesN) return(NULL)
    
    names(rList) <- names
    
    modifiedExperiments <- cdAnalyzer$experiments_to_absorbance_units(rList)
    append_record_to_logbook(paste0('Setting input units of ',modifiedExperiments))
  }
  
  # Convert everything from absorbance units to selected units
  cdAnalyzer$all_experiments_absorbance_units_to_other_units(input[["workingUnits"]])

  append_record_to_logbook(paste0('Converting CD units to: ',input[["workingUnits"]]))
  
  wlRange <- input$wavelengthRange
  cdAnalyzer$filter_data_by_wavelength(wlRange[1], wlRange[2])

}

observeEvent(input$loadExampleData,{

  if (length(cdAnalyzer$experimentNames) > 0) {

    popUpWarning('You already have data loaded. Please remove all the current experiments to load the example data.')

    return(NULL)
  }

  showModal(modalDialog(

    tags$h3('What kind of example would you like?'),

    # Include checkbox to track the previous processing step in the metadata
    tags$h4(selectInput('exampleDataType',NULL,c(
    'Basic sample and baseline'             = 'basic',
    'Spectrum in delta epsilon units'       = 'secondary',
    'Thermal ramp'                = 'thermal',
    'Chemical denaturation' = 'chemical'
    ))),

    footer=tagList(
      actionButton('submitExampleData', 'Submit'),
      modalButton('Cancel')
    )

  ))

})

observeEvent(input$submitExampleData,{

  removeModal()

  output$proccesingInfo <- NULL
  reactives$data_loaded <- NULL

  # Check that the working units are millidegrees to load the chemical denaturation example
  if (input$exampleDataType == 'chemical' & input$workingUnits != 'millidegrees') {

      popUpWarning('Please set the working units to millidegrees to load the chemical denaturation example.')
      return(NULL)
  }

  inputUnits <- 'millidegrees'

  exampleDataType <- input$exampleDataType

  if (exampleDataType == 'basic') {

    cd_data_files   <- c('./www/R78906.d01','./www/R78906.d02','./www/R78906.d03',
                        './www/R78907.d01','./www/R78907.d02','./www/R78907.d03')

    names           <- c('R78906.1 (buffer)','R78906.2 (buffer)','R78906.3 (buffer)','R78907.1 (sample)','R78907.2 (sample)','R78907.3 (sample)')

  }

  if (exampleDataType == 'thermal') {

    cd_data_files   <- c('./www/Lys-Tscan-Oct2024_20241101.csv')
    names           <- c('Lysozyme thermal ramp')

    updateSelectInput(session,
                      inputId="oligomeric_state_term",
                      label=NULL,
                      selected='monomer'
                      )

  }

  if (exampleDataType == 'chemical') {

    cd_data_files   <- c('./www/urea_chc.csv')
    names           <- c('Urea denaturation')

    updateSelectInput(session,
                      inputId="oligomeric_state_chem",
                      label=NULL,
                      selected='monomer'
                      )

  }

  if (exampleDataType == 'secondary') {

    cd_data_files   <- c('./www/Myoglobin.csv')
    names           <- c('Myoglobin')
    inputUnits      <- 'meanUnitMolarExtinction'

  }

  updateSelectInput(session,"workingUnits",NULL,choices = getChoices(inputUnits))
  Sys.sleep(1)

  i <- 0

  # iterate over the files
  for (name in names) {

    i              <- i + 1
    cd_data_file   <- cd_data_files[i]
    load_one_experiment(cd_data_file,name,inputUnits)

  }

  if (exampleDataType == 'basic') {

    popUpInfo('The basic sample and baseline data was loaded.
    If you are new to ChiraKit, we recommend following the Tutorial #1, starting from step 3.')

  }

  if (exampleDataType == 'secondary') {

    popUpInfo("The CD spectrum of myoglobin was loaded.
    We recommend navigating to the 'Secondary structure content' module (Section '2. Analysis').
    To quickly estimate the secondary structure content, press on 'Run estimation!'.")

  }

  if (exampleDataType == 'thermal') {

    popUpInfo("The thermal ramp data of lysozyme was loaded.
    We recommend navigating to the 'Thermal unfolding' module (Section '2. Analysis'),
    pressing on '2a. Create dataset', and changing the 'Analysis type'
    to 'Spectra decomposition (SVD)'.")

    temperature <- unlist(cdAnalyzer$get_experiment_properties('temperature'))
    id          <- unlist(cdAnalyzer$get_experiment_properties('spectraNames'))

    df           <- data.frame(id,temperature,'A')
    colnames(df) <- c('CD_curve','Temperature (°C/K)','Dataset')

  # Assign the created dataframe to the Table thermal_denaturation_data (available at the 2a. Thermal analysis Tab)
    output$thermal_denaturation_data <- renderRHandsontable({
      rhandsontable(df,rowHeaders=NULL)    %>%
        hot_col(col = c(1), readOnly=TRUE) %>%
        hot_table(stretchH='all')
    })

  }

  if (exampleDataType == 'chemical') {

    Sys.sleep(1)
    popUpInfo("The chemical denaturation data of the Clathrin heavy chain N-terminal domain was loaded.
    We recommend navigating to the 'Chemical unfolding' module (Section '2. Analysis'),
    and pressing on '2a. Create dataset'.
    To quickly fit the data, select the three-state model (N <-> I <-> U), allow fitting both slopes and use '2'
    as the initial guess for the parameter D50 (where ΔG1 equals 0).")

    reactives$urea_example_data_loaded <- TRUE
 } else {
    reactives$urea_example_data_loaded <- FALSE
 }

  Sys.sleep(0.2)

  # Copy experiments to allow modifying the wavelength range
  cdAnalyzer$initialize_experiment_modif()
  updateMaxVoltageValue()
  renderInputData()
  reactives$data_loaded <- TRUE

},priority = 10)

# React to changes in molecular weight, path length or concentration 
observeEvent(input$cdFilesInfo_cell_edit, {

  req(reactives$data_loaded)
  
  reactives$data_loaded <- NULL
  info <- input$cdFilesInfo_cell_edit

  columnIndex <- info$col[1] + 1 # Start index from 1
  
  rowIndex    <- info$row[1]
  value       <- as.numeric(info$value[1])

  if (cdAnalyzer$sharedParameters) {
    
    expNames          <- cdAnalyzer$experimentNames
    trueExperiments   <- !unlist(cdAnalyzer$get_experiment_properties('isFakeExperiment'))
    expNames          <- expNames[trueExperiments]
  } else {
    expNames <- c(cdAnalyzer$experimentNames[rowIndex])
  }
  
  for (expName in expNames) {
    
    valueNew <- value
    if (columnIndex == 2) {
      selectedVar <- "molecularWeight"
      selectedVar2print <- 'molecular weight (Da)'
    } else if (columnIndex == 3) {
      selectedVar <- "numberOfCroms"
      selectedVar2print <- 'chromophore units number'
    } else if (columnIndex == 4) {
      selectedVar <- "concentration"
      selectedVar2print <- 'concentration (mg/ml)'
    } else { 
      selectedVar  <- "pathLength"
      selectedVar2print <- 'path length (cm)'
      valueNew <- valueNew / 10 # To cm!
    } 
    
    cdAnalyzer$set_experiment_properties(expName,selectedVar,valueNew)
    
    append_record_to_logbook(paste0('Setting ',selectedVar2print,' of ',
                                    expName,' to ',valueNew))
    
  }
  
  reactives$data_loaded <- NULL
  convertExperimentToWorkingUnits()
  reactives$data_loaded <- TRUE
})

observeEvent(input$workingUnits,{
  
  req(reactives$data_loaded)
  reactives$data_loaded <- NULL
  # Convert everything from absorbance units to selected units
  cdAnalyzer$all_experiments_absorbance_units_to_other_units(input[["workingUnits"]])
  
  append_record_to_logbook(paste0('Converting CD units to: ',input[["workingUnits"]]))
  
  wlRange <- input$wavelengthRange
  cdAnalyzer$filter_data_by_wavelength(wlRange[1], wlRange[2])

  reactives$data_loaded <- TRUE
  
})

observeEvent(input$sharedExperimentParameters,{
  
  req(reactives$data_loaded)
  reactives$data_loaded <- FALSE
  
  cdAnalyzer$sharedParameters <- input$sharedExperimentParameters
  
  # Force CD samples to share the same molecular weight, concentration, path length, input units and concentration
  if (input$sharedExperimentParameters) {
    
    df   <- data.frame('All experiments',1,1,1,
                       as.character(selectInput('inputUnitsAll',label=NULL, choices = getChoices('Millidegrees'))),
                       1)
    
    colnames(df) <- c('File name',"Mol. weight (Dalton)","#Chromophore units","Conc. (mg/ml)",
                      "Input units","Path length (mm)")
    
    # Set default values
    exps              <- cdAnalyzer$experimentNames
    trueExperiments   <- !unlist(cdAnalyzer$get_experiment_properties('isFakeExperiment'))
    exps              <- exps[trueExperiments]
    
    for (exp in exps) {
      
      cdAnalyzer$set_experiment_properties(exp,'units'            , 'millidegrees')
      cdAnalyzer$set_experiment_properties(exp,'numberOfCroms'    , 1)
      cdAnalyzer$set_experiment_properties(exp,'concentration'    , 1)
      cdAnalyzer$set_experiment_properties(exp,'pathLength'       , 0.1)
      cdAnalyzer$set_experiment_properties(exp,'molecularWeight'  , 1)
      
      append_record_to_logbook(paste0('Setting number of chromophore units / concentration (mg/ml)/ molecular weight (Da) of ',
                                      exp,' to 1'))
      
      append_record_to_logbook(paste0('Setting path length (cm) of: ',
                                      exp,' to 0.1'))
      
    }
    
    cdAnalyzer$all_experiments_to_absorbance_units('millidegrees')
    
    append_record_to_logbook('Setting all experiment input units to millidegrees')
    
                       
  } else {
    
    df <- generateDTtable(cdAnalyzer)
  
  }
  
  cdAnalyzer$all_experiments_absorbance_units_to_other_units(input$workingUnits)

  append_record_to_logbook(paste0('Converting CD units to: ',input[["workingUnits"]]))
  
  # To allow the creation of selectInput inside DT table
  session$sendCustomMessage('unbind-DT', 'cdFilesInfo')
  output$cdFilesInfo <- renderDTtable(df,TRUE)
  reactives$data_loaded <- TRUE
  
})

observeEvent(input$legendInfo,{
  
  req(reactives$data_loaded)

  updateSelectInput(session,"mol2changeColor",NULL,
                    get_legend_from_rhandTable(input$legendInfo),input$mol2changeColor)
  

})

observeEvent(input$legendInfo$changes$changes, {
  
  # Get the changes made to the table
  changes <- input$legendInfo$changes$changes[[1]]
  
  # Check if any changes were made in the column with the legends
  if (any(changes[[2]] == 2)) {

    names          <- cdAnalyzer$experimentNames
    internalIDs    <- cdAnalyzer$get_experiment_properties('internalID')
    spectraID      <- get_ID_from_rhandTable(input$legendInfo)[changes[[1]]+1]
    newSpectraName <- get_legend_from_rhandTable(input$legendInfo)[changes[[1]]+1]
    
    id                 <- which(unlist(internalIDs) == spectraID)
    names(internalIDs) <- names
    id                 <- found_ids(internalIDs,c(id))
    
    append_record_to_logbook(paste0("Changing spectrum name from '",
                                    cdAnalyzer$experimentsOri[[names(id)[1]]]$spectraNames[id[1]],
                                    "' to '",newSpectraName,"'"))
    
    cdAnalyzer$experimentsOri[[names(id)[1]]]$spectraNames[id[1]]   <- newSpectraName
    cdAnalyzer$experimentsModif[[names(id)[1]]]$spectraNames[id[1]] <- newSpectraName

    legendDf <- getLegendDF(input$legendInfo)
    updateSESCA_ref(legendDf)

  }

})

observeEvent(input$colorForLegend,{
  
  req(reactives$data_loaded)
  legendDf <- getLegendDF(input$legendInfo)
  
  idx                  <- which(legendDf$Legend == input$mol2changeColor)
  legendDf$Color[idx]  <- input$colorForLegend
  output$legendInfo    <- helperRenderRHandsontable(legendDf)

})

observeEvent(input$triggerProcessing,{
  
  req(reactives$data_loaded)
  
  op <- input$operation
  
  fourth_col <- 'Spectrum 2: '
  if (op == 'Batch average') fourth_col <- ''
  if (op == 'Zero')          fourth_col <- 'Window range (nm): '
  if (op == 'Smooth')        fourth_col <- 'Window size (nm):  '
  if (op == 'Scale')         fourth_col <- 'Scaling factor:  '

  showModal(modalDialog(
    
    tags$h3('Please confirm the following processing step:'),
    tags$h4(paste0('Operation:          ',op)),
    tags$h4(paste0('Spectrum/a 1:       ',input$inputSpectra1)),
    tags$h4(HTML(paste0('Units of operation: ',workingUnits2ProperLabel(input$operationUnits)))),
    tags$h4(paste0(fourth_col,input$inputSpectra2)),
    tags$h4(paste0('Output name:        ',input$newSpectraName)),

    tags$h4(tags$br()),

    # Include checkbox to track the previous processing step in the metadata
    if (op %in% c('Sum','Subtract','Zero','Smooth','Scale')) {
      tags$h4(checkboxInput('trackPreviousChange','Include the previous processing step in the metadata',value = TRUE))
    },

    footer=tagList(
      actionButton('submitUnits2', 'Submit'),
      modalButton('Cancel')
    )
  ))
  
})

observeEvent(input$submitUnits2,{
  
  removeModal()
  # Retrieve the experiment names 
  names          <- cdAnalyzer$experimentNames
  
  # Retrieve the internal IDs
  internalIDs        <- cdAnalyzer$get_experiment_properties('internalID')
  names(internalIDs) <- names
  
  # Retrieve the user selected label for the new spectrum / spectra
  newSpectraName <- input$newSpectraName
  
  if (newSpectraName %in% names) {
    newSpectraName <- paste0(newSpectraName,' ',get_hour_minute_sec())
  }
  
  # Find the second spectrum label, the desired operation and the units of operation
  spectra2       <- input$inputSpectra2
  operation      <- input$operation
  operationUnits <- input$operationUnits
  
  legendDf <- getLegendDF(input$legendInfo)
  # Find the spectrum / spectra to be processed
  if (input$inputSpectra1 == "Selected in 'Show' Column") {
    
    # Sanity check - Stop if we haves less than 2 spectra and we selected the option "Selected in 'Show' Column"
    if (sum(legendDf$Show) < 2) {
      
      shinyalert(text = paste(
        "<b>Please select at least two spectra to use this option.</b>"),
        type = "info",closeOnEsc = T,closeOnClickOutside = T,
        html=T)
      
      updateCDFilesInfoTable()
      updateProcessingTable()
      
      return(NULL)
      
    }

    # Spectra from 'Show' legend
    spectra1 <- legendDf$Internal.ID[legendDf$Show]
    
  } else {
    
    spectra1 <- c(input$inputSpectra1)
    
  } 
  
  # Start of -- 1st step Find the first and second spectrum IDs
  
  id1s <- c()
  
  for (spectrum in spectra1) {
    id1s <- c(id1s,which(unlist(internalIDs) == spectrum))
  }
  
  # Find the second spectrum, if required
  
  if (operation %in% c('Sum','Subtract')) {
    id2                <- which(unlist(internalIDs) == spectra2)
    ids                <- found_ids(internalIDs,c(id1s,id2))
  } else {
    ids                <- found_ids(internalIDs,c(id1s))
  }
  
  # Find the experiment names that are required for the processing
  selected_id <- names %in% (unique(names(ids)))  

  # End of -- 1st step Find the first and second spectrum IDs

  # Decide if the generated experiments are of the type 'fake'
  # In other words, the calculated spectra will not be available for further processing
  operationIsInMolarUnits       <- grepl('molar', operationUnits, ignore.case = TRUE)
  operationIsInMeanUnits        <- grepl('mean' , operationUnits, ignore.case = TRUE)
  
  molWeights   <- (unlist(cdAnalyzer$get_experiment_properties('molecularWeight'  )))[selected_id]
  pathLengths  <- (unlist(cdAnalyzer$get_experiment_properties('pathLength'       )))[selected_id]
  concs        <- (unlist(cdAnalyzer$get_experiment_properties('concentration'    )))[selected_id]
  numberOfCU   <- (unlist(cdAnalyzer$get_experiment_properties('numberOfCroms' )))[selected_id]
  
  ##  Check that we can process normalized CD units
  
  mw_is_not_zero <- molWeights   != 0 & !is.na(molWeights)  
  pl_is_not_zero <- pathLengths  != 0 & !is.na(pathLengths) 
  cn_is_not_zero <- concs        != 0 & !is.na(concs)       
  nr_is_not_zero <- numberOfCU   != 0 & !is.na(numberOfCU) 
  
  check_list <- c()
  
  if (operationIsInMolarUnits) {
    
    check_list <- c(check_list,mw_is_not_zero,pl_is_not_zero,cn_is_not_zero)
    
  }
  
  if (operationIsInMeanUnits) {
    
    should_stop_process <- c(check_list,nr_is_not_zero)  
    
  }
  
  should_stop_process <- !all(check_list)

  # Check all the input units of the selected experiments
  all_input_units <- (unlist(cdAnalyzer$get_experiment_properties('units')))[selected_id]

  # Allow the process if all the input units match the operation units - because we don't need to convert to absorbance
  should_allow_process <- all(all_input_units == operationUnits)

  if (should_stop_process && !should_allow_process) {
    
    shinyalert(text = 
    paste("<b>Please add information about the experimental parameters.
          For example, the path length can not be zero.</b>"),
    type = "warning",closeOnEsc = T,closeOnClickOutside = T,
    html=T)
    
    updateCDFilesInfoTable()
    updateProcessingTable()
    
    return(NULL)
    
  }

  # Special check for fake experiments which can only be processed using the same units
  # that were used to generate them
  areSelSpecFake <- (unlist(cdAnalyzer$get_experiment_properties('isFakeExperiment')))[selected_id]
  fakeSignals     <- (unlist(cdAnalyzer$get_experiment_properties('fakeExperimentSignal')))[selected_id]
  
  if (any(areSelSpecFake)) {
    
    fakeSignals <- fakeSignals[fakeSignals != 'unknown']
    
    if (any(fakeSignals != operationUnits)) {
      
      shinyalert(text = 
      paste("<b>The desired operation units can't be used
      with the selected spectra. Please check the User Guide 
            ('Processing' Section).</b>"),
      type = "warning",closeOnEsc = T,closeOnClickOutside = T,
      html=T)
      
      updateCDFilesInfoTable()
      updateProcessingTable()
      
      return(NULL)
    }
  }
  
  append_record_to_logbook(paste0('A processing step was run with the following options: [',
                                  paste0(spectra1,collapse = ' ; '),'] [',
                                  operation,'] [',
                                  operationUnits,'] [',
                                  spectra2,']',
                                  '[Linear interpolation mode: ',input$allowLinearInterpolation,']'))
  
  # Check if they are equal and decide if we can further process 
  # the generated spectra or not ... 
  
  isFake <- length(c(unique(molWeights),unique(pathLengths),
                      unique(concs),unique(numberOfCU))) != 4    
  
  ## End of check that we can process normalized CD units
  
  # 2nd step - convert to the appropriate operation units
  rList2        <- as.list(rep(operationUnits,length(unique(names(ids)))))
  names(rList2) <- unique(names(ids))

  # Try to convert the units of the selected spectra
  cdAnalyzer$experiments_absorbance_units_to_other_units(rList2)

  # 3rd step - Retrieve the desired signals and process them

  if (operation %in% c('Sum','Subtract')) {
    
    id1s  <- head(ids, -1)        # Remove the last element (spectra #2), if required
    exp1s <- head(names(ids), -1) # Remove the last element (spectra #2), if required
    
  } else {
    id1s  <- ids
    exp1s <- names(ids)
  }
  
  # List of wavelength vectors - spectra 1
  wlSel     <- lapply(exp1s, function(exp) cdAnalyzer$experimentsOri[[exp]]$wavelength)

  # List of signal vectors - spectra 1
  signalSel <- lapply(1:length(exp1s), function(i) {
    cdAnalyzer$experimentsOri[[exp1s[i]]]$signalDesiredUnit[,id1s[i]]
  })

  # Check that we generated the signal in absorbance units and not NAs... because of missing metadata
  for (signal in signalSel) {

    # Produce an error if we only have NAs
    if (all(is.na(signal))) {

      shinyalert(text =
      paste("<b>Processing aborted: The selected spectra have missing metadata
      to use the selected operation units.
      Please add information about the experimental parameters.</b>"),
      type = "warning",closeOnEsc = T,closeOnClickOutside = T,
      html=T)

      #updateCDFilesInfoTable()
      #updateProcessingTable()

      return(NULL)
    }

  }

  htSel  <- lapply(1:length(exp1s), function(i) {
    cdAnalyzer$experimentsOri[[exp1s[i]]]$signalHT[,id1s[i]]
  })
  
  # List of temperature - spectra 1
  temperatureSel <- lapply(1:length(exp1s), function(i) {
    cdAnalyzer$experimentsOri[[exp1s[i]]]$temperature[id1s[i]]
  })
  
  # Check if all the wavelength vectors are equal
  if (!input$allowLinearInterpolation) {
    
    wlToCheck <- wlSel
    
    if (operation %in% c('Sum','Subtract')) {
      
      # Add  the second spectrum to the wavelength vectors list which is going to be checked
      wlToCheck[[length(wlToCheck)+1]] <- cdAnalyzer$experimentsOri[[last(names(ids))]]$wavelength
      
    }
    
    all_equal <- all(sapply(wlToCheck[-1], function(vec) identical(wlToCheck[[1]], vec)))
    
    # Abort processing if the wavelength vectors are not equal
    if (!all_equal) {
      
      shinyalert(text = 
      paste("<b>Processing aborted: The selected spectra were measured 
      at different wavelengths. 
      To continue you need to allow linear interpolation 
            (not recommended for baseline/reference subtraction).</b>"),
      type = "warning",closeOnEsc = T,closeOnClickOutside = T,
      html=T)
      
      updateCDFilesInfoTable()
      updateProcessingTable()
      return(NULL)
    }
      
  }
    
  # Change data_loaded to NULL to stop other reactives
  reactives$data_loaded <- NULL   
  
  processingUnitsString <- clean_html_sup_tag( workingUnits2ProperLabel(operationUnits) )
  
  # Case 1 - we generated an average spectrum
  
  if (operation == 'Average') {
    
    interpolatedAverage <- averageListOfVectors(wlSel,signalSel)
    signalNew           <- matrix(interpolatedAverage$average,ncol = 1)
    wlNew               <- interpolatedAverage$x
    
    htNew <- get_ht_average(wlSel,htSel)
    
    # 5th step (Case 1) - Create one new experiment per selected experiment
    newExperimentNames <- c(newSpectraName)
    legendDf <- appendRowsToLegendDataFrame(legendDf,c(newSpectraName))
    
    metadata_list <- list(
      'Processing'           = 'Generated by averaging many CD curves',
      'Original CD curve(s)'   = paste0(spectra1,collapse = ' '),
      'Units for processing' = processingUnitsString,
      'Linear interpolation mode' = input$allowLinearInterpolation
      )
    
    newTemp              <- mean(unlist(temperatureSel))
    
    append_cd_experiment(newSpectraName,newSpectraName,wlNew,signalNew,
                         operationUnits,metadata_list,newTemp,
                         htNew,isFake,operationUnits)
    
    # Case 2 - Batch average
  } else if (operation == 'Batch average') {
    
    # In this case, spectra 2 should be the string 'N = ' followed by a number
    n <- as.numeric(strsplit(spectra2,'=')[[1]][2])
    
    total_batches <- floor(length(spectra1) / n)
    newExperimentNames <- c()
    
    for (i in 1:total_batches) {
      
      newSpectraNameBatch <- paste0('Batch ',i,' ',newSpectraName)

      if (newSpectraNameBatch %in% names) {
        newSpectraNameBatch <- paste0(newSpectraNameBatch,' ',get_hour_minute_sec())
      }
      
      wlSelSubset          <- wlSel[(1+n*(i-1)):(n*i)]
      signalSelSubset      <- signalSel[(1+n*(i-1)):(n*i)]
      htSelSubset          <- htSel[(1+n*(i-1)):(n*i)]
      temperatureSelSubset <- temperatureSel[(1+n*(i-1)):(n*i)]
      newTemp              <- mean(unlist(temperatureSelSubset))
        
      interpolatedAverage  <- averageListOfVectors(wlSelSubset,signalSelSubset)
      signalNew            <- matrix(interpolatedAverage$average,ncol = 1)
      wlNew                <- interpolatedAverage$x
      
      htNew <- get_ht_average(wlSelSubset,htSelSubset)
      
      # 5th step (Case 1) - Create one new experiment per selected experiment
      metadata_list <- list(
        'Processing'           = 'Generated by averaging many CD curves',
        'Original CD curve(s)'   = paste0(spectra1[(1+n*(i-1)):(n*i)],collapse = ' '),
        'Units for processing' = processingUnitsString,
        'Linear interpolation mode' = input$allowLinearInterpolation
      )
      
      append_cd_experiment(newSpectraNameBatch,newSpectraNameBatch,wlNew,
                           signalNew,operationUnits,metadata_list,newTemp,
                           htNew,isFake,operationUnits)
      
      newExperimentNames <- c(newExperimentNames,newSpectraNameBatch)
      legendDf           <- appendRowsToLegendDataFrame(legendDf,c(newSpectraNameBatch))
      
    }
    
  } else if (operation == 'Split-pair subtraction  (1st - N/2+1, 2nd - N/2+2, ...)') {

    total_batches <- floor(length(signalSel) / 2)

    newExperimentNames <- c()

    wlSelSubset1          <- wlSel[1:(total_batches)]
    signalSelSubset1      <- signalSel[1:(total_batches)]

    wlSelSubset2          <- wlSel[(total_batches+1):(2*total_batches)]
    signalSelSubset2      <- signalSel[(total_batches+1):(2*total_batches)]

    for (i in 1:total_batches) {

      newSpectraNameBatch <- paste0(i,' ',newSpectraName)

      if (newSpectraNameBatch %in% names) {
        newSpectraNameBatch <- paste0(newSpectraNameBatch,' ',get_hour_minute_sec())
      }

      y1   <- signalSelSubset1[[i]]
      y2   <- signalSelSubset2[[i]]

      x1  <- wlSelSubset1[[i]]
      x2  <- wlSelSubset2[[i]]

      htNew                <- htSel[[i]]

      if (input$allowLinearInterpolation) {

        interpolatedSignals <- interpolateVectors(x1,y1,x2,y2)

        y1                  <- interpolatedSignals$y1
        y2                  <- interpolatedSignals$y2

        wlNew               <- interpolatedSignals$x1

        if (all(is.na(htNew))) {

          htNew <- rep(NA,length(interpolatedSignals$x1))

        } else {

            # We do not sum the HT signals!
            # We return the HT signal of the 'left' selected experiment
            # The HT signal data points are interpolated (if required)
            # to match the wavelength of the new generated experiment
          htNew  <- approx(x = x1,y = htNew,xout = wlNew , method = "linear")$y

        }

      } else {

        wlNew                <- wlSelSubset1[[i]]

      }

      signalNew            <- y1 - y2

      newTemp              <- temperatureSel[[i]]

      signalNew            <- matrix(signalNew,ncol = 1)
      htNew                <- matrix(htNew,ncol = 1)

      # 5th step (Case 1) - Create one new experiment per selected experiment
      metadata_list <- list(
        'Processing'  = 'Generated by subtraction ',
        'Curve 1'     = spectra1[i],
        'Curve 2 (to be subtracted)'     = spectra1[i+total_batches],
        'Units for processing' = processingUnitsString,
        'Linear interpolation mode' = input$allowLinearInterpolation
      )

      append_cd_experiment(newSpectraNameBatch,newSpectraNameBatch,wlNew,
                           signalNew,operationUnits,metadata_list,newTemp,
                           htNew,isFake,operationUnits)

      newExperimentNames <- c(newExperimentNames,newSpectraNameBatch)
      legendDf           <- appendRowsToLegendDataFrame(legendDf,c(newSpectraNameBatch))

    }

  } else if (operation == 'Consecutive subtraction (1st - 2nd, 3rd - 4th, ...)') {

    total_batches <- floor(length(signalSel) / 2)

    newExperimentNames <- c()

    wlSelSubset1          <- wlSel[seq(1,length(wlSel),2)]
    signalSelSubset1      <- signalSel[seq(1,length(wlSel),2)]

    wlSelSubset2          <- wlSel[seq(2,length(wlSel),2)]
    signalSelSubset2      <- signalSel[seq(2,length(wlSel),2)]

    spectraUneven         <- spectra1[seq(1,length(spectra1),2)]
    spectraEven           <- spectra1[seq(2,length(spectra1),2)]

    for (i in 1:total_batches) {

      newSpectraNameBatch <- paste0(i,' ',newSpectraName)

      if (newSpectraNameBatch %in% names) {
        newSpectraNameBatch <- paste0(newSpectraNameBatch,' ',get_hour_minute_sec())
      }

      y1   <- signalSelSubset1[[i]]
      y2   <- signalSelSubset2[[i]]

      x1  <- wlSelSubset1[[i]]
      x2  <- wlSelSubset2[[i]]

      htNew                <- htSel[[i]]

      if (input$allowLinearInterpolation) {

        interpolatedSignals <- interpolateVectors(x1,y1,x2,y2)

        y1                  <- interpolatedSignals$y1
        y2                  <- interpolatedSignals$y2

        wlNew               <- interpolatedSignals$x1

        if (all(is.na(htNew))) {

          htNew <- rep(NA,length(interpolatedSignals$x1))

        } else {

            # We do not sum the HT signals!
            # We return the HT signal of the 'left' selected experiment
            # The HT signal data points are interpolated (if required)
            # to match the wavelength of the new generated experiment
          htNew  <- approx(x = x1,y = htNew,xout = wlNew , method = "linear")$y

        }

      } else {

        wlNew                <- wlSelSubset1[[i]]

      }

      signalNew            <- y1 - y2

      newTemp              <- temperatureSel[[i]]

      signalNew            <- matrix(signalNew,ncol = 1)
      htNew                <- matrix(htNew,ncol = 1)

      # 5th step (Case 1) - Create one new experiment per selected experiment
      metadata_list <- list(
        'Processing'  = 'Generated by subtraction ',
        'Curve 1'     = spectraUneven[i],
        'Curve 2 (to be subtracted)'     = spectraEven[i],
        'Units for processing' = processingUnitsString,
        'Linear interpolation mode' = input$allowLinearInterpolation
      )

      append_cd_experiment(newSpectraNameBatch,newSpectraNameBatch,wlNew,
                           signalNew,operationUnits,metadata_list,newTemp,
                           htNew,isFake,operationUnits)

      newExperimentNames <- c(newExperimentNames,newSpectraNameBatch)
      legendDf           <- appendRowsToLegendDataFrame(legendDf,c(newSpectraNameBatch))

    }

  } else {
    
    signalNew  <- list()
    wlNew      <- list()
    htNew      <- list()

    if (operation == 'Smooth') {
      
      windowLength <- as.numeric(spectra2) 
      
      # smooth each of the spectra
      for (i in 1:length(exp1s)) {
        
        wlSel1     <- wlSel[[i]]
        signal1    <- signalSel[[i]]
        ht1        <- htSel[[i]]
        
        deltaWL      <- ( max(wlSel1) - min(wlSel1) ) / (length(wlSel1) - 1)
        odd_n_data_points_window_len <- ceiling(windowLength / deltaWL) %/% 2 * 2 + 1
        
        signalNew[[i]]  <- sgolayfilt(signal1, p = 2, n = odd_n_data_points_window_len)
        wlNew[[i]]      <- wlSel1
        htNew[[i]]      <- ht1
        
      }

    } else if (operation == 'Scale') {

      scalingFactor <- as.numeric(spectra2)

      # Apply a scaling factor to each of the spectra
      for (i in 1:length(exp1s)) {

        wlSel1     <- wlSel[[i]]
        signal1    <- signalSel[[i]]
        ht1        <- htSel[[i]]

        signalNew[[i]]  <- signal1 * scalingFactor
        wlNew[[i]]      <- wlSel1
        htNew[[i]]      <- ht1

      }

    } else if (operation == 'Zero') {
      
      wlRange <- as.numeric(spectra2) - 1

      # Apply an offset to each of the spectra. 
      # Mathematically, we subtract from each spectra the mean signal 
      # within a selected wavelength range
      for (i in 1:length(exp1s)) {
        
        wlSel1     <- wlSel[[i]]
        signal1    <- signalSel[[i]]
        ht1        <- htSel[[i]]
        
        wlSel_temp      <- wlSel1 >= (input$wavelengthRange[2] - wlRange) & wlSel1 <= input$wavelengthRange[2]
        
        offsetValue <- mean(signal1[wlSel_temp])

        signalNew[[i]]  <- signal1 - offsetValue
        
        wlNew[[i]]      <- wlSel1
        htNew[[i]]      <- ht1
        
      }
      
    } else if (operation %in% c('Subtract','Sum')) {

      # Retrieve the second spectra
      wlSel2      <- cdAnalyzer$experimentsOri[[last(names(ids))]]$wavelength
      signalSel2  <- cdAnalyzer$experimentsOri[[last(names(ids))]]$signalDesiredUnit[,last(ids)]
    
      for (i in 1:length(exp1s)) {
        
        wlSel1     <- wlSel[[i]]
        signal1    <- signalSel[[i]]
        ht1        <- htSel[[i]]
        
        if (input$allowLinearInterpolation) {
         
          interpolatedSignals <- interpolateVectors(wlSel1,signal1,wlSel2,signalSel2)
          y1                  <- interpolatedSignals$y1
          y2                  <- interpolatedSignals$y2
          wlNew[[i]]          <- interpolatedSignals$x1
          
          if (all(is.na(ht1))) {
            
            htNew[[i]] <- rep(NA,length(interpolatedSignals$x1))
            
          } else {
            
            # We do not sum the HT signals!
            # We return the HT signal of the 'left' selected experiment
            # The HT signal data points are interpolated (if required)
            # to match the wavelength of the new generated experiment
            htNew[[i]]      <- approx(x = wlSel1, 
                                      y = ht1, 
                                      xout = interpolatedSignals$x1 , method = "linear")$y
            
          }
          
        } else {
          
          y1                  <- signal1
          y2                  <- signalSel2
          wlNew[[i]]          <- wlSel1
          htNew[[i]]          <- ht1
        }
        
        if (operation == 'Subtract') {
          signalNew[[i]] <- y1 - y2
        } else {
          signalNew[[i]] <- y1 + y2
        }
      
      }
    } 
    
    newExperimentNames <- c()
    
    for (exp in unique(exp1s)) {
      
      correspondingIDs <- which(exp1s == exp)
      
      # All spectra from the same experiment should share the same wavelength data points
      wlExperiment     <- wlNew[[correspondingIDs[1]]]
      
      c1 <- length(unique(exp1s))    == 1
      c2 <- length(correspondingIDs) == 1
      
      # Assign the new experiment name
      if (c1) {
        newExperimentName <- newSpectraName
      } else {
        newExperimentName <- paste0(exp,' ',newSpectraName)
        
        if (newExperimentName %in% names) {
          newExperimentName <- paste0(newExperimentName,' ',get_hour_minute_sec())
        }
      }
      
      # Select the signal and assign spectra names 
      if (c2) {
        signalExperiment  <- matrix(signalNew[[correspondingIDs]],ncol = 1)
        htExperiment      <- matrix(htNew[[correspondingIDs]],ncol = 1)
        newSpectraNames   <- c(newExperimentName)
      } else {
        signalExperiment  <- as.matrix(do.call(cbind,signalNew[correspondingIDs]))
        htExperiment      <- as.matrix(do.call(cbind,htNew[correspondingIDs]))
        selSpectraNames   <- cdAnalyzer$experimentsOri[[exp]]$spectraNames[id1s[correspondingIDs]]
        newSpectraNames   <- paste0(selSpectraNames,' ',newSpectraName)
      }
     
      # Assign the temperature from the spectrum/a #1
      
      expTemperature    <- unlist(temperatureSel[correspondingIDs])

      # 5th step (Case 2) - Create one new experiment per selected experiment
      metadata_list <- list(
        'Processing'           = paste0('Generated by applying the ',operation,' operation'),
        'Original CD curve(s)' = paste0(cdAnalyzer$experimentsOri[[exp]]$spectraNames[id1s[correspondingIDs]],
                                        collapse = ' '),
        'Units for processing' = processingUnitsString
      )


      # Assign the temperature from the spectrum/a #1

      newExperimentNames <- c(newExperimentNames,newExperimentName)
      
      legendDf <- appendRowsToLegendDataFrame(legendDf,newSpectraNames)
      
      # Add the parameters for the smoothing
      if (operation == 'Smooth') {
        
        metadata_list[['Savitzky-Golay smoothing filter window']] = paste0(windowLength,' nm') 
        metadata_list[['Savitzky-Golay smoothing filter order']]  = '2'
        
        deltaWL      <- ( max(wlExperiment) - min(wlExperiment) ) / (length(wlExperiment) - 1)

        metadata_list[['Savitzky-Golay smoothing filter length']] = paste0(ceiling(windowLength / deltaWL) %/% 2 * 2 + 1)

      }
      
      # Add the parameters for the sum/subtract
      if (operation %in% c('Sum','Subtract')) {

        operationWord <- ifelse(operation == 'Subtract','subtraction','addition')

        second_curve_key <- paste0('Second CD curve (for ',operationWord,')')

        metadata_list[[second_curve_key]] = spectra2

        if (input$trackPreviousChange) {

          # Include the original CD curves of the right selected spectrum
          previous_exp2_metadata <- cdAnalyzer$experimentsOri[[last(names(ids))]]$metadata

          for (md in previous_exp2_metadata) {

            if (md == "Generated by averaging many CD curves") {

              metadata_list[[second_curve_key]] <- paste0(metadata_list[[second_curve_key]],
                                                          ' (average from: ',
                                                          previous_exp2_metadata[['Original CD curve(s)']],
                                                          ')')
              break
            }

            for (operationPrev in c('Smooth','Zero','Scale')) {

              if (grepl(paste0(operationPrev,' operation'),md)) {

                metadata_list[[second_curve_key]] <- paste0(metadata_list[[second_curve_key]],
                                                            ' (',operationPrev,' of: ',
                                                            previous_exp2_metadata[['Original CD curve(s)']],
                                                            ')')
                break
              }
            }
          }
        }
      }
      
      # Add the parameters for the smoothing
      if (operation == 'Zero') {
        
        wlRangeStr <- paste0(input$wavelengthRange[2] - wlRange,' - ',input$wavelengthRange[2], ' nm')
        metadata_list[['Wavelength range to calculate offset']] = wlRangeStr
        
      }

      # Add the parameters for the scaling
      if (operation == 'Scale') {

        metadata_list[['Scaling factor']] = paste0(spectra2)

      }

      if (input$trackPreviousChange) {

        # Include the original CD curves of the left selected spectrum
        previous_exp_metadata <- cdAnalyzer$experimentsOri[[exp]]$metadata

        for (md in previous_exp_metadata) {

          if (md == "Generated by averaging many CD curves") {

            metadata_list[['Original CD curve(s)']] <- paste0(metadata_list[['Original CD curve(s)']],
                                                            ' (average from: ',
                                                            previous_exp_metadata[['Original CD curve(s)']],
                                                            ')')
            break
          }

          for (operationPrev in c('Smooth','Zero','Scale')) {

            if (grepl(paste0(operationPrev,' operation'),md)) {

            metadata_list[['Original CD curve(s)']] <- paste0(metadata_list[['Original CD curve(s)']],
                                                              '. Origin: ',operationPrev,' of: ',
                                                              previous_exp_metadata[['Original CD curve(s)']]
                                                              )
              break
            }
          }

        }
      }

      # Include the right-selected spectrum metadata for the smooth zero and scale operations
      if (operation %in% c('Smooth','Zero','Scale') && input$trackPreviousChange) {

        for (md in previous_exp_metadata) {

          if (md %in% c("Generated by applying the Subtract operation",
                        "Generated by applying the Sum operation")) {

            operationPrev <- strsplit(md,' ')[[1]][5]

            operationWordPrev <- ifelse(operationPrev == 'Subtract','subtraction','addition')

            operationPrev <- ifelse(operationPrev == 'Subtract','minus','plus')

            second_curve_key <- paste0('Second CD curve (for ',operationWordPrev,')')

            metadata_list[['Original CD curve(s)']] <- paste0(metadata_list[['Original CD curve(s)']],
                                                                    '. Origin: ',
                                                                    previous_exp_metadata[['Original CD curve(s)']],
                                                                    ' ',operationPrev,' ',
                                                                    previous_exp_metadata[[second_curve_key]])

          }

        }
      }

      append_cd_experiment(newExperimentName,newSpectraNames,np_array(wlExperiment),
                           signalExperiment,operationUnits,metadata_list,
                           expTemperature,htExperiment,isFake,operationUnits)
      
    }
    
  } 
  
  ## If shared parameters, assign known values of MW, #AA, ...
  if (cdAnalyzer$sharedParameters & !isFake) {
    
    # Use first experiment as proxy for the parameters
    mol_weight  <- cdAnalyzer$experimentsOri[[1]]$molecularWeight
    path_length <- cdAnalyzer$experimentsOri[[1]]$pathLength
    cro_num     <- cdAnalyzer$experimentsOri[[1]]$numberOfCroms
    conc        <- cdAnalyzer$experimentsOri[[1]]$concentration
    
    # Assign the same values to the generated experiments
    for (exp in newExperimentNames) {
      
      cdAnalyzer$experimentsOri[[exp]]$molecularWeight  <- mol_weight
      cdAnalyzer$experimentsOri[[exp]]$pathLength       <- path_length
      cdAnalyzer$experimentsOri[[exp]]$numberOfCroms    <- cro_num
      cdAnalyzer$experimentsOri[[exp]]$concentration    <- conc
      
    }
    
  } else {
    
    # Assign the same value as the experiments used to generate the new CD curves
    if (operation %in% c('Smooth','Subtract','Sum','Zero','Scale')) {
      
      i <- 0
      for (expNew in newExperimentNames) {
        i <- i + 1
        
        mol_weight  <- cdAnalyzer$experimentsOri[[unique(exp1s)[i]]]$molecularWeight
        path_length <- cdAnalyzer$experimentsOri[[unique(exp1s)[i]]]$pathLength
        cro_num     <- cdAnalyzer$experimentsOri[[unique(exp1s)[i]]]$numberOfCroms
        conc        <- cdAnalyzer$experimentsOri[[unique(exp1s)[i]]]$concentration
        
        cdAnalyzer$experimentsOri[[expNew]]$molecularWeight  <- mol_weight
        cdAnalyzer$experimentsOri[[expNew]]$pathLength       <- path_length
        cdAnalyzer$experimentsOri[[expNew]]$numberOfCroms    <- cro_num
        cdAnalyzer$experimentsOri[[expNew]]$concentration    <- conc
        
      }
      
    }
    
    # Modify the input parameters Table
    updateCDFilesInfoTable()
  }

  # 6th step - convert back to the appropriate units
  rList2        <- as.list(rep(input$workingUnits,length(unique(names(ids)))+length(newExperimentNames)))
  names(rList2) <- c(unique(names(ids)),newExperimentNames)
  
  cdAnalyzer$experiments_absorbance_units_to_other_units(rList2)
  
  # Filter the new experiments
  wlRange <- input$wavelengthRange
  cdAnalyzer$filter_data_by_wavelength(wlRange[1], wlRange[2])
  
  updateSelectInput(session,"experiment2delete",selected = "ALL",
                    choices = c('ALL',cdAnalyzer$experimentNames))

  updateSelectInput(session,"selected_cd_exp",choices     = cdAnalyzer$experimentNames)
  
  output$legendInfo           <- helperRenderRHandsontable(legendDf)
  updateSESCA_ref(legendDf)

  updateProcessingTable(operation,operationUnits)
  
  # We need to use lapply instead of a for-loop to avoid R lazy evaluation!
  lapply(newExperimentNames, function (exp) {
    
    saneExpName <- gsub(':','',exp)
    
    tabP <- tabPanel(saneExpName,value=saneExpName,fluidRow(column(12,tableOutput(paste0('metadata_',saneExpName)))))
    
    appendTab("metadata",tabP,select=TRUE)
    
    metadata_info     <- cdAnalyzer$experimentsOri[[exp]]$metadata
    metadataFeature   <- names(metadata_info)
    metadataValue     <- unlist(metadata_info)
    
    metadata_df <- data.frame(metadataFeature,metadataValue)
    colnames(metadata_df) <- c('Metadata feature','Value') 
    
    output[[paste0('metadata_',saneExpName)]] <- renderTable({metadata_df})
    return(NULL)
    
  })
  
  Sys.sleep(2)
  reactives$data_loaded <- TRUE
  
})

# To update which the selectInput options
# For example, if we select 'Average' as the operation, there is no need to allow selecting
# a second spectrum
observeEvent(input$operation,{
  
  req(input$operation)
  updateProcessingTable(input$operation,input$operationUnits)
  
})

# Update the slider range based on the maximum accepted voltage value
# For this filter to work, the HT values of the experiments should have similar ranges
# E.g., it won't work if some CD experiments have HT values ranging from 400 to 600
# and other CD experiments have HT values from 2 to 6
observeEvent(input$maxHTvalue,{
  
  req(reactives$data_loaded)
  req(is.numeric(input$maxHTvalue))

  minWL <- cdAnalyzer$find_wavelength_limit(input$maxHTvalue)

  if (minWL != 0) {

    updateSliderInput(session,'wavelengthRange',NULL,value = c(minWL,input$wavelengthRange[2]))

  }

})

