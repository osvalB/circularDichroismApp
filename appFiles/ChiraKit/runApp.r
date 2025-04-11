options(browser='google-chrome-stable')
# Detect user name
user_name <- Sys.getenv("USER")
# Set the base directory for the app
base_dir <- paste0('/home/', user_name, '/circularDichroismApp/appFiles/ChiraKit/')
shiny::runApp(base_dir)
