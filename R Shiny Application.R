###################################### Development of Basic RShiny Application

library(shiny)
library(ggplot2)
library(dplyr)
library(DT)

# Define UI for the app
ui <- fluidPage(
  titlePanel("Pharmacokinetic (PK) Parameter Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = ".csv"),
      selectInput(
        "parameter", 
        "Select PK Parameter for Summary Statistics:", 
        choices = c("Cmax", "AUCt", "AUCi"),
        selected = "Cmax"
      ),
      actionButton("analyze", "Analyze")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Summary Statistics", DTOutput("summary_table")),
        tabPanel("Visualization", plotOutput("ratio_plot"))
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Reactive function to read uploaded data
  uploaded_data <- reactive({
    req(input$file)
    read.csv(input$file$datapath)
  })
  
  # Reactive function to calculate summary statistics
  summary_stats <- reactive({
    req(input$parameter, uploaded_data())
    
    data <- uploaded_data()
    data %>%
      summarize(
        Min = min(.data[[input$parameter]], na.rm = TRUE),
        Max = max(.data[[input$parameter]], na.rm = TRUE),
        Mean = mean(.data[[input$parameter]], na.rm = TRUE),
        Median = median(.data[[input$parameter]], na.rm = TRUE),
        SD = sd(.data[[input$parameter]], na.rm = TRUE),
        CV = (sd(.data[[input$parameter]], na.rm = TRUE) / mean(.data[[input$parameter]], na.rm = TRUE)) * 100
      )
  })
  
  # Reactive function to calculate Test/Reference ratios
  ratio_data <- reactive({
    req(uploaded_data())
    
    data <- uploaded_data()
    
    # Calculate ratios for each subject
    ratios <- data %>%
      group_by(Sub) %>%
      summarize(
        Cmax_Ratio = Cmax[Trmt == "T"] / Cmax[Trmt == "R"],
        AUCt_Ratio = AUCt[Trmt == "T"] / AUCt[Trmt == "R"],
        AUCi_Ratio = AUCi[Trmt == "T"] / AUCi[Trmt == "R"]
      ) %>%
      pivot_longer(cols = starts_with("Cmax_Ratio"), names_to = "Parameter", values_to = "Ratio")
    
    ratios
  })
  
  # Render summary statistics table
  output$summary_table <- renderDT({
    req(summary_stats())
    datatable(summary_stats(), options = list(pageLength = 5))
  })
  
  # Render ratio plot
  output$ratio_plot <- renderPlot({
    req(ratio_data())
    
    ggplot(ratio_data(), aes(x = Sub, y = Ratio, color = Parameter)) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 1.25, linetype = "dashed", color = "blue") +
      labs(
        title = "Test/Reference Ratios with Reference Lines",
        x = "Subject",
        y = "Ratio"
      ) +
      theme_minimal()
  })
}

# Run the application
shinyApp(ui = ui, server = server)

# Sample data to run the application
df = data.frame(Sub = result$Sub, Trmt = result$Trmt, Cmax = result$Cmax, AUCt = result$AUCt, AUCi = result$AUCi)
write.csv(df, "df.csv", row.names = FALSE)
