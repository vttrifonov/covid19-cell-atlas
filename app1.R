library(data.table)
library(magrittr)
library(shiny)
library(ggplot2)
source('helpers.R')

ui <- {
  fluidPage(
    fluidRow(
        selectizeInput(
            'genes',
            'Select cytokines',
            choices=NULL,
            multiple=TRUE
        )
    ),
    fluidRow(
        column(7,
            uiOutput('plot1.ui')
        )
    ),
    fluidRow(
        selectizeInput(
            'genes1',
            'Select gene',
            choices=NULL,
            multiple=FALSE
        )
    ),
    fluidRow(
        column(8,
            plotOutput('plot2', height=200*2)
        )
    ),
    fluidRow(
        column(8,
            plotOutput('plot3', height=200*24)
        )
    )
  )
}

server <- function(input, output, session) {
    {
        if (FALSE) {
            input <- list(genes=c('IL-6'), genes1=c('IL6'))
        }

        analysis <- reticulate::import('covid19_cell_atlas._analysis3')$analysis3

        plot1 <- reactive({
            genes <- req(input$genes)
            data <- analysis$app1_plot1_data(genes) %>% data.table

            ggplot(data)+
                aes(x=days_since_onset, y=log2(level))+
                geom_line(aes(group=donor), alpha=0.1)+
                geom_point(aes(fill=DSM_group), alpha=0.5, size=2, pch=21)+
                geom_smooth(aes(color=DSM_group), alpha=0.5)+
                facet_grid(cytokine~., scales='free')+
                labs(y='log2(pg/mL)')
        })

        plot2 <- reactive({
            genes <- req(input$genes1)
            data <- analysis$app1_plot2_data(genes) %>% data.table

            ggplot(data)+
                aes(days_since_onset, y=log2(X+1))+
                geom_point(aes(fill=dsm_severity_score_group), alpha=0.5, size=2, pch=21)+
                geom_line(aes(group=donor), alpha=0.1)+
                geom_smooth(aes(color=dsm_severity_score_group), alpha=0.5)+
                facet_grid(subset+gene~., scales='free')+
                labs(y='log2RPM')
        })

        plot3 <- reactive({
            genes <- req(input$genes1)
            data <- analysis$app1_plot3_data(genes) %>% data.table

            ggplot(data)+
                aes(days_since_onset, y=log(X+1))+
                geom_point(aes(fill=dsm_severity_score_group), alpha=0.5)+
                geom_line(aes(group=donor), alpha=0.1)+
                geom_smooth(aes(color=dsm_severity_score_group), alpha=0.5)+
                facet_grid(subset+cell_type+gene~., scales='free_x')+
                labs(y='log2RPM')
        })
    }

    {
        updateSelectizeInput(
            session, 'genes',
            choices = analysis$app1_genes,
            server = TRUE
        )

        output$plot1.ui <- renderUI({
            genes <- req(input$genes)
            plotOutput('plot1', height=200*length(genes))
        })

        output$plot1 <- renderPlot(plot1())

        updateSelectizeInput(
            session, 'genes1',
            choices = analysis$app1_genes1,
            server = TRUE
        )

        output$plot2 <- renderPlot(plot2())
        output$plot3 <- renderPlot(plot3())
    }
}

app <- shinyApp(
    ui=ui,
    server=server
)

runApp(app)
