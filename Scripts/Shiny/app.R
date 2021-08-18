library(shiny)
library(tidyverse)
library(ggmap)
library(viridis)
library(lubridate)
library(eDNAfuns)
library(shinycssloaders)
library(gifski)
library(here)

map.new.samples <-  read_rds(here("Data", "base_plot_eDNA.rds"))

Ab.table       <- read_csv(here("Data", "Abundance_table_BLAST.csv"))

Ab.table %>% 
    separate (Sample, into = c(NA, "Sample", NA)) %>% 
    group_by(Sample, Taxa) %>% 
    summarise(prevalence = n(), 
              lat = mean(lat),
              lon = mean (lon),
              Cruise = first(Cruise),
              Transect = mean (Transect)) -> Mid.step

Ab.table %>% 
    
    separate (Sample, into = c(NA, "Sample", NA))%>% 
eDNAindex(Sample_column = Sample, OTU_column = Taxa, Counts_column = nReads, Biological.replicate = Replicate ) %>% 
    left_join(Mid.step) -> Ab.table

Ab.table %>% 
    ungroup %>% 
    summarise_at(.vars = vars(matches("lon|lat")),
                 .funs = list(max = max, min = min)) -> limits
# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Geographical location of taxa in the Columbia River Plume"),
    
    fluidRow(
        column(5,
               selectInput(inputId =  "taxa",
                           label = "Taxa:", 
                           choices = Ab.table %>% distinct(Taxa),
                           # selected = HABs %>% distinct(label) %>% slice(1) %>% pull(),
                           multiple = T, width = '50%'),
               actionButton("make_plot", "Create Plot")
        )),
    fluidRow(
        column(9,
               
               plotOutput("distPlot",width = '100%') %>%  withSpinner(color="red", type = 6)))
        )
    

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    #TODO: start with a still image of the map
    output$distPlot <- renderPlot({
        
        map.new.samples +
            # geom_point(data = to.plot, aes(x = lon, y = lat, size = prevalence , color = Normalized.reads, group = grouping)) +
            scale_size(range = c(2,10)) 
        
    }, height = 800, width = 800)
    
    
    observeEvent(input$make_plot,{
        
        output$distPlot <- renderPlot({
            
            
            # generate bins based on input$bins from ui.R
            Ab.table %>% 
                filter(Taxa %in% input$taxa) %>% 
                rownames_to_column(var = "grouping") -> to.plot
            
            # start with a bar plot
            # ggplot (to.plot, aes(x = SiteDate, y = Normalized.reads, fill= taxa)) +
            #     geom_col()
            # Now the animation
            map.new.samples +
                geom_point(data = to.plot, aes(x = lon, y = lat, size = prevalence , color = Normalized.reads, group = grouping)) +
                scale_size(range = c(2,10)) +
                scale_color_viridis()-> p # The base plot
            
            p +
                coord_sf(xlim = c(-123, -125),
                         ylim = c(46, 47)) +
                facet_wrap(~Taxa, nrow = 2, ncol = 2) +
                theme(strip.text = element_text( size = 12, face = "bold"))#+
                # labs(title = 'Date: {closest_state}') +
                # transition_states(fct_reorder(paste(month(date,label = T), year(date), sep = "'"),date),
                #                   transition_length = 2,
                #                   state_length = 1) +
                # enter_fade()+
                # exit_shrink() -> anim
            
            # Wrap it with RenderImage?
            # anim_save(filename = paste0(input$taxa[1],"_ani.gif"),
            #           animation = anim,
            #           nframes = 50,
            #           fps = 5, 
            #           width = 800,
            #           height = 800)
            # 
            # output$distPlot <- renderImage({
            #     list(src = paste0(input$taxa[1],"_ani.gif"),
            #          contentType = 'image/gif')
            # }, deleteFile = FALSE)
        }, height = 1200, width = 1200) # size is ok but resolution is poor - need to change resolution on anim_save
        
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)