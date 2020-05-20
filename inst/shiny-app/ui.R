# UI for CMV resistance database
# Author - Charles, OJ
# Date - April 2020
################
# design idea https://community.rstudio.com/t/background-images-in-shiny/12261/3
library(shiny)


# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    # Application title
    titlePanel("Human Cytomegalovirus Drug Resistance Genotyping"),
    
    navlistPanel(widths = c(2,6),
       "Resistance Genotyping",
       tabPanel("Upload & Clinical",
                h3("Upload a file to return resistance data"),
                fileInput("vcf.file", label = "Accepts VCF, Varscan-tab, Fasta formats",multiple = F),
                div(style="display:inline-block",downloadButton("vcf.o.res", "download resistance data")),
                div(style="display:inline-block",downloadButton("vcf.o.all", "download all-mutants data")),
                br(),br(),h4("Clinical Overview"),
                DT::dataTableOutput("vcf.table_clin"),
                br(),br(),h4("Comprehensive resistance  - Values are EC50 fold change ratio to wild type virus"),
                tableOutput("vcf.table_res")
                #plotOutput("vcf.plot.res")

       ), #tabpanel
       

       tabPanel("Genetics",
                textOutput("vcf.title.UL54"),
                plotOutput("vcf.plot.lollipop.UL54"),
                textOutput("vcf.title.UL97"),
                plotOutput("vcf.plot.lollipop.UL97"),
                textOutput("vcf.title.UL89"),
                plotOutput("vcf.plot.lollipop.UL89"),
                plotOutput("vcf.plot.lollipop.UL56"),
                plotOutput("vcf.plot.lollipop.UL51"),
                plotOutput("vcf.plot.lollipop.UL27")
       ), #tabpanel
       
       
       "Database info",
       # tabPanel("Database-entries",
       #  # todo: make some fields mandatory? using shinyjs https://deanattali.com/2015/06/14/mimicking-google-form-shiny/
       #  # todo: handle submission, requires a data arhchitecture first.
       #   sidebarLayout(
       #     sidebarPanel(
       #       h3("Download resDB"),
       #       selectInput("db.db-version", "select db version",
       #                   c("v1.0", "v1.1-beta")),
       #       downloadButton("db.download", "Download Resistance CSV"),
       #       "",
       #       h3("Issues?"),
       #       textInput("db.gene", "Gene", placeholder = "e.g. UL54"),
       #       textInput("db.aamut", "Amino Acid Location & mutation",placeholder = "e.g. A571G"),
       #       selectInput("db.drugres", "Drug Resistance to assign / update",
       #                   c("", "GCV","CDV", "FOSarnate", "BCV", "LENtermovir", "TOMoglovir", "MARibavir", "other-in note"),multiple = T),
       #       textAreaInput("db.note", label = "Include in this note, what you would like to see updated",
       #                     placeholder = "e.g. UL54 Mutation A123S is assigned as resistant to GCV, it is missed as conferring resistance to FOS"),
       #       textInput("db.ref", "Online reference ",placeholder = "e.g. https://www.biorxiv.org/content/link/to/our/paper"),
       #       p(), #separete intput from output
       #       checkboxInput("db.contact", "I wish to be contacted abut this change", FALSE),
       #       textInput("db.contact_email", "email address for contact"),
       #       actionButton("db.submit", "Submit", class = "btn-primary")
       #     ),
       #     mainPanel(
       #       h3("Submit any noted issues with the CMV resistance database i.e. any missed mutations in literature,
       #          additional drugs a gene mutation confers resitance to"),
       #       h4("here we will show any current issues, noted conflicts in the literature for an AAmut")
       #     ) #main
       #   ) #sbl
       # ), #tabpanel
       
       
       tabPanel("Database-metrics",
                h2("gene - drug resistance heatmap, showing number of entries for each combination"),
                plotOutput("res.plot.dbheatmap")
                
       ), #tabpanel
       
       tabPanel("Terms of Use",
                p(""),
                br(),
                strong("Disclaimer"),
                div("This database has been created and is maintained by the Breuer Lab at UCL as a benefit to the research and education community. This tool is provided on an 'as is, best endeavours' basis only, and without guarantee to its accuracy or reliability.", tags$br(),
                  "We clearly state here that the tool is not validated for use in a diagnostic or clinical contexts.", tags$br(),
                  "We clearly state here that any file uploaded should be devoid of patient identifiable features, such as filename.", tags$br(),
                  "When you upload a file to use this web service, you are accepting these terms."),
                br(),
                
                strong("Privacy Policy"),
                p("Every time you access the hosted 'Human Cytomegalovirus Drug Resistance Genotyping' hosted UCL web address  we reserve the right to register the following data for statistical and troubleshooting purposes."),
                HTML("<ul><li>Date and time of access</li><li>IP address of request</li><li>Uploaded file metadata(size, file type). File contents are untouched.</li></ul>"),
                br(),
                br(),
                br(),
                br(),
                p("last updated: 07/05/2020")
                
                
                
                
       ), #tabpanel
       tabPanel("About",
                p(""),
                br(),
                
                strong("Data"),
                p("To capture new resistant mutations, articles in peer-reviewed journals are searched weekly using a PubMed search with key terms.
                  Regular expressions narrow down the list to those including mutations. The way the application functions, we can capture multiple data points per mutation, so as to minimise subjectivity. If you spot any mistakes or omissions in the database either contact oscar.charles.18@ucl.ac.uk or alter the github /inst/db files. "), 
                br(),
                
                strong("Method:"),
                p("Results generated are relative to the well characterised strain Merlin in GenBank.
When loading a fasta file, 'MAFFT --add --keeplength'  is used to create an alignment to MERLIN strain. Therefore, whole genome sequences as well as fasta files containing genetic fragments are processed seamlessly.
Alignments are converted to VCF format using snp-sites.\n
                                     Custom R functions then call resistance from the database."),
                br(),
                
                strong("Open Source"),
                p(HTML(paste0("The database, functions and application are available as an R package ",a("in Github here.", href ="https://github.com/ucl-pathgenomics/cmvdrg")))),
                br(),br(),
                
                strong("Ackowledgements"),
                p(HTML(paste0("We would like to thank the MRC and Wellcome Trust as sponsors of the Breuer lab, for enabling this work.\n
                  The Breuer Lab work in close co-operation with the ",a(href = 'https://www.ucl.ac.uk/infection-immunity/pathogen-genomics-unit', 'Pathogen Genomics Unit'), " "))),
                fluidRow( 
                  column(3, img(width = "90%", src = "ucl.png")),
                  column(3, img(width = "90%", src = "mrc.png")),
                  column(3, img(width = "90%", src = "wt.jpg")),
                  column(3, img(width = "90%", src = "923-pgu-logo.jpg"))
                )
                # forced, not nice
                # fluidRow( 
                #   column(4, img(height = 100, width = 300, src = "ucl.png")),
                #   column(4, img(height = 100, width = 100, src = "wt.jpg")),
                #   column(4, img(height = 100, width = 150, src = "923-pgu-logo.jpg"))
                # )
                #https://www.ucl.ac.uk/infection-immunity/pathogen-genomics-unit   add somehting to do with this
                
       ) #tabpanel
    ) #navlistpanel
  ) #fluidpage
) #shinyui
