
# # #debug
# inFile = list()
# inFile$datapath =  system.file("testdata", "A10.vcf", package = "cmvdrg")
# f.infile = as.character(inFile$datapath)
# session = list()
# session$token = "test"
# # end debug


shinyServer(function(input, output, session) {
  global = list()
  global$res_table = system.file("db", "cmvdrg-db1.csv", package = "cmvdrg")
  #create unique session folder
  global$date <- format(Sys.time(), "%Y-%m-%d")
  global$dir = ""
  global$genome = genome="NC_006273.2"
  global$path_gff3_file=system.file("ref", "NC_006273.2.gff3", package = "cmvdrg")
  global$path_fasta_file=system.file("ref", "NC_006273.2.fasta", package = "cmvdrg")
  global$path_txdb=system.file("ref", "NC_006273.2e.sqlite", package = "cmvdrg")

  
  
  
  ###************  processes  ****************###  
  
  
  vcf.d.all <- reactive({
    # Returns annotated variants dataframe
    # in  - any of the accepted file formats from file load
    # out - all mutants table
    
    inFile <- input$vcf.file
    if (is.null(inFile)){
      return(NULL)}
    
    
    # handle input .tab .vcf .fasta
    # returns minimum intermediate table
    dat1 <- read_input(inFile$datapath, global)
    
    ### annotate variants
    # used in optional processes, i.e. identifying syn / nonsynonymous mutations in resistance genes.
    dat2 <- annotate_variants(f.dat = dat1, global)
    
    
    return(dat2)
  })
  
  vcf.d.res <- reactive({
    # This is the main table used in the shinyapp, dataframe of resistant mutants.e ssentially call_resistance(..., all_mutations = F)
    # out - only resistant data. used by most plots.

    inFile <- input$vcf.file
    if (is.null(inFile)){
      return(NULL)}
    
    dat2 = vcf.d.all()

    ### add res info
    dat3 <- add_resistance_info(f.dat = dat2, resistance_table=global$res_table, all_muts = F)
    
    return(dat3)
  })
  

  
  
  
  
  ########################################### Outputs#######################################  
  output$vcf.table_clin <- DT::renderDataTable({
    if(is.null(vcf.d.res())){
      return(NULL)}
    out = make_clin_table(vcf.d.res())  
    return(out)
  })
  
  
  
  output$vcf.table_res <- renderTable({
    # resistance mutation table
    if(is.null(vcf.d.res())){
      return(NULL)
    }else if (nrow(vcf.d.res()) == 0){
      dat <- data.frame(notes = "there was no resistance identified in your sample")
      return
    }else{
    dat <- vcf.d.res()
    #full data is saved and user can DL
    #we really only need to present some of these columns
    dat = dat[,c(28,9,10,8,26,27,29,30,31,32:40,42,43,45,46,47,48)]
    dat[,19] =     paste0("<a href='",  dat[,20], "' target='_blank'>",substr(dat[,19],1,20),"</a>")
    dat[,20] = NULL
    #dat <- reshape2::dcast(dat, GENE + AA_CHANGE + freq ~ variable)
    }
    return(dat[,])
  },sanitize.text.function = function(x) x)
  

  output$res.plot.dbheatmap <- renderPlot({
    # plots as a heatmap, the amount of dat we have per gene, per drug.
    # //todo - currently records the number of data points, could record mean value etc
    resistance = read.csv(global$res_table, header = TRUE,na.strings = c("", "NA"), stringsAsFactors = F)
    resistance = melt(resistance, measure.vars = colnames(resistance[,6:14]))
    resistance = resistance[resistance$value != "",]
    resistance = resistance[!is.na(resistance$value),]
    resistance$value = str_replace(resistance$value, ">", "")
    resistance$value = str_replace(resistance$value, ",.{1,5}", "")
    resistance = resistance[resistance$value > 1,]
    
    resistance = reshape2::dcast(resistance, GENE ~ variable,fun.aggregate = length)
    resistance = reshape2::melt(resistance, id.vars = "GENE")
    resistance$Number_of_entries= resistance$value
    resistance$Drug = resistance$variable
    g = ggplot(resistance, aes(x = Drug, y = GENE, fill = Number_of_entries)) +
      geom_tile() +
      theme_classic() +
      scale_fill_gradient(low="white", high="red") +
      theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
    return(g)
  })
  

  output$vcf.plot.res <- renderPlot({
    validate(
     need(input$vcf.file != "", "Please load a file")
    )
    #if (is.null(input$vcf.file)){
    #  return(NULL)}
    if (nrow(vcf.d.res()) == 0){ # if no resistance data was identified
      return(NULL)
    }
    # plot resistance muts
    coding_df_res_resistance <- vcf.d.res()
    coding_df_res_resistance$freq <- readr::parse_number(coding_df_res_resistance$freq)
    
    if(is.null(coding_df_res_resistance$time_point)){
      g <- ggplot( coding_df_res_resistance ,aes(x = "mutations", y=freq,fill=factor(change)))
    }else{
      g <- ggplot( coding_df_res_resistance ,aes(x=time_point,y=freq,fill=factor(change)))
    }
    g <- g +
      geom_bar(position="dodge",stat="identity") +
      scale_fill_discrete(drop=FALSE) +
      scale_x_discrete(drop=FALSE) +
      theme_bw() +
      labs(fill = "Resistance")
    # ggsave(plot = g, filename = paste("session",global$date, session$token, "res_table.png", sep = "/"), device = "png")
    return(g)
  })
  
  
  ### lollipops
  
  # ## ul54 ##
  output$vcf.title.UL54 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(grep(unique(mut$Gene), pattern = "UL54", session)) > 0){
      paste("Resistance mutations in UL54 Gene")
    }
  })
    output$vcf.plot.lollipop.UL54 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- plot_lollipop(mut, f.gene = "UL54",global = global)
    #ggsave(plot = g, filename = paste("session/",global$date, "/",session$token, "/UL54-", "lollipop.png", sep = ""), device = "png")
    return(g)
  })
  ## ul97 ##
    output$vcf.title.UL97 <- renderText({
      mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
      if(length(grep(unique(mut$Gene), pattern = "UL97")) > 0){
        paste("Resistance mutations in UL97 Gene")
      }
    })
  output$vcf.plot.lollipop.UL97 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- plot_lollipop(mut, f.gene = "UL97",global = global)
    #ggsave(plot = g, filename = paste("session/",global$date, "/",session$token, "/UL97-", "lollipop.png", sep = ""), device = "png")
    return(g)
  })
  ## ul89 ##
  output$vcf.title.UL89 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(grep(unique(mut$GENE), pattern = "UL89")) > 0){
      paste("Resistance mutations in UL89 Gene")
    }
  })
  output$vcf.plot.lollipop.UL89 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- plot_lollipop(mut, f.gene = "UL89", global = global)
    #ggsave(plot = g, filename = paste("session/",global$date, "/",session$token, "/UL89-", "lollipop.png", sep = ""), device = "png")
    return(g)
  })

  # ## ul56 ##
  output$vcf.title.UL56 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(grep(unique(mut$Gene), pattern = "UL56")) > 0){
      paste("Resistance mutations in UL56 Gene")
    }
  })
  output$vcf.plot.lollipop.UL56 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- plot_lollipop(mut, f.gene = "UL56",global = global)
    #ggsave(plot = g, filename = paste("session/",global$date, "/",session$token, "/UL56-", "lollipop.png", sep = ""), device = "png")
    return(g)
  })
  #
  #
  # ## ul51 ##
  output$vcf.title.UL51 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(grep(unique(mut$Gene), pattern = "UL51")) > 0){
      paste("Resistance mutations in UL51 Gene")
    }
  })
  output$vcf.plot.lollipop.UL51 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- plot_lollipop(mut, f.gene = "UL51",global = global)
    #ggsave(plot = g, filename = paste("session/",global$date, "/",session$token, "/UL51-", "lollipop.png", sep = ""), device = "png")
    return(g)
  })
  #
  #
  # ## ul27 ##
  output$vcf.title.UL27 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(grep(unique(mut$Gene), pattern = "UL27")) > 0){
      paste("Resistance mutations in UL27 Gene")
    }
  })
  output$vcf.plot.lollipop.UL27 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- plot_lollipop(mut, f.gene = "UL27",global = global)
    #ggsave(plot = g, filename = paste("session/",global$date, "/",session$token, "/UL27-", "lollipop.png", sep = ""), device = "png")
    return(g)
  })

  
  
  
  
  ###************  Outputs  ****************###  
  output$vcf.ip <- renderText({
    print(session$request$REMOTE_ADDR)
  })
  
  # debug in brownser app mode
  output$vcf.o.res <- downloadHandler(
    filename = function(){paste(global$date, "_resmuts.csv", sep="")},
    content = function(filename){
      dat <- vcf.d.res()
      write.csv(x = dat, file = filename, row.names = F)
    }
  )
  
  # all mutations, synonymous and non synonymous
  output$vcf.o.all <- downloadHandler(
    filename = function(){paste(global$date, "_allmuts.csv", sep="")},
    content = function(filename){
      dat <- vcf.d.all()
      dat = add_resistance_info(f.dat = dat, resistance_table=global$res_table, all_muts = T)
      write.csv(x = dat, file = filename, row.names = F)
    }
  )

})
