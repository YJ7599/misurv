library(seqinr)
library(shiny)
library(shinydashboard)
library(dashboardthemes)
#library(tidyverse)
library(phyloseq)
library(plotly)
library(shinyWidgets)
library(shinyjs)
library(googleVis)
library(xtable)
library(DT)
library(htmltools)
library(biomformat)
library(phangorn)
library(bios2mds)
library(zip)
library(randomForestSRC)
library(rgl)
library(vegan3d)
#library(pca3d)
library(jpeg)
library(splitTools)

# COMMENTS
{
  TITLE = p("MiSurv: An Integrative Web Cloud Platform for Microbiome Data Analysis with Survival Responses", style = "font-size:18pt")
  HOME_COMMENT = p(strong("MiSurv", style = "font-size:15pt"), "is an integrative web cloud platform for processing, analyzing and visualizing microbiome data with survival responses. MiSurv consists of a data processing module and its following four data analytic modules: (1) Module 1: Comparative survival analysis between treatment groups, (2) Module 2: Comparative analysis in microbial composition between treatment groups, (3) Module 3: Association testing between microbial composition and survival responses, (4) Module 4: Prediction modeling using microbial taxa on survival responses. More details are as follows.", style = "font-size:13pt")
  
  HOME_COMMENT1 = p(strong("Data Processing : "), "Interactive procedures for (1) data inputs (.rdata, .rds, .biom, .txt, .csv, .tsv, .tre), (2) survival data and analytic plans (survival time, censored/event, follow-up period, subgroup analysis), (3) quality controls (kingdom, library size, mean proportion, taxonomic name), and (4) data transformations (alpha- and beta-diversity calculation, rarefaction, proportion, centered log-ratio, arcsine square root).", style = "font-size:13pt")
  HOME_COMMENT2 = p(strong("Module 1:"), "Comparative survival analysis between treatment groups, not involving microbiome data, with or without covariate adjustment(s). ", style = "font-size:13pt")
  HOME_COMMENT4 = p(strong("Module 3:"), "Association testing between microbial composition and survival responses with or without covariate adjustment(s).", style = "font-size:13pt")
  HOME_COMMENT3 = p(strong("Module 2:"), "Comparative analysis in microbial composition between treatment groups, not involving survival data, with or without covariate adjustment(s).", style = "font-size:13pt")
  HOME_COMMENT5 = p(strong("Module 4:"), "Prediction modeling using microbial taxa at different taxonomic ranks on survival responses", style = "font-size:13pt")
  
  HOME_COMMENT6 = p(strong("URLs:"), " Web server (online implementation):", tags$a(href = "http://misurv.micloud.kr", "http://misurv.micloud.kr"), 
                    "; GitHub repository (local implementation):", 
                    tags$a(href = "https://github.com/wg99526/misurvgit", "https://github.com/wg99526/misurvgit"), style = "font-size:13pt")
  HOME_COMMENT7 = p(strong("Maintainers:"), " Won Gu (", tags$a(href = "wpg5129@psu.edu", "wpg5129@psu.edu"), 
                    "); Hyojung Jang (", tags$a(href = "hyojung.jang@stonybrook.edu", "hyojung.jang@stonybrook.edu"), ")", style = "font-size:13pt")
  HOME_COMMENT8 = p(strong("Reference:"), "Gu W, Koh H, Jang HJ, Lee B, Kang B. MiSurv: An integrative web cloud platform for user-friendly microbiome data analysis with survival responses. (in review)", style = "font-size:13pt")
  
  INPUT_PHYLOSEQ_COMMENT1 = p("Description:", br(), br(), "This should be an '.Rdata' or '.rds' file, and the data should be in the 'phyloseq' format (see", tags$a(href = "https://bioconductor.org/packages/release/bioc/html/phyloseq.html", "https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), ") The phyloseq object should contain all the four necessary data, feature (OTU or ASV) table, taxonomic table, metadata/sample information, and phylogenetic tree.", 
                              br(), br(), "Details:", br(), br(), 
                              "1) The feature table should contain counts, where rows are features (OTUs or ASVs) and columns are subjects 
                              (row names are feature IDs and column names are subject IDs).", br(),
                              "2) The taxonomic table should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                              (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                              'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(),
                              "3) The metadata/sample information should contain variables for the subjects about host phenotypes, medical interventions, 
                              disease status or environmental/behavioral factors, where rows are subjects and columns are variables 
                              (row names are subject IDs, and column names are variable names).", br(), 
                              "4) The phylogenetic tree should be a rooted tree. Otherwise, MiSurv automatically roots the tree through midpoint rooting (phangorn::midpoint). 
                              The tip labels of the phylogenetic tree are feature IDs.", br(), br(), 
                              "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between feature table and metadata/sample information. 
                              MiSurv will analyze only the matched features and subjects."
                              , style = "font-size:11pt")
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download example microbiome data 'biom.Rdata' in the 'phyloseq' format.", br(), 
                              "For more details about 'phyloseq', see", 
                              tags$a(href = "https://bioconductor.org/packages/release/bioc/html/phyloseq.html", "https://bioconductor.org/packages/release/bioc/html/phyloseq.html"),
                              br(),
                              p(strong("Data description:"), "This example data are the public gut microbiome data (Zhang et al. 2018) we used in our real data applications (Gu et al., in review). The raw sequence data are deposited in QIITA (", tags$a(href = "https://qiita.ucsd.edu", "https://qiita.ucsd.edu"), ") with the ID number 10508 (", tags$a(href = "https://qiita.ucsd.edu/study/description/10508", "https://qiita.ucsd.edu/study/description/10508"), "). More detailed sample extraction and raw sequence data processing procedures are described in (Zhang et al. 2018; Gu et al., in review). "), br(), 
                              "> setwd('/yourdatadirectory/')", br(), br(), 
                              "> load(file = 'biom.Rdata')", br(), br(), 
                              "> library(phyloseq)", br(), br(), 
                              " > otu.tab <- otu_table(biom)", br(), 
                              " > tax.tab <- tax_table(biom)", br(), 
                              " > tree <- phy_tree(biom)", br(), 
                              " > sam.dat <- sample_data(biom)", br(), br(), 
                              " > sam.dat$T1Dweek",   HTML('&emsp;'),HTML('&nbsp;'),   "# Survival Time", br(), 
                              " > sam.dat$T1D ", HTML('&emsp;'),HTML('&emsp;'),HTML('&emsp;'), "# Event", br(), 
                              " > sam.dat$Antibiotics ", HTML('&emsp;'),  "# Treatment", br(), br(),
                              "You can check if the features are matched and identical across feature table, taxonomic table and 
                              phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information 
                              using following code.", br(), br(), 
                              " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                              " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT = p("Description:", br(), br(), 
                                    "1) The feature table (.txt or .csv) should contain counts, where rows are features (OTUs or ASVs) and columns are subjects 
                                    (row names are feature IDs and column names are subject IDs). Alternatively, you can upload .biom file processed by QIIME.", br(), 
                                    "2) The taxonomic table (.txt) should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                                    (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                                    'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'). Alternatively, you can upload .tsv file processed by QIIME.", br(), 
                                    "3) The metadata/sample information (.txt or .csv) should contain variables for the subjects about host phenotypes, medical interventions, 
                                    disease status or environmental/behavioral factors, where rows are subjects and columns are variables (row names are subject IDs, and 
                                    column names are variable names).", br(), 
                                    "4) The phylogenetic tree (.tre or .nwk) should be a rooted tree. Otherwise, MiSurv automatically roots the tree through midpoint 
                                    rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.", br(), br(), 
                                    "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                                    The subjects should be matched and identical between feature table and metadata/sample information. 
                                    MiSurv will analyze only the matched features and subjects.", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT2 = p("You can download example microbiome data 'biom.zip'. This zip file contains four necessary data, feature table (otu.tab.txt), 
                                     taxonomic table (tax.tab.txt), metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre).", br(), br(),
                                     "> setwd('/yourdatadirectory/')", br(), br(), 
                                     "> otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE) ", br(), 
                                     "> tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)", br(), 
                                     " > sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE) ", br(),
                                     "> tree <- read.tree(file = 'tree.tre')", br(), br(), 
                                     " > sam.dat$T1Dweek",   HTML('&emsp;'),HTML('&nbsp;'),   "# Survival Time", br(), 
                                     " > sam.dat$T1D ", HTML('&emsp;'),HTML('&emsp;'),HTML('&emsp;'), "# Event", br(), 
                                     " > sam.dat$Antibiotics ", HTML('&emsp;'),  "# Treatment", br(), br(),
                                     "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, 
                                     and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                                     " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                                     " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                                     " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt")
  
  QC_KINGDOM_COMMENT = p("A microbial kingdom to be analyzed. Default is 'Bacteria' for 16S data. Alternatively, you can type 'Fungi' for ITS data 
                         or any other kingdom of interest for shotgun metagenomic data.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT1 = p("Remove subjects that have low library sizes (total read counts). Default is 3,000.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT2 = p("Library size: The total read count per subject.", style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT1 = p("Remove features (OTUs or ASVs) that have low mean relative abundances (Unit: %). Default is 0.002%.",style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT2 = p("Mean proportion: The average of relative abundances (i.e., proportions) per feature.", style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT1 = p('Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings. 
                            Multiple character strings should be separated by a comma. Default is "", "metagenome", "gut metagenome", "mouse gut metagenome".',
                            style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT2 = p('Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings (i.e., taxonomic names that contain 
                            the specified character strings). Multiple character strings should be separated by a comma. Default is "uncultured", "incertae", "Incertae",
                            "unidentified", "unclassified", "unknown".',style = "font-size:11pt")
  
  ALPHA_COMMENT = p("Calculate alpha-diversity indices: Richness (Observed), Shannon (Shannon, 1948), Simpson (Simpson, 1949), Inverse Simpson (Simpson, 1949), 
                    Fisher (Fisher et al., 1943), Chao1 (Chao, 1984), ACE (Chao and Lee, 1992), ICE (Lee and Chao, 1994), PD (Faith, 1992).")
  ALPHA_REFERENCES = p("1. Chao A, Lee S. Estimating the number of classes via sample coverage. J Am Stat Assoc. 1992:87:210-217.", br(),
                       "2. Chao A. Non-parametric estimation of the number of classes in a population. Scand J Stat. 1984:11:265-270.", br(),
                       "3. Faith DP. Conservation evaluation and phylogenetic diversity. Biol Conserv. 1992:61:1-10.", br(),
                       "4. Fisher RA, Corbet AS, Williams CB. The relation between the number of species and the number of individuals 
                       in a random sample of an animal population. J Anim Ecol. 1943:12:42-58.", br(),
                       "5. Lee S, Chao A. Estimating population size via sample coverage for closed capture-recapture models. Biometrics. 1994:50:1:88-97.", br(),
                       "6. Shannon CE. A mathematical theory of communication. Bell Syst Tech J. 1948:27:379-423 & 623-656.", br(),
                       "7. Simpson EH. Measurement of diversity. Nature 1949:163:688.", br())
  BETA_COMMENT = p("Calculate beta-diversity indices: Jaccard dissimilarity (Jaccard, 1912), Bray-Curtis dissimilarity (Bray and Curtis, 1957), Unweighted UniFrac distance 
                   (Lozupone and Knight, 2005), Generalized UniFrac distance (Chen et al., 2012), Weighted UniFrac distance (Lozupone et al., 2007).")
  BETA_REFERENCES = p("1. Bray JR, Curtis JT. An ordination of the upland forest communities of Southern Wisconsin. Ecol Monogr. 1957;27(32549).", br(),
                      "2. Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD., et al. Associating microbiome composition with environmental 
                      covariates using generalized UniFrac distances. Bioinformatics. 2012;28(16):2106-13.", br(),
                      "3. Jaccard P. The distribution of the flora in the alpine zone. New Phytol. 1912;11(2):37-50.", br(),
                      "4. Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and qualitative β-diversity measures lead to 
                      different insights into factors that structure microbial communities. Appl Environ Microbiol. 2007;73(5):1576-85.", br(),
                      "5. Lozupone CA, Knight R. UniFrac: A new phylogenetic method for comparing microbial communities. Appl Environ Microbiol. 2005;71(12):8228-35.")
  DATA_TRANSFORM_COMMENT = p("Transform the count (original) data into four different formats 1) count (rarefied) (Sanders, 1968), 2) proportion, 3) centered log ratio (CLR) (Aitchison, 1982), 4) arcsine-root for each taxonomic rank (phylum, class, order, family, genus, species).")
  DATA_TRANSFORM_REFERENCE = p("1. Aitchison J. The statistical analysis of compositional data. J R Statist Soc B. 1982;44:2:139-77")
  
  #Surv_Treat_Comment = p("Used to compare two categories of the selected binary variable.",style = "font-size:11pt")
  Surv_Time__Comment = p( "A variable for follow-up time (required).", style = "font-size:11pt")
  Surv_CensorComment = p("A binary indicator for censored (0) or event (1) (required).", style = "font-size:11pt")
  Surv_FollowComment = p("You can adjust the follow-up period (start time and end time) of interest (or leave it as it is to survey the entier follow-up period just as given in the data), where the censored/event indicator is also automatically adjusted to be suited to the selected follow-up period (optional).", style = "font-size:11pt")
  Follow_Time_Comment = p("You can reset the follow-up period of interest or leave it as it is.", style = "font-size:11pt")
  Num_PredComment = p("The number of taxa at each taxonomic rank to be displayed on the importance plot (default: 10).", style = "font-size:11pt")
}

# UI
{
  ui = dashboardPage(
    title = "MiSurv",
    dashboardHeader(title = span(TITLE, style = "float:left;font-size: 20px"), titleWidth = "100%"),
    dashboardSidebar(
      tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
      sidebarMenu(id = "side_menu",
                  menuItem("Home", tabName = "home", icon = icon("home")),
                  menuItem("Data Processing",  icon = icon("file-text-o"),
                           menuSubItem("Data Input", tabName = "step1", icon = icon("mouse")),
                           menuSubItem("Quality Control", tabName = "step2", icon = icon("chart-bar")),
                           menuSubItem("Data Transformation", tabName = "step3", icon = icon("calculator"))), #or, icon th-large
                  menuItem("Module 1", tabName = "SurvAnalysis", icon = icon("bar-chart-o")
                  ),
                  menuItem("Module 2",  icon = icon("chart-pie"),
                           menuSubItem("Alpha Diversity", tabName = "alphaDivanalysis", icon = icon("font")),
                           menuSubItem("Beta Diversity", tabName = "betaDivanalysis", icon = icon("bold")),
                           menuSubItem("Taxonomic Abundance", tabName = "taxaAnalysis", icon = icon("align-left"))),
                  menuItem("Module 3",  icon = icon("disease"),
                           menuSubItem("Alpha Diversity", tabName = "alphaDivanalysisSurv", icon = icon("font")),
                           menuSubItem("Beta Diversity", tabName = "betaDivanalysisSurv", icon = icon("bold")),
                           menuSubItem("Taxonomic Abundance", tabName = "taxaAnalysisSurv", icon = icon("align-left"))),
                  menuItem("Module 4", tabName = "RandomForest", icon = icon("tree")))),
    dashboardBody(
      tags$head(tags$style(HTML(".content { padding-top: 2px;}"))),
      tags$script(src = "fileInput_text.js"),
      useShinyjs(),
      shinyDashboardThemes(theme = "poor_mans_flatly"),
      uiOutput("themes"),
      tabItems(
        
        ##### HOME ####
        tabItem(tabName = "home",
                div(id = "homepage", br(), HOME_COMMENT, 
                    p(" ", style = "margin-bottom: 10px;"),
                    div(tags$img(src='MiSurv_workflow.png', height = 656, width = 1200), style="text-align: center;"),
                    br(),
                    tags$ul(
                      tags$li(HOME_COMMENT1), tags$li(HOME_COMMENT2), tags$li(HOME_COMMENT3), tags$li(HOME_COMMENT4), tags$li(HOME_COMMENT5),
                      style = "font-size:13pt"),
                    br(),
                    HOME_COMMENT6,
                    HOME_COMMENT7,
                    HOME_COMMENT8,
                    br()
                )),
        
        ##### DATA INPUT ####
        tabItem(tabName = "step1", br(),
                fluidRow(column(width = 6, style='padding-left:+20px',
                                box(
                                  width = NULL, status = "primary", solidHeader = TRUE,
                                  title = strong("Data Input"),
                                  selectInput("inputOption", h4(strong("Data Type?")), c("Choose one" = "", "Phyloseq", "Individual Data"), width = '30%'),
                                  div(id = "optionsInfo", tags$p("You can choose phyloseq or individual data.", style = "font-size:11pt"), style = "margin-top: -15px"),
                                  uiOutput("moreOptions"))),
                         column(width = 6, style='padding-left:0px', 
                                uiOutput("addDownloadinfo"), uiOutput("ref_micloud")))
        ),
        
        ##### QC ####
        tabItem(tabName = "step2", 
                br(), 
                fluidRow(column(width = 3, style = "padding-left:+20px", 
                                box(
                                  width = NULL, status = "primary", solidHeader = TRUE, 
                                  title = strong("Survival Data"), 
                                  h4(strong("Survival Time?", style = "color:black")),
                                  Surv_Time__Comment,
                                  p(" ", style = "margin-bottom: -20px;"),
                                  selectInput("surv.Time.select", label = "",
                                              c("Choose one" = "", ""), selected = FALSE, multiple = FALSE, selectize = FALSE, width = '70%'), 
                                  
                                  h4(strong("Censored/Event?", style = "color:black")),
                                  Surv_CensorComment,
                                  p(" ", style = "margin-bottom: -20px;"),
                                  selectInput("censor.select", label = "",
                                              c("Choose one" = "", ""), selected = FALSE, multiple = FALSE, selectize = FALSE, width = '70%'),
                                  
                                  p(" ", style = "margin-bottom: +20px;"),
                                  h4(strong("Adjust Follow-up Period?")),
                                  Surv_FollowComment,
                                  p(" ", style = "margin-bottom: -30px;"),
                                  sliderInput("follow.time", label = "", min=0, max=20, value = c(0,15), step = 1)
                                ),
                                box(
                                  width = NULL, status = "primary", solidHeader = TRUE, 
                                  title = strong("Microbiome Data"), 
                                  textInput("kingdom", h4(strong("Bacteria?")), value = "Bacteria"),
                                  QC_KINGDOM_COMMENT,
                                  tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                                  tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'),
                                  
                                  
                                  sliderInput("slider1", h4(strong("Library Size?")), min=0, max=10000, value = 3000, step = 1000),
                                  QC_LIBRARY_SIZE_COMMENT1,
                                  QC_LIBRARY_SIZE_COMMENT2,
                                  
                                  sliderInput("slider2", h4(strong("Mean Proportion?")), min = 0, max = 0.1, value = 0.002, step = 0.001,  post  = " %"),
                                  QC_MEAN_PROP_COMMENT1,
                                  QC_MEAN_PROP_COMMENT2,
                                  
                                  h4(strong("Erroneous Taxonomic Names?")),
                                  textInput("rem.str", label = "Complete Match", value = ""),
                                  QC_TAXA_NAME_COMMENT1,
                                  
                                  textInput("part.rem.str", label = "Partial Match", value = ""),
                                  QC_TAXA_NAME_COMMENT2,
                                  
                                  actionButton("run", (strong("Run!")), class = "btn-info"), br(), br(),
                                  uiOutput("moreControls")
                                )),
                         column(width = 9, style = "padding-left:+20px",
                                box(width = NULL, status = NULL, solidHeader = FALSE, 
                                    fluidRow(width = 12,
                                             status = "primary", solidHeader = TRUE, 
                                             valueBoxOutput("sample_Size", width = 3),
                                             valueBoxOutput("OTUs_Size", width = 3),
                                             valueBoxOutput("phyla", width = 3),
                                             valueBoxOutput("classes", width = 3)),
                                    fluidRow(width = 12, 
                                             status = "primary", solidHeader = TRUE,  #여깅
                                             valueBoxOutput("orders", width = 3),
                                             valueBoxOutput("families", width = 3),
                                             valueBoxOutput("genera", width = 3),
                                             valueBoxOutput("species", width = 3)),
                                    fluidRow(style = "position:relative",
                                             tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                                                    tabPanel("Histogram",
                                                             plotlyOutput("hist"),
                                                             sliderInput("binwidth", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                             chooseSliderSkin("Round", color = "#112446")),
                                                    tabPanel("Box Plot", 
                                                             plotlyOutput("boxplot"))),
                                             tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                                                    tabPanel("Histogram",
                                                             plotlyOutput("hist2"),
                                                             sliderInput("binwidth2", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                             chooseSliderSkin("Round", color = "#112446")),
                                                    tabPanel("Box Plot",
                                                             plotlyOutput("boxplot2"))))
                                ))),
                
                # mainPanel(width = 9,
                #           fluidRow(width = 12,
                #                    status = "primary", solidHeader = TRUE, 
                #                    valueBoxOutput("sample_Size", width = 3),
                #                    valueBoxOutput("OTUs_Size", width = 3),
                #                    valueBoxOutput("phyla", width = 3),
                #                    valueBoxOutput("classes", width = 3)),
                #           fluidRow(width = 12, 
                #                    status = "primary", solidHeader = TRUE,
                #                    valueBoxOutput("orders", width = 3),
                #                    valueBoxOutput("families", width = 3),
                #                    valueBoxOutput("genera", width = 3),
                #                    valueBoxOutput("species", width = 3)),
                #           fluidRow(style = "position:relative",
                #                    tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                #                           tabPanel("Histogram",
                #                                    plotlyOutput("hist"),
                #                                    sliderInput("binwidth", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                #                                    chooseSliderSkin("Round", color = "#112446")),
                #                           tabPanel("Box Plot", 
                #                                    plotlyOutput("boxplot"))),
                #                    tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                #                           tabPanel("Histogram",
                #                                    plotlyOutput("hist2"),
                #                                    sliderInput("binwidth2", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                #                                    chooseSliderSkin("Round", color = "#112446")),
                #                           tabPanel("Box Plot",
                #                                    plotlyOutput("boxplot2")))))
                
        ),
        
        ##### Data Transformation ####
        tabItem(tabName = "step3", br(),
                fluidRow(column(width = 6, style = 'padding-left:+20px',
                                box(title = strong("Data Transformation"), width = NULL, status = "primary", solidHeader = TRUE,
                                    ALPHA_COMMENT, 
                                    BETA_COMMENT, 
                                    DATA_TRANSFORM_COMMENT,
                                    actionButton("divCalcRun", (strong("Run!")), class = "btn-info"), 
                                    p(" ", style = "margin-bottom: +10px;"),
                                    p(strong("Attention:"), "Once you changed any setting in the preceding Data Input or Quality Control, you have to click this Run button again to use any of the following data analytic modules."),
                                ),
                                uiOutput("divCalcDownload")),
                         column(width = 6, style='padding-left:0px',
                                box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                                    p("Alpha Diversity", style = "font-size:12pt"),
                                    ALPHA_REFERENCES,
                                    p("Beta Diversity", style = "font-size:12pt"),
                                    BETA_REFERENCES,
                                    p("Taxonomic Abundance", style = "font-size:12pt"),
                                    DATA_TRANSFORM_REFERENCE)))
                
        ),
        
        ##### Model 1: Survival Analysis ####
        tabItem(tabName = "SurvAnalysis", br(),
                sidebarLayout(
                  position = "left", 
                  div(style="width: 98%;",
                      sidebarPanel(width = 3,
                                   uiOutput("primvarsSurv"),
                                   uiOutput("prim_vars_types.Surv"),
                                   uiOutput("surviveTime"),
                                   uiOutput("censoring"),
                                   uiOutput("subgroup.sel"),
                                   uiOutput("covariates.Surv"), br(),
                                   uiOutput("referencesM1"))),
                  mainPanel(width = 9,
                            fluidRow(width = 8,
                                     status = "primary", solidHeader = TRUE, 
                                     div(style="margin-left: -1.2%; margin-right:-1.2%", valueBoxOutput("control_size", width = 3 ),
                                         valueBoxOutput("treatment_size", width = 3 ),
                                         valueBoxOutput("censored_size", width = 3 ),
                                         valueBoxOutput("uncensored_size", width = 3 ))),
                            #valueBoxOutput("significant_test", width = 3)),
                            fluidRow(
                              #div(style='width: 95%',
                              tabPanel(width = 8,"Survival Curve", #width = 12,
                                       uiOutput("surv_display_results")))))),#),
        
        
        ##### Model 2: ALPHA DIVERSITY ####
        tabItem(tabName = "alphaDivanalysis", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("primvars"),
                                                       uiOutput("prim_vars_types"),
                                                       uiOutput("covariates"), 
                                                       p(" ", style = "margin-bottom: +25px;"),
                                                       actionButton("runbtn_bin", (strong("Run!")), class = "btn-info"),
                                                       #br(), 
                                                       uiOutput("alpha_downloadTable"),
                                                       uiOutput("referencesM2.alpha"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("alpha_display_results"))))),
        
        ##### Model 2: Taxa Analysis ####
        tabItem(tabName = "taxaAnalysis", br(),
                sidebarLayout( 
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("primvars_taxa"),
                                                       uiOutput("morePrimvar_opt_taxa"),
                                                       uiOutput("covariates_taxa"), 
                                                       actionButton("taxa_runbtn_bin", (strong("Run!")), class = "btn-info"),
                                                       uiOutput("dttttfragile"),
                                                       uiOutput("rtttt.anti"))),
                  
                  #mainPanel(width = 9, uiOutput("airpods"))
                  mainPanel(width = 9, 
                            fluidRow(width = 8,
                                     uiOutput("airpods")
                                     ))
                  
                  #mainPanel("mainpanel")
                  #mainPanel(width = 9,
                  #          fluidRow(width = 8, 
                  #                   div(style='height:800px;overflow-y: scroll;', uiOutput("airpods")), br(),br(),
                  #                   uiOutput("nine_five_four"))
                  )),
        
        ##### Model 2: BETA DIVERSITY ####
        tabItem(tabName = "betaDivanalysis", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("beta_primvar_cross"),
                                                       uiOutput("beta_prim_vars_types_cross"),
                                                       uiOutput("beta_covariates_cross"), 
                                                       actionButton("beta_runbtn_cross_bin", (strong("Run!")), class = "btn-info"),
                                                       #br(), 
                                                       uiOutput("beta_downloadTable"),
                                                       uiOutput("referencesM2.beta"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("beta_display_results_cross"))))),
        
        
        
        
        ##### Model 3: ALPHA DIVERSITY ####
        tabItem(tabName = "alphaDivanalysisSurv", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("subgroupSurv.A"),
                                                       uiOutput("subgroup.sel.A"),
                                                       uiOutput("surviveTimeCoxA"),
                                                       uiOutput("covariatesCoxA"), 
                                                       actionButton("runbtn_CoxA", (strong("Run!")), class = "btn-info"),
                                                       #br(),
                                                       uiOutput("alpha_surv_downloadTable"),
                                                       uiOutput("referencesM3.alpha"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("alpha_surv_display_resultsCoxA"))))),
        
        ##### Model 3: BETA  DIVERSITY ####
        tabItem(tabName = "betaDivanalysisSurv", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("subgroupSurv.B"),
                                                       uiOutput("subgroup.sel.B"),
                                                       uiOutput("surviveTimeCoxB"),
                                                       uiOutput("covariatesCoxB"), 
                                                       actionButton("runbtn_CoxB", (strong("Run!")), class = "btn-info"),
                                                       #br(),
                                                       uiOutput("beta_surv_downloadTable"),
                                                       uiOutput("referencesM3.beta"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("betaS_display_results_cross"))))),
        ##### Model 3: Taxonomic Abundance ####
        tabItem(tabName = "taxaAnalysisSurv", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("censor.t"),
                                                       uiOutput("subgroupSurv.T"),
                                                       uiOutput("group.t.surv"),
                                                       uiOutput("covariatesCoxT"), 
                                                       actionButton("runbtn_CoxT", (strong("Run!")), class = "btn-info"),
                                                       #br(),
                                                       uiOutput("taxa_surv_downloadTable"),
                                                       uiOutput("referencesM3.taxa"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     div(style='height:700px;overflow-y: scroll;', uiOutput("ha_shi")), br(),br(), #s_taxa_display_result 
                                     
                                     uiOutput("ple_surv_dend"))
                            
                  ))),
        
        ##### Model 4: Prediction ####
        tabItem(tabName = "RandomForest", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("model4_data_format"),
                                                       uiOutput("subgroupSurv.4"),
                                                       uiOutput("method.4"),
                                                       actionButton("runbtn_model4", (strong("Run!")), class = "btn-info"),
                                                       #br(),
                                                       uiOutput("downloadTable_m4"),
                                                       uiOutput("referencesM4"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12,
                                     div(style='height:800px;overflow-y: scroll;', uiOutput("surv4_display_results"))))))
      )
    )
  )
}

