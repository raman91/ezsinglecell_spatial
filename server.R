## max data size
options(shiny.maxRequestSize = 1024^10)
options(shiny.launch.browser = T)

shinyServer(function(input, output, session) {
    v <- reactiveValues(scData = NULL,
                        scRNAData = NULL,
                        idents = NULL,
                        isUMAPdone = NULL,
                        isDeconvdone = NULL,
                        isClusterdone = NULL,
                        pcGenes = NULL,
                        plotlySelection = NULL,
                        ips.markers = NULL,
                        ips.spatialfeatures = NULL,
                        ip.markers1 = NULL)
    celltypes <- NULL
    prePlot <- function(){
      while(names(dev.cur()) != "null device"){
        dev.off()
      }
    }
    observe({
        #s <- event_data("plotly_selected")
        #cells <- s[["key"]]
        v$plotlySelection <- event_data("plotly_selected")[["key"]]
    })
    ##-------------------Side Panel-------------------

    normMethod <- NULL

    output$name.field <- renderUI({
        if(is.null(input$cellAnnoFiles)){
            numericInput(inputId = "field",
                         label = "Field",
                         value = 1,
                         min = 1)
        }else{
            annoFile <- input$cellAnnoFiles
            anno.data <- read.table(annoFile$datapath[1], header = T,
                                    sep = "\t", stringsAsFactors = FALSE)[1:5,]
            groupings <- colnames(anno.data)
            selectInput("groupby",
                        label = "Group by:",
                        choices = groupings)
        }
    })

    observeEvent(input$loadButton, {
        tpmFiles <- input$tpmFiles
        annoFile <- input$cellAnnoFiles
        names.field <- input$field
        if(input$norm){
          normMethod <- "LogNormalize"
        }
        if (is.null(tpmFiles)){
            v$scData <- NULL
        }else{
            withProgress(message="Loading and Processing Data...", value=0, {
                print(tpmFiles$datapath)
                print(tpmFiles$name)
                print(file.exists(paste(tpmFiles$datapath[1], "/", tpmFiles$name[1], sep="")))
                path = "/acrc_raman/jinmiao/CJM_lab/Raman/Projects/hyperion_cytofkit2/spatial_shiny/stxBrain/"
                exp.data <- Load10X_Spatial(path, filename = "V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5", assay = "Spatial")
                additional.ident <- NULL
                if(!is.null(annoFile)){
                    anno.data <- read.table(annoFile$datapath[1], header = T,
                                            sep = "\t", stringsAsFactors = FALSE)
                    to.append <- apply(anno.data, 1, paste, collapse = "_")
                    colnames(exp.data) <- to.append
                    names.field <- match(input$groupby, colnames(anno.data))
                    additional.ident <- data.frame(data.frame(anno.data[,-1], row.names = to.append))
                    additional.ident[] <- lapply(additional.ident, factor)
                }
                incProgress(0.5, "Creating Seurat Object")
                #sObj <- CreateSeuratObject(exp.data,
                #              project = input$projName,
                #              names.field = names.field,
                #              names.delim = input$delim,
                #              is.expr = input$expThres,
                #              normalization.method = normMethod,
                #              min.genes = input$min.genes,
                #              min.cells = input$min.cells)
                #mito.genes <- grep("^MT-", rownames(sObj@assays$RNA@data), ignore.case = TRUE, value = TRUE)
                #sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
                #incProgress(0.5, "Adding metadata")
                #sObj <- AddMetaData(sObj, percent.mito, "percent.mito")
                #if(!is.null(additional.ident)){
                #    sObj <- AddMetaData(sObj, additional.ident)
                #}
                v$scData <- exp.data
            })
        }
        dir.create("Seurat_results")
    })

    observeEvent(input$reset, {
      session$reload()
      print("Reset done")
    })

    observeEvent(input$saveButton, {
        if(!is.null(input$tpmFiles)){
            withProgress(message="Saving Results...", value=0, {
                print(getwd())
                dir.create("Seurat_results")
                resultDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results")
                filename <- paste0(resultDir, .Platform$file.sep, v$scData@project.name, "_", Sys.Date())
                sObj <- v$scData
                save(sObj, file= paste0(resultDir, .Platform$file.sep, sObj@project.name, "_", Sys.Date(), ".Robj"))
            })
          ## open the results directory
          opendir(resultDir)
        }
    })

    output$ident.swap <- renderUI({
        if(is.null(v$scData)){
            return(NULL)
        }else{
            groupings <- names(v$scData@meta.data[,!names(v$scData@meta.data) %in% c("nFeature_RNA", "nCount_RNA", "percent.mt")])
            tagList(
                h4("Set current identity:"),
                fluidRow(
                    column(6,
                           selectInput("active.ident", label = NULL,
                                       choices = groupings)
                           ),
                    column(6,
                           actionButton("swap.ident",label = NULL, icon = icon("arrow-right"))
                           )
                )

            )
        }
    })

    observeEvent(input$swap.ident, {
        v$scData <- SetIdent(v$scData, value = as.character(v$scData@meta.data[,input$active.ident]))
    })

    output$logo <- renderImage({
      return(list(
        src = "inst/extdata/logo.png",
        contentType = "image/png",
        alt = "Singapore Immunology Network"
      ))
    }, deleteFile = FALSE)

    opendir <- function(dir = getwd()){
      if (.Platform['OS.type'] == "windows"){
        shell.exec(dir)
      } else {
        system(paste(Sys.getenv("R_BROWSER"), dir))
      }
    }

    observeEvent(input$OpenDir, {
      resultDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results")
      if(!dir.exists(resultDir)){
        dir.create("Seurat_results")
      }
      pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
      if(dir.exists(pdfDir)){
        opendir(pdfDir)
      }else{
        warning("No reports created yet!")
        dir.create(pdfDir)
      }
    })

    ##---------------QC tabset-------------------

    output$nCount_SpatialPlot <- renderPlot({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            VlnPlot(v$scData, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
            }
    })

    output$SpatialFeaturePlot <- renderPlot({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
           SpatialFeaturePlot(v$scData, features = "nCount_Spatial") + theme(legend.position = "right")
           }
    })

    #output$nCount_RNAPlot <- renderPlotly({
    #    if(is.null(v$scData)){
    #        plotly_empty()
    #    }else{
    #        VlnPlot(v$scData, "nCount_RNA")
    #    }
    #})

    output$name <- renderPrint({
        s <- event_data("plotly_selected")
        c(s[["key"]], class(s[["key"]]))
    })

    observeEvent(input$PDFa, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"QC_violin_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "QC_spatial_feature_plot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                nG <- VlnPlot(v$scData, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
                pM <- SpatialFeaturePlot(v$scData, features = "nCount_Spatial") + theme(legend.position = "right")
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(nG)
                print(pM)
                #print(nU)
                dev.off()
            })
        }
    })

   ##---------------PCA tabset-------------------
    # PCA plot
    observeEvent(input$doPCA, {
      withProgress(message = "Scaling Data...", value = 0,{
        incProgress(0.5, message = "Running Dimension Reduction...")
        v$scData <- SCTransform(v$scData, assay = "Spatial", verbose = FALSE)

        output$deg1.gene.select <- renderUI({
          if(is.null(v$scData)){
            return(NULL)
          }else{
            selectInput("deg1.gene", label = "Gene to visualise",
                        choices = rownames(v$scData))
          }
        })

        output$Deg1.plot <- renderPlot({
          if(is.null(v$scData)){
            return(NULL)
          }else{
            withProgress(message="Generating DEG Plot...", value=0, {
              SpatialFeaturePlot(object = v$scData, features = input$deg1.gene, alpha = c(0.1, 1))
            })
          }
        })

        v$scData <- RunPCA(v$scData, assay = "SCT", verbose = FALSE)
        v$scData <- FindNeighbors(v$scData, reduction = "pca", dims = 1:30)
        v$scData <- FindClusters(v$scData, verbose = FALSE)
        v$scData <- RunUMAP(v$scData, reduction = "pca", dims = 1:30)
        print(v$scData[["pca"]], dims = 1:5, nfeatures = 5)
        v$isUMAPdone <- TRUE
        Dim_plot1 <- DimPlot(v$scData, reduction = "umap", label = TRUE)
        Dim_plot2 <- SpatialDimPlot(v$scData, label = TRUE, label.size = 3)
        print(Dim_plot1)
        print(Dim_plot2)
        incProgress(0.4, message = "Getting list of PC genes...")
        pc.table <- list()
        for(i in 1:20){
            pcg <- TopCells(v$scData)
            pc.table[[i]] <- pcg
        }
        pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
        v$pcGenes <- pc.table
      })
    })

    output$clustUI <- renderUI({
      if(is.null(v$isUMAPdone)){
        return(NULL)
      }else{
        tagList(
          fluidRow(
              column(6,
                     numericInput("clus.res",
                                  label = "Cluster Resolution",
                                  value = 0.6,
                                  min = 0.1,
                                  step = 0.1)
                     ),
              column(6,
                     actionButton("findCluster", "Find Clusters", icon = icon("hand-pointer-o")),
                     textOutput("cluster.done")
                     )
          )
        )
      }
    })

    observeEvent(input$findCluster, {
      withProgress(message = "Finding clusters...", value = 0.3, {
        v$scData <- FindNeighbors(v$scData, dims = 1:10)
        v$scData <- FindClusters(v$scData, resolution = input$clus.res)
        output$cluster.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone <- TRUE
      })
    })

    output$DimPlot <- renderPlotly({
        if(is.null(v$isUMAPdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating Dim Plot...", value=0, {
                DimPlot(v$scData, label = TRUE)
            })
        }
    })

    output$SpatialDimPlot <- renderPlot({
        if(is.null(v$isUMAPdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating SpatialDim Plot...", value=0, {
               SpatialDimPlot(v$scData, label = TRUE, label.size = 3)
            })
        }
    })

    observeEvent(input$PDFd, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"UMAP_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "UMAP_plot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                dimplot1 <- DimPlot(v$scData, reduction = "umap", label = T)
                dimplot2 <- SpatialDimPlot(v$scData, label = TRUE, label.size = 3)
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(dimplot1)
                print(dimplot2)
                dev.off()
            })
            withProgress(message="Downloading UMAP coordinates...", value=0.5, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"umap_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "umap_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                    i = i + 1;
                }
                write.csv(v$scData@reductions$umap@cell.embeddings, file = filename2)
            })
            withProgress(message="Downloading cluster IDs...", value=0.9, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                    i = i + 1;
                }
                write.csv(v$scData@active.ident, file = filename2)
            })
        }
    })

    # Viz plot

    output$vizPlot <- renderPlot({
        if(is.null(v$scData)){
            return(NULL)
        }else{
          VizDimLoadings(v$scData, dims = as.numeric(input$select.pc))
        }
    })

    output$PCHeatmap <- renderPlot({
        if(is.null(v$scData)){
            return(NULL)
        }else{
            DimHeatmap(v$scData, dims = as.numeric(input$select.pc))
        }
    })

    output$PCtable <- DT::renderDataTable({
        if(is.null(v$scData) ){
            return(NULL)
        }else{
            v$pcGenes
        }
    }, options = list(scrollX = TRUE))

    observeEvent(input$PDFe, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"Viz_Heatmap_plots_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "Viz_Heatmap_plots_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                #isolate({
                    vizdim <- VizDimLoadings(v$scData, dims = as.numeric(input$select.pc))
                    dimheat <- DimHeatmap(v$scData)
                    print (vizdim)
                    print (dimheat)
                #})
                dev.off()
                pcGenes <- v$pcGenes
                write.csv(v$pcGenes, file = paste0(pdfDir, .Platform$file.sep,"PC_genes_", Sys.Date(), ".csv"))
            })
        }
    })

    ##---------------DEGs tabset-------------------

    observeEvent(input$doDeg, {
      if(is.null(v$scData)){
        return(NULL)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
            ips.markers <- FindAllMarkers(v$scData, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc)
            v$ips.markers <- ips.markers
        })
      }
    })

    output$deg.gene.select <- renderUI({
      if(is.null(v$ips.markers)){
        return(NULL)
      }else{
        selectInput("deg.gene", label = "Gene to visualise",
                    choices = rownames(v$ips.markers))
      }
    })

    output$Deg.plot <- renderPlot({
      if(is.null(v$scData)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          SpatialFeaturePlot(object = v$scData, features = input$deg.gene, alpha = c(0.1, 1))
        })
      }
    })

    # output$Deg1.plot <- renderPlotly({
    #  if(is.null(v$ips.markers)){
    #    return(NULL)
    #  }else{
    #      withProgress(message="Generating DEG Plot...", value=0, {
    #       FeaturePlot(v$scData, input$deg.gene)
    #      })
    #    }
    #  })

    output$Deg.table <- DT::renderDataTable(
      v$ips.markers, options = list(scrollX = TRUE, scrollY = "400px"))

    observeEvent(input$doDegn, {
      if(is.null(v$scData)){
        return(NULL)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
          v$scData <- FindSpatiallyVariableFeatures(v$scData, assay = "SCT", features = VariableFeatures(v$scData)[1:1000], selection.method = "markvariogram")
          ips.spatialfeatures <- SpatiallyVariableFeatures(v$scData, selection.method = "markvariogram")
          v$ips.spatialfeatures <- ips.spatialfeatures
        })
      }
    })

    output$degn.gene.select <- renderUI({
      if(is.null(v$ips.spatialfeatures)){
        return(NULL)
      }else{
        selectInput("degn.gene", label = "Gene to visualise",
                    choices = v$ips.spatialfeatures)
      }
    })

    output$Degn.plot <- renderPlot({
      if(is.null(v$scData)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          SpatialFeaturePlot(object = v$scData, features = input$degn.gene, alpha = c(0.1, 1))
        })
      }
    })

    observeEvent(input$PDFk, {
      if(!is.null(v$scData)){
        withProgress(message="Downloading plot PDF files...", value=0, {
          print(getwd())
          pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
          if(!dir.exists(pdfDir)){
            dir.create(pdfDir)
          }
          filename2 <- paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), ".pdf")
          i = 0
          while(file.exists(filename2)){
            filename2 <- paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
            i = i + 1;
          }
          degSpatial <- SpatialFeaturePlot(object = v$scData, features = input$deg.gene, alpha = c(0.1, 1))
          #degFeature <- FeaturePlot(v$scData, input$deg.gene)
          prePlot()
          pdf(filename2,
              width=as.numeric(input$pdf_w),
              height=as.numeric(input$pdf_h))
          print(degSpatial)
          #print(degFeature)
          dev.off()
          write.csv(v$ips.markers, file = paste0(pdfDir, .Platform$file.sep,"DEG_table_", input$c1, "vs", input$c2, "_", Sys.Date(), ".csv"))
        })
      }
    })

    ##---------------Deconvolution-------------------

      observeEvent(input$loadButton1, {
      tpmFiles1 <- input$tpmFiles1
      if (is.null(tpmFiles1)){
        v$scRNAData <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(tpmFiles1$datapath)
          print(tpmFiles1$name)
          print(file.exists(paste(tpmFiles1$datapath[1], "/", tpmFiles1$name[1], sep="")))
          exp.data1 <- readRDS(tpmFiles1$datapath)
          incProgress(0.5, "Creating Seurat Object")
          v$scRNAData <- exp.data1
        })
      }
    })

    observeEvent(input$process_scRNA, {
      withProgress(message = "Processing scRNA-seq dataset...", value = 0,{
        incProgress(0.5, message = "Processing...")
        v$scRNAData <- Seurat::SCTransform(v$scRNAData, verbose = FALSE)
        v$scRNAData <- Seurat::RunPCA(v$scRNAData, verbose = FALSE)
        v$scRNAData <- Seurat::RunUMAP(v$scRNAData, dims = 1:30, verbose = FALSE)
        v$scRNAData <- Seurat::FindNeighbors(v$scRNAData, dims = 1:30, verbose = FALSE)
        v$scRNAData <- Seurat::FindClusters(v$scRNAData, verbose = FALSE)
        scRNA_plot <- Seurat::DimPlot(v$scRNAData, group.by = "subclass")
        print (scRNA_plot)
        #output$process_sc.done <- renderText(paste0("Processing of scRNA-seq data done!"))
        v$isProcess_scRNAdone <- TRUE
      })
    })

    output$scRNAPlot <- renderPlot({
      if(is.null(v$isProcess_scRNAdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating scRNA-seq Plot...", value=0, {
          DimPlot(v$scRNAData, group.by = "subclass", label = T)
        })
      }
    })

    observeEvent(input$vis_spRNA, {
      withProgress(message = "Processing scRNA-seq dataset...", value = 0,{
        incProgress(0.5, message = "Processing...")
        spRNA_plot <- SpatialDimPlot(v$scData, label = TRUE, label.size = 3)
        print (spRNA_plot)
        #output$vis_sp.done <- renderText(paste0("Processing of scRNA-seq data done!"))
        v$isVis_spRNAdone <- TRUE
      })
    })

    output$spRNAPlot <- renderPlot({
      if(is.null(v$isVis_spRNAdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating SpatialDim Plot...", value=0, {
          SpatialDimPlot(v$scData, label = TRUE, label.size = 3)
        })
      }
    })

    observeEvent(input$get_markers, {
      withProgress(message = "Getting marker genes...", value = 0,{
        incProgress(0.5, message = "Processing...")
        Seurat::Idents(object = v$scRNAData) <- v$scRNAData@meta.data$subclass
        ips.markers1 <- Seurat::FindAllMarkers(object = v$scRNAData,
                                                      assay = "SCT",
                                                      slot = "data",
                                                      verbose = TRUE,
                                                      only.pos = TRUE,
                                                      logfc.threshold = input$logfc,
                                                      min.pct = input$min_pct)
        v$ips.markers1 <- ips.markers1
        #output$vis_sp.done <- renderText(paste0("Processing of scRNA-seq data done!"))
      })
    })

    observeEvent(input$doDeconv, {
      withProgress(message = "Performing deconvolution...", value = 0,{
        incProgress(0.5, message = "Deconvoluting...")
        spotlight_ls <- spotlight_deconvolution(v$scRNAData,
                                                counts_spatial = v$scData@assays$Spatial@counts,
                                                clust_vr = "subclass",
                                                cluster_markers = v$ips.markers1,
                                                cl_n = 50,
                                                hvg = 3000,
                                                ntop = NULL,
                                                transf = "uv",
                                                method = "nsNMF",
                                                min_cont = 0.09)
        v$spotlight_ls <- spotlight_ls
        v$isDeconvdone <- TRUE
        decon_mtrx <- v$spotlight_ls[[2]]
        v$decon_mtrx <- decon_mtrx
        cell_types_all <- colnames(v$decon_mtrx )[which(colnames(v$decon_mtrx ) != "res_ss")]
        v$cell_types_all <- cell_types_all
        v$scData@meta.data <- cbind(v$scData@meta.data, v$decon_mtrx)
        img_path = "/acrc_raman/jinmiao/CJM_lab/Raman/Projects/hyperion_cytofkit2/spatial_shiny/tissue_lowres_image.png"
        v$img_path <- img_path
        deconv_plot <- SPOTlight::spatial_scatterpie(se_obj = v$scData,
                                      cell_types_all = v$cell_types_all,
                                      img_path = v$img_path,
                                      pie_scale = 0.4)
        print(deconv_plot)
        #output$vis_sp.done <- renderText(paste0("Processing of scRNA-seq data done!"))
      })
    })

    output$DeconvPlot <- renderPlot({
      if(is.null(v$isDeconvdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating SpatialDim Plot...", value=0, {
          SPOTlight::spatial_scatterpie(se_obj = v$scData,
                                        cell_types_all = v$cell_types_all,
                                        img_path = v$img_path,
                                        pie_scale = 0.4)
        })
      }
    })


    output$PCAplot_c <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData3.combined, reduction = "pca", label = T, group.by = 'orig.ident')
        })
      }
    })

    observeEvent(input$runUMAP, {
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
        v$scData3.combined <- FindNeighbors(v$scData3.combined, reduction = "pca", dims = 1:30)
        v$scData3.combined <- FindClusters(v$scData3.combined, resolution = 0.5)
        v$scData3.combined <- RunUMAP(v$scData3.combined, reduction = "pca", dims = 1:30)
        v$isUMAPdone1 <- TRUE
        UMAP_plot1a <- DimPlot(v$scData3.combined, reduction = "umap", label = T)
        UMAP_plot1b <- DimPlot(v$scData3.combined, reduction = "umap", label = T, group.by = 'batch')
        UMAP_plot1c <- DimPlot(v$scData3.combined, reduction = "umap", label = T, group.by = 'orig.ident')
        print(UMAP_plot1a)
        print(UMAP_plot1b)
        print(UMAP_plot1c)
      })
    })

    output$UMAPplot_a <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData3.combined, reduction = "umap", label = T)
        })
      }
    })

    output$UMAPplot_b <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData3.combined, reduction = "umap", label = T, group.by = 'batch')
        })
      }
    })

    output$UMAPplot_c <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData3.combined, reduction = "umap", label = T, group.by = 'orig.ident')
        })
      }
    })

    observeEvent(input$runTSNE, {
      withProgress(message = "Running TSNE...", value = 0,{
        incProgress(0.5, message = "Running TSNE...")
        v$scData3.combined <- RunTSNE(v$scData3.combined, reduction = "pca", dims = 1:30)
        v$isTSNEdone1 <- TRUE
        TSNE_plot1a <- DimPlot(v$scData3.combined, reduction = "tsne", label = T)
        TSNE_plot1b <- DimPlot(v$scData3.combined, reduction = "tsne", label = T, group.by = 'batch')
        TSNE_plot1c <- DimPlot(v$scData3.combined, reduction = "tsne", label = T, group.by = 'orig.ident')
        print(TSNE_plot1a)
        print(TSNE_plot1b)
        print(TSNE_plot1c)
      })
    })

    output$TSNEplot_a <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData3.combined, reduction = "tsne", label = T)
        })
      }
    })

    output$TSNEplot_b <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData3.combined, reduction = "tsne", label = T, group.by = 'batch')
        })
      }
    })

    output$TSNEplot_c <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData3.combined, reduction = "tsne", label = T, group.by = 'orig.ident')
        })
      }
    })

    output$ident.swap1 <- renderUI({
      if(is.null(v$scData3.combined)){
        return(NULL)
      }else{
        groupings1 <- names(v$scData3.combined@meta.data[,!names(v$scData3.combined@meta.data) %in% c("nGene", "nUMI", "percent.mito")])
        tagList(
          h4("Set current identity:"),
          fluidRow(
            column(3,
                   selectInput("active.ident1", label = NULL,
                               choices = groupings1)
            ),
            column(3,
                   actionButton("swap.ident1",label = NULL, icon = icon("arrow-right"))
            )
          )

        )
      }
    })

    observeEvent(input$PDFl, {
      if(!is.null(v$scData3.combined) ){
        withProgress(message="Downloading plot PDF files...", value=0, {
          print(getwd())
          pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
          if(!dir.exists(pdfDir)){
            dir.create(pdfDir)
          }
          filename2 <- paste0(pdfDir, .Platform$file.sep,"Data_integration_Seurat_", Sys.Date(), ".pdf")
          i = 0
          while(file.exists(filename2)){
            filename2 <- paste0(pdfDir, .Platform$file.sep,"Dimension_reduction_plots_Seurat_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
            i = i + 1;
          }
          PCA_plot1a <- DimPlot(v$scData3.combined, reduction = "pca", label = T)
          PCA_plot1b <- DimPlot(v$scData3.combined, reduction = "pca", label = T, group.by = 'type')
          PCA_plot1c <- DimPlot(v$scData3.combined, reduction = "pca", label = T, group.by = 'orig.ident')
          UMAP_plot1a <- DimPlot(v$scData3.combined, reduction = "umap", label = T)
          UMAP_plot1b <- DimPlot(v$scData3.combined, reduction = "umap", label = T, group.by = 'type')
          UMAP_plot1c <- DimPlot(v$scData3.combined, reduction = "umap", label = T, group.by = 'orig.ident')
          TSNE_plot1a <- DimPlot(v$scData3.combined, reduction = "tsne", label = T)
          TSNE_plot1b <- DimPlot(v$scData3.combined, reduction = "tsne", label = T, group.by = 'type')
          TSNE_plot1c <- DimPlot(v$scData3.combined, reduction = "tsne", label = T, group.by = 'orig.ident')
          prePlot()
          pdf(filename2,
              width=as.numeric(input$pdf_w),
              height=as.numeric(input$pdf_h))
          print(PCA_plot1a)
          print(PCA_plot1b)
          print(PCA_plot1c)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          print(UMAP_plot1c)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          dev.off()
        })
        withProgress(message="Downloading Dimension_reduction coordinates...", value=0.6, {
          print(getwd())
          pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
          if(!dir.exists(pdfDir)){
            dir.create(pdfDir)
          }
          filename2a <- paste0(pdfDir, .Platform$file.sep,"seurat_pca_", Sys.Date(), ".txt")
          filename2b <- paste0(pdfDir, .Platform$file.sep,"seurat_umap_", Sys.Date(), ".txt")
          filename2c <- paste0(pdfDir, .Platform$file.sep,"seurat_tsne_", Sys.Date(), ".txt")
          i = 0
          while(file.exists(filename2a)){
            filename2a <- paste0(pdfDir, .Platform$file.sep,"seurat_pca_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
            i = i + 1;
          }
          while(file.exists(filename2b)){
            filename2b <- paste0(pdfDir, .Platform$file.sep,"seurat_umap_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
            i = i + 1;
          }
          while(file.exists(filename2c)){
            filename2c <- paste0(pdfDir, .Platform$file.sep,"seurat_tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
            i = i + 1;
          }
          write.csv(v$scData3.combined@reductions$pca@cell.embeddings, file = filename2a)
          write.csv(v$scData3.combined@reductions$umap@cell.embeddings, file = filename2b)
          write.csv(v$scData3.combined@reductions$tsne@cell.embeddings, file = filename2c)
        })
      }
    })

    ##---------------Summary tab

    ##------Clean up when ending session----
    session$onSessionEnded(function(){
      prePlot()
    })
})


