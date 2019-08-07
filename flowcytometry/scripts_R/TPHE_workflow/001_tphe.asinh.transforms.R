## asinh.transform object building

library(flowCore)
## read in 'target/representative' .fcs file; compensate
dat <- fset.compensate(read.flowSet("./data_source/TPHE2_with.spillover/NGL091/HD0189_TPHE_NGL091 TPHE.fcs", transformation = FALSE))
## fluor channels
channels <- colnames(dat)
channels <- channels[-grep("FSC|SSC|Time", channels)]

##initialize named-list to store transformation parameters
asinh.transforms.TPHE <- vector(mode = "list", length = length(channels))
names(asinh.transforms.TPHE) <- channels

##build transforms per channel; interactively test/visualize using 'asinh.viz' shiny function
library(shiny)
library(ggplot2)

asinh.viz(data.frame(exprs(dat[[1]]))[, colnames(dat) %in% channels])
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Blue A"
channels[grep(channel.name, channels)]

a = 0
b = 0.2/700
c = 1
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = 0,
                                                   b = 0.2/700,
                                                   c = 1
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Blue B"
channels[grep(channel.name, channels)]

a = 0
b = 1/150
c = 2
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = 0,
                                                   b = 1/150,
                                                   c = 2
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Green A"
channels[grep(channel.name, channels)]

a = 0
b = 1/1000
c = 2
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = 0,
                                                   b = 1/1000,
                                                   c = 2
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Green C"
channels[grep(channel.name, channels)]

a = -1.5
b = 1/400
c = 3
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = -1.5,
                                                   b = 1/400,
                                                   c = 3
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Green D"
channels[grep(channel.name, channels)]

a = -1
b = 1/450
c = 2
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = -1,
                                                   b = 1/450,
                                                   c = 2
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Green E"
channels[grep(channel.name, channels)]

a = -1
b = 1/400
c = 2
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = -1,
                                                   b = 1/400,
                                                   c = 2
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Red A"
channels[grep(channel.name, channels)]

a = 0
b = 1/400
c = 2
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = 0,
                                                   b = 1/400,
                                                   c = 2
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Red B"
channels[grep(channel.name, channels)]

a = -1
b = 1/400
c = 3
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = -1,
                                                   b = 1/400,
                                                   c = 3
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Red C"
channels[grep(channel.name, channels)]

a = -1
b = 1/300
c = 2
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = -1,
                                                   b = 1/300,
                                                   c = 2
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Violet A"
channels[grep(channel.name, channels)]

a = 0
b = 1/450
c = 1.5
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = 0,
                                                   b = 1/450,
                                                   c = 1.5
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Violet B"
channels[grep(channel.name, channels)]

a = -1
b = 1/400
c = 3
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = -1,
                                                   b = 1/400,
                                                   c = 3
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Violet C"
channels[grep(channel.name, channels)]

a = 0
b = 1/400
c = 1
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = 0,
                                                   b = 1/400,
                                                   c = 1
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Violet D"
channels[grep(channel.name, channels)]

a = -2
b = 1/400
c = 2.5
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = -2,
                                                   b = 1/400,
                                                   c = 2.5
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Violet E"
channels[grep(channel.name, channels)]

a = -2
b = 1/400
c = 2.5
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = -2,
                                                   b = 1/400,
                                                   c = 2.5
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Violet G"
channels[grep(channel.name, channels)]

a = -3
b = 1/400
c = 2.5
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = -3,
                                                   b = 1/400,
                                                   c = 2.5
                                  )
  ),
  asinh.vals = asinh.vals
)
##################################################################################################################################################################################
##################################################################################################################################################################################
channel.name <- "Violet H"
channels[grep(channel.name, channels)]

a = -0.5
b = 1/250
c = 1
asinh.vals <- list(a = a, b = b, c = c)

asinh.transforms.TPHE[[channels[grep(channel.name, channels)]]] <- list(
  asinh.translist = transformList(channels[grep(channel.name, channels)],
                                  arcsinhTransform(transformationId = "arcsinh-transformation",
                                                   a = -0.5,
                                                   b = 1/250,
                                                   c = 1
                                  )
  ),
  asinh.vals = asinh.vals
)

saveRDS(asinh.transforms.TPHE, "./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")