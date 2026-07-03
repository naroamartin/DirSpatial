################################################################################
##  RHEA -- 2 x 3 PANELS               
##                                                                            
##  Rows:    leading hemisphere  /  trailing hemisphere                        
##  Columns: (1) craters, coloured and sized by diameter (log scale)         
##           (2) NW surface estimate with h_cv                                
##           (3) NW surface estimate with h_best (MCV)                                                                                             ##
################################################################################
rm(list = ls ())
library(DirStatsOld)                        # loc.directional.linear()
if (!requireNamespace("png", quietly = TRUE)) install.packages("png")
read_png <- function(path) png::readPNG(path)
load("rhea_results_log10.RData")

## Colour palette and a tiny value -> colour mapper.
colsbv <- colorRampPalette(c("royalblue", "lightblue", "orange", "red"))
col_cuts <- function(x, breaks) {
  colsbv(length(breaks) - 1)[cut(x, breaks, include.lowest = TRUE, 
                                 labels = FALSE)]
}

outdir <- "Rhea/Plots/"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

#-------------------------------------------------------------------------------
# CRATER TABLE (column 1)  
#-------------------------------------------------------------------------------
craters <- read.table("jgre20485-sup-0002-TableS1.txt",
                      skip = 2, sep = "\t", dec = ".")
names(craters) <- c("Center_Longitude", "Center_Latitude", "diam", "name")

craters$Center_Longitude <- 360 - craters$Center_Longitude
craters$phi <-  craters$Center_Latitude  / 180 * pi              # latitude
craters$theta <- (craters$Center_Longitude / 180 * pi) %% (2 * pi) # longitude

## Colour and size on the log-diameter scale (shared by both hemispheres).
logd <- log(craters$diam)
sca.diam <- range(logd)
brks.diam <- seq(sca.diam[1], sca.diam[2], length.out = 101)
cex_rng <- c(0.40, 2.2)    # min / max dot size

craters$col <- col_cuts(logd, brks.diam)
craters$cex <- cex_rng[1] + diff(cex_rng) * (logd - sca.diam[1]) / diff(sca.diam)

#-------------------------------------------------------------------------------
# SETTINGS FOR THE NW SURFACES (columns 2-3)
#-------------------------------------------------------------------------------
ng <- 180            # surface grid resolution (per side)
white_amt <- 0.30           # fade terrain toward white (0 = none, 1 = white)
surf_alpha <- 0.55           # surface opacity (0 = invisible, 1 = opaque)
cov_cex <- 0.5            # size of the white data dots

# range(nw_cv)
# [1] 1.172696 1.444811
# range(nw_best)
# 1.201928 1.349688

sca.surf  <- c(1.1, 1.5)      #for log10(Y) 
#sca.surf  <- c(2.5, 3.5).    #for log(Y)
brks.surf <- seq(sca.surf[1], sca.surf[2], length.out = 101)

bandwidths <- c(cv = h_cv, best = h_best)   #columns 2 and 3

#-------------------------------------------------------------------------------
## DRAW AND SAVE ONE FILE PER PANEL
#-------------------------------------------------------------------------------
for (type in c("leading", "trailing")) {
  
  shift <- if (type == "trailing") pi else 0
  
  ## -- read the hemisphere image and locate the moon's disk in it -----------
  img <- read_png(paste0(type, ".png"))            # array [H, W, channels]
  H <- dim(img)[1]; W <- dim(img)[2]
  
  lum  <- if (length(dim(img)) == 3) (img[, , 1] + img[, , 2] + img[, , 3]) / 3 else img
  mask <- lum < 0.98                               # non-white = the disk
  rr <- range(which(rowSums(mask) > 0))
  cc <- range(which(colSums(mask) > 0))
  cx <- mean(cc)                                 # disk centre (plot x)
  cy <- H - mean(rr) + 1                         # disk centre (plot y, up)
  R <- (diff(cc) + diff(rr)) / 4                # disk radius
  
  ## -- project the craters onto this hemisphere -----------------------------
  u <- cos(craters$theta - shift) * cos(craters$phi)
  v <- sin(craters$phi)
  vis <- sin(craters$theta - shift) * cos(craters$phi) < 0   # front hemisphere
  ord <- order(craters$diam[vis])                  # small first, big on top
  
  ## ======================= COLUMN 1: craters ===============================
  png(file.path(outdir, sprintf("rhea_%s_Y.png", type)),
      width = 5, height = 5, units = "in", res = 300)
  par(mar = rep(0.2, 4))
  plot(NA, xlim = c(1, W), ylim = c(1, H), asp = 1, axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  rasterImage(as.raster(img), 1, 1, W, H)
  points((cx + u[vis] * R)[ord], (cy + v[vis] * R)[ord],
         pch = 21, bg = craters$col[vis][ord], col = "black",
         lwd = 0.3, cex = craters$cex[vis][ord])
  dev.off()
  cat("Saved:", sprintf("rhea_%s_Y.png", type), "\n")
  
  ## ================== COLUMNS 2-3: NW surfaces =============================
  for (b in seq_along(bandwidths)) {
    
    h <- bandwidths[b]
    name <- names(bandwidths)[b]                   
    
    ## grid over the projected unit disk
    g <- seq(-1, 1, length.out = ng)
    U <- matrix(g,      ng, ng, byrow = TRUE)
    Vg <- matrix(rev(g), ng, ng, byrow = FALSE)
    ins <- as.vector(U^2 + Vg^2) <= 1
    d <- -sqrt(pmax(0, 1 - U^2 - Vg^2))          # front hemisphere (depth<0)
    
    ## un-rotate the screen frame back to the data frame
    if (type == "trailing") { gX1 <- -U; gX2 <- -d } else { gX1 <- U; gX2 <- d }
    gX3 <- Vg
    dir_grid <- cbind(as.vector(gX1)[ins],
                      as.vector(gX2)[ins],
                      as.vector(gX3)[ins])
    
    ## Nadaraya-Watson estimate on the grid, mapped to colours
    zhat <- as.vector(loc.directional.linear(x = dir_grid, data.dir = X,
                                             data.lin = Y, h = h, p = 0)$Yhat)
    zcol <- rep(NA_character_, ng * ng)
    zcol[ins] <- adjustcolor(
      col_cuts(pmin(pmax(zhat, sca.surf[1]), sca.surf[2]), brks.surf),
      alpha.f = surf_alpha)
    ras <- as.raster(matrix(zcol, ng, ng))
    
    ## draw: faded terrain, then the surface, then the data locations
    png(file.path(outdir, sprintf("rhea_%s_nw_%s.png", type, name)),
        width = 5, height = 5, units = "in", res = 300)
    par(mar = rep(0.2, 4))
    plot(NA, xlim = c(1, W), ylim = c(1, H), asp = 1, axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    rasterImage(as.raster(img * (1 - white_amt) + white_amt), 1, 1, W, H)
    rasterImage(ras, cx - R, cy - R, cx + R, cy + R, interpolate = TRUE)
    points(cx + u[vis] * R, cy + v[vis] * R, pch = 16, cex = cov_cex,
           col = adjustcolor("white", alpha.f = 0.75))
    dev.off()
    cat("Saved:", sprintf("rhea_%s_nw_%s.png", type, name), "\n")
  }
}


#-------------------------------------------------------------------------------
# COLOUR BARS  
#-------------------------------------------------------------------------------
color_bar <- function(breaks, ticks, labels = ticks,
                      pal = colorRampPalette(c("blue", "white", "red")), ...) {
  par(mar = c(2, 3, 2, 3))
  n_cols <- length(breaks) - 1
  cols   <- pal(n_cols)
  image(x = seq_along(cols), y = 1, z = matrix(seq_along(cols), ncol = 1),
        col = cols, axes = FALSE, xlab = "", ylab = "", ...)
  tick_positions <- (ticks - min(breaks)) / (max(breaks) - min(breaks)) *
    (n_cols - 1) + 1
  axis(1, at = tick_positions, labels = labels, lwd = 4, cex.axis = 2)
  box(lwd = 4)
}

## -- diameter (column 1): natural-log colors -----------
dticks <- c(10, 20, 50, 100, 200, 400)
dticks <- dticks[log(dticks) >= sca.diam[1] & log(dticks) <= sca.diam[2]]
png(file.path(outdir, "rhea_colorbar_diam.png"),
    width = 10, height = 1.9, units = "in", res = 200, pointsize = 18)
color_bar(breaks = brks.diam, ticks = log(dticks), labels = dticks, pal = colsbv)
dev.off()

## -- NW surface (columns 2-3): log10 colors------------
nwticks <- c(13, 15, 20, 25, 30)
nwticks <- nwticks[log10(nwticks) >= sca.surf[1] & log10(nwticks) <= sca.surf[2]]
png(file.path(outdir, "rhea_colorbar_nw.png"),
    width = 10, height = 1.9, units = "in", res = 200, pointsize = 18)
color_bar(breaks = brks.surf, ticks = log10(nwticks), labels = nwticks, pal = colsbv)
dev.off()