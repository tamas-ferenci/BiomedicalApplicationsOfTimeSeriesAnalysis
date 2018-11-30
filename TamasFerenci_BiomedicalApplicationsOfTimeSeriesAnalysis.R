## ----setup, include=FALSE------------------------------------------------
library( lattice )
library( data.table )
library( rms )
knitr::opts_chunk$set( cache = TRUE )

## ------------------------------------------------------------------------
ts( rnorm( 20 ), frequency = 4, start = c( 2010, 2 ) )
ts( rnorm( 30 ), frequency = 12, start = c( 2010, 2 ) )

## ---- fig.height=5-------------------------------------------------------
ldeaths
plot( ldeaths )

## ---- fig.height = 5.5---------------------------------------------------
SimDataFourier <- data.frame( t = 1:1000 )
SimDataFourier <- transform( SimDataFourier, y = 0.5*sin( t*2 ) + sin( t/10*2 ) +
                               rnorm( length( t ), 0, 0.1 ) )
xyplot( y ~ t, data = SimDataFourier, type = "l", xlim = c( 0, 200 ) )

## ------------------------------------------------------------------------
xyplot( spec ~ freq, data = spectrum( SimDataFourier$y, plot = FALSE ), type = "l",
        scales = list( y = list( log = 10 ) ) )

## ------------------------------------------------------------------------
locmaxpanel <- function( x, y, width, maxmeddiff = 1, rounddigit = 2, ... ) {
  if( width%%2==0 )
    width <- width+1
  panel.xyplot( x, y, ... )
  maxs <- zoo::rollapply( y, width, function(x) (which.max(x)==(width+1)/2)&
                            (max(x)-median(x)>maxmeddiff),
                          align = "center", fill = NA )
  panel.abline( v = x[ maxs ], col = "gray" )
  panel.points( x[ maxs ], y[ maxs ] )
  panel.text( x[ maxs ], y[ maxs ], round( x[ maxs ], rounddigit ), pos = 4 )
}

## ------------------------------------------------------------------------
xyplot( spec ~ freq, data = spectrum( SimDataFourier$y, plot = FALSE ), type = "l",
        scales = list( y = list( log = 10 ) ), panel = locmaxpanel, width = 21, maxmeddiff = 2 )

## ----message=FALSE,warning=FALSE,fig.height=5.5--------------------------
## require( tuneR ) ## require( pastecs ) ## devtools::install_github( "mkfs/r-physionet-ptb" )
## https://www.physionet.org/physiobank/database/ptbdb/
## system2( system.file( "exec", "download_ptb.sh", package = "r.physionet.ptb" ) )
## system2( system.file( "exec", "ptb_patient_to_json.rb", package = "r.physionet.ptb" ),
## args="patient001" )
library( r.physionet.ptb )
ptb <- r.physionet.ptb::ptb.from.file( "patient001.json" )
ptbecg <- r.physionet.ptb::ptb.extract.lead( ptb, "i" )$`1-10010`
xyplot( ptbecg~seq_along( ptbecg ), type = "l", xlim = c( 0, 5000 ), xlab = "Time", ylab = "" )

## ---- message=FALSE,warning=FALSE----------------------------------------
xyplot( spec ~ freq, data = spectrum( ptbecg, plot = FALSE, span = rep( 201, 3 ) ), type = "l",
        scales = list( y = list( log = 10 ) ), panel = locmaxpanel, width = 21,
        maxmeddiff = 2e-4 )

## ---- fig.height = 6-----------------------------------------------------
SimDataWavelet <- data.frame( t = 1:2000 )
SimDataWavelet <- transform( SimDataWavelet,
                             y = ifelse( t<=1000, sin( t*2 ), sin( t/10*2 ) ) +
                               rnorm( length( t ), 0, 0.1 ) )
xyplot( y ~ t, data = SimDataWavelet, type = "l" )

## ------------------------------------------------------------------------
xyplot( spec ~ freq, data = spectrum( SimDataWavelet$y, plot = FALSE ), type = "l",
        scales = list( y = list( log = 10 ) ) )

## ------------------------------------------------------------------------
WaveletComp::wt.image( WaveletComp::analyze.wavelet( SimDataWavelet, "y",
                                                     verbose = FALSE, make.pval = FALSE ) )

## ---- fig.height = 6-----------------------------------------------------
SimDataWavelet <- data.frame( t = 1:2000 )
SimDataWavelet <- transform( SimDataWavelet,
                             y = WaveletComp::periodic.series( start.period = 20,
                                                               end.period = 200,
                                                  length = length( t ) ) +
                               0.1*rnorm( length( t ) ) )
xyplot( y ~ t, data = SimDataWavelet, type = "l" )

## ------------------------------------------------------------------------
xyplot( spec ~ freq, data = spectrum( SimDataWavelet$y, plot = FALSE ), type = "l",
        scales = list( y = list( log = 10 ) ) )

## ------------------------------------------------------------------------
WaveletComp::wt.image( WaveletComp::analyze.wavelet( SimDataWavelet, "y",
                                                     verbose = FALSE, make.pval = FALSE ) )

## ------------------------------------------------------------------------
tmpfile <- tempfile( fileext = ".xlsx" )
download.file( url = paste0( "https://www.gov.uk/government/uploads/system/uploads/",
                             "attachment_data/file/339410/NoidsHistoricAnnualTotals.xlsx" ),
               destfile = tmpfile, mode = "wb" )
res1 <- XLConnect::loadWorkbook( tmpfile )
XLConnect::setMissingValue( res1, value = c( "*" ) )
res1 <- do.call( plyr::rbind.fill, lapply( XLConnect::getSheets( res1 ), function( s ) {
  temp <- XLConnect::readWorksheet( res1, sheet = s, startRow = 4 )
  temp <- temp[ , grep( "Disease", colnames( temp ) ):ncol( temp ) ]
  temp <- temp[ 1:( if( sum( is.na( temp$Disease ) )==0 ) nrow( temp ) else
    which( is.na( temp$Disease ) )[ 1 ]-1 ), ]
  for( i in 2:ncol( temp ) )
    temp[ , i ] <- as.numeric( gsub( "[[:space:].,‡†]", "", temp[ , i ] ) )
  temp2 <- as.data.frame( t( temp[ , - 1 ] ) )
  colnames( temp2 ) <- temp[ , 1 ]
  temp2$Year <- as.numeric( substring( rownames( temp2 ), 2, 5 ) )
  temp2
} ) )
unlink( tmpfile )

## ------------------------------------------------------------------------
tmpfile <- tempfile( fileext = ".xlsx" )
download.file( url = paste0( "https://www.gov.uk/government/uploads/system/uploads/",
                             "attachment_data/file/664864/",
                             "Annual_totals_from_1982_to_2016.xlsx" ),
               destfile = tmpfile, mode = "wb" )
res2 <- XLConnect::loadWorkbook( tmpfile )
XLConnect::setMissingValue( res2, value = c( "--" ) )
res2 <- do.call( plyr::rbind.fill, lapply( XLConnect::getSheets( res2 )[ -1 ], function( s ) {
  temp <- XLConnect::readWorksheet( res2, sheet = s, startRow = 5 )
  temp <- temp[ 1:( nrow( temp )-1 ), ]
  temp2 <- as.data.frame( t( temp[ , - 1 ] ) )
  colnames( temp2 ) <- temp[ , 1 ]
  temp2$Year <- as.numeric( substring( rownames( temp2 ), 2, 5 ) )
  temp2
} ) )
unlink( tmpfile )

## ------------------------------------------------------------------------
tmpfile <- tempfile( fileext = ".xls" )
download.file( url = paste0( "https://www.ons.gov.uk/file?uri=/",
                             "peoplepopulationandcommunity/populationandmigration/",
                             "populationestimates/adhocs/",
                             "004358englandandwalespopulationestimates1838to2014/",
                             "englandandwalespopulationestimates18382014tcm77409914.xls" ),
               destfile = tmpfile, mode = "wb" )
res3 <- XLConnect::readWorksheetFromFile( tmpfile, sheet = "EW Total Pop 1838-2014", startRow = 2,
                               endRow = 179 )
unlink( tmpfile )
names( res3 )[ 1 ] <- "Year"
res3$Persons <- ifelse( res3$Persons < 100000, res3$Persons*1000, res3$Persons )
res3 <- res3[ , c( "Year", "Persons" ) ]
res4 <- read.csv( paste0( "https://www.ons.gov.uk/generator?format=csv&uri=/",
                          "peoplepopulationandcommunity/populationandmigration/",
                          "populationestimates/timeseries/ewpop/pop" ), skip = 7 )
names( res4 ) <- c( "Year", "Persons" )
res4 <- res4[ res4$Year>=2015, ]
UKEpid <- merge( plyr::rbind.fill( res1, res2 ), rbind( res3, res4 ) )
UKPertussis <- UKEpid[ , c( "Year", "Whooping cough", "Persons" ) ]
UKPertussis$Inc <- UKPertussis$`Whooping cough`/UKPertussis$Persons*100000
UKPertussis <- UKPertussis[ !is.na( UKPertussis$`Whooping cough` ), ]

## ------------------------------------------------------------------------
xyplot( Inc ~ Year, data = UKPertussis, type = "l", ylab = "Incidence [/100 000/year]" )

## ------------------------------------------------------------------------
WaveletComp::wt.image( WaveletComp::analyze.wavelet( UKPertussis, "Inc",
                                                     verbose = FALSE, make.pval = FALSE ) )

## ---- fig.height = 6-----------------------------------------------------
data( "CVDdaily", package = "season" )
rownames( CVDdaily ) <- NULL
xyplot( cvd ~ date, data = CVDdaily, type = "l", xlab = "Time", ylab = "Number of deaths" )

## ------------------------------------------------------------------------
CVDdaily$year <- lubridate::year( CVDdaily$date )
CVDdaily$wday <- as.factor( lubridate::wday( CVDdaily$date, week_start = 1 ) )
CVDdaily$yday <- lubridate::yday( CVDdaily$date )/yearDays( CVDdaily$date )
head( CVDdaily[ , c( "date", "year", "wday", "yday", "cvd" ) ] )

## ---- message=FALSE------------------------------------------------------
library( mgcv )
fit <- gam( cvd ~ s( as.numeric( date ) ) + wday + s( yday, bs = "cc" ), data = CVDdaily,
            family = nb( link = log ) )
summary( fit )

## ---- fig.height = 6-----------------------------------------------------
plot( fit, select = 1, scale = 0, rug = FALSE, trans = exp, shade = TRUE,
      col = trellis.par.get()$superpose.line$col[1], xaxt = "n", xlab = "Year", ylab = "IRR" )
axis( 1, at = seq( CVDdaily$date[1], tail( CVDdaily$date, 1 )+1, by = "year" ),
      labels = year( CVDdaily$date[1] ):year( tail( CVDdaily$date, 1 )+1 ) )

## ---- fig.height = 6-----------------------------------------------------
plot( fit, select = 2, scale = 0, rug = FALSE, trans = exp, shade = TRUE,
      col = trellis.par.get()$superpose.line$col[1], xaxt = "n", xlab = "Day of year",
      ylab = "IRR" )
axis( 1, at = seq( 0, 1, 1/12 ), labels = seq( 0, 1, 1/12 )*30*12 )

## ---- fig.height = 6-----------------------------------------------------
source( "https://pastebin.com/raw/hBmStX4Y" )
termplot2( fit, terms = "wday", se = TRUE, yscale = "exponential",
           col.term = trellis.par.get()$superpose.line$col[1],
           col.se = "gray80", se.type = "polygon", xlab = "Day of week" )

## ---- fig.height = 4-----------------------------------------------------
do.call( gridExtra::grid.arrange, lapply( c( 2, 6, 12, 24 ), function( o ) {
  xyplot( y ~ t, groups = grp, data = rbind( data.frame( grp = "data", SimDataFourier ),
                                             data.frame( grp = "smooth", t = SimDataFourier$t,
                                                         y = forecast::ma( SimDataFourier$y,
                                                                           o ) ) ),
          type = "l", xlim = c( 0, 200 ), main = paste0( "Order: ", o ) )
} ) )

## ------------------------------------------------------------------------
plot( decompose( ldeaths ) )

## ------------------------------------------------------------------------
plot( stl( ldeaths, s.window = "periodic" ) )

## ------------------------------------------------------------------------
plot( forecast::seasadj( stl( ldeaths, s.window = "periodic" ) ) )

## ---- fig.height = 5.5---------------------------------------------------
TLCData <- read.table( "https://content.sph.harvard.edu/fitzmaur/ala2e/tlc-data.txt",
                       col.names = c( "ID", "Trt", paste0( "Wk", c( 0, 1, 4, 6 ) ) ) )
TLCData <- reshape( TLCData, varying = paste0( "Wk", c( 0, 1, 4, 6 ) ), v.names = "LeadLevel",
                    timevar = "Week", times = c( 0, 1, 4, 6 ), idvar = "ID", direction = "long" )
TLCData$Trt <- relevel( TLCData$Trt, ref = "P" )
TLCData$Week.f <- as.factor( TLCData$Week )
xyplot( LeadLevel ~ Week | Trt, groups = ID, data = TLCData, type = "b" )

## ---- fig.height = 5.5---------------------------------------------------
TLCData <- data.table( TLCData )
TLCData$Time <- as.numeric( TLCData$Week.f )
dd <- datadist( TLCData )
options( datadist = "dd" )
xYplot( Cbind( mean, lwr, upr ) ~ Week, groups = Trt, type = "b",
        data = TLCData[ , .( mean = mean( LeadLevel ),
                             lwr = t.test( LeadLevel )$conf.int[1],
                             upr = t.test( LeadLevel )$conf.int[2] ) , .( Trt, Week ) ],
        ylim = c( 10, 30 ), ylab = "Mean lead level" )

## ------------------------------------------------------------------------
ols( LeadLevel ~ Week.f*Trt, data = TLCData )

## ------------------------------------------------------------------------
fit <- Gls( LeadLevel ~ Week.f*Trt, data = TLCData, corr = nlme::corSymm( form = ~ Time | ID ),
     weights = nlme::varIdent( form = ~ 1 | Week.f ) )
fit

## ------------------------------------------------------------------------
temp <- Predict( fit, Trt, Week.f )
temp$Week.f <- as.numeric( levels( temp$Week.f ) )[ temp$Week.f ]
xYplot( Cbind( yhat, lower, upper ) ~ Week.f, groups = Trt, data = temp, type = "b",
        ylim = c( 10, 30 ) )

## ------------------------------------------------------------------------
data( "Orthodont", package = "nlme" )
OrthoFem <- Orthodont[ Orthodont$Sex=="Female", ]
plot( OrthoFem )

## ------------------------------------------------------------------------
fit1 <- nlme::lmList( distance ~ I( age - 11 ), data = OrthoFem )
plot( nlme::intervals( fit1 ) )

## ---- fig.height = 5.5---------------------------------------------------
fit2 <- nlme::lme( distance ~ age, data = OrthoFem, random = ~1|Subject )
xyplot( distance + fitted( fit1 ) + fitted( fit2, level = 1 ) +
          fitted( fit2, level = 0 ) ~ age | Subject, data = OrthoFem,
        type = c( "p", "l", "l", "l" ), distribute.type = TRUE, ylab = "Distance", grid = TRUE,
        auto.key = list( text = c( "Measured", "Fitted (individual models)",
                                   "Fitted (mixed effects)", "Fitted (fixed effects only)" ),
                         columns = 4, points = FALSE, lines = TRUE ) )

