BiocManager::install(c("Fragman"))
library("Fragman")
library(ggplot2)


# storing.inds() = read fsa files in folder and store them as list
# ladder.info.attach() = match sizes according to ladder
# overview2 = plotting
# score.markers = find peaks; can be automated

#ladder.corrector = allows to manually click on peaks.

rm(list=ls());
try(dev.off(),silent=TRUE);

folder <- "/Users/angueyraaristjm/Documents/Python/PythonCode/fatools-master"
testdata <- storing.inds(folder)
### here we just load our sample data and use the first 2 plants
?testdata
# this doesn't work: data(testdata)
# get a sublist of fsa objects
exData <- testdata[1]
class(exData) <- "fsa_stored"
# plot(exData) # to visualize

#overview plot of all fsa files
# plot(x, lay=c(2,1), channel=NULL, cex.legend=.5, ncol.legend=4,lims=NULL, color=NULL, ...)
plot.fsa_stored(exData, channel = c(1,2), cex.legend=1)

# figure out threshold 
plot(exData$test.fsa[0:2000,5], typ='l')
plot(exData$test.fsa[0:nrow(exData$test.fsa),5], typ='l')
ilim01 = 1550;
plot(exData$test.fsa[ilim01:3000,5], typ='l')
ilim02 = 4200;
plot(exData$test.fsa[ilim01:ilim02,5], typ='l') # 16 peaks for liz500

# for fragment analysis using FAM13 and LIZ500, relevant channels are 1 and 5 respectively
exData$test.fsa = exData$test.fsa[ilim01:ilim02,c(1,5)]

# match ladder
#LIZ500 (https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_042491.pdf)
liz500 <- c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500) 
ladderData = ladder.info.attach(stored=exData, ladder=liz500, method='iter2', draw=TRUE, ladd.init.thresh=250)

#run correction if needed
ladderData = ladder.corrector(stored=exData, to.correct = "test.fsa", ladder=liz500)

# ladder is stored as an environment (hidden) variable
list.data.covarrubias$test.fsa

# try automatic peak detection or build a range of expected peaks (a "panel")
expectedPeaks = overview2(my.inds=exData, channel = 1, ladder=liz500, init.thresh=1100)
# if not click on peaks
expectedPeaks = locator(type="p", pch=20, col="red")$channel_1
# or just define manually
expectedPeaks$channel_1 = 335



#analyze peaks
results <- score.markers(my.inds=exData, channel = 1, panel=expectedPeaks$channel_1, ladder=liz500, electro=FALSE)
# create summary
finalOut = get.scores(results)
finalOut


# # fitting ladder with a line (seems somewhat unreliable)
# plot(list.data.covarrubias$test.fsa$pos,list.data.covarrubias$test.fsa$wei)
# linearModel = lm(list.data.covarrubias$test.fsa$wei ~ list.data.covarrubias$test.fsa$pos)
# abline(linearModel,cex = 1.3,pch = 16)
# weightX = (linearModel$coefficients[2] * list.data.covarrubias$test.fsa$pos) + linearModel$coefficients[1]
# points(list.data.covarrubias$test.fsa$pos,weightX, typ='l')
# MW = (linearModel$coefficients[2] * 1:length(exData$test.fsa[,1])) + linearModel$coefficients[1]
# plot(MW, exData$test.fsa[,1] ,typ='l',  xlim=c(300, 400))



# fitting ladder with a 5th degree polynomial (this is what fragman uses)

ladder = data.frame("p" = list.data.covarrubias$test.fsa$pos, "w"=list.data.covarrubias$test.fsa$wei)
polyModel = lm(w ~ poly(p,5), data = ladder)

predicted.intervals <- predict(polyModel,ladder,interval='confidence',level=0.99)
plot(list.data.covarrubias$test.fsa$pos,list.data.covarrubias$test.fsa$wei)
points(list.data.covarrubias$test.fsa$pos,predicted.intervals[,1], typ='l')
points(list.data.covarrubias$test.fsa$pos,predicted.intervals[,1], pch=20)

full_ladder = data.frame("p"=1:length(exData$test.fsa[,1]))
predicted.intervals <- predict(polyModel,full_ladder)
plot(predicted.intervals, exData$test.fsa[,1], typ='l', xlim=c(330, 340))
# plot(predicted.intervals, exData$test.fsa[,2], typ='l', xlim=c(0, 100))





