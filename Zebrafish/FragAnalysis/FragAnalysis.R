## FRAGMAN functions
# storing.inds() = read fsa files in folder and store them as list
# ladder.info.attach() = match sizes according to ladder
# overview2 = plotting
# score.markers = find peaks; can be automated

#ladder.corrector = allows to manually click on peaks.
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
BiocManager::install(c("Fragman"))
library("Fragman")
library(ggplot2)
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# Clear environment
rm(list=ls()); + try(dev.off(),silent=TRUE);
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# define folder
dir.fsa = "/Users/angueyraaristjm/Documents/LiMolec/zfGenotyping/20200629_fragAnalysis/gnat2"
# load all fsa files in folder
fsaData = storing.inds(dir.fsa)
fsaNames = names(fsaData)
cat(fsaNames, sep="\n") 
# define ladder
# LIZ500 (https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_042491.pdf)
# liz500 <- c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500) 
# finding that first marker is usually contaminated, so decided to remove it
liz500 <- c(50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500) 

# # plot all (takes a while) (remember this is not trimmed yet)
# plot.fsa_stored(fsaData[1], lay=c(ceiling(length(fsaNames)/2),3), channel = c(1), cex.legend=1)
# # plot one
# i= 1;
# plot.fsa_stored(fsaData[i], channel = c(1,5), cex.legend=1)
# try(dev.off(),silent=TRUE);
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# get single file
i=7
tempName = fsaNames[i] 
message(paste("Analyzing:",fsaNames[i]))
tempData <- fsaData[tempName]
class(tempData) <- "fsa_stored"
# check if this has been done
if (file.exists(paste(paste(dir.fsa,gsub('.{0,4}$', '', tempName),sep = "/"),".csv", sep=""))) {
   message("Already calibrated and exported; no need to redo")
}
# plto data to assess if it's worth it remapping
# plot(tempData[[tempName]][,1], typ='l',xlim=c(1700,length(tempData[[tempName]][,1])),ylim=c(0,30000))
plot(tempData[[tempName]][,1], typ='l',xlim=c(1800,3000),ylim=c(min(tempData[[tempName]][1800:3000,1]),max(tempData[[tempName]][1800:3000,1])), main = tempName)
# figure out threshold by checking liz500 channel
plot(tempData[[tempName]][1500:2500,5], typ='l')
ilim01 = 1650;
plot(tempData[[tempName]][ilim01:5000,5], typ='l')
ilim02 = 5000;
plot(tempData[[tempName]][ilim01:ilim02,5], typ='l') # 16 peaks for liz500 (15 peaks if removed '35' marker)

# for fragment analysis using FAM13 and LIZ500, relevant channels are 1 and 5 respectively
tempData[[tempName]] = tempData[[tempName]][ilim01:ilim02,c(1,5)]



guessThreshold = 400
# match ladder
# works better if higher values when noise is high
ladderData = ladder.info.attach(stored=tempData, ladder=liz500, method='iter2', draw=TRUE, ladd.init.thresh=guessThreshold)
# replot ladder if needed to play with threshold
# plot(tempData[[tempName]][,2], typ='l') # 16 peaks for liz500
#run correction if needed
# ladderData = ladder.corrector(stored=tempData, to.correct = tempName, ladder=liz500)

# ladder is stored as an environment (hidden) variable ->      list.data.covarrubias[[tempName]]
ladder = data.frame("p" = list.data.covarrubias[[tempName]]$pos, "w"=list.data.covarrubias[[tempName]]$wei)
# if all fails, do completely manual mapping
# plot(tempData[[tempName]][,2], typ='l', xlim = c(500,1000))
# ladder = data.frame("p" = c(20,170,310,550,600,650,880,1150,1450,1680,1740,2040,2310,2550,2590), "w"=liz500)
ladder


# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# fit ladder with a 5th degree polynomial (this is what fragman uses)
polyModel = lm(w ~ poly(p,5), data = ladder)

fitWeights<- predict(polyModel,ladder,interval='confidence',level=0.99)
plot(list.data.covarrubias[[tempName]]$pos,list.data.covarrubias[[tempName]]$wei) + 
   points(list.data.covarrubias[[tempName]]$pos,fitWeights[,1], typ='l') + 
   points(list.data.covarrubias[[tempName]]$pos,fitWeights[,1], pch=20)


# export data as csv into same folder and with same name
# write.csv(data.frame("size"=fitWeights,"fluo"=tempData[[tempName]][,1],"ladder"=tempData[[tempName]][,1]),paste(paste(dir.fsa,gsub('.{0,4}$', '', tempName),sep = "/"),".csv", sep=""),row.names = FALSE)


# check remapped data
full_ladder = data.frame("p"=1:length(tempData[[tempName]][,1]))
fitWeights <- predict(polyModel,full_ladder)
# plot the data
plot(fitWeights, tempData[[tempName]][,1], typ='l', xlim=c(0, 600), main=tempName)
# zoom into ROI
p_lo = 100;
p_hi = 550;
# plot(fitWeights, tempData[[tempName]][,1], typ='l', xlim=c(p_lo, p_hi), ylim=c(0,1000), main=tempName)
plot(fitWeights, tempData[[tempName]][,1], typ='l', xlim=c(p_lo, p_hi), 
     ylim=c(min(tempData[[tempName]][fitWeights>p_lo&fitWeights<p_hi,1]),1.1*max(tempData[[tempName]][fitWeights>p_lo&fitWeights<p_hi,1])), main=tempName)
# add ladder channel
points(fitWeights, tempData[[tempName]][,2], typ='l', col=rgb(1, 0, 0, .5))

# or plot whole ladder by itself
# plot(fitWeights, tempData[[tempName]][,2], typ='l', xlim=c(0, 100))




# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# # This would be the rest of the automated peak analysis, but it's not what I ultimately want to do
# # try automatic peak detection or build a range of expected peaks (a "panel")
# expectedPeaks = overview2(my.inds=tempData, channel = 1, ladder=liz500, init.thresh=1100)
# # if not click on peaks
# expectedPeaks = locator(type="p", pch=20, col="red")$channel_1
# # or just define manually
# expectedPeaks$channel_1 = 335
# #analyze peaks
# results <- score.markers(my.inds=tempData, channel = 1, panel=expectedPeaks$channel_1, ladder=liz500, electro=FALSE)
# # create summary
# finalOut = get.scores(results)
# finalOut
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------









