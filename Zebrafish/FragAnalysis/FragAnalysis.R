## FRAGMAN functions
# storing.inds() = read fsa files in folder and store them as list
# ladder.info.attach() = match sizes according to ladder
# overview2 = plotting
# score.markers = find peaks; can be automated

#ladder.corrector = allows to manually click on peaks.
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# BiocManager::install(c("Fragman"))
library("Fragman")
library(ggplot2)
library(plotly)
library(ggpubr) 
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# Clear environment
rm(list=ls()); 
try(dev.off(),silent=TRUE);
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# define folder
dir.fsa = "/Users/angueyraaristjm/Documents/LiMolec/zfGenotyping/20210804_tbx2bfoxq2F0s/foxq2"
# load all fsa files in folder
fsaData = storing.inds(dir.fsa)
fsaNames = names(fsaData)
cat(gsub(".fsa","",fsaNames), sep="\n") 
# define ladder
# ROX400 = c(50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380, 400 );
# LIZ500 (https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_042491.pdf)
# liz500 <- c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500)
# finding that first marker is usually contaminated, so decided to remove it
# liz500 <- c(50, 75, 100, 139, 150, 160, 200, 300, 340, 350, 400, 450, 490, 500)
# liz500 <- c(50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500)
liz500 <- c(75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500)

# liz500 <- c(50, 75, 100, 340, 350, 400, 450, 490, 500)

# # plot all (takes a while) (remember this is not trimmed yet)
# plot.fsa_stored(fsaData[1], lay=c(ceiling(length(fsaNames)/2),3), channel = c(1), cex.legend=1)
# # plot one
i= 0;
# plot.fsa_stored(fsaData[i], channel = c(1,5), cex.legend=1)
# try(dev.off(),silent=TRUE);
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------- -----------------------------------------------
# expected peaks: gnat2 = 295bp;     syt5a (FiiRii) = 477bp;     efna1b = 495bp; 
#                 eml1 = 227bp (mut = -11bp); 
#                 sema7a = 358 bp;   tbx2a = 487 bp;   tbx2b = 332 bp
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# get single file (run section with command+alt+T)
i=i-1
tempName = fsaNames[i]; message(paste("Analyzing:",fsaNames[i]))
tempData <- fsaData[tempName] 
class(tempData) <- "fsa_stored"
# check if this has been done
if (file.exists(paste(paste(dir.fsa,gsub('.{0,4}$', '', tempName),sep = "/"),".csv", sep=""))) {message("Already calibrated and exported; no need to redo")}
# for fragment analysis using FAM13 and LIZ500, relevant channels are 1 and 5 respectively
chDNA = 1; chLadder = 10;
# plot data to assess if it's worth it remapping
# plot(tempData[[tempName]][,1], typ='l',xlim=c(1700,length(tempData[[tempName]][,1])),ylim=c(0,30000))
dlo=1800;
dhi=5900;
plot(tempData[[tempName]][,chDNA], typ='l',xlim=c(dlo,dhi),ylim=c(min(tempData[[tempName]][dlo:dhi,chDNA]),max(tempData[[tempName]][dlo:dhi,chDNA])), main = tempName)
# figure out threshold by checking liz500 channel
plot(tempData[[tempName]][,chLadder], xlim=c(100,2200), typ='l')
ilim01 = 1350;
plot(tempData[[tempName]][,chLadder], xlim =c(ilim01,6000), typ='l')
ilim02 = 6500;
plot(tempData[[tempName]][ilim01:ilim02,chLadder], typ='l') # 16 peaks for liz500 (15 peaks if removed '35' marker)


tempData[[tempName]] = tempData[[tempName]][ilim01:ilim02,c(chDNA,chLadder)]
# plot(tempData[[tempName]][,2], typ='l') # 16 peaks for liz500 (15 peaks if removed '35' marker)

guessThreshold = quantile(tempData[[tempName]][,2],.992);
# guessThreshold = 50;
# match ladder (works better if higher values when noise is high)
ladderData = ladder.info.attach(stored=tempData, ladder=liz500, method='iter2', draw=TRUE, ladd.init.thresh=guessThreshold)
# replot ladder if needed to play with threshold
# plot(tempData[[tempName]][,2], typ='l') # 16 peaks for liz500
# run manual correction if needed
# ladderData = ladder.corrector(stored=tempData, to.correct = tempName, ladder=liz500)

# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ladder is stored as an environment (hidden) variable -> list.data.covarrubias[[tempName]]
ladder = data.frame("p" = list.data.covarrubias[[tempName]]$pos, "w"=list.data.covarrubias[[tempName]]$wei)
# if all fails, do completely manual mapping
# plot(tempData[[tempName]][,2], typ='l', xlim = c(500,1000))
# ladder = data.frame("p" = c(118,270,410,628,685,738,958,1215,1508,1718,1776,2058,2318,2532,2578), "w"=liz500)
# ladder
# fit ladder with a 5th degree polynomial (this is what fragman uses)
polyModel = lm(w ~ poly(p,5), data = ladder)
fitWeights<- predict(polyModel,ladder,interval='confidence',level=0.99);
plot(list.data.covarrubias[[tempName]]$pos,list.data.covarrubias[[tempName]]$wei) + 
   points(list.data.covarrubias[[tempName]]$pos,fitWeights[,1], typ='l') + 
   points(list.data.covarrubias[[tempName]]$pos,fitWeights[,1], pch=20)
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# check remapped data
full_ladder = data.frame("p"=1:length(tempData[[tempName]][,1]))
fitWeights <- predict(polyModel,full_ladder)
# plot the data
plot(fitWeights, tempData[[tempName]][,1], typ='l', xlim=c(0, 600), main=tempName)
# zoom into ROI
# p_lo = 350; p_hi =  600; #syt5a | tbx2a
# p_lo = 250; p_hi = 420; #sema7a | tbx2b
p_lo = 300; p_hi = 500; # foxq2
# p_lo = 200; p_hi = 275; #eml1
# p_lo = 100; p_hi = 300; #ntng2b
# p_lo = 100; p_hi =  600; # whole range
# p_lo = 350; p_hi =  420; # temp
# plot(fitWeights, tempData[[tempName]][,1], typ='l', xlim=c(p_lo, p_hi), ylim=c(0,1000), main=tempName)
tempPeak = max(tempData[[tempName]][fitWeights>p_lo&fitWeights<p_hi,1]); tempBP =  fitWeights[which(tempData[[tempName]]==tempPeak)];
p = ggplot() + geom_line(aes(x = fitWeights, y = tempData[[tempName]][,1]), size=.5) + # frag Data
   geom_line(aes(x = fitWeights, y = tempData[[tempName]][,2]), color=rgb(.8, 0, 0, .25)) + #ladder Data
   geom_line(aes(x = c(tempBP,tempBP), y = c(0,tempPeak)), color=rgb(0, .5, .8, .8)) + # identified peak
   annotate("text", x=tempBP, y=tempPeak*1.05, label=paste(toString(round(tempBP,digits=2)),'bp',sep=' '), size=5) +
   annotate("text", x=tempBP+15, y=tempPeak, label=paste(toString(round(tempPeak,digits=0)),'au',sep=' '), size=5, hjust = 0) +
   ylab("Fluo (a.u.)") +
   xlim(p_lo, p_hi) +
   ylim(min(tempData[[tempName]][fitWeights>p_lo&fitWeights<p_hi,1]),1.1*tempPeak) +
   ggtitle(tempName) +
   theme_classic(base_size = 16, base_rect_size = 5) +
   theme(axis.line = element_line(size = 2), axis.text = element_text(size=14))

# p2 = ggplot() + geom_line(aes(x = fitWeights, y = tempData[[tempName]][,2]), color=rgb(1, 0, 0, .5)) +
#    geom_line(aes(x = c(tempBP,tempBP), y = c(0,max(tempData[[tempName]][,2]))), color=rgb(0, .5, .8, .8)) +
#    ylab("Fluo (a.u.)") +
#    xlim(p_lo, p_hi) + 
#    theme_classic()
# 
# ggarrange(p, p2, heights = c(2, 0.7),ncol = 1, nrow = 2, align = "v")

ggplotly(p)

# or plot whole ladder by itself
# plot(fitWeights, tempData[[tempName]][,2], typ='l', xlim=c(0, 100))
# export data as csv into same folder and with same name
write.csv(data.frame("size"=fitWeights,"fluo"=tempData[[tempName]][,1],"ladder"=tempData[[tempName]][,2]),paste(paste(dir.fsa,gsub('.{0,4}$', '', tempName),sep = "/"),".csv", sep=""),row.names = FALSE) + 
   message(paste("Saved analysis for:",tempName)) +
dev.off()
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------- -----------------------------------------------
# expected peaks: gnat2 = 295bp;     syt5a = 477bp;     efna1b = 495bp; 
#                 eml1 = 227bp (mut = -11bp); 
#                 sema7a = 358 bp;   tbx2a = 487 bp;   tbx2b = 332 bp
#                 foxq2 = 431 bp; slc45a2 = AA115,AB161,AC124,
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# p = ggplot() + geom_line(aes(x = fitWeights, y = tempData[[tempName]][,1]), size=.5) + # frag Data
#    geom_line(aes(x = fitWeights, y = tempData[[tempName]][,2]), color=rgb(.8, 0, 0, .8)) + #ladder Data
#    geom_line(aes(x = c(tempBP,tempBP), y = c(0,tempPeak)), color=rgb(0, .5, .8, .8)) + # identified peak
#    annotate("text", x=tempBP, y=tempPeak*1.05, label=paste(toString(round(tempBP,digits=0)),'bp',sep=' '), size=5) +
#    annotate("text", x=tempBP+15, y=tempPeak, label=paste(toString(round(tempPeak,digits=0)),'au',sep=' '), size=5, hjust = 0) +
#    ylab("Fluo (a.u.)") +
#    xlim(p_lo, p_hi) +
#    ylim(min(tempData[[tempName]][fitWeights>p_lo&fitWeights<p_hi,1]),1.1*tempPeak) +
#    ggtitle(tempName) +
#    theme_classic(base_size = 16, base_rect_size = 5) +
#    theme(axis.line = element_line(size = 2), axis.text = element_text(size=14))


p
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









