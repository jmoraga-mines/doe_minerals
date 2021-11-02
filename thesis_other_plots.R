### Other charts

# Draw polygon for MTMF valid values
fun1 <- curve(300*x+15, col="red", main="MTMF Valid pixels area", ylab="y", add=T)
fun1_x <- fun1$x[fun1$x>=0.15]
fun1_y <- fun1$y[fun1$x>=0.15]
# polygon(c(0.15,fun1_x,1), c(0,fun1_y, 0), col="gray")
polygon(c(0.07,0.07, 1, 1), 
        c(0,300*0.07+15, 315, 0), col="gray")
abline(v=0.07, col="blue")
fun1 <- curve(300*x+15, col="red", add=T, main="MTMF Valid pixels area", ylab="y")

mtmf <- stack("d:/mm_final/Brady/HyMapBrady_mtmf")
chalcedony_mtmf_x <- values(mtmf[[1]])
chalcedony_mtmf_y <- values(mtmf[[8]])

### Target Minerals Resample Library
if (!require(pacman)){
  install.packages("pacman")
  library(pacman)
  }
pacman::p_load(RStoolbox)
TargetMinerals <- readSLI("E:/mm_final/TargetMineralsResampled.sli")
f1 <- write.csv(TargetMinerals, "e:/mm_final/TargetMineralsResampled.csv")
NonTargetMinerals <- readSLI("E:/mm_final/Non_Target_Spectra_2.sli")
f1 <- write.csv(NonTargetMinerals, "e:/mm_final/Non_Target_Spectra_2.csv")
max(ceiling(NonTargetMinerals/1000)*1000)
NonTargetMinerals
n <- NonTargetMinerals[c(1:61,68:126),]
plot(n$wavelength, n$Asphalt_01, type="l", ylim= range(n$Asphalt_01, n$Asphalt_02+2000), col="blue")
par(new=TRUE) # will add to current plot
plot(n$wavelength, n$Asphalt_02+2000, type="l", ylim= range(n$Asphalt_01, n$Asphalt_02+2000), col="red", add=T, axes=F, xlab="", ylab="")
par(new=TRUE) # will add to current plot
plot(n$wavelength, n$Asphalt_01, col="blue", ylim= range(n$Asphalt_01, n$Asphalt_02+2000), add=T)

op <- par(mar=c(5,4,4,13))
plot(n$wavelength, n$Asphalt_01, type="b", lwd=2, pch=16, 
     ylim= range(n$Asphalt_01, n$Asphalt_4+3000), col="blue", 
     xlab="Wavelength (micrometers)", ylab="Reflectance")
lines(n$wavelength, n$Asphalt_02+1000, type="b", lwd=2, pch=16, 
      ylim= range(n$Asphalt_01, n$Asphalt_4+3000), col="red")
lines(n$wavelength, n$Asphalt_03+2000, type="b", lwd=2, pch=16, 
      ylim= range(n$Asphalt_01, n$Asphalt_4+3000), col="green")
lines(n$wavelength, n$Asphalt_4+3000, type="b", lwd=2, pch=16, 
      ylim= range(n$Asphalt_01, n$Asphalt_4+3000), col="orange")
title("Non target pixels spectral library\n(displaced for legibility)")
par(op)
legend("topright", legend = c("Asphalt 1", "Asphalt 2", "Asphalt 3", "Asphalt 4"), 
       lwd=2, pch = 16, col=c("blue", "red", "green", "orange"), inset = c(-0.34,0.1))

