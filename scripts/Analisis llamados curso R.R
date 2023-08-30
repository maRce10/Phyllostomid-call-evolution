
library(warbleR)
library(nlme)
library(ape)

warbleR_options(wav.path = "~/Dropbox/Projects/Phyllostomid call evolution/grabaciones/", bp = c(1, 193), parallel = 3, wl = 124)


#Verificar que los archivos puedan ser leidos
checkwavs()

#importar selecciones 
selinfin <- read.csv("selinfin.csv", header = T, sep = ",")
cs checksels(X=selinfin, parallel = 3, path = "C:/Users/Nazareth/Documents/Datos acusticos/Tesis")
read_wave(X=selinfin, index = 1, path="C:/Users/Nazareth/Documents/Datos acusticos/Tesis")
setwd("C:/Users/Nazareth/Desktop/R_avanzado/pr0yect0")

#Filtrar se?ales de mala calidad
fsel<-sig2noise(X=selinfin[seq(1,nrow(selinfin),2),],mar=0.01, parallel=3, path="C:/Users/Nazareth/Documents/Datos acusticos/Tesis")

#Medir parametros acusticos
params<-specan(fsel,bp=c(1,193), parallel=3, wl=124, path="C:/Users/Nazareth/Documents/Datos acusticos/Tesis")
write.csv(params, "params.csv")
head(params, 3)

#agregar columna con nombre de spp
pys<-cbind(fsel$spp, params)
str(pys)
names(pys)[1]<-"spp"
levels(pys$spp)

#promediar parametros por spp con un loop
#objeto vacio 
mdat<-NULL

#asignar c/spp como objeto
cspp <- pys$spp
CS <- split(pys, cspp)
l <- length(CS) 

#loop
for (i in 1:l) {
  dts<- CS[[i]]
  head(dts, 3)
  sp<-as.character(unique(dts$spp))
  
  duration<-mean(dts$duration)
  meanfreq<-mean(dts$meanfreq)
  sd<-mean(dts$sd)
  freq.median<-mean(dts$freq.median)
  freq.Q25<-mean(dts$freq.Q25)
  freq.Q75<-mean(dts$freq.Q75)
  freq.IQR<-mean(dts$freq.IQR)
  time.median<-mean(dts$time.median)
  time.Q25<-mean(dts$time.Q25)
  time.Q75<-mean(dts$time.Q75)
  time.IQR<-mean(dts$time.IQR)
  skew<-mean(dts$skew)
  kurt<-mean(dts$kurt)
  sp.ent<-mean(dts$sp.ent)
  time.ent<-mean(dts$time.ent)
  entropy<-mean(dts$entropy)
  sfm<-mean(dts$sfm)
  meandom<-mean(dts$meandom)
  mindom<-mean(dts$mindom)
  maxdom<-mean(dts$maxdom)
  dfrange<-mean(dts$dfrange)
  modindx<-mean(dts$modindx)
  startdom<-mean(dts$startdom)
  enddom<-mean(dts$enddom)
  dfslope<-mean(dts$dfslope)
  meanpeakf<-mean(dts$meanpeakf)
  
  if (is.null(mdat)) mdat <- data.frame(sp, duration,meanfreq, sd, freq.median, freq.Q25, freq.Q75, freq.IQR, time.median, time.Q25, time.Q75, time.IQR, skew, kurt, sp.ent, time.ent, entropy, sfm, meandom, mindom, maxdom, dfrange, modindx, startdom, enddom, dfslope, meanpeakf)  else
    mdat <- rbind(mdat, data.frame(sp, duration,meanfreq, sd, freq.median, freq.Q25, freq.Q75, freq.IQR, time.median, time.Q25, time.Q75, time.IQR, skew, kurt, sp.ent, time.ent, entropy, sfm, meandom, mindom, maxdom, dfrange, modindx, startdom, enddom, dfslope, meanpeakf))
}

head(mdat, 2)

#Agregar columna de gremio a los datos
grem<-read.csv("gremios_x_spp.csv", header = T, sep = ",")
head(grem, 1)
levels(grem$spp)
levels(mdat$sp)
pyg <- cbind(grem$guild, mdat)
str(pyg)
names(pyg)[1]<- "guild"
levels(pyg$guild)
pyg<-subset(pyg, select=c(sp, guild, duration:meanpeakf))
head(pyg, 2)
write.csv(pyg, "C:/Users/Nazareth/Desktop/R_avanzado/pr0yect0/pyg.csv")

############################################################################################################################################################
#datos filogeneticos
setwd("C:/Users/Nazareth/Desktop/R_avanzado/pr0yect0")
arb<-read.nexus("S17613_arblmurcis.nex")
plot(arb)

#eliminar del arbol las spp que no interesan para el analisis 
spp <- paste(sapply(arb$tip.label, function(x) strsplit(x, "_")[[1]][1]),sapply(arb$tip.label, function(x) strsplit(x, "_")[[1]][2]), sep = "_")

arb2 <- drop.tip(arb, tip = setdiff(arb$tip.label, pyg$sp))

arb2$tip.label
plot(arb2)

supf <- pyg[pyg$sp %in% arb2$tip.label,]

Ntip(arb2)
nrow(supf)
head(supf, 1)

#######################################################################################################################################

#analisis PGLS
#duration
mod1 <- gls(duration ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod1)
boxplot(duration ~ guild, data = supf)

#meanfreq
mod2 <- gls(meanfreq ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod2)
boxplot(meanfreq ~ guild, data = supf)

#sd
mod3 <- gls(sd ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod3)
boxplot(sd ~ guild, data = supf)

#freq.median
mod4 <- gls(freq.median ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod4)
boxplot(freq.median ~ guild, data = supf)

#freq.Q25
mod5 <- gls(freq.Q25 ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod5)
boxplot(freq.Q25 ~ guild, data = supf)

#freq.Q75
mod6 <- gls(freq.Q75 ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod6)
boxplot(freq.Q75 ~ guild, data = supf)

#freq.IQR
mod7 <- gls(freq.IQR ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod7)
boxplot(freq.IQR ~ guild, data = supf)

#time.median
mod8 <- gls(time.median ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod8)
boxplot(time.median ~ guild, data = supf)

#time.Q25
mod9 <- gls(time.Q25 ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod9)
boxplot(time.Q25 ~ guild, data = supf)

#time.Q75
mod10 <- gls(time.Q75 ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod10)
boxplot(time.Q75 ~ guild, data = supf)

#time.IQR
mod11 <- gls(time.IQR ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod11)
boxplot(time.IQR ~ guild, data = supf)

#skew
mod12 <- gls(skew ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod12) 
boxplot(skew ~ guild, data = supf)

#kurt
mod13 <- gls(kurt ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod13)
boxplot(kurt ~ guild, data = supf)

#sp.ent
mod14 <- gls(sp.ent ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod14)
boxplot(sp.ent ~ guild, data = supf)

#time.ent
mod15 <- gls(time.ent ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod15)
boxplot(time.ent ~ guild, data = supf)

#entropy
mod16 <- gls(entropy ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod16)
boxplot(entropy ~ guild, data = supf)

#sfm
mod17 <- gls(sfm ~ guild, correlation = corBrownian(phy = arb2),  data = supf, method = "ML")
anova(mod17)
boxplot(sfm ~ guild, data = supf)

#meandom
mod18 <- gls(meandom ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod18)
boxplot(meandom ~ guild, data = supf)

#mindom
mod19 <- gls(mindom ~ guild, correlation = corBrownian(phy = arb2), 
data = supf, method = "ML")
anova(mod19)
boxplot(mindom ~ guild, data = supf)

#maxdom
mod20 <- gls(maxdom ~ guild, correlation = corBrownian(phy = arb2), 
data = supf, method = "ML")
anova(mod20)
boxplot(maxdom ~ guild, data = supf)

#dfrange
mod21 <- gls(dfrange ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod21)
boxplot(dfrange ~ guild, data = supf)

# modindx
mod22 <- gls(modindx ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod22)
boxplot(modindx ~ guild, data = supf)

#startdom
mod23 <- gls(startdom ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod23)
boxplot(startdom ~ guild, data = supf)

#enddom
mod24 <- gls(enddom ~ guild, correlation = corBrownian(phy = arb2), 
data = supf, method = "ML")
anova(mod24)
boxplot(enddom ~ guild, data = supf)

#dfslope
mod25 <- gls(dfslope ~ guild, correlation = corBrownian(phy = arb2),  data = supf, method = "ML")
anova(mod25)
boxplot(dfslope ~ guild, data = supf)

# meanpeakf
mod26 <- gls(meanpeakf ~ guild, correlation = corBrownian(phy = arb2), data = supf, method = "ML")
anova(mod26)
boxplot(meanpeakf ~ guild, data = supf)



#Graficos de pruebas significativas
par(mfcol = c(3, 2))
par(mar = c(2, 3, 2, 1), oma = c(2, 3, 0, 0))

#meanfreq
boxplot(meanfreq ~ guild, data = supf, col="grey90", type="n", xaxt = "n", yaxt = "n", frame.plot = F, main="frecuencia promedio ")
axis(1, at = seq(1, 6, by=1), labels = c("APEC", "DAPC", "FAEE", "FBEA", "FPEE", "HEMA"), cex.axis = 1)
axis(2, las=1)
box(bty = "l", lwd = 1.1)

#sd
boxplot(sd ~ guild, data = supf, col="grey90", type="n", xaxt = "n", yaxt = "n", frame.plot = F, main="sd de la frecuencia promedio")
axis(1, at = seq(1, 6, by=1), labels = c("APEC", "DAPC", "FAEE", "FBEA", "FPEE", "HEMA"), cex.axis = 1)
axis(2, las=1)
box(bty = "l", lwd = 1.1)

#freq.Q75
boxplot(freq.Q75 ~ guild, data = supf, col="grey90", type="n", xaxt = "n", yaxt = "n", frame.plot = F, main="frecuencia Q75")
axis(1, at = seq(1, 6, by=1), labels = c("APEC", "DAPC", "FAEE", "FBEA", "FPEE", "HEMA"), cex.axis = 1)
axis(2, las=1)
box(bty = "l", lwd = 1.1)

#freq.IQR
boxplot(freq.IQR ~ guild, data = supf, col="grey90", type="n", xaxt = "n", yaxt = "n", frame.plot = F, main="frecuencia IQR")
axis(1, at = seq(1, 6, by=1), labels = c("APEC", "DAPC", "FAEE", "FBEA", "FPEE", "HEMA"), cex.axis = 1)
axis(2, las=1)
box(bty = "l", lwd = 1.1)

#modindx
boxplot(modindx ~ guild, data = supf, col="grey90", type="n", xaxt = "n", yaxt = "n", frame.plot = F, main="?ndice de modulaci?n ")
axis(1, at = seq(1, 6, by=1), labels = c("APEC", "DAPC", "FAEE", "FBEA", "FPEE", "HEMA"), cex.axis = 1)
axis(2, las=1)
box(bty = "l", lwd = 1.1)

#sfm
boxplot(sfm ~ guild, data = supf, col="grey90", type="n", xaxt = "n", yaxt = "n", frame.plot = F, main="SFM")
axis(1, at = seq(1, 6, by=1), labels = c("APEC", "DAPC", "FAEE", "FBEA", "FPEE", "HEMA"), cex.axis = 1)
axis(2, las=1)
box(bty = "l", lwd = 1.1)

mtext("Gremio", side = 1, line = 0.6, cex = 1.35, outer = T)
mtext("Frecuencia", side = 2, line = 0.6, cex = 1.35, outer = T)
