
# Load package for biodiversity analysis
library(vegan)

# In case user is running windows and can't use quartz 
if (.Platform$OS.type=="windows") {quartz<-function() windows()}

##############################################################################################################################
# Preparations
##############################################################################################################################

# Set working directory (you need to adjust this depending on where you put the project)
setwd("~/Desktop/GitHub/Mosquito_ecology/Analysis/")

# Read data into R
data = read.csv("FinalData.csv", header=T, stringsAsFactors=T)

# Subset species data 
species = data[,16:38]
row.names(species) = data$AbrSiteName

# Subset landscape data, adding categorical habitat type as first column
Habitat = rep(c("Field", "NearField", "Edge", "NearForest", "Forest"), 9)
envir = data.frame(Habitat, data[,c("BarrenLand", "Building", "CultivatedCrops", "DeciduousForest", "EvergreenForest", "Grassland", "MixedForest", "Pavement", "ShrubScrub")])
row.names(envir) = data$AbrSiteName

##############################################################################################################################
# Species richness, Shannon-Wiener, and rarefied diversity across habitat types (Fig 2A-C)
##############################################################################################################################

# Define vector of habitats as factors (needed for boxplots to display correctly)
Habitat_f = factor(envir$Habitat, c("Field", "NearField", "Edge", "NearForest", "Forest"))

# Compute species richness across habitat types
evenness = data.frame(Habitat = Habitat_f, Transect = data$Transect, Abundance = rowSums(species), Richness = rowSums(species>0))

# Compute Shannon-Wiener diversity across habitat types
shannondiv = c()
for (i in 1:nrow(diversity.data)) {shannondiv[i] = diversity(diversity.data[i,1:23])}
evenness = data.frame(evenness, ShannonDiversity=shannondiv)

# Estimate rarefied richness across habitat types
rarefied = c()
for (i in 1:nrow(diversity.data)) {rarefied[i] = rarefy(diversity.data[i,1:23], sample=15)}
evenness = data.frame(evenness, Rarefied=rarefied)

# ANOVAs to compare diversity across habitat types
anova(lm(Richness~Habitat, data=evenness))	# significant
anova(lm(ShannonDiversity~Habitat, data=evenness))	# nearly significant
anova(lm(Rarefied~Habitat, data=evenness))	# significant

# Make boxplots of diversity measures across habitat types (Fig 2A-C in manuscript) 
quartz()
colors = c("darkgoldenrod4", "darkgoldenrod1", "lightgoldenrod2", "green", "darkgreen")
box.spacing = 1:5/7
x.lims = c(0,box.spacing[5] + box.spacing[1])
layout(matrix(1:3,1,3))
par(mar=c(3,2,3,1))
boxplot(Richness~Habitat, data=evenness, at=box.spacing, boxwex=0.1, pch=20, cex=0.5, medlwd = 1, col=colors, xlim=x.lims, cex.axis=0.5, main="Richness")
boxplot(ShannonDiversity~Habitat, data=evenness, at=box.spacing, boxwex=0.1, pch=20, cex=0.5, medlwd = 1, col=colors, xlim=x.lims, cex.axis=0.5, main="Shannon-Weiner diversity")
boxplot(Rarefied~Habitat, data=evenness, at=box.spacing, boxwex=0.1, pch=20, cex=0.5, medlwd = 1, col=colors, xlim=x.lims, cex.axis=0.5, main="Rarefied richness")

########################################################################################################################################
# Partial canonincal correspondence analysis (Fig 3A-B)
########################################################################################################################################

# with scaling=1, sites are centroids of species, they are plotted close to the species that tend to occur there
# with scaling=2, species are centroids of sites, they are plotted close to the sites where they occur

# Compute average square root transformation to dampen effects of large data points
species.sqrt = sqrt(species)/data$TrapNights

# pcca with environmental variables as predictors, conditioned on non forest/field related environmental variables
factors = c("DeciduousForest", "EvergreenForest", "Grassland", "MixedForest", "ShrubScrub")
covariables = c("BarrenLand", "Building", "CultivatedCrops", "Pavement")
pcca = cca(X=species.sqrt, Y=envir[factors], Z=envir[covariables])

# Significance test for variation explained by pcca
anova.cca(pcca)
# Percent variation explained by environmental variables
round(pcca$CCA$tot.chi/pcca$tot.chi, 4)*100
# Percent of variation explained by first axis
round(sum(pcca$CCA$eig[1])/pcca$CCA$tot.chi, 4)*100
# Percent of variation explained by second axis
round(sum(pcca$CCA$eig[2])/pcca$CCA$tot.chi, 4)*100
# Percent of variation explained by first 2 axes combined
round(sum(pcca$CCA$eig[1:2])/pcca$CCA$tot.chi, 4)*100

# pCCA plot (sites and species separately)
quartz()
layout(matrix(1:2,1,2))
site.cols = Habitat
site.cols[site.cols=="Field"] = "darkgoldenrod4"
site.cols[site.cols=="NearField"] = "darkgoldenrod1"
site.cols[site.cols=="Edge"] = "lightgoldenrod2"
site.cols[site.cols=="NearForest"] = "green"
site.cols[site.cols=="Forest"] = "darkgreen"
site.pch = c(rep(6,15), rep(5,15), rep(0,15))
plot(pcca, display=c("sites", "bp"), type="n", xlim=c(-2.5, 2.5), ylim=c(-1, 1), main="Sites")
abline(v=0, col="white", lwd=2)
abline(h=0, col="white", lwd=2)
abline(v=0, col="lightgray", lwd=0.5, lty=3)
abline(h=0, col="lightgray", lwd=0.5, lty=3)
text(pcca, display="bp", col="darkgray", cex=0.5, arrow.mul=2.2)
points(pcca, scaling=1, display="sites", col=site.cols, pch=site.pch, cex=0.75, font=2)
points(x=rep(-2.5, 8), y=(c(2.1,2,1.9,1.8,1.7,1.6,1.5,1.4)), pch=c(rep(19,5),6,5,0), cex=0.4, col=c(site.cols[1:5], rep("black", 3)))
legendlabs = c("Field (0 m)", "Near field (90 m)", "Edge (100 m)", "Near forest (110 m)", "Forest (200 m)", "Durham", "Prairie Ridge", "Lake Wheeler")
text(x=rep(-2.5, 8), y=(c(2.1,2,1.9,1.8,1.7,1.6,1.5,1.4)-0.02), labels=legendlabs, cex=0.4, pos=4)
plot(pcca, display=c("species", "bp"), col="red", type="n", ylab="", xlim=c(-2, 2), ylim=c(-1, 1), main="Species")
abline(v=0, col="white", lwd=2)
abline(h=0, col="white", lwd=2)
abline(v=0, col="lightgray", lwd=0.5, lty=3)
abline(h=0, col="lightgray", lwd=0.5, lty=3)
text(pcca, display="bp", col="darkgray", cex=0.5, arrow.mul=1.7)
text(pcca, scaling=2, display="species", col="black", cex=0.4, font=4)

##############################################################################################################################
# Fit GLMMS
##############################################################################################################################

# Load package to do GLMMs with negative binomial errors
library(glmmADMB)

# Store model output, AIC differences, and predicted values of the best model
linear = list()
quadratic = list()
predicted = list()
AICdiff = c()

# Define predictor variables 
field.dist = rep(c(0, 90, 100, 110, 200), 9)
transect = data$Transect

# Define values for generating predictions
newdata = data.frame(field.dist=0:20*10)

# RUN MODELS, COMPARE MODELS, SAVE RESULTS
  # Diagnostics were done for individual models to determine whether transformation or
  # fitting a zero-inflated model improved model fit, and to determine whether a species
  # should be included (many rare species could not be adequately modeled)

# Cx.sal
i = 1
linear[[i]] = glmmadmb(log1p(species[,"Cx.sal"]) ~ field.dist + (1|transect), family="nbinom")
quadratic[[i]] = glmmadmb(log1p(species[,"Cx.sal"]) ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(linear[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "Cx.sal"
# Ae.albo
i = 2
linear[[i]] = glmmadmb(log1p(species[,"Ae.albo"]) ~ field.dist + (1|transect), family="nbinom")
quadratic[[i]] = glmmadmb(log1p(species[,"Ae.albo"]) ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(quadratic[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "Ae.albo"
# Ae.cin
i = 3
linear[[i]] = glmmadmb(species[,"Ae.cin"] ~ field.dist + (1|transect), family="nbinom")
quadratic[[i]] = glmmadmb(species[,"Ae.cin"] ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(linear[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "Ae.cin"
# Ae.vex
i = 4
linear[[i]] = glmmadmb(log1p(species[,"Ae.vex"]) ~ field.dist + (1|transect), family="nbinom")
quadratic[[i]] = glmmadmb(log1p(species[,"Ae.vex"] )~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(linear[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "Ae.vex"
# Ps.fer
i = 5
linear[[i]] = glmmadmb(log1p(species[,"Ps.fer"]) ~ field.dist + (1|transect), family="nbinom", zeroInflation = TRUE)
quadratic[[i]] = glmmadmb(log1p(species[,"Ps.fer"]) ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom", zeroInflation = TRUE)
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(linear[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "Ps.fer"
# Cx.err
i = 6
linear[[i]] = glmmadmb(species[,"Cx.err"] ~ field.dist + (1|transect), family="nbinom")
quadratic[[i]] = glmmadmb(species[,"Cx.err"] ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(linear[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "Cx.err"
# Ps.col
i = 7
linear[[i]] = glmmadmb(log1p(species[,"Ps.col"]) ~ field.dist + (1|transect), family="nbinom")
quadratic[[i]] = glmmadmb(log1p(species[,"Ps.col"]) ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(linear[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "Ps.col"
# Oc.tris
i = 8
linear[[i]] = glmmadmb(log1p(species[,"Oc.tris"]) ~ field.dist + (1|transect), family="nbinom")
quadratic[[i]] = glmmadmb(log1p(species[,"Oc.tris"]) ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(linear[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "Oc.tris"
# Cx.pip.q
i = 9
linear[[i]] = glmmadmb(log1p(species[,"Cx.pip.q"]) ~ field.dist + (1|transect), family="nbinom")
quadratic[[i]] = glmmadmb(log1p(species[,"Cx.pip.q"]) ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(linear[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "Cx.pip.q"
# An.quad
i = 10
linear[[i]] = glmmadmb(log1p(species[,"An.quad"]) ~ field.dist + (1|transect), family="nbinom")
quadratic[[i]] = glmmadmb(log1p(species[,"An.quad"]) ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(linear[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "An.quad"
# An.punc
i = 11
linear[[i]] = glmmadmb(log1p(species[,"An.punc"]) ~ field.dist + (1|transect), family="nbinom")
quadratic[[i]] = glmmadmb(log1p(species[,"An.punc"]) ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
predicted[[i]] = predict(linear[[i]], newdata)
names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- "An.punc"

# Make table of results
BestModel = c()
LeadingCoef = c()
S.E. = c()
Interpretation = c()
for (i in 1:length(AICdiff)) {
  if (AICdiff[i] > -2) {
    BestModel[i] = "Linear"
    coefs = coefficients(summary(linear[[i]]))
    LeadingCoef[i] = coefs[2,"Estimate"]
    S.E.[i] = coefs[2,"Std. Error"]
    if (LeadingCoef[i] > 0) {dir = "positive"} else {dir = "negative"}
    if (coefs[2,"Pr(>|z|)"] < 0.05) {
      if (LeadingCoef[i] > 0) {
        Interpretation[i] = "Significant forest preference."
      } else {Interpretation[i] = "Significant field preference."}
    } else {Interpretation[i] = "No significant preference."}
  } else {
    BestModel[i] = "Quadratic"
    coefs = coefficients(summary(quadratic[[i]]))
    LeadingCoef[i] = coefs[3,"Estimate"]
    S.E.[i] = coefs[3,"Std. Error"]
    if (LeadingCoef[i] < 0) {
      Interpretation[i] = "Significant edge preference."
    } else {Interpretation[i] = "Significant edge avoidance"}
  }
}
glmm.results = data.frame(AICdiff, LeadingCoef, S.E., Interpretation)
glmm.results

##############################################################################################################################
# Plot best fitting GLMMs for sufficiently abundant species (Fig 4A-K)
##############################################################################################################################

quartz()
layout(matrix(1:12, 4, 3))
par(mar=c(2.9,3.2,1.5,1.5))
long.names = c("Culex salinarius", "Aedes albopictus", "Aedes cinereus", "Aedes vexans", "Psorophora ferox", "Culex erraticus", "Psorophora columbiae", "Ochlerotatus triseriatus", "Culex pipiens/quinquefasciatus", "Anopheles quadrimaculatus", "Anopheles punctipennis")
for (i in names(linear)) {
  if (which(names(linear)==i) > 4) {ylab = ""} else {ylab = "Adbundance"}
  plot(species[,i] ~ field.dist, ylab=ylab, xlab="", pch=20, cex=0.5, cex.axis=0.5, cex.lab=0.7, tck=-0.03, mgp=c(2, 0.2, 0), xaxt="n")
  axis(side=1, at=c(0, 90, 100, 110, 200), labels = c("-100", "-10", "0", "10", "100"), cex.axis=0.5, tck=-0.03, mgp=c(2, 0.2, 0))
  axis(side=1, at=c(0, 100, 200), labels=c("(Field)", "(Edge)", "(Forest)"), tick=F, cex.axis=0.5, mgp=c(2,0.7,0))
  if (i %in% c("Ae.cin", "Cx.err")) {lines(exp(predicted[[i]]) ~ newdata$field.dist)} else {lines(expm1(exp(predicted[[i]])) ~ newdata$field.dist)}
  mtext(long.names[i], adj=0, cex=0.5, font=3)
  p = coefficients(summary(linear[[i]]))[2,"Pr(>|z|)"]
  if (p > 0.05) {p = "p > 0.05"}
  if (p < 0.05) {p = "p < 0.05*"}
  if (p < 0.01) {p = "p < 0.01**"}
  if (p < 0.001) {p = "p < 0.001***"}
  mtext(paste("n = 45;", p), cex=0.5, adj=1)
}


##############################################################################################################################
# Mood's median test for rare species
##############################################################################################################################

# Load package for processing data
library(plyr)

# Prepare data frame
species.rare = species[,setdiff(names(species), names(linear))]

# Function to perform Mood's median test (takes dataframe with 2 cols: col1=continuous data, col2=categories)
moodmedian.test = function (data) {
  med = median(data[,1])
  tab = ddply(data, 2, function (x) {return(data.frame(greater=sum(x[,1]>med), lessequal=sum(x[,1]<=med)))})
  return(fisher.test(matrix(c(tab$greater, tab$lessequal), 2, length(unique(data[,2])), byrow=T))$p.value)
}

# Run Mood's median test for all species
n <- p <- s <- c()
for (i in 1:ncol(species.rare)) {
  n[i] = names(species.rare)[i]
  pval = moodmedian.test(data.frame(species.rare[,i], Habitat_f))
  p[i] = pval
  if (pval > 0.05) {s = c(s, "p > 0.05")} 
  if (pval < 0.05 & pval > 0.01) {s = c(s, "p < 0.05*")} 
  if (pval < 0.01 & pval > 0.001) {s = c(s, "p < 0.01**")} 
  if (pval < 0.001) {s = c(s, "p < 0.001***")} 
}
mood.results = data.frame(species = n, p.value=p, significance = s)
mood.results

##############################################################################################################################
# Box-plots with Mood's median test results (Fig S1)
##############################################################################################################################

quartz()
colors = c("darkgoldenrod4", "darkgoldenrod1", "lightgoldenrod2", "green", "darkgreen")
box.labs = c("Field", "Field edge", "Edge", "Forest edge", "Forest")
layout(matrix(1:12, 4, 3))
par(mar=c(2,2,3,1), mgp=c(2,0.2,0), tck=-0.02)
box.spacing = 1:5/7
x.lims = c(0,box.spacing[5] + box.spacing[1])
long.names2 = c("Ochlerotatus canadensis", "Coquillettidea perturbans", "Anopheles crucians", "Ochlerotatus hendersoni", "Ochlerotatus atlanticus", "Ochlerotatus dupreei", "Psorophora cyanescens", "Psorophora howardi","Aedes japonicus", "Aedes japonicus", "Ochlerotatus infirmatus", "Psorophora ciliata", "Urotanaenia sapphirina")

for (i in names(species.rare)) {
  if (which(names(species.rare)==i) > 4) {ylab = ""} else {ylab = "Adbundance"}
  boxplot(species.rare[,i]~Habitat_f, ylab=ylab, at=box.spacing, boxwex=0.1, pch=20, cex=0.5, medlwd = 1, cex.axis=0.5, xlim=x.lims, col=colors, names=box.labs)
  mtext(long.names2[which(names(species.rare)==i)], line=0.4, cex=0.6, adj=0, font=4)
  mtext(paste("n = ", sum(species.rare[,i]), "; ", mood.results$significance[which(mood.results$species == i)], sep=""), line=0.4, cex=0.5, adj=1)
}

##############################################################################################################################
# END
##############################################################################################################################







