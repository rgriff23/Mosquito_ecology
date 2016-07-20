
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
for (i in 1:nrow(species)) {shannondiv[i] = diversity(species[i,])}
evenness = data.frame(evenness, ShannonDiversity=shannondiv)

# Estimate rarefied richness across habitat types
rarefied = c()
for (i in 1:nrow(species)) {rarefied[i] = rarefy(species[i,], sample=15)}
evenness = data.frame(evenness, Rarefied=rarefied)

# ANOVAs to compare diversity across habitat types
anova(lm(Richness~Habitat, data=evenness))	# significant
anova(lm(ShannonDiversity~Habitat, data=evenness))	# nearly significant
anova(lm(Rarefied~Habitat, data=evenness))	# significant

# Make boxplots of diversity measures across habitat types (Fig 2A-C in manuscript) 
quartz()
par(mar=c(3,4,3,2), mfrow=c(2,2))
colors = c("floralwhite", "lightgoldenrod1","chartreuse2", "chartreuse4", "darkgreen")
box.spacing = 1:5/7
x.lims = c(0,box.spacing[5] + box.spacing[1])
boxplot(Richness~Habitat, data=evenness, at=box.spacing, boxwex=0.1, pch=20, cex=0.5, medlwd = 1, col=colors, xlim=x.lims, cex.axis=0.5, ylab="Richness")
mtext("A", adj=0, line=1)
boxplot(ShannonDiversity~Habitat, data=evenness, at=box.spacing, boxwex=0.1, pch=20, cex=0.5, medlwd = 1, col=colors, xlim=x.lims, cex.axis=0.5, ylab="Shannon-Weiner diversity")
mtext("B", adj=0, line=1)
boxplot(Rarefied~Habitat, data=evenness, at=box.spacing, boxwex=0.1, pch=20, cex=0.5, medlwd = 1, col=colors, xlim=x.lims, cex.axis=0.5, ylab="Rarefied richness")
mtext("C", adj=0, line=1)

########################################################################################################################################
# Partial canonincal correspondence analysis (Fig 3A-B)
########################################################################################################################################

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
site.cols[site.cols=="Field"] = "floralwhite"
site.cols[site.cols=="NearField"] = "lightgoldenrod1"
site.cols[site.cols=="Edge"] = "chartreuse2"
site.cols[site.cols=="NearForest"] = "chartreuse4"
site.cols[site.cols=="Forest"] = "darkgreen"
site.pch = c(rep(6,15), rep(5,15), rep(0,15))
plot(pcca, display=c("sites", "bp"), type="n", xlim=c(-2.5, 2.5), ylim=c(-1, 1))
mtext("A", adj=0, line=1)
abline(v=0, col="white", lwd=2)
abline(h=0, col="white", lwd=2)
abline(v=0, col="lightgray", lwd=0.5, lty=3)
abline(h=0, col="lightgray", lwd=0.5, lty=3)
text(pcca, display="bp", col="darkgray", cex=0.5, arrow.mul=2.2)
points(pcca, scaling=1, display="sites", col=site.cols, pch=site.pch, cex=0.75, font=2)
points(x=rep(-2.5, 8), y=(c(2.1,2,1.9,1.8,1.7,1.6,1.5,1.4)+0.3), pch=c(rep(19,5),6,5,0), cex=0.4, col=c(site.cols[1:5], rep("black", 3)))
legendlabs = c("Field (0 m)", "Near field (90 m)", "Edge (100 m)", "Near forest (110 m)", "Forest (200 m)", "Durham", "Prairie Ridge", "Lake Wheeler")
text(x=rep(-2.5, 8), y=(c(2.1,2,1.9,1.8,1.7,1.6,1.5,1.4)+0.3), labels=legendlabs, cex=0.4, pos=4)
plot(pcca, display=c("species", "bp"), col="red", type="n", ylab="", xlim=c(-2, 2), ylim=c(-1, 1))
mtext("B", adj=0, line=1)
abline(v=0, col="white", lwd=2)
abline(h=0, col="white", lwd=2)
abline(v=0, col="lightgray", lwd=0.5, lty=3)
abline(h=0, col="lightgray", lwd=0.5, lty=3)
text(pcca, display="bp", col="darkgray", cex=0.5, arrow.mul=1.7)
text(pcca, scaling=2, display="species", col="black", cex=0.4, font=1)

##############################################################################################################################
# Fit GLMMS
##############################################################################################################################

# Load package to do GLMMs with negative binomial errors
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")), type="source")
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

# Prepare data frame (species with >18 non-zero abundance values)
species.common = species[,which(colSums(species>0)>18)]

# Loop through common species
# Log1p transform for all but Ae.cin and Cx.err
# Add zero-inflation for Ps.fer
for (i in 1:ncol(species.common)) {
  if (names(species.common)[i] %in% c("Ae.cin", "Cx.err")) {
    linear[[i]] = glmmadmb(species.common[,i] ~ field.dist + (1|transect), family="nbinom")
    quadratic[[i]] = glmmadmb(species.common[,i] ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom")
  } else {
    if (names(species.common)[i] != "Ps.fer") {z = F} else {z = T}
    linear[[i]] = glmmadmb(log1p(species.common[,i]) ~ field.dist + (1|transect), family="nbinom", zeroInflation=z)
    quadratic[[i]] = glmmadmb(log1p(species.common[,i]) ~ field.dist + I(field.dist^2) + (1|transect), family="nbinom", zeroInflation=z)
    }
  AICdiff[i] = AIC(quadratic[[i]]) - AIC(linear[[i]])
  predicted[[i]] = predict(linear[[i]], newdata)
  names(linear)[i] <- names(quadratic)[i] <- names(predicted)[i] <- names(AICdiff)[i] <- names(species.common)[i]
}

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
long.names = c("Culex salinarius", "Aedes albopictus", "Aedes cinereus", "Aedes vexans", "Psorophora ferox", "Culex erraticus", "Psorophora columbiae", "Aedes triseriatus", "Culex pipiens/quinquefasciatus", "Anopheles quadrimaculatus", "Anopheles punctipennis")
for (i in names(linear)) {
  if (which(names(linear)==i) > 4) {ylab = ""} else {ylab = "log(adbundance + 1)"}
  plot(log1p(species[,i]) ~ field.dist, ylab=ylab, xlab="", pch=20, cex=0.5, cex.axis=0.5, cex.lab=0.7, tck=-0.03, mgp=c(2, 0.2, 0), xaxt="n")
  axis(side=1, at=c(0, 90, 100, 110, 200), labels = c("-100", "-10", "0", "10", "100"), cex.axis=0.5, tck=-0.03, mgp=c(2, 0.2, 0))
  axis(side=1, at=c(0, 100, 200), labels=c("(Field)", "(Edge)", "(Forest)"), tick=F, cex.axis=0.5, mgp=c(2,0.7,0))
  if (i %in% c("Ae.cin", "Cx.err")) {lines(log1p(exp(predicted[[i]])) ~ newdata$field.dist)} else {lines(exp(predicted[[i]]) ~ newdata$field.dist)}
  mtext(long.names[which(names(linear)==i)], adj=0, cex=0.5, font=4)
  p = coefficients(summary(linear[[i]]))[2,"Pr(>|z|)"]
  if (p > 0.05) {p = "p > 0.05"}
  if (p < 0.05) {p = "p < 0.05*"}
  if (p < 0.01) {p = "p < 0.01**"}
  if (p < 0.001) {p = "p < 0.001***"}
  mtext(paste("n = 45;", p), cex=0.5, adj=1)
}

##############################################################################################################################
# Paired t-test for difference in abundance between -10m and 10m (not included in published paper)
##############################################################################################################################

trap2 = species.common[c(c(1:9)*5)-3,]
trap4 = species.common[c(c(1:9)*5)-1,]
ttest = list()
alt = c("g","t","l","t","l","t","g","l","g","g","t")
for (i in 1:11) {ttest[[i]] = t.test(trap2[,i], trap4[,i], alternative=alt[i], paired=TRUE)}
names(ttest) = names(species.common)

##############################################################################################################################
# Mood's median test for rare species
##############################################################################################################################

# Load package for processing data
library(plyr)

# Prepare data frame
species.rare = species[,setdiff(names(species), names(species.common))]

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
colors = c("floralwhite", "lightgoldenrod1","chartreuse2", "chartreuse4", "darkgreen")
box.labs = c("Field", "Field edge", "Edge", "Forest edge", "Forest")
layout(matrix(1:12, 4, 3))
par(mar=c(2,2,3,1), mgp=c(2,0.2,0), tck=-0.02)
box.spacing = 1:5/7
x.lims = c(0,box.spacing[5] + box.spacing[1])
long.names2 = c("Aedes canadensis", "Coquillettidea perturbans", "Anopheles crucians", "Aedes hendersoni", "Aedes atlanticus", "Aedes dupreei", "Psorophora cyanescens", "Psorophora howardi","Aedes japonicus", "Aedes infirmatus", "Psorophora ciliata", "Urotanaenia sapphirina")

for (i in names(species.rare)) {
  if (which(names(species.rare)==i) > 4) {ylab = ""} else {ylab = "Adbundance"}
  boxplot(species.rare[,i]~Habitat_f, ylab=ylab, at=box.spacing, boxwex=0.1, pch=20, cex=0.5, medlwd = 1, cex.axis=0.5, xlim=x.lims, col=colors, names=box.labs)
  mtext(long.names2[which(names(species.rare)==i)], line=0.4, cex=0.6, adj=0, font=4)
  mtext(paste("n = ", sum(species.rare[,i]), "; ", mood.results$significance[which(mood.results$species == i)], sep=""), line=0.4, cex=0.5, adj=1)
}

##############################################################################################################################
# END
##############################################################################################################################







