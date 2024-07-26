setwd("~/Carter/Research/Papers/In Progress/Biodiversity and Carbon Tradeoffs/Data")
library(stringr)
library(vegan)
library(dplyr)
library(FD)
library(rtry)


#### data import and cleaning ####
understory <- read.csv("AK_VEG_SUBPLOT_SPP.csv") ## understory data
plotdata <- read.csv("plot_data.csv") ## plot data with associated environmentla information

understory <- understory[which(understory$PLOT %in% plotdata$Plot),] ## pulling out only the understory data in the study area

# understory$plot_yr <- paste(understory$PLOT, understory$INVYR, sep = "-") ## making a plot year column
length(unique(understory$PLOT)) ## 27 plots, exactly what we should have
unique(understory$PLOT) ## these are the plots we have information on

understory <- understory[, colnames(understory) == "PLOT" | 
                           colnames(understory) == "SUBP" | 
                           colnames(understory) == "VEG_SPCD" | 
                           colnames(understory) == "SP_CANOPY_COVER_TOTAL"] ## pulling out: PLOT, SUBP, VEG_SPCD, SP_CANOPY_COVER_TOTAL
table(understory$SUBP) ## looking at the number of subplots

colnames(understory) ## checking that I got it right
understory <- understory %>% 
  group_by(PLOT, VEG_SPCD) %>% 
  summarise(plotavg = mean(SP_CANOPY_COVER_TOTAL)) ## summarizing by transect then by plot
understory <- as.data.frame(understory) ## turning the tibble back into a dataframe

understory <- reshape(understory, idvar = "PLOT", timevar = "VEG_SPCD", direction = "wide") ## reshaping the data
rownames(understory) <- understory$PLOT ## renaming rows to be PlotID
understory$PLOT <- NULL ## removing nonspecies column
colnames(understory) <- sub("plotavg.", "", colnames(understory)) ## changing column names to be species codes

plotdata <- plotdata[order(plotdata$Plot),] ## ordering the plots so the rows line up
plotdata$TotalPlotCover <- apply(understory, 1, sum, na.rm = TRUE) ## getting the total cover per plot
plotdata$Rich <- apply(as.data.frame(ifelse(understory[,] >= 0,1,0)), 1, sum, na.rm = TRUE) ## getting the count of each species

understory[understory == 0] <- 0.000001 ## assigning 0 values a very small value 
understory[is.na(understory)] <- 0 ## did previous step so NAs can be accounted for
understory <- as.matrix(understory) ## turning community data.frame into a matrix

plotdata$simpson <- diversity(understory, index = "simpson") ## adding another response variable
plotdata$aspect <- abs(plotdata$aspect - 180) ## making the aspect value a linear predictor

#### bringing in the trait data ####
Traits <- read.csv("CleanedTraitData.csv") ## pulling in Trait Data
USDA <- read.csv("USDAcodeKey.csv")
SEAK <- read.csv("SEAKspeciesUSDAcodes.csv")

SEAK$latin <- USDA$Scientific.Name.with.Author[match(SEAK$species, USDA$Symbol)] ## matching species codes to give latin names to the species in the data
SEAK$family <- USDA$Family[match(SEAK$species, USDA$Symbol)] ## making a family col to get FG later
SEAK$latin <- word(SEAK$latin, 1,2, sep=" ") ## making the SEAK latin names 2 words

rm(USDA) ## keeping global env clean

table(word(SEAK$latin, 1)) ## double checking these look okay
length(unique(word(SEAK$latin, 1))) ## 138 genera

FG <- read.csv("FunctionalGroups.csv") ## reading in the functional group data
SEAK$FG <- FG$FunctionalGroup[match(SEAK$family, FG$family)] ## finding the codes for the SEAK species
SEAK$family <- NULL;rm(FG) ## cleaning data frames
SEAK <- SEAK[SEAK$species %in% colnames(understory),] ## SEAK represents all species found in the plots

rownames(Traits) <- Traits$code ## turning rownames into the species codes
Traits <- Traits[,c(4:6)] ## cleaning the df to only be the traits 
sum(understory[,colnames(understory) %in% rownames(Traits)], na.rm = TRUE)/sum(understory, na.rm = TRUE) ## only 24% of species cover has trait data
sum(understory[,colnames(understory) %in% SEAK$species[SEAK$FG == "Angiosperm"]], na.rm = TRUE)/sum(understory, na.rm = TRUE) ## angiosperms are 56% of species data
sum(understory[,colnames(understory) %in% rownames(Traits)], na.rm = TRUE)/sum(understory[,colnames(understory) %in% SEAK$species[SEAK$FG == "Angiosperm"]], na.rm = TRUE) ## 42% of angiosperms have trait data

## sensitivity of different species 
understoryAngio <- understory[,colnames(understory) %in% SEAK$species[SEAK$FG == "Angiosperm"]] ## making df for subsets of species included
understoryTrait <- understory[,colnames(understory) %in% rownames(Traits)]
RichAngio <- as.data.frame(ifelse(understoryAngio[,] >= 0.000001,1,0)) ## creating a P/A matrix for species present in subsets of data
RichTrait <- as.data.frame(ifelse(understoryTrait[,] >= 0.000001,1,0))

plotdata$RichAngio <- apply(RichAngio, 1, sum, na.rm = TRUE) ## getting the count of each species
plotdata$RichTrait <- apply(RichTrait, 1, sum, na.rm = TRUE) ## getting the count of each species

plotdata$simpAngio <- diversity(understoryAngio, index = "simpson") ## adding another response variable
plotdata$simpTrait <- diversity(understoryTrait, index = "simpson") ## adding another response variable

plotdata <- plotdata[,c(1,16,18,19,17,20,21,13,2,3,7:10)]
rm(RichAngio);rm(RichTrait);rm(SEAK)

#### modelling ####
## Checking for Collinearity
## Custum function from Zuur for collineariy
# Library files for courses provided by: Highland Statistics Ltd.
# To cite these functions, use:
# Mixed effects models and extensions in ecology with R. (2009).
# Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.
# Copyright Highland Statistics LTD.
# VIF FUNCTION.
# To use:  corvif(YourDataFile)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}


#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}
#END VIF FUNCTIONS

## models
colnames(plotdata)
MyVar <- c("slope", "elev", "aspect", "AMT", "AP", "treecover", "bio50")
corvif(plotdata[,MyVar]) ## cannot include both slides and things that determine slides

rm(corvif);rm(MyVar);rm(myvif)

pseudo.R.squared <- function(model){
  1 - (model$deviance/model$null.deviance)
} ## writing an R squared function 


#### Richness Models ####
## m1 - richness of all species ~ bio50 + other covariates
m1 <- glm(Rich ~ elev + aspect + AMT + AP + bio50,
          data = plotdata,
          family = poisson(link = "log"))
summary(m1) ## 266.72
pseudo.R.squared(m1) ## 0.325
plot(residuals(m1)); abline(h = 0, col = "red", lwd = 2, lty = 2)
qqnorm(residuals(m1), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(m1), col = "red", lwd = 2, lty = 2)

## m2 - richness of angiosperms ~ bio50 + other covariates
m2 <- glm(RichAngio ~ elev + aspect + AMT + AP + bio50,
          data = plotdata,
          family = poisson(link = "log"))
summary(m2) ## 259.42
pseudo.R.squared(m2) ## 0.344
plot(residuals(m2)); abline(h = 0, col = "red", lwd = 2, lty = 2)
qqnorm(residuals(m2), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(m2), col = "red", lwd = 2, lty = 2)

## m3 - richness of trait species ~ bio50 + other covariates
m3 <- glm(RichTrait ~ elev + aspect + AMT + AP + bio50,
          data = plotdata,
          family = poisson(link = "log"))
summary(m3) ## 173.06
pseudo.R.squared(m3) ## 0.205
plot(residuals(m3)); abline(h = 0, col = "red", lwd = 2, lty = 2)
qqnorm(residuals(m3), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(m3), col = "red", lwd = 2, lty = 2)


#### Diversity Models ####
## m4 - simpson ~ normally important covariates + bio 50
m4 <- glm(simpson ~ elev + aspect + AMT + AP + bio50, 
          data = plotdata)
summary(m4) ## AIC -28.7
pseudo.R.squared(m4) ## 0.240
plot(residuals(m4)); abline(h = 0, col = "red", lwd = 2, lty = 2)
qqnorm(residuals(m4), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(m4), col = "red", lwd = 2, lty = 2)

m5 <- glm(simpAngio ~ elev + aspect + AMT + AP + bio50, 
          data = plotdata)
summary(m5) ## AIC -20.471
pseudo.R.squared(m5) ## 0.187
plot(residuals(m5)); abline(h = 0, col = "red", lwd = 2, lty = 2)
qqnorm(residuals(m5), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(m5), col = "red", lwd = 2, lty = 2)

m6 <- glm(simpTrait ~ elev + aspect + AMT + AP + bio50, 
          data = plotdata)
summary(m6) ## AIC -2.187
pseudo.R.squared(m6) ## 0.283
plot(residuals(m6)); abline(h = 0, col = "red", lwd = 2, lty = 2)
qqnorm(residuals(m6), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(m6), col = "red", lwd = 2, lty = 2)


#### Figure 2 ####
newdata <-data.frame(bio50 = seq(min(plotdata$bio50), max(plotdata$bio50), length.out = 30),
                     elev = mean(plotdata$elev),
                     aspect = mean(plotdata$aspect),
                     AMT = mean(plotdata$AMT),
                     AP = mean(plotdata$AP))
preds1 <- predict(m1, newdata, type="response", se.fit = TRUE)
preds2 <- predict(m2, newdata, type="response", se.fit = TRUE)
preds3 <- predict(m3, newdata, type="response", se.fit = TRUE)
preds4 <- predict(m4, newdata, type="response", se.fit = TRUE)
preds5 <- predict(m5, newdata, type="response", se.fit = TRUE)
preds6 <- predict(m6, newdata, type="response", se.fit = TRUE)
u.ci.1 <- preds1$fit + (qt(0.975, df = 21)*preds1$se.fit)
u.ci.2 <- preds2$fit + (qt(0.975, df = 21)*preds2$se.fit)
u.ci.3 <- preds3$fit + (qt(0.975, df = 21)*preds3$se.fit)
u.ci.4 <- preds4$fit + (qt(0.975, df = 21)*preds4$se.fit)
u.ci.5 <- preds5$fit + (qt(0.975, df = 21)*preds5$se.fit)
u.ci.6 <- preds6$fit + (qt(0.975, df = 21)*preds6$se.fit)
l.ci.1 <- preds1$fit - (qt(0.975, df = 21)*preds1$se.fit)
l.ci.2 <- preds2$fit - (qt(0.975, df = 21)*preds2$se.fit)
l.ci.3 <- preds3$fit - (qt(0.975, df = 21)*preds3$se.fit)
l.ci.4 <- preds4$fit - (qt(0.975, df = 21)*preds4$se.fit)
l.ci.5 <- preds5$fit - (qt(0.975, df = 21)*preds5$se.fit)
l.ci.6 <- preds6$fit - (qt(0.975, df = 21)*preds6$se.fit)

par(mfrow = c(1,2))

## Panel A
## all species
plot(Rich ~ bio50, 
     data = plotdata,
     type = "n",
     xlab = "Biomass (Mg ha-1)",
     ylab = "Number of Species",
     ylim = c(0,100),
     las = 1)
polygon(x = c(newdata$bio50, rev(newdata$bio50)),
        y = c(l.ci.1, rev(u.ci.1)),
        col = adjustcolor("black", alpha.f = 0.1),
        border = NA)
lines(newdata$bio50, preds1$fit, lty=c(1), col = "black")
points(Rich ~ bio50, data = plotdata,
       col = "black",
       pch = 16,
       cex = .75)
## only angiosperms
polygon(x = c(newdata$bio50, rev(newdata$bio50)),
        y = c(l.ci.2, rev(u.ci.2)),
        col = adjustcolor("darkgoldenrod", alpha.f = 0.1),
        border = NA)
lines(newdata$bio50, preds2$fit, lty=c(1), col = "darkgoldenrod")
points(RichAngio ~ bio50, data = plotdata,
       col = "darkgoldenrod",
       pch = 16,
       cex = .75)
## only species w/ trait data
polygon(x = c(newdata$bio50, rev(newdata$bio50)),
        y = c(l.ci.3, rev(u.ci.3)),
        col = adjustcolor("dodgerblue4", alpha.f = 0.1),
        border = NA)
lines(newdata$bio50, preds3$fit, lty=c(1), col  = "dodgerblue4")
points(RichTrait ~ bio50, data = plotdata,
       pch = 16,
       col = "dodgerblue4",
       cex = .75)
legend("topright", 
       legend = c("All Species", "Species with Trait Data"), 
       col = c("black", "dodgerblue4"), 
       pch = 19, 
       cex = 1, 
       ncol = 1, 
       bty = "n")
# legend("topright", 
#        legend = c("All Species", "Angiosperms", "Species with Trait Data"), 
#        col = c("black", "darkgoldenrod", "dodgerblue4"), 
#        pch = 19, 
#        cex = 1, 
#        ncol = 1, 
#        bty = "n")
mtext(side = 3, line = .75, cex = 1, "a)", adj = -0.05)

## Panel B
plot(simpson ~ bio50, data = plotdata,
     type = "n",
     ylim = c(0,1),
     xlab = "Biomass (Mg ha-1)",
     ylab = "Simpson Diversity",
     las = 1)
polygon(x = c(newdata$bio50, rev(newdata$bio50)),
        y = c(l.ci.4, rev(u.ci.4)),
        col = adjustcolor("black", alpha.f = 0.1),
        border = NA)
lines(newdata$bio50, preds4$fit, lty=c(2))
points(simpson ~ bio50, data = plotdata,
       pch = 16,
       cex = .75)
# my.slope <- summary(m2)$coef[6,1]
# my.int <- summary(m2)$coef[1,1]
# leg = vector('expression', 1)
# leg[1] = substitute(expression(y == m*x + b), 
#                     list(m = format(my.slope,dig=3),
#                          b = format(my.int,dig=3)))[2]
# legend("bottomright", legend = leg, pch = NA, bty = "n")
## only angiosperms
polygon(x = c(newdata$bio50, rev(newdata$bio50)),
        y = c(l.ci.5, rev(u.ci.5)),
        col = adjustcolor("darkgoldenrod", alpha.f = 0.1),
        border = NA)
lines(newdata$bio50, preds5$fit, lty=c(2), col = "darkgoldenrod")
points(simpAngio ~ bio50, data = plotdata,
       col = "darkgoldenrod",
       pch = 16,
       cex = .75)
## only species w/ trait data
polygon(x = c(newdata$bio50, rev(newdata$bio50)),
        y = c(l.ci.6, rev(u.ci.6)),
        col = adjustcolor("dodgerblue4", alpha.f = 0.1),
        border = NA)
lines(newdata$bio50, preds6$fit, lty=c(2), col  = "dodgerblue4")
points(simpTrait ~ bio50, data = plotdata,
       pch = 16,
       col = "dodgerblue4",
       cex = .75)
# legend("topright", 
#        legend = c("All Species", "Angiosperms", "Species with Trait Data"), 
#        col = c("black", "darkgoldenrod", "dodgerblue4"), 
#        pch = 19, 
#        cex = 1, 
#        ncol = 1, 
#        bty = "n")
mtext(side = 3, line = .75, cex = 1, "b)", adj = -0.05)


coef(m1)-confint(m1)[,1]
coef(m2)-confint(m2)[,1]
coef(m3)-confint(m3)[,1]
coef(m4)-confint(m4)[,1]
coef(m5)-confint(m5)[,1]
coef(m6)-confint(m6)[,1]

rm(m1);rm(m2);rm(m3);rm(m4)
rm(preds1);rm(preds2);rm(preds3);rm(preds4)
rm(l.ci.1);rm(l.ci.2);rm(l.ci.3);rm(l.ci.4)
rm(u.ci.1);rm(u.ci.2);rm(u.ci.3);rm(u.ci.4)
rm(leg);rm(my.int);rm(my.slope)

#### understory NMDS for supplement ####
set.seed(1)
NMDSord <- metaMDS(understory, try = 1000, distance = "bray") ## convergence
env <- data.frame(bio = plotdata$bio50) ## making an enviroment dataframe
rownames(env) <- plotdata$Plot ## making sure plots align
env$BioFactor <- factor(ifelse(env$bio >= quantile(plotdata$bio50, probs = seq(from = 0, to = 1, by = 0.05))[14], ">= 65 percentile Biomass", "< 65 percentile Biomass"))

par(mfrow = c(3,1))
plot(NMDSord, type = "n",display = "sites") 
points(NMDSord, display = "sites", pch = 19, cex = .75, col = factor(env[,2]))
ordiellipse(NMDSord, display = "sites", factor(env[,2]), draw = "lines",
            col = c(1,2), label = FALSE) ## ellipses based on biomass
leg <- c(as.expression(bquote("" >= 65*"th percentile biomass")),
         as.expression(bquote("" < 65*"th percentile biomass")))
legend("topleft", legend = leg, col = c(2,1), pch = 19, cex = 1, ncol = 1, bty = "n")
mtext(side = 3, line = .75, cex = 1, "a)", adj = -0.05)

env <- data.frame(bio = plotdata$bio50) ## making an enviroment dataframe
rownames(env) <- plotdata$Plot ## making sure plots align
env$BioFactor <- factor(ifelse(env$bio >= quantile(plotdata$bio50)[4], ">= 75 percentile Biomass", "< 75 percentile Biomass"))
plot(NMDSord, type = "n",display = "sites") 
points(NMDSord, display = "sites", pch = 19, cex = .75, col = factor(env[,2]))
ordiellipse(NMDSord, display = "sites", factor(env[,2]), draw = "lines",
            col = c(1,2), label = FALSE) ## ellipses based on biomass
leg <- c(as.expression(bquote("" >= 75*"th percentile biomass")),
         as.expression(bquote("" < 75*"th percentile biomass")))
legend("topleft", legend = leg, col = c(2,1), pch = 19, cex = 1, ncol = 1, bty = "n")
mtext(side = 3, line = .75, cex = 1, "b)", adj = -0.05)


env <- data.frame(bio = plotdata$bio50) ## making an enviroment dataframe
rownames(env) <- plotdata$Plot ## making sure plots align
env$BioFactor <- factor(ifelse(env$bio >= quantile(plotdata$bio50, probs = seq(from = 0, to = 1, by = 0.05))[18], ">= 85 percentile Biomass", "< 85 percentile Biomass"))
plot(NMDSord, type = "n",display = "sites") 
points(NMDSord, display = "sites", pch = 19, cex = .75, col = factor(env[,2]))
ordiellipse(NMDSord, display = "sites", factor(env[,2]), draw = "lines",
            col = c(1,2), label = FALSE) ## ellipses based on biomass
leg <- c(as.expression(bquote("" >= 85*"th percentile biomass")),
         as.expression(bquote("" < 85*"th percentile biomass")))
legend("topleft", legend = leg, col = c(2,1), pch = 19, cex = 1, ncol = 1, bty = "n")
mtext(side = 3, line = .75, cex = 1, "c)", adj = -0.05)


stressplot(NMDSord)
rm(env);rm(newdata);rm(NMDSord);rm(leg)

#### CWM trait analysis ####
understoryTrait <- understoryTrait[,order(colnames(understoryTrait))]
Traits <- Traits[rownames(Traits) %in% colnames(understoryTrait),]
Traits <- Traits[order(rownames(Traits)),]

CWM <- functcomp(Traits, understoryTrait)
CWM$bio <- plotdata$bio50
CWM$treecover <- plotdata$treecover

summary(m1 <- glm(LDMC ~ bio, data= CWM))
summary(m2 <- glm(SLA ~ bio, data = CWM))
summary(m3 <- glm(height ~ bio, data = CWM))


#### Figure 3 ####
newdata <- data.frame(bio = seq(min(CWM$bio), max(CWM$bio), length.out = 30))
preds1 <- predict(m1, newdata, type = "response", se.fit = TRUE)
preds2 <- predict(m2, newdata, type = "response", se.fit = TRUE)
preds3 <- predict(m3, newdata, type = "response", se.fit = TRUE)

par(mfrow = c(3,1), oma = c(0,3,0,0))
# par(mfrow = c(1,1))
plot(LDMC ~ bio, data = CWM,
     ylim = c(0.05,0.4),
     xlab = "Biomass (Mg ha-1)",
     ylab = "",
     las = 1,
     pch = 16,
     cex = 0.75,
     cex.lab = 1.5,
     cex.axis = 1.5)
# mtext("LDMC (g/g)", side = 2, line = 4, cex = 1.5)
lines(newdata$bio, preds1$fit,
      lty = 1,
      lwd = 2,
      col = "black")
polygon(x = c(newdata$bio, rev(newdata$bio)),
        y = c(preds1$fit-(preds1$se.fit*qt(0.975, df = 26)), rev(preds1$fit+(preds1$se.fit*qt(0.975, df = 26)))),
        col = adjustcolor("black", alpha.f = 0.1),
        border = NA)
my.slope <- summary(m1)$coef[2,1]
my.int <- summary(m1)$coef[1,1]
leg = vector('expression', 1)
leg[1] = substitute(expression(y == m*x + b), 
                    list(m = format(my.slope,dig=3),
                         b = format(my.int,dig=3)))[2]
legend("topright", legend = leg, pch = NA, bty = "n", cex = 1.5)
mtext(side = 3, line = .75, cex = 1, "a)", adj = -0.05)

plot(SLA ~ bio, data = CWM,
     xlab = "Biomass (Mg ha-1)",
     ylab = "",
     las = 1,
     ylim = c(-10, 200),
     pch = 16,
     cex = 0.75,
     cex.lab = 1.5,
     cex.axis = 1.5)
mtext(side = 2, line = 4, cex = 1.5, as.expression(bquote("SLA" ~ (mm^2/mg))), adj = 0.5)
lines(newdata$bio, preds2$fit,
      lty = 2,
      lwd = 2,
      col = "black")
polygon(x = c(newdata$bio, rev(newdata$bio)),
        y = c(preds2$fit-(preds2$se.fit*qt(0.975, df = 26)), rev(preds2$fit+(preds2$se.fit*qt(0.975, df = 26)))),
        col = adjustcolor("black", alpha.f = 0.1),
        border = NA)
my.slope <- summary(m2)$coef[2,1]
my.int <- summary(m2)$coef[1,1]
leg = vector('expression', 1)
leg[1] = substitute(expression(y == m*x + b), 
                    list(m = format(my.slope,dig=3),
                         b = format(my.int,dig=3)))[2]
legend("topright", legend = leg, pch = NA, bty = "n", cex = 1.5)
mtext(side = 3, line = .75, cex = 1, "b)", adj = -0.05)

plot(height ~ bio, data = CWM,
     xlab = "Biomass (Mg ha-1)",
     ylab = "",
     las = 1,
     pch = 16,
     ylim = c(0, 30),
     cex = 0.75,
     cex.lab = 1.5,
     cex.axis = 1.5)
# mtext("Height (mm)", side = 2, line = 4, cex = 1.5)
lines(newdata$bio, preds3$fit,
      lty = 2,
      lwd = 2,
      col = "black")
polygon(x = c(newdata$bio, rev(newdata$bio)),
        y = c(preds3$fit-(preds3$se.fit*qt(0.975, df = 26)), rev(preds3$fit+(preds3$se.fit*qt(0.975, df = 26)))),
        col = adjustcolor("black", alpha.f = 0.1),
        border = NA)
my.slope <- summary(m3)$coef[2,1]
my.int <- summary(m3)$coef[1,1]
leg = vector('expression', 1)
leg[1] = substitute(expression(y == m*x + b), 
                    list(m = format(my.slope,dig=3),
                         b = format(my.int,dig=3)))[2]
legend("topright", legend = leg, pch = NA, bty = "n", cex = 1.5)
mtext(side = 3, line = .75, cex = 1, "c)", adj = -0.05)

rm(m1);rm(m2);rm(m3);rm(newdata);rm(preds1);rm(preds2);rm(preds3)
rm(leg);rm(my.int);rm(my.slope)



