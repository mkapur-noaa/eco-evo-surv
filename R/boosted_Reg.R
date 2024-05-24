## analyze outputs using BRTs via gbm package implemented in dismo
## https://rspatial.org/raster/sdm/9_sdm_brt.html#google_vignette
## https://blog.devgenius.io/r-tutorial-boosted-regression-trees-9f243d88a921

data(Anguilla_train)
head(Anguilla_train)


## x is predictors, y is response
# angaus.tc5.lr01 <- gbm.step(data=Anguilla_train, gbm.x = 3:13, gbm.y = 2,
#                             family = "bernoulli", tree.complexity = 5,
#                             learning.rate = 0.01, bag.fraction = 0.5)

mre_all.01 <- gbm.step(data=mre_all_mod, gbm.x = c(1,7,9,26,27), gbm.y = 18,
                            family = "bernoulli", tree.complexity = 5,
                            learning.rate = 0.01, bag.fraction = 0.5)

names(angaus.tc5.lr01)
summary(angaus.tc5.lr01)

angaus.tc5.lr005 <- gbm.step(data=Anguilla_train, gbm.x = 3:13, gbm.y = 2,
                             family = "bernoulli", tree.complexity = 5,
                             learning.rate = 0.005, bag.fraction = 0.5)
angaus.simp <- gbm.simplify(angaus.tc5.lr005, n.drops = 5)

gbm.plot(angaus.tc5.lr005)

plot(gbm::pretty.gbm.tree(angaus.tc5.lr005))
gbm.perspec(angaus.tc5.lr005, 7, 1, y.range=c(15,20), z.range=c(0,0.6), col = 'blue')






library(lme4); library(sjPlot)

mod1 <- lmer(abs(MRE) ~
               (1|year) +
               # year +
               trophic +
               # scenario+
               # cv_median+
               climate_change +
               evolution +
               factor(replicate),
             data = mre_all_mod)
summary(mod1)
plot(mod1)
sjPlot::plot_model(mod1, type = 'pred')
sjPlot::tab_model(mod1,
                  show.re.var= TRUE,
                  dv.labels= "Effects of Predictors on MRE_SSB")
sjPlot::plot_model(mod1,
                   # axis.labels=c("Urchin", "Depth", "Fish"),
                   show.values=TRUE,
                   show.p=TRUE)

ranef(mod1)
lattice::dotplot(fixef(mod1))
lattice::dotplot(ranef(mod1))



mod2 <- lmer(abs(MRE) ~  species + climate_change + factor(replicate) +
               (1|year),
             data = mre_all_mod)

lattice::dotplot(fixef(mod2))
plot_model(mod2, type = 'pred')

mod3 <- lmer(abs(MRE) ~  species + scenario + factor(replicate) + (1|year),
             data = mre_all_mod)
lattice::dotplot(fixef(mod3))
plot_model(mod3, type = 'pred')


mod4 <- lmer(abs(MRE) ~  species + evolution + factor(replicate) + (1|year),
             data = mre_all_mod)
lattice::dotplot(fixef(mod4))
plot_model(mod4, type = 'pred')
