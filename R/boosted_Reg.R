## analyze outputs using BRTs via gbm package implemented in dismo
## https://rspatial.org/raster/sdm/9_sdm_brt.html#google_vignette
## https://blog.devgenius.io/r-tutorial-boosted-regression-trees-9f243d88a921

data(Anguilla_train)
head(Anguilla_train)
angaus.tc5.lr01 <- gbm.step(data=Anguilla_train, gbm.x = 3:13, gbm.y = 2,
                            family = "bernoulli", tree.complexity = 5,
                            learning.rate = 0.01, bag.fraction = 0.5)

names(angaus.tc5.lr01)
summary(angaus.tc5.lr01)

angaus.tc5.lr005 <- gbm.step(data=Anguilla_train, gbm.x = 3:13, gbm.y = 2,
                             family = "bernoulli", tree.complexity = 5,
                             learning.rate = 0.005, bag.fraction = 0.5)
angaus.simp <- gbm.simplify(angaus.tc5.lr005, n.drops = 5)
