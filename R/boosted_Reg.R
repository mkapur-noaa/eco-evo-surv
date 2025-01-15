## analyze outputs using BRTs via gbm package implemented in dismo
## https://rspatial.org/raster/sdm/9_sdm_brt.html#google_vignette
## https://blog.devgenius.io/r-tutorial-boosted-regression-trees-9f243d88a921
#
# data(Anguilla_train)
# head(Anguilla_train)
#
#
# ## x is predictors, y is response
# # angaus.tc5.lr01 <- gbm.step(data=Anguilla_train, gbm.x = 3:13, gbm.y = 2,
# #                             family = "bernoulli", tree.complexity = 5,
# #                             learning.rate = 0.01, bag.fraction = 0.5)
#
# mre_all.01 <- gbm.step(data=full_df, gbm.x = c(1,7,9,26,27), gbm.y = 18,
#                             family = "bernoulli", tree.complexity = 5,
#                             learning.rate = 0.01, bag.fraction = 0.5)
#
# names(angaus.tc5.lr01)
# summary(angaus.tc5.lr01)
#
# angaus.tc5.lr005 <- gbm.step(data=Anguilla_train, gbm.x = 3:13, gbm.y = 2,
#                              family = "bernoulli", tree.complexity = 5,
#                              learning.rate = 0.005, bag.fraction = 0.5)
# angaus.simp <- gbm.simplify(angaus.tc5.lr005, n.drops = 5)
#
# gbm.plot(angaus.tc5.lr005)
#
# plot(gbm::pretty.gbm.tree(angaus.tc5.lr005))
# gbm.perspec(angaus.tc5.lr005, 7, 1, y.range=c(15,20), z.range=c(0,0.6), col = 'blue')
#
#
#
#


library(lme4); library(sjPlot)

# full_df$species <- factor(full_df$species, levels = sppLabs2$Var5[sppLabs2$Var4])
full_df <- merge(full_df, sppLabs2, by.x = 'species', by.y = 'Var2')
full_df <- merge(full_df, scenLabs2, by.x = 'scenario', by.y = 'Var2')
full_df$replicate <- factor(full_df$replicate)
full_df$climate_change <- factor(grepl('\\bCC',full_df$scenario))
full_df$evolution <- factor(grepl('_Evo',full_df$scenario))

mod1 <- lmer(abs(MRE_totbio) ~
               (1|year) +
               # year +
               Var5 + ## renamed spp
               Var1.y+ ## renamed scen
               # cv_median+
               # climate_change +
               # evolution +
               factor(replicate),
             data = full_df)

## pull stuff to plot accordingly
mod_vals <- data.frame(param = c(sppLabs2$Var5[sppLabs2$Var4],
                     scenLabs2$Var1),

           coef = as.vector(c(fixef(mod1)[1], fixef(mod1)[1]+fixef(mod1)[2:3],
                    fixef(mod1)[1], fixef(mod1)[1]+fixef(mod1)[4:6])),
           se = as.vector(summary(mod1)$coefficients[1:7,2])) %>%
  mutate(lwr = coef - 1.96*se, upr = coef+1.96*se,
         categ = c(rep('spp',3),rep('scen',4))) %>%
  dplyr::select(param, categ, coef, lwr, upr)

mod_vals$param<- factor(mod_vals$param, levels = mod_vals$param)

ggplot(mod_vals, aes(x = param, y = coef, color = param)) +
  geom_point() +
  coord_flip()+
  scale_color_manual(values = c('grey20','grey30','grey40',
                                rev(scenPal))) +
  theme(legend.position = 'none')+
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0) +
  # facet_grid(categ~., scales = 'free')+
  labs(x = '', y = 'Predicted Values of Relative Error')

ggsave(last_plot(),
       file = here::here('figs','manuscript','lmer.png'),
       width = 6, height = 4, dpi = 500, unit = 'in')


summary(mod1)
plot(mod1)
sjPlot::plot_model(mod1, type = 'pred')
sjPlot::tab_model(mod1,
                  show.re.var= TRUE,
                  dv.labels= "Effects of Predictors on Relative Error")
sjPlot::tab_model(mod1, file = "tables/model_coef.xls")

sjPlot::plot_model(mod1,
                   # axis.labels=c("Urchin", "Depth", "Fish"),
                   show.values=TRUE,
                   show.p=TRUE)

ranef(mod1)
lattice::dotplot(fixef(mod1))
lattice::dotplot(ranef(mod1))



mod2 <- lmer(abs(MRE_totbio) ~  species + climate_change + factor(replicate) + (1|year),
             data = full_df)

lattice::dotplot(fixef(mod2))
plot_model(mod2, type = 'pred')

mod3 <- lmer(abs(MRE_totbio) ~  species + scenario + factor(replicate) + (1|year),
             data = full_df)
lattice::dotplot(fixef(mod3))
plot_model(mod3, type = 'pred')

mod4 <- lmer(abs(MRE_totbio) ~  interaction(species,scenario) + factor(replicate) + (1|year),
             data = full_df)
# mod4 <- lmer(abs(MRE) ~  species + evolution + factor(replicate) + (1|year),
#              data = full_df)
lattice::dotplot(fixef(mod4))
plot_model(mod4, type = 'pred')
