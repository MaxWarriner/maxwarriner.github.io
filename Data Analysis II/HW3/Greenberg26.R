# https://www.nature.com/articles/s41893-025-01736-0#Fig2
# https://springernature.figshare.com/articles/dataset/County-Level_Environmental_and_Social_Influences_on_EPA_Criminal_Prosecutions_in_the_United_States/26180530?file=55989278
library(tidyverse)
dat.greenberg <- read_csv("Greenberg26.csv") |>
  mutate(EPAregion = factor(EPAregion))

mod.p <- glm(
  FinalEC ~ eqi_2jan2018_vc +
    pctnonwhite10 +
    metro +
    gdp2017b +
    fac_penalty_count +
    CIDDist +
    EPAregion,
  data = dat.greenberg,
  family=poisson(),
  control = glm.control(maxit = 100)
)

mod.nb <- MASS::glm.nb(
  FinalEC ~ eqi_2jan2018_vc +
    pctnonwhite10 +
    metro +
    gdp2017b +
    fac_penalty_count +
    CIDDist +
    EPAregion,
  data = dat.greenberg,
  control = glm.control(maxit = 100)
)

library(lmtest) # note this is a conservative test
lrtest(mod.p, mod.nb)

library(emmeans) #redo this TERRIBLE plot!
plot(emmeans(mod.nb, ~eqi_2jan2018_vc, 
        at=list(eqi_2jan2018_vc = seq(-3,3,1))))



mod.p2 <- glm(
  FinalEC ~ water_eqi_2jan2018_vc + 
          land_eqi_2jan2018_vc +
          air_eqi_2jan2018_vc +
          sociod_eqi_2jan2018_vc +
          pctnonwhite10 +
          metro +
          gdp2017b +
          fac_penalty_count +
          CIDDist + 
          EPAregion,
  data = dat.greenberg,
  family=poisson(),
  control = glm.control(maxit = 100)
)

mod.nb2 <- MASS::glm.nb(
  FinalEC ~ water_eqi_2jan2018_vc + 
          land_eqi_2jan2018_vc +
          air_eqi_2jan2018_vc +
          sociod_eqi_2jan2018_vc +
          pctnonwhite10 +
          metro +
          gdp2017b +
          fac_penalty_count +
          CIDDist + 
          EPAregion,
  data = dat.greenberg,
  control = glm.control(maxit = 100)
)

library(lmtest) # note this is a conservative test
lrtest(mod.p2, mod.nb2)

library(nonnest2)
# Run the Clarke test
vuongtest(mod.nb, mod.nb2)

library(emmeans)
plot(emmeans(mod.nb, ~eqi_2jan2018_vc, 
        at=list(eqi_2jan2018_vc = seq(-3,3,by=1)),
      type = "response"))


library(ggeffects)
plot(ggeffect(mod.nb, ~eqi_2jan2018_vc))
