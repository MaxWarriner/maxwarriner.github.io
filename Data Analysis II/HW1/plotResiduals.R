plotResiduals<-function(mod){
  require(tidyverse)
  require(patchwork)
  ggdat <- data.frame(y.hat = fitted(mod),
                      e     = residuals(mod),
                      sr    = rstudent(mod))
  
  ggdat.gaussian<-data.frame(x=seq(min(ggdat$e)-0.25*diff(range(ggdat$e)),
                                   max(ggdat$e)+0.25*diff(range(ggdat$e)),
                                   length.out = 1000)) %>%
    mutate(f=dnorm(x,
                   #ei should have mean zero
                   mean=0,
                   #ei should have common variance  (regression estimated)
                   sd=summary(mod)$sigma))
  ####################################
  # Create PDF of Residuals Plot"
  ####################################
  default.bins <- round(log2(nrow(ggdat)+1))
  p1<-ggplot(ggdat,aes(x=e))+
    geom_histogram(aes(y=after_stat(density)),
                   bins=default.bins,
                   fill="lightgrey", color="black")+
    geom_hline(yintercept=0)+
    stat_density(aes(x=e, color="Empirical"),
                 geom="line",position="identity")+
    geom_line(data=ggdat.gaussian,aes(x=x,y=f,color="Gaussian-Assumed"))+
    theme_bw()+
    xlab("Residual")+
    ylab("Density")+
    labs(color = "")+
    theme(legend.position="bottom")
  
  ####################################
  # Create Data for CDF of Residuals Plot"
  ####################################
  e.cdf.func<-ecdf(residuals(mod))
  e.cdf<-e.cdf.func(sort(residuals(mod)))

  ggdat<-data.frame(e=sort(residuals(mod)),
                    e.cdf=e.cdf)
  ####################################
  # Create QQ Plot"
  ####################################
  p2<-ggplot(data=ggdat, aes(sample=e))+ #standardize e
    geom_qq() +
    geom_qq_line() +
    theme_bw()+
    xlab("Gaussian Quantiles")+
    ylab("Sample Quantiles")
  
  ####################################
  # Create Data for Fitted vs Residual Plot"
  ####################################
  ggdat<-data.frame(x=fitted(mod),
                    e=residuals(mod)) %>%
    mutate(id=1:n())
    
  ggdat.out3 <- ggdat %>% 
    filter(abs(e)>3*summary(mod)$sigma)
  
  ####################################
  # Create Fitted vs Residual Plot"
  ####################################
  p3<-ggplot(data=ggdat,aes(x=id,y=e))+
    geom_point(shape=1)+
    geom_hline(yintercept = 0,color="red",linetype="dashed")+
    xlab("Observation Number")+
    ylab("Residual")+
    theme_bw()+
    geom_hline(yintercept = c(-3,3)*summary(mod)$sigma, 
               color="red", linetype="dotted",linewidth=0.75)+
    geom_point(data=ggdat.out3, aes(x=x,y=e), fill="red", shape=21)
  
  ####################################
  # Create Fitted vs Residual Plot"
  ####################################
  p4<-ggplot(data=ggdat,aes(x=x,y=e))+
    geom_point(shape=1)+
    geom_hline(yintercept = 0,color="red",linetype="dashed")+
    xlab(bquote("Fitted Values"~(hat(Y))))+
    ylab("Residual")+
    theme_bw()+
    geom_hline(yintercept = c(-3,3)*summary(mod)$sigma, 
               color="red", linetype="dotted",linewidth=0.75)+
    geom_point(data=ggdat.out3, aes(x=x,y=e), fill="red", shape=21)
  
  ####################################
  # Print Plots"
  ####################################
  (p1|p2)/(p3|p4) +
    plot_layout(guides = 'collect') &
    theme(legend.position = "bottom")
}
