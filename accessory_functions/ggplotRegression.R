#### Linear model  ####
ggplotRegression <- function (fit, x, Title) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste(Title, "//", "Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)), y = "Relative abundance (%)", x= x)
}
