####----- EDA functions ------------------------

# For creating bar graphs of the categorical variables
mybarplot <- function(var, palette="Set2",
                      xlab="", ylab="",title="",dat='', flip_axis=F) 
{
  # palette=c("Set2","Set3","Accent","Set1","Paired","Dark2")
  xlab <- ifelse(xlab=="", deparse(substitute(var)), xlab)
  ylab <- ifelse(ylab=="", "Number of observations", ylab)
  
  # prepare data frame for plotting
  plotdf <- dat %>%
    group_by({{ var }}) %>%
    summarise(n=n()) %>%
    mutate(pct = n/sum(n),
           lbl = percent(pct)) 
  
  
  # set extra layers to include
  other_layers <- list(
    if (flip_axis) coord_flip()
  )
  
  # create the plot
  return(
    ggplot(plotdf, aes({{ var }},n, fill={{ var }})) +
      geom_bar(stat = "identity",position = "dodge")  +
      geom_text(aes(label = n), size=3, position = position_stack(vjust = 0.5)) +
      scale_fill_brewer(palette=palette) +
      labs(x=xlab, y=ylab, title=title) +
      theme_classic() +
      theme(legend.position = "none") +
      other_layers
  )
}

#===============================
# get distribution of one continuous and one categorical

# my_boxplot <- function(x_var, fill_var, dat=NULL)
# {
#   ggplot(dat, aes({{x_var}}, fill={{fill_var}})) +
#     geom_boxplot() +
#     scale_fill_brewer(palette="Set2") +
#     facet_grid(rows = vars({{fill_var}}), labeller = label_both) +
#     theme_bw() +
#     theme(legend.position = "none",
#           # remove axis ticks on y
#           axis.text.y = element_blank(),
#           axis.ticks.y = element_blank()
#     )
# }

my_boxplot <- function(x_var, y_var, dat=NULL) 
{
  ggplot(dat, aes({{x_var}}, {{y_var}}, fill={{x_var}})) +
    geom_boxplot() +
    scale_fill_brewer(palette="Set2") +
    theme_bw() +
    theme(legend.position = "none") 
}


####----- Modeling functions ------------------------
fit_mod <- function(response, data, trans_y = c("none", "log", "sqrt", "cubert"))
{
  # get the formula that matches any transformation of the response
  transform_y <- match.arg(trans_y)
  if(transform_y == "log") {
    model_formula <- formula(paste0("log(", quo_name(enquo(response)), ")", "~", "pred_value"))
  } else if(transform_y == "sqrt") {
    model_formula <- formula(paste0("sqrt(", quo_name(enquo(response)), ")", "~", "pred_value"))
  } else if(transform_y == "cubert") {
    model_formula <- formula(paste0("(", quo_name(enquo(response)), "^(1/3))", "~", "pred_value"))
  } else {
    # default formula without any transformation (when trans_y is set to 'none' or not specified)
    model_formula <- formula(paste0(quo_name(enquo(response)), "~", "pred_value")) 
  }
  
  # Create a nested dataset for convenient modeling of all predictors.
  suppressWarnings(models <- data %>%
    # Select variables of interest
    select(ends_with("_avg"), Tribal_Community, Rural_Score, is_metro_micro,
           HHA_Score) %>%
    # Transform data to long format to aid nesting
    pivot_longer(-ends_with("_avg"), names_to="predictor", values_to="pred_value") %>%
    nest(-predictor) %>%
    # Nest the continuous predictor separately and attach to the previously nested data:
    # this step is necessary since LIA_CT_PP is numeric and can't be added to factors in
    # same column
    bind_rows(merged_TX_sub %>%
                select(ends_with("_avg"), LIA_CS_PP) %>%
                mutate(predictor = rep("LIA_CS_PP", nrow(merged_TX_sub))) %>%
                rename(pred_value = LIA_CS_PP) %>%
                nest(-predictor))%>%
    # fit individual models all at once
    mutate(model = map(data, ~ lm(model_formula, data = .))))
  
  # Extract metrics (e.g. R-squared, MSE) for model comparison
  metrics <- models %>%
    mutate(metric = map(model, glance)) %>%
    unnest(metric) %>%
    select(predictor, r.squared, adj.r.squared, sigma, F.statistic=statistic, p.value, AIC) %>%
    mutate(MSE=sigma^2, .keep="unused", .after="adj.r.squared")
  
  
  return(list(models=models$model, metrics=metrics))
}

# A function for displaying table results for initial models that saves us the trouble of 
# repeating the same extra kable parameters for better display.
# 
init_mdl_tbl <- function(mdl_object, tf="", var="") {
  # tf: transformation
  # var: dependent variable 
  end_text <- c("-"="untransformed", 
                "log"="with log transform",
                "sqrt"="with square root transform",
                "cbrt"="with cube root transform"
                )
  cap_text <- paste("Response:", var, end_text[tf])
  
  kable(mdl_object$metrics, format="latex",booktabs = T, caption = cap_text) %>%
    kable_styling(font_size = 10,  latex_options = c("HOLD_position"))
}

#----------- This function generates tables of model results
mdl_result <- function(model, conf_level = 0.95, type="param")
{ # param: parameter level metrics
  # overall: for overall model metrics
  
  if(type == "param") {
    # extracting model results at the parameter level
    out <- model %>%
      tidy(conf.int=TRUE, conf.level = conf_level) %>%
      kable(format="latex", booktabs = T,
            caption = "Parameter estimates") %>%
      kable_styling(font_size = 10,  latex_options = c("HOLD_position"))
  } else if (type=="overall"){
    # extracting results at the model level
    out  <- model %>%
      glance() %>%
      select(r.squared,adj.r.squared, sigma, F.statistic=statistic, p.value) %>%
      mutate(MSE=sigma^2, .keep="unused", .after="adj.r.squared") %>%
      kable(format="latex", booktabs = T, caption="Overall Model Performance Results") %>%
      kable_styling(font_size = 10,  latex_options = c("HOLD_position"))
  }
  return(out)
}


BP_tbl <- function(model, alpha = 0.05, df=1) 
{
  # Get BP stastistic
  bp <- bptest(model, student = FALSE)
  # Compute the chi-square acceptance region
  chqtest <- qchisq(1 - alpha, df=df) 
  
  out <- data.frame(BP.statistic=bp$statistic, p.value=bp$p.value, chisq=chqtest)
  
  kable(out, format="latex", booktabs = T, caption = "Breusch-Pagan test") %>%
    kable_styling(font_size = 10,  latex_options = c("HOLD_position"))
}

## Normal probability plot of the residuals 
my_NPP <- function(model)
{
  qqnorm(resid(model), xlab = "Expected", ylab = "Residual", main = "")
  qqline(resid(model))
  title("Normal Probability Plot")
}


## Coefficient of correlation between the ordered residuals and their expected values under ## normality
my_CC <- function(model, dat)
{
  corr <- cbind(
    "Residual"        = round(resid(model), 2),
    "Rank"            = rank(resid(model)),
    "Exp. Value under Normality" = round(sqrt(deviance(model) / df.residual(model))*                     qnorm((rank(resid(model)) - 0.375) / (nrow(dat) + 0.25)), 2)
  )
  cor_test <- cor.test(corr[,3], corr[,1])
  
  return(cor_test)
}

### Lowess curve and Regression confidence band
my_LCCBS <- function(x,y,dat, mod)
{
  plot( y~ x, dat, xlab = "Predictor", ylab = "Response")
  title("Lowess Curve and Linear Regression Confidence Bands")
  with(dat, lines(loess.smooth(x, y)), col = "red")
  
  # Gather confidence bands, ordered by x, and add lines to plot
  order_x <- order(as.vector(x))
  ci <- cbind(model.frame(mod), predict(mod, int = "c"))[order_x, ]
  lines(lwr ~ x, ci, col = "blue", lty = "dashed" )
  lines(upr ~ x, ci, col = "blue", lty = "dashed" )
}









