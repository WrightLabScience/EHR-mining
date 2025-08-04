library(dplyr)
library(survival)

extractCoxCoefs <- function(model_summary, row_name='treatment', get_coefficients=TRUE, print_conf_vars) {
   if (!get_coefficients) {
      t <- model_summary
   } else {
      t <- model_summary$coefficients
   }
   t <- t[row_name,]
   if ('robust se' %in% names(t)) {
      se <- t['robust se']
   } else {
      se <- t['se(coef)']
   }
   me <- qnorm(0.975) * unname(se)
   est <- unname(t['coef'])
   return(c(
      'est' = exp(est),
      'lower' = exp(est - me),
      'upper' = exp(est + me),
      'pval' = unname(t['Pr(>|z|)'])
   ))
}

modelCoxIPTWcoefs <- function(df, 
                              time_var = 'time', 
                              status_var = NULL, 
                              treatment_var, 
                              censor_time=30L, 
                              conf_vars=NULL, 
                              print_conf_vars=FALSE,
                              normalize_weights=TRUE,
                              one_hot_encode_vars=NULL) {
   stopifnot('Time variable > length 1.' = length(time_var) == 1L)
   if (!is.null(status_var)) stopifnot('Status variable > length 1.' = length(status_var) == 1L)
   stopifnot('Treatment variable > length 1.' = length(treatment_var) == 1L)
   stopifnot('Provided time variable not present in dataset.' = time_var %in% names(df))
   stopifnot('Provided status variable not present in dataset.' = status_var %in% names(df))
   stopifnot('Provided treatment variable not present in dataset.' = treatment_var %in% names(df))
   stopifnot('Provided confounding variable(s) not present in dataset.' = all(conf_vars %in% names(df)))
   stopifnot('Invalid censoring time' = is.integer(censor_time) & censor_time > 0L)
   
   # clean up variable names
   if (time_var != 'time') 
      df <- df %>% rename(time = !!time_var)
   
   if (is.null(status_var)) {
      df <- df %>% 
         mutate(
            status = ifelse(is.na(time) | time > censor_time, 0L, 1L), # 0 is right-censored, 1 is event at time
            time = ifelse(status == 0L, censor_time, time)
         )
      status_var <- 'status'
   }
   
   if (status_var != 'status')
      df <- df %>% rename(status = !!status_var)
   
   if (treatment_var != 'treatment')
      df <- df %>% rename(treatment = !!treatment_var)
   
   if (is.null(conf_vars))
      conf_vars <- df %>% select(-time, -status, -treatment) %>% names()
   
   if (print_conf_vars)
      cat(conf_vars, sep=', ')
   
   # set.seed(123L)
   # var_order <- sample(seq_along(df), length(df), replace=FALSE)
   # df <- df[var_order]
   
   # get unadjusted treatment effect + multivariate reg. as baseline
   uni <- df %>%
      coxph(formula = Surv(time, status) ~ treatment, data = .) %>%
      summary() %>%
      extractCoxCoefs()
   mult <- df %>%
      coxph(formula = Surv(time, status) ~ treatment + ., data = .) %>%
      summary() %>%
      extractCoxCoefs()
   
   
   # fit propensity score model
   prop_model <- df %>%
      select(-time, -status) %>%
      glm(formula = as.formula(paste0('treatment ~ ', paste(conf_vars, collapse=' + '))),
          family = binomial(link='logit'),
          data = .)
   df$prop_score <- predict(prop_model, df, type='response')
   
   # hist(df$prop_score, breaks=40)
   
   df <- df %>% 
      mutate(extreme_prop_score = prop_score < 0.01 | prop_score > 0.99)
   
   print(df %>% 
            select(treatment, extreme_prop_score) %>% 
            table())
   
   cat('Removing', sum(df$extreme_prop_score), 'patients with near 0 or 1 propensity scores.\n')
   
   df_og <- df
   df <- df_og %>%
      filter(!extreme_prop_score) %>%
      select(-extreme_prop_score)
   
   # compute IPT weights - overlap and inverse
   df <- df %>%
      mutate(
         w_SW = case_when(
            treatment ~ (sum(df$treatment) / nrow(df)) / prop_score,
            !treatment ~ (1 - (sum(df$treatment) / nrow(df))) / (1 - prop_score)
         ),
         w_ATE = case_when(
            treatment ~ 1 / prop_score,
            !treatment ~ 1 / (1 - prop_score)
         ),
         w_ATO = case_when( # overlap weights
            treatment ~ 1 - prop_score,
            !treatment ~ prop_score
         )
      )
   
   # normalize weights to sum to original sample size
   if (normalize_weights) {
      df$w_SW <- df$w_SW / (sum(df$w_SW) / nrow(df))
      df$w_ATE <- df$w_ATE / (sum(df$w_ATE) / nrow(df))
      df$w_ATO <- df$w_ATO / (sum(df$w_ATO) / nrow(df))  
   }
   
   # data.frame(SMD = round(cobalt::col_w_smd(mat = as.matrix(df %>% select(-treatment, -time, -status, -prop_score, -w_ATE, -w_ATO, -w_SW)),
   #                   treat = df$treatment,
   #                   weights = df$w_SW),
   #       2))
   
   # get propensity weighted treatment effects
   mult_SW <- df %>%
      select(-w_ATE, -w_ATO, -prop_score) %>%
      coxph(formula = Surv(time, status) ~ treatment + ., data = .,
            weights = w_SW) %>%
      summary() %>%
      extractCoxCoefs()
   mult_ATE <- df %>%
      select(-w_SW, -w_ATO, -prop_score) %>%
      coxph(formula = Surv(time, status) ~ treatment + ., data = .,
            weights = w_ATE) %>%
      summary() %>%
      extractCoxCoefs()
   mult_ATO <- df %>%
      select(-w_SW, -w_ATE, -prop_score) %>%
      coxph(formula = Surv(time, status) ~ treatment + ., data = .,
            weights = w_ATO) %>%
      summary() %>%
      extractCoxCoefs()
   
   return(list(
      'og_data' = df_og,
      'trimmed_data' = df,
      'coefs' = do.call(rbind, 
                        list(
                           'uni' = uni, 
                           'mult' = mult, 
                           'mult_ATE' = mult_ATE,
                           'mult_ATO' = mult_ATO, 
                           'mult_SW' = mult_SW
                        ))
   ))
}


plotCoxCoefs <- function(coefs, treatment_var, par_given=FALSE) {
   stopifnot('coefs should be matrix' = is.matrix(coefs))
   stopifnot('Column names wrong. Should be `est`, `lower`, `upper`' = all(colnames(coefs) == c('est', 'lower', 'upper', 'pval')))
   stopifnot('Values do not look right. `est` should be > `lower` and < `upper`' = all(coefs[,'est'] > coefs[,'lower']) & all(coefs[,'est'] < coefs[,'upper']))
   
   model_dict <- c('uni' = 'univar.', 
                   'mult' = 'multivar.', 
                   'mult_ATE' = 'multivar.\ninverse\nweights',
                   'mult_ATO' = 'multivar.\noverlap\nweights', 
                   'mult_SW' = 'multivar.\nstabilized\nweights')
   rownames(coefs) <- unname(model_dict[rownames(coefs)])
   
   num_vals <- nrow(coefs)
   xvals <- 1:num_vals
   ymax <- max(c(1 / coefs[,colnames(coefs) != 'pval'], coefs[,colnames(coefs) != 'pval']))
   ymin <- 1 / ymax
   
   xaxshift <- 0.25
   if (!par_given) {
      par(mar=c(4,3,1,1), mgp=c(2, 0.4, 0), tck=-0.01)
   }
   plot(NA, xlim=c(1-xaxshift, num_vals + xaxshift), ylim=c(ymin, ymax), log='y', yaxt='n', xaxt='n', xlab='', 
        ylab='Hazards ratio', main='ATE estimates', cex.main=1)
   abline(h=1, lty=2)
   text(x=num_vals + 2.1*xaxshift, srt=270, xpd=NA, y=c(sqrt(prod(1, ymin)), sqrt(prod(1, ymax))), col='#aaaaaa',
        labels=paste0(treatment_var, ' is ', c('better', 'worse')))
   
   axis(side=2, las=1)
   par(mgp=c(2, 0, 0))
   axis(side=1, at=xvals, labels=rownames(coefs), gap.axis=-100, padj=1)
   arrows(x0=xvals, y0=coefs[,'lower'], y1=coefs[,'upper'], code=3, angle=90, length=0.1)
   points(x=xvals, y=coefs[,'est'], pch=16)
   text(x=xvals, y=coefs[,'upper']*1.1, labels=paste0('p=', formatC(coefs[,'pval'], digits=2, format='g')))
}
