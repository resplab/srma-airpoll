
# Part 1: Data cleaning

# Converts OR and HR into RR, and leaves RR unconverted
# unit: unit of effect size measure, measured as odds ratio (OR), hazard ratio (HR), or relative risk (RR).
# ee: effect size, as column in table
# i: incidence rate (person per year)
rrConvert <- function(unit, ee, i) {case_when(
    unit == "OR" ~ ee/((1-i)+(i*ee)),
    unit == "HR" ~ (1 - exp(ee*log(1-i)))/i,
    unit == "RR" ~ ee
    )
}


# Adjusts to ug/m^3 increment for effect size
# lee: log effect size 
# ie: original increment of exposure (ug/m^3)
# i: target increment (ug/m^3), in our case it will be 5 for all.
adjInc <- function(lee, ie, i){
    lee * i/ie
}

# Adjusts to PM2.5 (per 5 ug/m3)
# lee: log effect size 
# ie: original increment of exposure (ug/m^3)
adjIncPM <- function(lee, ie){
  lee * 5/ie
}

# Adjusts to NO2 (per 10 ug/m3)
# lee: log effect size 
# ie: original increment of exposure (ug/m^3)
adjIncNO <- function(lee, ie){
  lee * 10/ie
}

# Adjusts to O3 (per 60 ug/m3)
# lee: log effect size 
# ie: original increment of exposure (ug/m^3)
adjIncO <- function(lee, ie){
  lee * 60/ie
}


# Meta analysis conversions, adds 8 columns for 
# df: meta-analysis data source
# u: unit of effect size measure (column unit, OR, HR, RR)
# e: effect size (column ee)
# u95: effect size upper 95% CI (column up95)
# l95: effect size lower 95% CI (column low95)
# io: study incidence rate (column incidence)
# ic: study increment (column increment)
# ti: target increment (define as number, in ug/m^3)

metaPrep <- function(df, u, e, u95, l95, io, ic, ti){
    mutate(log.EE = log(e)) %>%
    mutate(log.SE = (log(u95)-log(l95))/3.92) %>%
    mutate(alog.EE = adjInc(log.EE, ic, ti)) %>%
    mutate(alog.SE = adjInc(log.SE, ic, ti)) %>%
    mutate(log.rEE = log(rrConvert(u, e, io))) %>%
    mutate(log.rSE = (log(rrConvert(u, u95, io)) - log(rrConvert(u, l95, io)))/3.92) %>%
    mutate(alog.rEE = adjInc(log.rEE, ic, ti)) %>%
    mutate(alog.rSE = adjInc(log.rSE, ic, ti))
  }

# Part 2: Meta-analysis functions
# For step 1 of meta-analysis, pooling of effect sizes (fully adjusted), not standardized to exposure increment.
# df: data
# ti: title for meta-analysis, input as string
meta_analyze_1 <- function(df, ti){ 
  metagen(
    data = df,
    TE = log.EE,
    seTE = log.SE,
    sm = "OR",
    fixed = TRUE,
    random = TRUE,
    method.tau = "DL",
    studlab = study,
    title = ti
  )
}


# For step 2 of meta-analysis, pooling of effect sizes (fully adjusted), standardized to 5 ug/m^3 exposure increment
meta_analyze_2 <- function(df, ti){ 
  metagen(
    data = df,
    TE = alog.EE,
    seTE = alog.SE,
    sm = "OR",
    fixed = TRUE,
    random = TRUE,
    method.tau = "DL",
    studlab = study,
    title = ti
  )
}

# For step 3 of meta-analysis, pooling of effect sizes (fully adjusted), standardized to 5 ug/m^3 exposure increment and converted to relative risk.
meta_analyze_3 <- function(df, ti){ 
  metagen(
    data = df,
    TE = alog.rEE,
    seTE = alog.rSE,
    sm = "RR",
    fixed = FALSE,
    random = TRUE,
    method.tau = "DL",
    studlab = study,
    title = ti
  )
}

# Funnel plot
# ma: input meta-analysis object (from metagen)
# sl: include study lable (true/false)
fun <- function(ma, sl){
  funnel(ma,
           studlab = sl)
}

# Begg's test for funnel plot bias
# ma: input meta-analysis object (from metagen)
funBias <- function(ma){
      metabias(ma, method.bias = "Begg",
               k.min = 5)
      }


# Print forest plot for manuscript
# ma: input meta-analysis object (from metagen)
# fix: include pooled fixed effect size estimate (true/false)
# rand: include pooled random effect size estimate (true/false)
forest_basic <- function(ma, fix, rand) {
  forest(ma,
          prediction = F,
          print.tau2 = T,
          print.chi2 = T,
          sortvar = avg.exp,
          comb.fixed = fix,
          comb.random = rand,
          leftcols = c("study"),
          leftlabs = c ("Study"),
          label.left = "No effect",
          label.right = "Effect"
         )
}

# Print forest plot with RevMan5 format for manuscript
# ma: input meta-analysis object (from metagen)
# xl: lower bound for scale bar on forest plot
# xu: upper bound for scale bar on forest plot
forest_rev5 <- function(ma, xl, xu) {
  forest(ma,
         prediction = F,
         print.tau2 = T,
         xlim = c(xl, xu),
         sortvar = avg.exp,
         layout = "RevMan5",
         label.left = "Decreased risk",
         label.right = "Increased risk" 
  )
}

# Meta-regression for average exposure
# ma: input meta-analysis object (from metagen)
mregExp <- function(ma){
  metareg(ma, ~avg.exp, rm.na = T)
}
# Meta-regression for all other variables
# ma: input meta-analysis object (from metagen)
mregAll <- function(ma){
  metareg(ma, ~avg.exp+ region + sex.pf + case.def + adj.pol, rm.na = T)
}

