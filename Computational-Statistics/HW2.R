library(dplyr)
library(leaps)
library(ggplot2)

  baseball <- read.table("D:\\2학기\\통계계산특론1\\HW2\\baseball.txt", header = T, sep = " ")
  str(baseball)
  summary(baseball)

# step
  
  full <- lm(log(salary) ~ ., baseball) ; summary(full) ; anova(full)
  null <- lm(log(salary) ~ 1, baseball) ; summary(null)

  extractAIC(full) 
  extractAIC(null)
  
  bothstep <- step(full, direction = "both") ; summary(bothstep) ; anova(bothstep)
  backstep <- step(full, direction = "backward") ; summary(backstep) 
  forstep <- step(null, scope = list(lower=null, upper=full), direction = "forward") ; summary(forstep) ; anova(forstep)
  
  extractAIC(bothstep)
  extractAIC(backstep)
  extractAIC(forstep)

  
# all possible
  
  allsubs <- regsubsets(log(salary) ~ ., data=baseball, nvmax=27)
  summary(allsubs)
  
  plot(allsubs, scale = "Cp")
  
  allsubsinfo <- data.frame(p=2:28, summary(allsubs)$outmat, Cp=summary(allsubs)$cp)
  allsubsinfo %>% filter(rank(Cp) %in% 1:3)

  sub12 = lm(log(salary) ~ average + runs + rbis + sos + freeagent + 
                           arbitration + runsperso + hitsperso + soserrors + sbsobp + sbsruns, baseball)
  sub13 = lm(log(salary) ~ obp + runs + triples + rbis + sos + freeagent + 
                           arbitration + runsperso + hitsperso + soserrors + sbsobp + sbsruns, baseball)
  sub14 = lm(log(salary) ~ average + runs + triples + rbis + sos + freeagent + 
                           arbitration + runsperso + hitsperso + hrsperso + soserrors + sbsobp + sbsruns, baseball)
  
  extractAIC(sub12)
  extractAIC(sub13)
  extractAIC(sub14)
  
  ggplot(filter(allsubsinfo, p>=4), aes(p, Cp)) + geom_line(size=1, color="grey") + geom_point(size=2) + 
    geom_vline(aes(xintercept=13), linetype="dashed") + geom_hline(aes(yintercept=5.196692), linetype="dashed") +
    theme_test() + theme(aspect.ratio=3/4, axis.title = element_text(size=12))
   