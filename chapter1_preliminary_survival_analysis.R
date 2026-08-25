### Chapter 1: Preliminary survival analysis

library(survival)

## Data setup

rm(list=ls())
load("/Users/geruxi/Downloads/The breast cancer survival data/Data to use for project/br_can_data_subset_of_vars.RData")
brcan.data <- data.frame(cbind(stime, cens, age))
brcan.data

## Missing values and censoring


summary(brcan.data)
colSums(is.na(brcan.data))
table(brcan.data$cens)

brcan.data <- brcan.data[complete.cases(brcan.data), ]
brcan.data <- subset(brcan.data, stime > 0 & cens %in% c(0, 1))
summary(brcan.data)
table(brcan.data$cens)

## Kaplan-Meier analysis


km.fit <- survfit(Surv(stime,cens)~1, data=brcan.data,
                  conf.int=0.95, conf.type="log-log")
km.fit
summary(km.fit)
plot(km.fit, xlab="time in days", ylab="survival probability",
     main="Kaplan-Meier estimate for Manchester breast cancer data")

## Cox proportional hazards model


cox.fit <- coxph(Surv(stime,cens)~age, data=brcan.data,
                 ties="exact")
cox.fit
summary(cox.fit)
cox.zph(cox.fit)

## Exponential survival model


exp.fit <- survreg(Surv(stime,cens)~age, data=brcan.data,
                   dist="exponential")
exp.fit
summary(exp.fit)
exp(-coef(exp.fit))

## Comparison of the three models


age0 <- median(brcan.data$age)
new.data <- data.frame(age=age0)

cox.surv <- survfit(cox.fit, newdata=new.data)

t.grid <- seq(0, max(brcan.data$stime), length=500)
lp.exp <- predict(exp.fit, newdata=new.data, type="lp")
s.exp <- exp(-t.grid/exp(lp.exp))


plot(km.fit, conf.int=FALSE,
     xlab="time in days", ylab="survival probability",
     main=paste("Estimated survival functions, age =", age0),
     lwd=2, col="black")
lines(cox.surv, conf.int=FALSE, col="red", lwd=2)
lines(t.grid, s.exp, col="blue", lwd=2, lty=2)
legend("bottomleft",
       legend=c("Kaplan-Meier unconditional",
                paste("Cox PH, age =", age0),
                paste("Exponential, age =", age0)),
       col=c("black", "red", "blue"),
       lty=c(1, 1, 2),
       lwd=2,
       bty="o")

## Model summaries and diagnostics

exp(coef(cox.fit))
exp(confint(cox.fit))

exp(-coef(exp.fit)["age"])

cox.zph(cox.fit)

AIC(cox.fit)
AIC(exp.fit)
