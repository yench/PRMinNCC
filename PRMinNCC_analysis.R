library(data.table)
library(survival)
library(ggplot2)
n = 5000 # Cohort size
mctrl = c(1, 1)  # Number of ctrls/case for each failure type

##################################################################
#                       Proportional risks                       #
##################################################################
# Generate full cohort and NCC data
datncc = rbindlist(attr(datful, "ph2data")) 

# Replicate data and set up cell-mean coding
replincc = repli.dat(datncc)

# Fit model (1)
fitncc1 = coxph(Surv(eventime, status) ~ z1.1 + z2.1 + z1.2 + z2.2 + strata(matched.id, cause),  data = replincc)

## Equivalently, fit model (1) separately to each cause using the original dataset
# fitncc1.1 = coxph(Surv(eventime, ind.fail==1) ~ z1 + z2 + strata(matched.id),  data = datncc)
# fitncc1.2 = coxph(Surv(eventime, ind.fail==2) ~ z1 + z2 + strata(matched.id),  data = datncc)

# Fit model (2)
fitncc2 = coxph(Surv(eventime, status) ~ z1.1 + z2.1 + z1.2 + z2.2 + strata(matched.id) + cause, data = replincc)

# Estimate baseline hazards
bhazncc1 = bhaz.wgt(fitncc1, cause = replincc$cause, wgt = replincc$wgt, matched.id = replincc$matched.id)
bhazncc2 = bhaz.wgt(fitncc2, cause = replincc$cause, wgt = replincc$wgt, matched.id = replincc$matched.id)

# Graphical evaluation of proportionality
ggplot(bhazncc1, aes(x = failtime, y = bHaz))       + geom_step(aes(group = cause, color = cause))
ggplot(bhazncc2, aes(x = failtime, y = bHaz.cause)) + geom_step(aes(group = cause, color = cause))

# True baselines for comparison
length.out = 3000
b10 = log(-log(0.9)/10)
b20 = log(-log(0.95)/10)
t = rep(seq(0, 10, length.out=length.out))
Lam10 = exp(b10)*t
Lam20 = exp(b20)*t
Lam.plot = data.table(t = t, cause = rep(c("1", "2"), each=length(t)), bHaz =c(Lam10, Lam20))
ggplot(data=Lam.plot, aes(x = t, y = bHaz)) + geom_line(aes(group = interaction(cause), color = cause))

##################################################################
#                     Non-proportional risks                     #
##################################################################
# Generate full cohort and NCC data
datful = simul.data(n, mctrl, alpha2 = 5)   
datncc = rbindlist(attr(datful, "ph2data"))

# Replicate data and set up cell-mean coding
replincc = repli.dat(datncc)

# Create a proportionality factor linear in case time
replincc[ , t := 0]
replincc[cause=="2" , t := min(eventime/10), by = .(matched.id)]

# Fit model (1) and (2)
fitncc1 = coxph(Surv(eventime, status) ~ z1.1 + z2.1 + z1.2 + z2.2 + strata(matched.id, cause),  data = replincc)
fitncc2 = coxph(Surv(eventime, status) ~ z1.1 + z2.1 + z1.2 + z2.2 + strata(matched.id) + cause + t, data = replincc)

# Estimate baseline hazards
bhazncc1 = bhaz.wgt(fitncc1, cause = replincc$cause, wgt = replincc$wgt, matched.id = replincc$matched.id)
bhazncc2 = bhaz.wgt(fitncc2, cause = replincc$cause, wgt = replincc$wgt, matched.id = replincc$matched.id)


# graphical evaluation of proportionality
ggplot(bhazncc1, aes(x = failtime, y = bHaz))       + geom_step(aes(group = cause, color = cause))
ggplot(bhazncc2, aes(x = failtime, y = bHaz.cause)) + geom_step(aes(group = cause, color = cause))

##################################################################
#                         Matching on age                        #
##################################################################
# Generate full cohort and NCC data
datful = simul.data(n, mctrl, age.match = TRUE)   
datncc = rbindlist(attr(datful, "ph2data"))      

# Replicate NCC data and set up cell-mean coding
replincc = repli.dat(datncc)
replincc[ , z3.match := ifelse(cause=="2", z3, 0)] 

# Fit model (1) and (2) 
fitncc1 = coxph(Surv(eventime, status) ~ z1.1 + z2.1 + z1.2 + z2.2 + strata(matched.id, cause),             data = replincc)
fitncc2 = coxph(Surv(eventime, status) ~ z1.1 + z2.1 + z1.2 + z2.2 + strata(matched.id) + cause + z3.match, data = replincc)
