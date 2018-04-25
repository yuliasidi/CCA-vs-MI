library(MASS)
library(norm)
library(dplyr)
library(purrr)
library(mice)
library(ggplot2)

setwd()

set.seed(987586)

#Specify sample size
ssf <- 40

#Specify mean vector and covariance matrix
M <- matrix( c(10,13)  , nrow = 1, ncol = 2, byrow = T)
S <- matrix( c(3,0,0,3), nrow = 2, ncol = 2, byrow = T)
y <- data.frame(id = 1:ssf, mvrnorm(n = ssf, mu = M, Sigma = S))

y$X1[y$X1 > 20] <- 20 #if value >max replace with max
y$X1[y$X1 < 0 ] <-  0 #if value <min replace with min
y$X2[y$X2 > 20] <- 20 #if value >max replace with max
y$X2[y$X2 < 0 ] <-  0 #if value <min replace with min

Y <- y%>%
  reshape2::melt(.,id = c('id'), variable.name = 'group')%>%
  mutate(group2 = ifelse(group=='X1',1,2))

# MCAR 25% of the values 
Y.mcar.25 <- Y%>%
  do(.,cbind(.,r.mcar = as.numeric(1:(2*ssf)%in%sample(2*ssf, 2*ssf*0.25, replace = FALSE))))%>%
  mutate(aval = ifelse(r.mcar==0,value,NA))%>%
  select(c('aval', 'group2', 'id'))

# MAR 25% of the values 
Y.mar.25 <- Y%>%
  arrange(group)%>%
  group_by(group)%>%
  mutate(n_check = ifelse(group=='X1',floor(ssf*(1-.15)),floor(ssf*(1-.35))),
         n=1:n(),r.mar = as.numeric(n>n_check))%>%
  mutate(aval = ifelse(r.mar==0,value,NA))%>%ungroup()%>%
  select(c('aval','group2','id'))

# MNAR 25% of the values 
Y.mnar.25 <- Y%>%
  arrange(group,value)%>%
  group_by(group)%>%
  mutate(n_check = ifelse(group=='X1',floor(ssf*(1-.15)),floor(ssf*(1-.35))),
         n=1:n(),r.mnar = as.numeric(n>n_check))%>%
  mutate(aval = ifelse(r.mnar==0,value,NA))%>%ungroup()%>%
  select(c('aval','group2','id'))


# Complete case analysis for each missing data
m1 <- lm(aval ~ group2, Y.mcar.25)
m2 <- lm(aval ~ group2, Y.mar.25)
m3 <- lm(aval ~ group2, Y.mnar.25)

# MI  
rngseed(8687427)

#Specify number of imputaitons
m <- 5

#MI fpr MCAR
data1 <- as.matrix(Y.mcar.25)
s <- prelim.norm(data1)
thetahat <- em.norm(s)
mylm <- NULL

theta <- da.norm(s, thetahat, steps = 20, showits = TRUE)
getparam.norm(s, theta) # look at result

#Data augmentation
da.list <- list()
for (i in 1:1000){
  print(i)
  out <- da.norm(s, thetahat, steps=1)
  theta.hat <- out
  da.list[[i]] <- getparam.norm(s, out)
}

temp<-list()
for (i in 1:1000){temp[[i]] <- da.list[[i]][[1]][2]}

#Output and safe the ACF graphs (autocorrelation graphs)
pdf("DA_MCAR25.pdf")
par(mfrow = c(1,2))
acf(unlist(temp), main = "All Draws - No Burn-in")
acf(unlist(temp[-c(1:100)]), main = "Burn-in 100 Draws")
dev.off()

for(i in 1:m){
  a <- imp.norm(s, thetahat, data1)
  a <- as.data.frame(a)
  mylm[[i]] <- lm(aval ~ group2, data = a)
}
mylm <- as.mira(mylm)
analysis1 <- summary(pool(mylm, method = "smallsample"))


# MI for MAR
data2 <- as.matrix(Y.mar.25)
s2 <- prelim.norm(data2)
thetahat2 <- em.norm(s2)
mylm2 <- NULL
theta2 <- da.norm(s2, thetahat2, steps=20, showits=TRUE)
getparam.norm(s2, theta2) # look at result

##Data augmentation
da.list <- list()
for (i in 1:1000){
  print(i)
  out <- da.norm(s2, thetahat2, steps=1)
  theta.hat <- out
  da.list[[i]] <- getparam.norm(s2, out)
}

temp <- list()
for (i in 1:1000){temp[[i]] <- da.list[[i]][[1]][2]}

#Output and safe the ACF graphs (autocorrelation graphs)
pdf("DA_MAR25.pdf")
par(mfrow = c(1,2))
acf(unlist(temp), main = "All Draws - No Burn-in")
acf(unlist(temp[-c(1:100)]), main = "Burn-in 100 Draws")
dev.off()

for(i in 1:m){
  a2 <- imp.norm(s2, thetahat2, data2)
  a2 <- as.data.frame(a2)
  mylm2[[i]] <- lm(aval ~ group2, data = a2)
}
mylm2 <- as.mira(mylm2)
analysis2 <- summary(pool(mylm2, method = "smallsample"))

#MI for MNAR
data3 <- as.matrix(Y.mnar.25)
s3 <- prelim.norm(data3)
thetahat3 <- em.norm(s3)
mylm3 <- NULL
theta3 <- da.norm(s3, thetahat3, steps = 20, showits = TRUE)
getparam.norm(s3, theta3) # look at result

##Data augmentation
da.list <- list()
for (i in 1:1000){
  print(i)
  out <- da.norm(s3, thetahat3, steps = 1)
  theta.hat <- out
  da.list[[i]] <- getparam.norm(s3, out)
}

temp<-list()
for (i in 1:1000){temp[[i]] <- da.list[[i]][[1]][2]}

#Output and safe the ACF graphs (autocorrelation graphs)
pdf("DA_MNAR25.pdf")
par(mfrow = c(1,2))
acf(unlist(temp), main = "All Draws - No Burn-in")
acf(unlist(temp[-c(1:100)]), main = "Burn-in 100 Draws")
dev.off()

for(i in 1:m){
  a3 <- imp.norm(s3, thetahat3, data3)
  a3 <- as.data.frame(a3)
  mylm3[[i]] <- lm(aval ~ group2, data = a3)
}
mylm3 <- as.mira(mylm3)
analysis3 <- summary(pool(mylm3, method = "smallsample"))

#Combine results
ci.mcar.25.mi <- data.frame(low.95 = analysis1[2,6], 
                            high.95 = analysis1[2,7], model = "MI MCAR 25%")
ci.mar.25.mi  <- data.frame(low.95 = analysis2[2,6]  , 
                            high.95 = analysis2[2,7], model = "MI MAR 25%")
ci.mnar.25.mi <- data.frame(low.95 = analysis3[2,6], 
                            high.95 = analysis3[2,7], model = "MI MNAR 25%")
ci.mar.25.cca <- data.frame(low.95 = confint(m2)[2,1], 
                            high.95 = confint(m2)[2,2], model = "CCA MAR 25%")
ci.mnar.25.cca<- data.frame(low.95 = confint(m3)[2,1], 
                            high.95 = confint(m3)[2,2], model = "CCA MNAR 25%")


conf95.comp <- data.frame(low.95 = confint(m1)[2,1], high.95
                          = confint(m1)[2,2], model = "CCA MCAR 25%")
conf95.comp <- rbind(conf95.comp, ci.mcar.25.mi, ci.mar.25.cca, ci.mar.25.mi,
                   ci.mnar.25.cca, ci.mnar.25.mi)

conf95.comp$len <- conf95.comp$high.95 - conf95.comp$low.95

