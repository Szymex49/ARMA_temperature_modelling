library(dplyr)
library(ggplot2)
library(zeallot)
library(tidyr)
library(reshape2)
library(cowplot)
library(latex2exp)
options(repr.plot.width = 12, repr.plot.height = 8)

analysis <- function(X, Y, display_plot=TRUE){
  spearman <- cor(X, Y, use="pairwise.complete.obs", method="spearman")
  pearson <- cor(X, Y, use="pairwise.complete.obs", method="pearson")
  print(c("Pearson", pearson))
  print(c("Spearman", spearman))
  if (display_plot) {
    ggplot() + geom_point(aes(x=X, y=Y))
  }
}


regression <- function(X, Y, display_plot=TRUE){
  r <- cor(X, Y, use="pairwise.complete.obs")
  Sx <- sd(X)
  Sy <- sd(Y)
  a <- r * Sy / Sx
  b <- mean(Y) - a * mean(X)
  
  if(display_plot){
    xs <- seq(min(X), max(X), 0.01)
    plt <- ggplot() +
      geom_point(aes(x=X, y=Y), alpha=0.5) +
      geom_line(aes(x = xs, y = a * xs + b), linewidth=1, col="red")
    show(plt)
  }
  
  return(c(a, b))
}


cor_plot <- function(Xt, type="acf", h=30){
  cors <- acf(Xt, plot=FALSE, lag.max=h)
  pcors <- pacf(Xt, plot=FALSE, lag.max=h)
  
  if (type == "acf"){
    plt <- ggplot() + 
      geom_point(aes(cors$lag, cors$acf), size=3, col="blue", alpha=0.5) + ggtitle("ACF") +
      geom_line(aes(cors$lag, 0), col="#006bc398", size=1) +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20)) +
      xlab("h") +
      ylab("")
  }
  if (type == "pacf"){
    plt <- ggplot() +
      geom_point(aes(pcors$lag, pcors$acf), size=3, col="blue", alpha=0.5) + ggtitle("PACF") +
      geom_point(aes(0, 1), size=3, col="blue", alpha=0.5) +
      geom_line(aes(cors$lag, 0), col="#006bc398", size=1) +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20)) +
      xlab("h") +
      ylab("")
  }
  return(plt)
}


pvar <- function(X, k, plot=FALSE){
  vars <- c()
  for (i in 1:(length(X) - k)){
    vars <- append(vars, var(X[i:(i+k)]))
  }
  if (plot){
    plt <- ggplot() + geom_line(aes(1:(length(X) - k), vars)) + xlab("Dzień") + ylab("Wariancja cząstkowa")
    return(plt)
  }
}

data <- read.csv("data/data.csv")
data
ggsave("images/raw_data.pdf", plt, width=14, height=6)

plt <- ggplot() + geom_line(aes(x=data$day_number, y=data$temp))  + xlab("Dzień") + ylab("Temperatura")
ggsave("images/raw_data.pdf", plt, width=14, height=6)

new_data <- data[2730:4622,]
day <- new_data$day_number
temp <- new_data$temp

plt <- ggplot() + geom_line(aes(x=day, y=temp))  + xlab("Dzień") + ylab("Temperatura")
ggsave("images/range_data.pdf", plt, width=14, height=6)

length(temp)



ggplot() + geom_line(aes(x=day, y=temp))

plt1 <- cor_plot(temp, h=80)
plt2 <- cor_plot(temp, "pacf", 80)
plt <- plot_grid(plt1, plt2, align="h")
plt
ggsave("images/raw_acf.pdf", plt, width=14, height=6)

temp_diff <- diff(temp)
plt <- ggplot() +
  geom_line(aes(1:length(temp_diff), temp_diff), col="#0065d9b1") +
  xlab("t") +
  ylab(TeX("$\\tilde{X}_t$")) +
  geom_line(aes(1:length(temp_diff), 0), col="#00000073", size=1) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
plt
ggsave("images/diff_temp.pdf", plt, width=12, height=8)

p_max <- 8
q_max <- 8

AIC_values <- matrix(nrow=p_max+1, ncol=q_max+1)

for (p in 0:p_max){
  for (q in 0:q_max){
    model <- arima(temp_diff, order = c(p, 0, q))
    AIC_values[p+1, q+1] <- model$aic
  }
}

indices <- which(AIC_values == min(AIC_values), arr.ind = TRUE)
p_found <- indices[1]-1
q_found <- indices[2]-1
c(p_found, q_found)

model <- arima(temp_diff, order = c(3, 0, 4), method = 'ML')

print(model)

coefficients <- model$coef
coef_sigma2 <- model$sigma2
AIC_value <- model$aic
residuals <- model$residuals

M <- 1000
h <- 10
alpha <- 0.05

cor_matrix <- matrix(0, nrow = M, ncol = h)
pcor_matrix <- matrix(0, nrow = M, ncol = h)

for (j in 1:M){
  
  Xt <- arima.sim(model = list(order = c(3, 0, 4), ar = c(0.7874, 0.8668, -0.664), ma = c(-0.747, -1.1119, 0.655, 0.2111)), n = 1000, sd = sqrt(16.39))
  
  cors <- acf(Xt, lag.max=h, plot=FALSE)$acf
  cors <- cors[2:(length(cors))]
  cor_matrix[j,] <- cors
  
  pcors <- pacf(Xt, lag.max=h, plot=FALSE)$acf
  pcor_matrix[j,] <- pcors
}

q_d <- c()
q_g <- c()
pq_d <- c()
pq_g <- c()
for (i in 1:h){
  q_d <- append(q_d, quantile(cor_matrix[,i], alpha/2))
  q_g <- append(q_g, quantile(cor_matrix[,i], 1 - alpha/2))
  pq_d <- append(pq_d, quantile(pcor_matrix[,i], alpha/2))
  pq_g <- append(pq_g, quantile(pcor_matrix[,i], 1 - alpha/2))
}

cors <- acf(temp_diff, plot=FALSE, lag.max=h)$acf
cors <- cors[2:(length(cors))]
ggplot() + 
  geom_point(aes(1:h, q_d, col="Przedział ufności"), size=3, alpha=0.5) +
  geom_point(aes(1:h, q_g, col="Przedział ufności"), size=3, alpha=0.5) +
  geom_point(aes(1:h, cors, col="Dane"), size=3, alpha=0.5)

pcors <- pacf(temp_diff, plot=FALSE, lag.max=h)$acf
ggplot() + 
  geom_point(aes(1:h, pq_d, col="Przedział ufności"), size=3, alpha=0.5) +
  geom_point(aes(1:h, pq_g, col="Przedział ufności"), size=3, alpha=0.5) +
  geom_point(aes(1:h, pcors, col="Dane"), size=3, alpha=0.5)

plot(residuals)
length(day[2:length(day)])
length(residuals)
plt <- ggplot() + geom_line(aes(x=day[2:length(day)], y=residuals))  + xlab("Dzień") + ylab("Residua")
ggsave("images/res_data.pdf", plt, width=14, height=6)

vars <- pvar(residuals, 200, TRUE)
vars <- vars + xlab("Dzień") + ylab("Wariancja cząstkowa")
ggsave("images/var_res_data.pdf", vars, width=14, height=6)

plt1 <- cor_plot(residuals)
plt2<- cor_plot(residuals, "pacf")
plt <- plot_grid(plt1, plt2, align="h")
ggsave("images/acf_res.pdf", plt, width=14, height=6)


xs <- seq(-15, 15, 0.01)
plt <- ggplot() + 
  geom_histogram(aes(residuals, after_stat(density)), col="black") +
  geom_line(aes(xs, dnorm(xs, 0, sd(residuals)), col="Gęstość teoretyczna"), linewidth=1) + 
  scale_color_manual(values=c("red")) +
  scale_fill_manual(values=c("grey")) + xlab("Residua") + ylab("Gęstość")
ggsave("images/des_res.pdf", plt, width=14, height=6)

xs <- seq(-15, 15, 0.01)
F <- ecdf(residuals)
plt <- ggplot() + 
  geom_line(aes(xs, F(xs), col="Dystrybuanta empir."), linewidth=1) +
  geom_line(aes(xs, pnorm(xs, 0, sd(residuals)), col="Dystrybuanta teoret."), linewidth=1) + xlab("") + ylab("Dystrybuanta")
ggsave("images/dys_res.pdf", plt, width=14, height=6)

