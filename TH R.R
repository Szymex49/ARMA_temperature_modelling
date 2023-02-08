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
    plt <- ggplot() + geom_line(aes(1:(length(X) - k), vars))
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

length(temp)



ggplot() + geom_line(aes(x=day, y=temp))




plt1 <- cor_plot(temp, h=80)
plt2 <- cor_plot(temp, "pacf", 80)
plt <- plot_grid(plt1, plt2, align="h")
plt
ggsave("images/raw_acf.pdf", plt, width=14, height=6)





