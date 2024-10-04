library(tidyverse)

# 95% confidence intervals
popmean <- 30
popsd <- 6
popsize <- 100000

set.seed(31947)
popscores <- rnorm(n = popsize, mean = popmean, sd = popsd)

# take samples
nsamp <- 20
smplsize <- 100

samplist <- list()

for(i in 1:nsamp) {
  samplist[[i]] <- sample(x = popscores, size = smplsize)
}

# make a histogram
sample3_df <- data.frame(weights = samplist[[3]])

ggplot(sample3_df, aes(x = weights)) +
  geom_histogram(fill = "brown", colour = "black", bins = 20) +
  labs(y = "Count") +
  theme_classic(base_size = 14) +
  ggtitle("Distribution of Weights from Sample 3") +
  theme(plot.title = element_text(hjust = .5))

# calculate descriptives
means <- sapply(samplist, mean)
SD <- sapply(samplist, sd)
SE <- SD / sqrt(smplsize)

# 95% CI
t_value <- qt(p = .975, df = smplsize - 1)
MoE <- t_value * SE
l_bound <- means - MoE
u_bound <- means + MoE
samples <- (seq(1, nsamp))

# put into a dataframe visualize
fake1 <- data.frame(samples, means, SE, l_bound, u_bound)
head(fake1)

ggplot(fake1, aes(x = samples, y = means)) +
  geom_hline(yintercept = popmean, colour = "red", linetype = "dashed") +
  geom_point(shape = 17, size = 3) +
  geom_errorbar(ymax = u_bound, ymin = l_bound, width = .5) +
  theme_classic(base_size = 14) +
  xlab("Sample No.") +
  ylab("Sample Mean") +
  ggtitle("Sample Means and 95% CIs") +
  theme(plot.title = element_text(hjust = .5)) +
  ylim(26, 34)

# new samples with samller size
set.seed(53191)
nsamp <- 20
smplsize <- 40

samplist <- list()

for(i in 1:nsamp) {
  samplist[[i]] <- sample(x = popscores, size = smplsize)
}

# descriptives
means <- sapply(X = samplist, FUN = mean)
SD <- sapply(X = samplist, FUN = sd)
SE <- SD / sqrt(smplsize)

# calculate confidence intervals
t_value <- qt(p = .975, df = smplsize - 1)
MoE <- t_value * SE
l_bound <- means - MoE
u_bound <- means + MoE
samples <- (seq(1, nsamp))

# make another dataframe and plot
fake2 <- data.frame(samples, means, SE, l_bound, u_bound)
head(fake2)

ggplot(fake2, aes(x = samples, y = means)) +
  geom_hline(yintercept = popmean, colour = "red", linetype = "dashed") +
  geom_point(shape = 17, size = 3) +
  geom_errorbar(ymax = u_bound, ymin = l_bound, width = .5) +
  theme_classic(base_size = 14) +
  xlab("Sample No.") +
  ylab("Sample Mean") +
  ggtitle("Sample Means and 95% CIs") +
  theme(plot.title = element_text(hjust = .5)) +
  ylim(26, 34)

# using CIs to test hypotheses about the population mean
Xmean <- 35
Yweights1 <- read.table(file = "breederY1.txt", header = T)

# histogram of breeder Y
ggplot(Yweights1, aes(x = weights)) +
  geom_histogram(fill = "cyan", colour = "black", bins = 12) +
  labs(y = "Count") +
  theme_classic(base_size = 14) +
  ggtitle("Distribution of Weights from Breeder Y") +
  theme(plot.title = element_text(hjust = .5))

# null: mean X = mean Y
# testing the null
Ymean <- mean(Yweights1$weights)
Ysd <- sd(Yweights1$weights)
Yse <- Ysd / sqrt(20)

print(round(Ymean, digits = 2))
print(round(Ysd, digits = 2))
print(round(Yse, digits = 2))

# 95% CI
t_value <- qt(p = .975, df = 19)
MoE <- t_value * Yse
lowerCI <- Ymean - MoE
upperCI <- Ymean + MoE
Y_CI <- c(lowerCI, upperCI)

print(Y_CI)


