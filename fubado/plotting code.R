ggplot(fake1, aes(x = samples, y = means)) +
  geom_hline(yintercept = popmean, colour = "red", linetype = "dashed") +
  geom_point(shape = 17, size = 3) +
  geom_errorbar(ymax = u_bound, ymin = l_bound, width = .5) +
  theme_classic(14) +
  xlab("Sample No.") +
  ylab("Sample Mean") +
  ggtitle("Sample Means with 95% CIs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(26, 34)
