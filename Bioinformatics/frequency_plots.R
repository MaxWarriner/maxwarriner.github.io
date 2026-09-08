library(tidyverse)
dat <- tibble(region = rep(c("Entire mRNA", "5'UTR", "CDs", "3'UTR"), each = 4), 
              nucleotide = rep(c("A", "C", "G", "T"), times = 4), 
              freq = c(0.24, 0.23, 0.28, 0.25, 
                       0.22, 0.19, 0.36, 0.23, 
                       0.23, 0.25, 0.31, 0.22,
                       0.27, 0.22, 0.20, 0.30))

dat$region <- factor(dat$region, levels = c("Entire mRNA", "5'UTR", "CDs", "3'UTR"))

plot <- ggplot(data = dat, aes(x = region, y = freq, fill = nucleotide)) + 
  geom_col() + 
  theme_minimal() + 
  geom_hline(yintercept = 0) + 
  xlab("Region") + 
  ylab("Frequency")

ggsave(plot, filename = "nucleotide_frequencies.png", height = 3, width = 5)


dimer_dat <- read_csv('dimer_frequency.csv')
plot <- ggplot(data = dimer_dat, aes(x = Dimer, y = Frequency)) + 
  geom_col() + 
  geom_hline(yintercept = 0) + 
  theme_minimal() + 
  ggtitle("Entire mRNA Dimer Frequencies")


ggsave(plot, filename = "dimer_frequencies.png", height = 3, width = 5)


codon_dat <- tibble(Codon = c("AGA", "AGG", 'CGA', 'CGC', 'CGG', 'CGT'), 
                    Freq = c(0.044, 0.289, 0.2, 0.111, 0.311, 0.044))
plot <- ggplot(data = codon_dat, aes(x = Codon, y = Freq)) + 
  geom_col() + 
  geom_hline(yintercept = 0) + 
  theme_minimal() + 
  ylab("Relative Frequency") + 
  xlab('Arginine Codon')


ggsave(plot, filename = "arg_frequencies.png", height = 3, width = 5)


aa_dat <- read_csv('aa_frequencies.csv')

plot <- ggplot(data = aa_dat, aes(x = `Amino Acid`, y = Frequency)) + 
  geom_col() + 
  geom_hline(yintercept = 0) + 
  theme_minimal() + 
  ylab("Relative Frequency") + 
  xlab('Amino Acid')


ggsave(plot, filename = "aa_frequencies.png", height = 3, width = 5)

