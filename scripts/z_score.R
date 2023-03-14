### DATA PREPARATION

## INDUCED
# Load data
i <- read.table("distance_AsiSI_perm_hg19_induced", col.names = c("id", "distance"))
# Subset real pattern
smp_i <- i[i$id == "GCGATCGC",2]
# Subset permutations
pop_i <- i[,2]
# Absolute values real pattern
smp_i_abs <- abs(smp_i)
# Absolute values permutations
pop_i_abs <- abs(pop_i)
# Logarithmic values real pattern
smp_i_log <- log(smp_i_abs)
# Logarithmic values permutations
pop_i_log <- log(pop_i_abs)

## UNINDUCED
# Load data
u <- read.table("distance_AsiSI_perm_hg19_uninduced", col.names = c("id", "distance"))
# Subset real pattern
smp_u <- u[u$id == "GCGATCGC",2]
# Subset permutations
pop_u <- u[,2]
# Absolute values real pattern
smp_u_abs <- abs(smp_u)
# Absolute values permutations
pop_u_abs <- abs(pop_u)
# Logarithmic values real pattern
smp_u_log <- log(smp_u_abs)
# Logarithmic values permutations
pop_u_log <- log(pop_u_abs)

## REMOVE RAW DATA
rm(i,u)


### MEANS AND STANDARD DEVIATIONS

## INDUCED
# Mean distance real pattern
m_smp_i <- mean(smp_i)
# Standard deviation real pattern
sd_orig_i <- sd(smp_i)
# Mean absolute distance real pattern
m_smp_i_abs <- mean(smp_i_abs)
# Mean distance permutations
m_pop_i <- mean(pop_i)
# Standard deviation permutations
sd_pop_i <- sd(pop_i)
# Mean absolute distance permutations
m_pop_i_abs <- mean(pop_i_abs)

## UNINDUCED
# Mean distance real pattern
m_smp_u <- mean(smp_u)
# Standard deviation real pattern
sd_orig_u <- sd(smp_u)
# Mean absolute distance real pattern
m_smp_u_abs <- mean(smp_u_abs)
# Mean distance permutations
m_pop_u <- mean(pop_u)
# Standard deviation permutations
sd_pop_u <- sd(pop_u)
# Mean absolute distance permutations
m_pop_u_abs <- mean(pop_u_abs)


### STATISTICAL TESTS

## STUDENT'S T-TESTS
# Real pattern vs permutations, induced
t_sp_i <- t.test(smp_i, pop_i)
# Real pattern vs permutations, uninduced
t_sp_u <- t.test(smp_u, pop_u)
# Induced vs uninduced, real pattern
t_iu_s <- t.test(smp_i, smp_u)
# Induced vs uninduced, permutations
t_iu_p <- t.test(pop_i, pop_u)

## STUDENT'S T-TESTS (ABSOLUTE MEANS)
# Real pattern vs permutations, induced
t_sp_i_abs <- t.test(smp_i_abs, pop_i_abs)
# Real pattern vs permutations, uninduced
t_sp_u_abs <- t.test(smp_u_abs, pop_u_abs)
# Induced vs uninduced, real pattern
t_iu_s_abs <- t.test(smp_i_abs, smp_u_abs)
# Induced vs uninduced, permutations
t_iu_p_abs <- t.test(pop_i_abs, pop_u_abs)

## Z-SCORES
# Real pattern vs permutations, induced
z_sp_i <- (m_smp_i - m_pop_i) / sd_pop_i
# Real pattern vs permutations, uninduced
z_sp_u <- (m_smp_u - m_pop_u) / sd_pop_u

## Z-SCORES (ABSOLUTE MEANS)
# Real pattern vs permutations, induced
z_sp_i_abs <- (m_smp_i_abs - m_pop_i_abs) / sd_pop_i
# Real pattern vs permutations, uninduced
z_sp_u_abs <- (m_smp_u_abs - m_pop_u_abs) / sd_pop_u

## WILCOXON RANK SUM TEST
# Real pattern vs permutations, induced
w_sp_i <- wilcox.test(smp_i, pop_i)
# Real pattern vs permutations, uninduced
w_sp_u <- wilcox.test(smp_u, pop_u)
# Induced vs uninduced, real pattern
w_iu_s <- wilcox.test(smp_i, smp_u)
# Induced vs uninduced, permutations
w_iu_p <- wilcox.test(pop_i, pop_u)

## WILCOXON RANK SUM TEST (ABSOLUTE MEANS)
# Real pattern vs permutations, induced
w_sp_i_abs <- wilcox.test(smp_i_abs, pop_i_abs)
# Real pattern vs permutations, uninduced
w_sp_u_abs <- wilcox.test(smp_u_abs, pop_u_abs)
# Induced vs uninduced, real pattern
w_iu_s_abs <- wilcox.test(smp_i_abs, smp_u_abs)
# Induced vs uninduced, permutations
w_iu_p_abs <- wilcox.test(pop_i_abs, pop_u_abs)


### PLOTS

## HISTOGRAMS
pdf("smp_i.pdf", width = 24, height = 12)
hist(smp_i, breaks=100, freq=FALSE, main="Real pattern, induced", xlab="Distance")
dev.off()
# Permutations induced
pdf("pop_i.pdf", width = 24, height = 12)
hist(pop_i, breaks=100, freq=FALSE, main="permutations, induced", xlab="Distance")
dev.off()
# Real pattern uninduced
pdf("smp_u.pdf", width = 24, height = 12)
hist(smp_u, breaks=100, freq=FALSE, main="Real pattern, uninduced", xlab="Distance")
dev.off()
# Permutations uninduced
pdf("pop_u.pdf", width = 24, height = 12)
hist(pop_u, breaks=100, freq=FALSE, main="permutations, uninduced", xlab="Distance")
dev.off()

## HISTOGRAMS ABSOLUTE VALUES
pdf("smp_i_abs.pdf", width = 24, height = 12)
hist(smp_i_abs, breaks=100, freq=FALSE, main="Real pattern, induced", xlab="Distance")
dev.off()
# Permutations induced
pdf("pop_i_abs.pdf", width = 24, height = 12)
hist(pop_i_abs, breaks=100, freq=FALSE, main="permutations, induced", xlab="Distance")
dev.off()
# Real pattern uninduced
pdf("smp_u_abs.pdf", width = 24, height = 12)
hist(smp_u_abs, breaks=100, freq=FALSE, main="Real pattern, uninduced", xlab="Distance")
dev.off()
# Permutations uninduced
pdf("pop_u_abs.pdf", width = 24, height = 12)
hist(pop_u_abs, breaks=100, freq=FALSE, main="permutations, uninduced", xlab="Distance")
dev.off()

## HISTOGRAMS LOG SPACE
# Real pattern induced
pdf("smp_i_ln.pdf", width = 24, height = 12)
hist(smp_i_log, breaks=100, freq=FALSE, main="Real pattern, induced", xlab="Distance")
dev.off()
# Permutations induced
pdf("pop_i_ln.pdf", width = 24, height = 12)
hist(pop_i_log, breaks=100, freq=FALSE, main="permutations, induced", xlab="Distance")
dev.off()
# Real pattern uninduced
pdf("smp_u_ln.pdf", width = 24, height = 12)
hist(smp_u_log, breaks=100, freq=FALSE, main="Real pattern, uninduced", xlab="Distance")
dev.off()
# Permutations uninduced
pdf("pop_u_ln.pdf", width = 24, height = 12)
hist(pop_u_log, breaks=100, freq=FALSE, main="permutations, uninduced", xlab="Distance")
dev.off()

## HISTOGRAMS SHORT DISTANCES (1x10^6)
# Subset data
smp_i_1000000 <- subset(smp_i, smp_i > -500000 & smp_i < 500000)
pop_i_1000000 <- subset(pop_i, pop_i > -500000 & pop_i < 500000)
smp_u_1000000 <- subset(smp_u, smp_u > -500000 & smp_u < 500000)
pop_u_1000000 <- subset(pop_u, pop_u > -500000 & pop_u < 500000)
# Create bins
bins <- seq(-500000, 500000, by=10000)
# Real pattern induced
pdf("smp_i_1000000.pdf", width = 24, height = 12)
hist(smp_i_1000000, breaks=bins, freq=FALSE, main="Real pattern, induced", xlab="Distance")
dev.off()
# Permutations induced
pdf("pop_i_1000000.pdf", width = 24, height = 12)
hist(pop_i_1000000, breaks=bins, freq=FALSE, main="permutations, induced", xlab="Distance")
dev.off()
# Real pattern uninduced
pdf("smp_u_1000000.pdf", width = 24, height = 12)
hist(smp_u_1000000, breaks=bins, freq=FALSE, main="Real pattern, uninduced", xlab="Distance")
dev.off()
# Permutations uninduced
pdf("pop_u_1000000.pdf", width = 24, height = 12)
hist(pop_u_1000000, breaks=bins, freq=FALSE, main="permutations, uninduced", xlab="Distance")
dev.off()


### SAVE

## Image
save.image(file = "z_scores.Rdata")
