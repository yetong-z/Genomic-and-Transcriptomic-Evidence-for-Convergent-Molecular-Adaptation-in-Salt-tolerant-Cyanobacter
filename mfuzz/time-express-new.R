BiocManager::install('Mfuzz')
library('Mfuzz')


# read data
dat_yt <- read.csv('20250621/all_DEGs_tpm_counts_matrix.csv', row.names = 1)
dat_yt <- as.matrix(dat_yt)

mfuzz_class <- new('ExpressionSet',exprs = dat_yt)

# solve missing data
mfuzz_class1 <- filter.NA(mfuzz_class,thres = 0.25)
mfuzz_class2 <- fill.NA(mfuzz_class1,mode = 'mean')
mfuzz_class3 <- filter.std(mfuzz_class2,min.std = 0)

# normalize
mfuzz_class4 <- standardise(mfuzz_class3)


# FCM cluster
n_yt <- 6

# Evaluate the optimal m value (fuzzification parameter)
m_yt <- mestimate(mfuzz_class4)

# Set a random number seed
set.seed(123)

# cluster
cl_yt <- mfuzz(mfuzz_class4, c = n_yt, m = m_yt)

mfuzz.plot2(mfuzz_class4, cl = cl_yt, mfrow = c(2,3), time.labels = c(0,1,2,4))

# merge data
protein_cluster <- cl_yt$cluster
protein_cluster1 <- cbind(dat_yt[names(protein_cluster),], protein_cluster)
write.table(protein_cluster1,'clusters_for_normalized_DEGs_20250621.csv', sep = ',', col.names = NA, quote = F)







