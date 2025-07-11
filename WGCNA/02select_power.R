############################### 2.select best β ###################################
# rm(list = ls())  
# load("step1_input.Rdata")

R.sq_cutoff = 0.8
if(T){
  # set range
  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2))
  # Call the network topology analysis function
  sft <- pickSoftThreshold(datExpr, 
                           powerVector = powers, 
                           verbose = 5)
  # SFT.R.sq > 0.8 , slope ≈ -1
  pdf("step2_power-value.pdf",width = 16,height = 12)
  # Plot the results
  par(mfrow = c(1,2));
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
       type = "n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels = powers,cex = 0.9,col = "red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h = R.sq_cutoff, col = "red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)",ylab = "Mean Connectivity", type = "n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex=0.9,col = "red")
  abline(h = 100,col = "red")
  dev.off()
}

sft$powerEstimate  # check best power
power = sft$powerEstimate

save(sft, power_, file = "step2_power_value.Rdata")
