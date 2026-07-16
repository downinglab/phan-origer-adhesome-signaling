# Figure 5i Mean-CV Plots
# By Andrew Phan

#calculate CV and mean of normalized expression for each gene in the LacZ condition
lacZ_norm = grepl("_LacZ",colnames(expr))
lacZ_norm = expr[,lacZ_norm]
lacZ_cvmean = calcCV(lacZ_norm)
lacZ_cvmean$gene = rownames(expr)
lacZ_cvmean$cond = "LacZ"
lacZ_cvmean$gene_sig = expr$gene_sig

#calculate CV and mean of normalized expression for each gene in the SHROOM3 condition
S3_norm = grepl("_Shroom3",colnames(expr))
S3_norm = expr[,S3_norm]
S3_cvmean = calcCV(S3_norm)
S3_cvmean$gene = rownames(expr)
S3_cvmean$cond = "SHROOM3"
S3_cvmean$gene_sig = expr$gene_sig

#merge the two conditions to one dataframe
comb_cvmean = rbind(lacZ_cvmean, S3_cvmean)

#create plots of CV vs mean
cvL = ggplot(lacZ_cvmean, aes(x=mean, y=cv)) + geom_point() + theme_classic()
cvS = ggplot(S3_cvmean, aes(x=mean, y=cv)) + geom_point() + theme_classic()
cv_comb = ggplot(comb_cvmean[sample(nrow(comb_cvmean)),], aes(x=mean, y=cv, col=cond)) + geom_point(alpha=0.75) + theme_classic() + ylim(-3.25,3.6) + xlim(-4.5,1) #shuffles rows to randomize order of points
#cv_comb

#look at CV vs mean for only subsets of genes
#sub = pcp
sub = "pluri"
#sub = FN
#sub = EMT

cv_sub = grepl(paste(sub, collapse = "|"), comb_cvmean$gene)
cv_sub = comb_cvmean[cv_sub,]
cv_pl = ggplot(cv_sub[sample(nrow(cv_sub)),], aes(x=mean, y=cv, col=cond)) + geom_point(alpha=0.75) + theme_classic()
#cv_pl

#plot the difference of mean and CV of the same genes
dif_cv <- merge(lacZ_cvmean, S3_cvmean, by=0, all=TRUE) #if you want to see comparisons of conditions for each gene
#te_pluri = te[grepl(paste(pluri, collapse = "|"), te$gene.x),] #only look at genes from a specified signature
dif_cv = dif_cv[complete.cases(dif_cv), ] #keep only rows without NA

#dif_cv$dMean = (dif_cv$mean.y - dif_cv$mean.x)/dif_cv$mean.x
dif_cv$dMean = (dif_cv$mean.y - dif_cv$mean.x)#/dif_cv$mean.x
dif_cv$dCV = dif_cv$cv.y - dif_cv$cv.x
dif_cv$dFano = dif_cv$fano.y - dif_cv$fano.x


dif_cv = dif_cv[order(dif_cv$gene_sig.x, decreasing = F),]

#plot all genes, colored by gene signature
lmFit_pluri = lm(dif_cv[dif_cv$gene_sig.x=="pluri",15] ~ dif_cv[dif_cv$gene_sig.x=="pluri",14]) #find linear regression of pluripotency genes
lmFit = lm(dif_cv[,15] ~ dif_cv[,14]) #find linear regression of all genes

ggplot(dif_cv, aes(x=dMean, y=dCV, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#D3D3D3", "#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_abline(intercept = coefficients(lmFit)[1], slope = coefficients(lmFit)[2]) +
  geom_abline(intercept = coefficients(lmFit_pluri)[1], slope = coefficients(lmFit_pluri)[2], col = "red") +
  ylim(-1.5, 1.5) + xlim(-1.6, 1.6) +
  theme(legend.position='none')

lmFit
lmFit_pluri

#plot only genes in gene signatures
ggplot(dif_cv[dif_cv$gene_sig.x!=0,], aes(x=dMean, y=dCV, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")
ggplot(dif_cv[dif_cv$gene_sig.x=="pluri",], aes(x=dMean, y=dCV, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")

#separate genes by change in CV and mean (increase or decrease in each) in S3 vs LacZ
cv_gene_clusters = list(mposcpos = subset(dif_cv, dMean > 0 & dCV > 0), mnegcpos = subset(dif_cv, dMean < 0 & dCV > 0), mposcneg = subset(dif_cv, dMean > 0 & dCV < 0), mnegcneg = subset(dif_cv, dMean < 0 & dCV < 0))
