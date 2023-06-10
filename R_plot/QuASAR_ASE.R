library(QuASAR)
ase.dat <- UnionExtractFields(c("stomachadult2.input.quasar.in.gz", "brainadult3.input.quasar.in.gz", "brainadult1.input.quasar.in.gz", "brainadult2.input.quasar.in.gz", "brainadult4.input.quasar.in.gz"), combine=TRUE)
ase.dat.gt <- PrepForGenotyping(ase.dat, min.coverage=5)
sample.names <- colnames(ase.dat.gt$ref)
sample.names
ase.joint <- fitAseNullMulti(ase.dat.gt$ref, ase.dat.gt$alt, log.gmat=log(ase.dat.gt$gmat))
ourInferenceData <- aseInference(gts=ase.joint$gt, eps.vect=ase.joint$eps, priors=ase.dat.gt$gmat, ref.mat=ase.dat.gt$ref, alt.mat=ase.dat.gt$alt, min.cov=5, sample.names=sample.names, annos=ase.dat.gt$annotations)
##
write.table(ourInferenceData[[1]]$dat, file='inference_results/stomachadult2.txt', row.names=FALSE, quote=FALSE,sep="\t")
