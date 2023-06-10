library(seqLogo)
m <- read.table("/home/galaxy/project/alleleSpecific_analysis/plots/e_fig1d.txt")
p <- makePWM(m, alphabet="DNA") #
seqLogo(p@pwm, ic.scale=TRUE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15)

mFile <- system.file("extdata/pwm1", package="seqLogo")
m <- read.table(mFile)
m
p <- makePWM(m)
seqLogo(p)
seqLogo(p, ic.scale=FALSE)
seqLogo(p, fill=c(A="#4daf4a", C="#377eb8", G="#ffd92f", T="#e41a1c"),
        ic.scale=FALSE)

r <- makePWM(m, alphabet="RNA")
seqLogo(r, ic.scale=FALSE)