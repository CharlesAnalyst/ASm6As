# https://cran.r-project.org/web/packages/eulerr/vignettes/introduction.html
# library(eulerr)

# library(VennDiagram)
# library(ggvenn)

library(venneuler)
vd <- venneuler(c(A=85, B=44, "A&B"=44))
plot(vd)


library(eulerr)

vd <- euler(c(A = 0.3, B = 0.1,"A&B" = 0.1))
plot(vd)
