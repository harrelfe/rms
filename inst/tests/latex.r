require(rms)
set.seed(1)
x <- runif(100)
y <- abs(x - 0.5) + runif(100)
f <- ols(y ~ rcs(x, 5))
latex(f, file='')



require(rms)
x1 <- runif(200); x2 <- runif(200)
y <- sample(0:1, 200, TRUE)
f <- lrm(y ~ rcs(x1) + rcs(x2))
cat('\\documentclass{article}\\begin{document}\\usepackage{longtable}\n', file='/tmp/e.tex')
lat <- latex(f, file='/tmp/e.tex', append=TRUE)
sink('/tmp/e.tex', append=TRUE)
print(f, latex=TRUE)
sink()
cat('\\end{document}\n', file='/tmp/e.tex', append=TRUE)

