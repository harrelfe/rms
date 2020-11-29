if(FALSE) {
require(rms)
set.seed(1)
n <- 5000
x1 <- runif(n); x2 <- runif(n); x3 <- runif(n); x4 <- runif(n); x5 <- runif(n)
y <- round(400 * runif(n))
fm <- y ~ x1 + x2 + x3 + x4 + x5
print(system.time(f <- lrm(fm)))
print(system.time(f <- orm(fm)))
ti <- numeric(0)
rs <- c(5,10,20,40,80,120,160,200,250,seq(300, 1000, by=100),1500,2000,2500,3000)
for(r in rs) {
  cat(r, '\n')
  y <- round(r * runif(n))
  ti <- c(ti, system.time(orm(fm))['elapsed'])
}
plot(rs, ti)   # linear in no. of intercepts!

y <- round(1000 * runif(n))
print(system.time(f <- orm(fm, x=TRUE, y=TRUE)))
print(system.time(validate(f, B=10)))   # 15x longer vs. 10x
Rprof()
# for(i in 1 : 10) f <- orm(fm)
print(validate(f, B=10))
Rprof(NULL)
# s <- summaryRprof()


if(require(proftools)) {
  tmp.dot <- tempfile()
  tmp.pdf <- tempfile()
  pd <- readProfileData()
  profileCallGraph2Dot(pd, filename = tmp.dot)
  system(sprintf("dot -Tpdf -o %s %s", tmp.pdf, tmp.dot))
  browseURL(sprintf("file://%s", tmp.pdf))

  unlink(tmp.dot)
  unlink(tmp.pdf)
}

}
