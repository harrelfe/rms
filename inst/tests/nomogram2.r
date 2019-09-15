require(rms)
create <- FALSE
if(create) {
d <- csv.get("hongLiData.csv", lowernames=TRUE, charfactor=TRUE)
d <-
  upData(d,
         labels=c(surg.discharge.days="Discharge Days",
                  primarysite="Primary Site",
                  facilitytype="Facility Type",
                  insurance1="Insurance",
                  race1="Race",
                  region="Region",
                  samefacility="Same Facility",
                  port.nccn.fail="outcome"))
d <- d[Cs(port.nccn.fail, surg.discharge.days, primarysite, facilitytype,
          insurance1, race1, samefacility, region)]
save(d, file='nomogram2.rda', compress=TRUE)
} else load('nomogram2.rda')

ddist  <-  datadist(d);  options(datadist='ddist')

f <- lrm(port.nccn.fail ~ surg.discharge.days + primarysite +
           facilitytype + insurance1 + race1 + samefacility + region,
         data=d)

for(abbrev in c(FALSE, TRUE)) {

  n <- nomogram(f, lp.at=seq(-2, 5, by=0.5), fun=plogis,
                fun.at=c(seq(.1, .9, by=.1), .95, .99, .999),
                funlabel="Risk of Delayed PORT Initiation",
                abbrev=abbrev, minlength=1, lp=FALSE)
  if(! abbrev) n1 <- n else n2 <- n
  }

plot(n1)
plot(n2)
attr(n2, 'info')$Abbrev

# Hong Li <liho@musc.edu>
# The variable samefacility has two categories and region has 4 categories. But in the nomogram, the variable samefacility and region are switched, i.e. samefacility has 4 categories and region has 2 categories. All other variables are correct. 
