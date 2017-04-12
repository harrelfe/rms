rms
=====

Regression Modeling Strategies

Current Goals
=============
* A non-downward compatible change will occur in the next release of the package
* The survfit.formula function (seen by the user as just survfit) for obtaining nonparametric survival estimates will be replaced by the npsurv function
* The purpose is to avoid conflicts with the survival package
* survfit.coxph has a new id option that generalizes individual=TRUE; need to change survfit.cph and survest.cph to use that

Web Sites
=============
* Overall: http://biostat.mc.vanderbilt.edu/Rrms
* Book: http://biostat.mc.vanderbilt.edu/rms
* CRAN: http://cran.r-project.org/web/packages/rms
* Changelog: https://github.com/harrelfe/rms/commits/master

To Do
=====
* Fix survplot so that explicitly named adjust-to values are still in subtitles.  See tests/cph2.s.
* Fix fit.mult.impute to average sigma^2 and then take square root, instead of averaging sigma
* Implement user-added distributions in psm - see https://github.com/harrelfe/rms/issues/41
