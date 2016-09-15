## From John.Stickley@bch.nhs.uk

require(rms)

test_data <-
structure(
    list(group = structure( c(1L, 2L, 2L, 1L, 2L, 1L, 2L,
        1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L,
        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L,
        2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L),
    .Label = c("a", "b"), class = "factor"),
    survival = c(17.46772065, 98.92209528,78.29864249, 9.822413669, 79.55050186, 82.36595474, 1.42766417,
                71.48805748, 61.33571345, 84.62631825, 93.03022837, 44.04354499,
                81.06711649, 26.19891261, 68.64477557, 52.2160246, 17.780942,
                4.515968877, 95.46066172, 73.63010059, 40.13833451, 20.39467002,
                50.80529216, 70.23087236, 23.89309088, 53.86527662, 3.422234859,
                35.30675488, 50.07307746, 4.68602929, 86.04636345, 72.98976535,
                33.18048902, 37.94566436, 83.17678398, 16.95356411, 80.5844794,
                8.599290846, 46.06581857, 1.644574571, 34.81582745, 49.96017595,
                11.74200883, 60.07697075, 80.40946019, 55.00705828, 17.75483404,
                98.69523629, 68.15668013, 4.959304343),
    outcome = c(0L, 0L, 0L,0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L,
                    1L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L,
                    1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L)),
    .Names = c("group", "survival", "outcome"), class = "data.frame",
    row.names = c(NA, -50L))

table(test_data$group)
s <- npsurv(Surv(survival, outcome) ~ group, data = test_data)
s$n.risk[c(1, s$strata['group=a'] + 1)]

survplot(s, conf='none', lwd=2, lty=c(1,1,1), n.risk=TRUE, time.inc=5,
         label.curves=FALSE, cex.lab=1.75)

