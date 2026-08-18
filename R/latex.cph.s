## -----------------------------------------------------------------------
## latex.cph, with the survival table's non-html (LaTeX/Typst) path
## rebuilt on tinytable, replacing the old Hmisc::latex()/tempfile-based
## code. Everything else (the caption fix already applied, the html/md
## kable-based htmlTab() path, the equation header construction) is
## unchanged.
##
## New helper: surv_tiny_table (underscore-named, per convention -- not
## an S3 method). Builds the survival table for either output="latex" or
## output="typst", handling alignment (right-justified data, centered
## header, confirmed via i=0 targeting the header row), centering the
## whole table (typst: theme_typst(figure=TRUE, align_figure="c");
## latex: theme_latex(placement="H"), which was confirmed by direct
## compile test to center automatically via an inserted \centering, no
## separate mechanism needed), and extracts the rendered string via
## save_tt(x, output) -- the documented public function for this,
## confirmed reliable, unlike as.character() (a real internal bug) and
## build_tt() (works but unexported).
##
## Label construction is genuinely prType()-dependent here (LaTeX's
## S_{0}(t) vs Typst's S_(0)(t) subscript grouping) -- unlike the
## equation content elsewhere, tinytable's cell/header content does not
## go through texmath, confirmed directly by the braces-showing-literally
## bug a few turns back. For the stratified case, stratum-derived label
## text (n) is run through typstTranslate() for the typst branch, for
## the same defensive reasons typstTranslate() gets used on arbitrary
## text elsewhere -- the existing latex branch does not currently run
## the equivalent latexTranslate() on n, a pre-existing gap left
## untouched here to keep this diff focused.
## -----------------------------------------------------------------------

surv_tiny_table <- function(s, output) {
  ## Local, non-exported helper -- builds just the raw-Typst fence
  ## string (matching typstAsis()'s own fence, four backticks for the
  ## same nesting-safety reason) WITHOUT the asis_output()/console-print
  ## behavior typstAsis() adds -- this result is one intermediate piece
  ## of latex.cph's larger `w` vector, which itself gets asis-emitted
  ## once, later, via rendHTML(). Defined locally rather than shared
  ## from Hmisc to avoid exporting a function this small.
  typstFence <- function(x) paste0("````{=typst}\n", x, "\n````\n")

  d  <- as.data.frame(unclass(s))
  d  <- cbind(t = as.numeric(rownames(s)), d)
  nc <- ncol(d)

  x <- tt(d, output = output)
  x <- style_tt(x, j = 1 : nc, align = 'r')
  x <- style_tt(x, i = 0, align = 'c')

  if(output == 'latex')
    ## Confirmed by direct compile test: theme_latex(placement=...)
    ## already inserts \centering automatically into the table float --
    ## no separate centering mechanism needed for LaTeX.
    x <- theme_latex(x, placement = 'H')

  ## save_tt(x, output): documented public function -- returns the
  ## rendered string directly when its second argument is a format name
  ## rather than a file path. Confirmed reliable, unlike as.character()
  ## (which appears to have a real internal bug) and build_tt() (which
  ## works but is unexported).
  result <- save_tt(x, output)

  ## typst output must be fence-marked (typstFence) so Pandoc's markdown
  ## reader treats it as raw pass-through Typst rather than plain text
  ## needing its own escaping -- confirmed necessary, since Typst syntax
  ## (unlike LaTeX) has no "looks like a raw command, preserve it"
  ## heuristic Pandoc applies to untagged text. latex output stays bare:
  ## already confirmed working end-to-end without any fence, since
  ## Pandoc does pass bare raw-LaTeX-looking text through for a LaTeX
  ## target.
  ##
  ## Centering: theme_typst(figure=TRUE, align_figure="c") -- the
  ## previous approach here -- was confirmed NOT reliable (real cph/Gls/
  ## psm output showed an outer #block[ of unrelated origin surrounding
  ## the table with no definite width, leaving align(center) nothing to
  ## center within). Same fix now applied here as in
  ## coef_tiny_table/stats_tiny_table: an explicit
  ## #block(width:100%)[#align(center)[...]] wrap, robust regardless of
  ## what, if anything, surrounds it.
  if(output == 'typst')
    result <- typstFence(paste0('#block(width: 100%)[#align(center)[\n',
                                result, '\n]]'))
  result
}

latex.cph <-
  function(object, title,
           file='',
           append=FALSE, surv=TRUE, maxt=FALSE, which=NULL, varnames,
           columns=65, inline=FALSE,
           before=if(inline)"" else "&", after="",
           dec=3, pretrans=TRUE,
           caption=NULL, digits=.Options$digits, size='', ...)
{
  md <- prType() %in% c('html', 'md', 'markdown')
  filegiven <- file != ''
  f <- object
  whichThere <- length(which)

  atr <- f$Design

  lev <- names(f$freq)
  Intercept <- -f$center
  strata <- levels(f$strata)    ## was f$strata
  w <- if(length(caption)) {
         if(md) paste('<div align=center><strong>', caption,
                      '</strong></div>')
         else
           paste0("$$\\textbf{", caption, "}$$")
         }
  if(! length(which) & !inline) {
    if(length(strata)==0)
      w <- c(w, paste("$$P(T\\geq t~|~X) = S_{0}(t)^{\\mathrm{e}^{X\\beta}},~~ \\mathrm{where}$$",sep=""))
    else {
      sname <- atr$name[atr$assume.code==8]
      strata.sub <- letters[8 + (1 : length(sname))]
      s <- paste("\\mathrm{",sname,"}=",strata.sub,sep="")
      s <- paste(s, collapse=",")
      w <- c(w,paste("$$P(T\\geq t~|~X,",s,")=S_{",
                     paste(strata.sub,collapse=""),
                     "}(t)^{\\mathrm{e}^{X\\beta}},~~\\mathrm{where}$$", sep=""))
    }
  }
  if(!length(which)) which <- 1:length(atr$name)
  if(missing(varnames)) varnames <- atr$name[atr$assume.code!=9]
  ltx <-
    latexrms(f, file='', append=TRUE, which=which, varnames=varnames,
             columns=columns,
             before=before, after=after,
             prefix=if(!whichThere)"X\\hat{\\beta}" else NULL,
             intercept=Intercept, inline=inline,
             pretrans=pretrans, digits=digits, size=size)
  if(inline) return(ltx)
  w <- c(w, ltx)

  htmlTab <- function(s) {
    s <- cbind('$t$'= as.numeric(rownames(s)), s)
    for(j in 1 : ncol(s)) s[, j] <- round(s[, j], dec)
    if (requireNamespace("kableExtra", quietly=TRUE)) {
      as.character(
        knitr::kable(s, format='html',
                     align='r', row.names=FALSE) |>
        kableExtra::kable_styling(full_width=FALSE)   )
    } else {
      as.character(
        knitr::kable(s, format='html',
                     align='r', row.names=FALSE) )
    }
  }

  ss <- f$surv.summary
  if(surv && length(ss)) {
    fs <- levels(f$strata)
    nstrat <- 0; if(length(fs)) nstrat <- length(fs)
    times <- as.numeric(dimnames(ss)[[1]])
    maxtime <- f$maxtime
    if(max(times) >= maxtime) maxt <- FALSE
    output <- if(prType() == 'typst') 'typst' else 'latex'   ## NEW

    if(nstrat == 0) {
      s <- matrix(ss[, , 1], ncol=1)
      if(maxt) {
        s <- cbind(s, f$surv[L <- length(f$surv)])
        times <- c(times, f$time[L])
      }
      lab0 <- if(prType() == 'typst') "$S_(0)(t)$" else "$S_{0}(t)$"   ## NEW
      dimnames(s) <- list(t=format(times), lab0)
      if(md) {
        w <- c(w, htmlTab(s))
      } else {
        w <- c(w, surv_tiny_table(s, output))   ## NEW: was Hmisc::latex()/tempfile
      }
    } else {

      ## Change . to ,blank
      n <- sedit(paste(fs,',',sep=''), '.', ', ')
      ## Change sname=*, to *,
      n <- sedit(n, paste(sname,'=*,',sep=''), rep('*, ', length(sname)))
      n <- substring(n, 1, nchar(n) - sum(atr$assume.code == 8) - 1)
      s <- ss[, , 1]
      if(maxt) {
        smax <- rep(NA, nstrat)
        for(i in 1 : nstrat)
          smax[i] <- f$surv[[i]][abs(f$time[[i]]-maxtime) < 0.001]
        s <- rbind(s, smax)
        times <- c(times, maxtime)
      }

      labn <- if(prType() == 'typst')                            ## NEW
                ## Quoted as a Typst string literal -- bare multi-letter
                ## text inside Typst math mode (e.g. S_(female)) is read
                ## as an undefined variable identifier, not literal
                ## text; confirmed directly by compile error ("unknown
                ## variable: female"), with Typst's own error message
                ## pointing at quoting as the fix. Not needed for the
                ## no-strata case (lab0, above) since that embeds a bare
                ## numeric literal ("0"), which Typst math handles
                ## correctly without quoting.
                paste0("$S_(\"", typstTranslate(n), "\")(t)$")
              else
                paste0("$S_{", n, "}(t)$")
      dimnames(s) <- list(t=format(times), labn)
      if(md) {
        w <- c(w, htmlTab(s))
      }
      else {
        w <- c(w, surv_tiny_table(s, output))   ## NEW: was Hmisc::latex()/tempfile
      }
    }
  }
  if(filegiven || prType() == 'plain') {
    cat('\n', w, sep='\n', file=file, append=append)
    return(invisible())
  }
  rendHTML(w, html=FALSE)
}
