#' Bry-Boshan or Harding and Pagan
#'
#' Compute Bry-Boschan or Harding and Pagan turning points dating rules.
#'
#' Purpose: compute BB or HP turning points dating rules
#'
#' @section References:
#'     D. Harding & A. Pagan (2002), "Dissecting the Cycle: a Methodological
#'          Investigation", Journal of Monetary Economics, n 49, pp. 365-381
#'     G. Bry & C. Boschan (1971), "Cyclical Analysis of Time Series:
#'    		 Selected Procedures and Computer Programs", NBER Technical Paper 20
#'
#' @param x A time series (seasonally adjusted).
#' @param ilog 1 for logs of data (multiplicative model)
#' @param ioutput 1 for full output, 0 otherwise
#' @param imcd Short-term MA for step IV:
#'             0 = Let program decide based on MCD
#'             3 = MA(3), 4 = MA(4), 5 = MA(5), 6 = MA(6)
#' @param nma number of terms in MA filter for initial TP
#'            (BB use 12-term asymmetric)
#' @return Dates for 'peaks' and 'troughs'  in the yearmon format.
#'
brybos <- function(x, ilog, imcd, nma) {
  nd <- length(x)

  if (ilog == 1) {
    if (min(x) > 0) {
      x <- log(x)
    } else {
      warning("Log option chosen with negative values -- levels used")
    }
  }

  # Execute Step I of Bry-Boschan procedure
  # Check for outliers and replace with fitted values
  xspa <- spencer(x)
  xo <- bbout(x, xspa, ilog)

  # Calculate spencer curve using outlier adjusted series
  xspb <- spencer(xo)

  # Calculate 12 month moving average
  # The original work uses a non-symmetric filter when order is even
  # here the procedure always result centered
  xma <- MA(xo, nma)

  # Calculate MCD smoothed series
  # Bry-Boschan pages 24-25
  mcdnum <- mcd(x, 12, ilog)

  if (imcd == 0) {
    mcdfilt <- mcdnum
    if (mcdnum < 3)
      mcdfilt <- 3
    if (mcdnum > 6)
      mcdfilt <- 6
    mcdx <- MA(x, mcdfilt)
    # Set initial and terminal values as in Bry-Boschen
    n2 <- mcdfilt / 2
    n2t <- trunc(n2)
    # if ((n2 - n2t) > 0.01) {
    #   nm <- n2t
    #   np <- n2t
    # } else {
    #   nm <- n2t
    #   np <- n2t - 1
    # }
    nm <-
      np <- n2t # no need for the above code as MA is always centered
    # First and last observations (As in BB -- see program TP2)
    mcdx[1:nm] <- mcdx[nm + 1]
    mcdx[(length(mcdx) - np + 1):length(mcdx)] <-
      mcdx[length(mcdx) - np]
  } else {
    mcdfilt <- imcd
    mcdx <- MA(x, mcdfilt)
    n2 <- mcdfilt / 2
    n2t <- trunc(n2)
    nm <- np <- n2t
    mcdx[1:nm] <- mcdx[nm + 1]
    mcdx[(length(mcdx) - np + 1):length(mcdx)] <-
      mcdx[length(mcdx) - np]
  }

  # Execute step II of Bry-Boschan procedure
  pt1 <- dates1(xma)

  if ((sum(pt1$peak) == 0) | (sum(pt1$trough) == 0)) {
    warning("No peaks or no troughs for this series. Processing stops")
    bcp5 <- NA
    bct5 <- NA
    # GOTO bottom
  } else {
    tt <- as.matrix(1:nd)
    bcp2 <- tt[pt1$peak == 1,,drop = F]
    bct2 <- tt[pt1$trough == 1,,drop = F]

    # Make sure dates alternate, peak then trough
    # if not choose highest peak and lowest trough
    bc2 <- alter2(bcp2, bct2, xma)
    bcp2 <- bc2$peak
    bct2 <- bc2$trough
    rm(bc2)

    # Execute step III of BRy-Boschan procedure
    # A. refine previous dates to n6 using xspb
    #    note: (+/-)6 used instead of BB-Book (+/-)5 -- see documentation
    #    of BB program TP2 0570 and MDP 0220, 02230 and 0610

    bc3 <- refine(bcp2, bct2, xspb, 6)
    bcp3 <- bc3$peak
    bct3 <- bc3$trough
    rm(bc3)
    # B. enforce 15 month minimum P-P and T-T cycles
    bc3 <- enfvd(bcp3, bct3, xspb)
    bcp3 <- bc3$peak
    bct3 <- bc3$trough
    rm(bc3)

    # Execute Step IV of Bry-BOschan Procedure
    # A. refine previous dates to n6 using of mcdx
    #    note: (+/-) 6 used instead of BB-Book (+/-)5 -- see documentation
    #    of BB program TP2 0570 and MDP 0220, 02230 and 0610
    bc4 <- refine(bcp3, bct3, mcdx, 6)
    bcp4 <- bc4$peak
    bct4 <- bc4$trough
    rm(bc4)

    span <- max(c(mcdnum, 4))

    # A. Refine previous dates to nmax(4,mcd) using X
    bc5 <- refine(bcp4, bct4, x, span)
    bcp5 <- bc5$peak
    bct5 <- bc5$trough
    rm(bc5)

    # B. Eliminate turns with 6 months of beginning and end of series
    bc5 <- enfvb(bcp5, bct5, x)
    bcp5 <- bc5$peak
    bct5 <- bc5$trough
    rm(bc5)

    # C. Elimination of peaks and troughs at both ends which are lower
    #    or higher than values closer to end
    bc5 <- enfvc(bcp5, bct5, x)
    bcp5 <- bc5$peak
    bct5 <- bc5$trough
    rm(bc5)

    # D. Elimination of P-P and T-T cycles less than 15 months
    bc5 <- enfvd(bcp5, bct5, x)
    bcp5 <- bc5$peak
    bct5 <- bc5$trough
    rm(bc5)

    # E. Eliminate phases whose duration is less than 5 months
    bc5 <- enfve(bcp5, bct5, x)
    bcp5 <- bc5$peak
    bct5 <- bc5$trough
    rm(bc5)

  } # a partir daqui eh o label 'bottom'



  #' cal_ind_p = zeros(rows(cal),1);
  #' cal_ind_t = zeros(rows(cal),1);
  #' tmp_p = zeros(rows(cal),1);
  #' tmp_t = zeros(rows(cal),1);
  #' if sumc(bcp5) .> 0;
  #' cal_ind_p[bcp5]=ones(rows(bcp5),1);
  #' tmp_p[bcp5+1]=ones(rows(bcp5),1);
  #' endif;
  #' if sumc(bct5) .> 0;
  #' cal_ind_t[bct5]=ones(rows(bct5),1);
  #' tmp_t[bct5+1]=ones(rows(bct5),1);
  #' endif;
  #' rec_ind=zeros(rows(cal),1);
  #'
  #' rec_ind=cumsumc(tmp_p-tmp_t);
  #' if bct5[1] .< bcp5[1];
  #' rec_ind=1+cumsumc(tmp_p-tmp_t);
  #' endif;
  #'
  #' @ Correct for Missing Values @
  #'   bcp5=bcp5+i1;
  #' bct5=bct5+i1;
  #'
  #' y_adj=xo;                  @ X adjusted for outliers @
  #'   if i1 .> 0;
  #' cal_ind_p = miss(zeros(i1,1),0)|cal_ind_p;
  #' cal_ind_t = miss(zeros(i1,1),0)|cal_ind_t;
  #' rec_ind = miss(zeros(i1,1),0)|rec_ind;
  #' y_adj=miss(zeros(i1,1),0)|y_adj;
  #' endif;
  #' if i2 .> 0;
  #' cal_ind_p = cal_ind_p|miss(zeros(i2,1),0);
  #' cal_ind_t = cal_ind_t|miss(zeros(i2,1),0);
  #' rec_ind = rec_ind|miss(zeros(i2,1),0);
  #' y_adj = y_adj|miss(zeros(i2,1),0);
  #' endif;
  #'
  #' retp(rec_ind,cal_ind_p,cal_ind_t,bcp5,bct5);
  #'
return(list(peak = zoo::as.yearmon(time(x)[bcp5]),
            trough = zoo::as.yearmon(time(x)[bct5])))
}


# Function alter2 ---------------------------------------------------------

alter2 <- function(bcpn, bctn, x) {
  # Enforces alternation between peaks and troughs

  # skip the procedure if either peak or trough contain missing values
  if (nrow(bcpn) == 0 | nrow(bctn) == 0) {
    warning("Missing values in the peak or trough")
    return(list(peak = bcpn, trough = bctn))
  }

  bcpn <- cbind(bcpn, 1)
  bctn <- cbind(bctn, 0)
  anmat <- rbind(bcpn, bctn)
  anmat <- anmat[order(anmat[, 1]), ]

  nv <- 0
  j  <- 1
  while (j <= (nrow(anmat) - 1)) {
    cval <- anmat[j, 2]
    nval <- anmat[j + 1, 2]
    if (nval == cval) {
      if (cval == 1) {
        # eliminate lesser of two peaks in case of eliminate last peak
        s <- as.matrix(rep(FALSE, nrow(anmat)))
        d1 <- anmat[j, 1]
        d2 <- anmat[j + 1, 1]
        if (x[d1] >= x[d2])
          s[j + 1,1] <- TRUE
        else
          s[j,1] <- TRUE
        anmat <- anmat[!s,,drop = F]
        nv <- 1
        j  <- j - 1
      } else {
        # eliminate greater of two troughs in case of eliminate first trough
        s <- as.matrix(rep(FALSE, nrow(anmat)))
        d1 <- anmat[j, 1]
        d2 <- anmat[j + 1, 1]
        if (x[d1] >= x[d2])
          s[j,1] <- TRUE
        else
          s[j + 1,1] <- TRUE
        anmat <- anmat[!s,,drop = F]
        nv <- 1
        j  <- j - 1
      }
    }
    j <- j + 1
  }

  bcpn <- anmat[anmat[, 2] == 1, 1, drop = F]
  bctn <- anmat[anmat[, 2] == 0, 1, drop = F]
  return(list(peak = bcpn, trough = bctn))
}


# Function bbout1, bbout2 (joined in bbout) -------------------------------

bbout <- function(x, xsp, ilog) {
  # Checks fo outliers in series and replaces with fitted values
  # from spencer curve.
  # Additive deviations from series used to check for outliers
  # (appropriate for data in logarithms)
  # Ratio-deviations from series used to check for outliers
  # (appropriate for data in levels)

  nd <- length(x)
  if (ilog == 1) {
    d <- x - xsp
  } else {
    d <- x / xsp
  }
  ds <- scale(d) # returns a matrix (no ts attributes)
  dsi <- abs(ds) >= 3.5
  # as in BB TP1 line 1340-1350 no changes to first and last obs
  dsi[1] <- FALSE
  dsi[nd] <- FALSE

  xo <- x
  tt <- 1:nd
  if (sum(dsi) > 0) {
  xo[dsi] <- xsp[dsi]
  }

  return(xo)
}


# Function dates1 ---------------------------------------------------------

dates1 <- function(x12) {
  nd <- length(x12)
  peak1 <- as.matrix(rep(0, nd))
  trough1 <- as.matrix(rep(0, nd))

  # Step 1
  # Determine dates that are higher than any values n6 months
  # Determine dates that are lower than any values n6 months
  # These are initial peaks and troughs
  # Note:
  #  BB documentation and books says that checks are made for
  #  n5 months. However, from their subroutine WTP, they look
  #  for local min/mas in n6 range.

  # Handling of initial and terminal dates follows Bry-Bosch
  #  subroutine WTP

  for (i in 2:6) {
    if (x12[i] == max(x12[1:(i + 6)]))
      peak1[i,1]   <- 1
    if (x12[i] == min(x12[1:(i + 6)]))
      trough1[i,1] <- 1
  }
  for (i in 7:(nd - 1)) {
    np <- min(c(i + 6, nd))
    if (x12[i] == max(x12[(i - 6):np]))
      peak1[i,1]   <- 1
    if (x12[i] == min(x12[(i - 6):np]))
      trough1[i,1] <- 1
  }
  return(list(peak = peak1, trough = trough1))
}


# Function enfvb ----------------------------------------------------------

enfvb <- function(bcp, bct, x) {
  # Eliminate turns with 6 months of beginning and end of series

  # skip the procedure if either peak or trough contain missing values
  if (nrow(bcp) == 0 | nrow(bct) == 0) {
    warning("Missing values in the peak or trough")
    return(list(peak = bcp, trough = bct))
  }

  nv <- 0
  # check peaks
  s <- as.matrix(rep(FALSE, nrow(bcp)))
  if (bcp[1,1] <= 6) {
    s[1,1] <- TRUE
    nv <- 1
  }
  if (bcp[nrow(bcp),1] > (length(x) - 6)) {
    s[nrow(bcp),1] <- TRUE
    nv <- 1
  }
  bcp <- bcp[!s,,drop = F]
  if (nrow(bcp) == 0 | nrow(bct) == 0) {
     return(list(peak = bcp, trough = bct))
  }

  # check troughs
  s <- as.matrix(rep(FALSE, nrow(bct)))
  if (bct[1,1] <= 6) {
    s[1,1] <- TRUE
    nv <- 1
  }
  if (bct[nrow(bct),1] > (length(x) - 6)) {
    s[nrow(bct),1] <- TRUE
    nv <- 1
  }
  bct <- bct[!s,,drop = F]
  if (nrow(bcp) == 0 | nrow(bct) == 0) {
    return(list(peak = bcp, trough = bct))
  }

  if (nv != 0) {
    bc <- alter2(bcp, bct, x)
    bcp <- bc$peak
    bct <- bc$trough
    rm(bc)
  }
  return(list(peak = bcp, trough = bct))
}



# Function enfvc ----------------------------------------------------------

enfvc <- function(bcp, bct, x) {
  # Elimination of peaks (and troughs) at both ends which are lower (or higher)
  # than values closer to end
  nv <- 0
  #TOP
  while (1) {
    # A. check first date

    # Skip this procedure if either bcp or bct contain missing values
    if (nrow(bcp) == 0 | nrow(bct) == 0) {
      warning("Missing values in the peak or trough")
      return(list(peak = bcp, trough = bct))
    }

    m <- min(c(bcp, bct))
    if (m == min(bcp)) {
      # first TP is a peak
      s <- as.matrix(rep(FALSE, nrow(bcp)))
      d1 <- bcp[1, 1]
      if (x[1] > x[d1]) {
        s[1, 1] <- TRUE
        nv <- 1
        bcp <- bcp[!s,,drop = F]
        next # GOTO TOP
      }
    } else {
      # first TP is a trough
      s <- as.matrix(rep(FALSE, nrow(bct)))
      d1 <- bct[1, 1]
      if (x[1] < x[d1]) {
        s[1, 1] <- TRUE
        nv <- 1
        bcp <- bcp[!s,,drop = F]
        next # GOTO TOP
      }
    }

    # B. check last date
    m <- max(c(bcp, bct))
    if (m == max(bcp)) {
      # last TP is a peak
      s <- as.matrix(rep(FALSE, nrow(bcp)))
      d1 <- bcp[nrow(bcp), 1]
      if (x[length(x)] > x[d1]) {
        s[nrow(bcp), 1] <- TRUE
        nv <- 1
        bcp <- bcp[!s,,drop = F]
        next # GOTO TOP
      } else {
        # last TP is a trough
        s <- as.matrix(rep(FALSE, nrow(bct)))
        d1 <- bct[nrow(bct), 1]
        if (x[length(x)] < x[d1]) {
          s[nrow(bct), 1] <- 1
          nv <- 1
          bct <- bct[!s, , drop = F]
          next # GOTO TOP
        }
      }
    }
    break
  } # end TOP look

  if (nv != 0) {
    bc <- alter2(bcp, bct, x)
    bcp <- bc$peak
    bct <- bc$trough
    rm(bc)
  }
  return(list(peak = bcp, trough = bct))
}


# Function enfvd ----------------------------------------------------------

enfvd <- function(bcp, bct, x) {
  # Enforce 15 month P-P and T-T cycles.
  # If violated, choose highest peak and lowest trough

  # Skip this procedure if either bcp or bct contain missing values
  if (nrow(bcp) == 0 | nrow(bct) == 0) {
    warning("Missing values in the peak or trough")
    return(list(peak = bcp, trough = bct))
  }

  nv <- 0
  # Check peak to peak
  i <- 2
  while (i <= nrow(bcp)) {
    if ((bcp[i,1] - bcp[i - 1,1]) < 15) {
      s <- as.matrix(rep(FALSE, nrow(bcp)))
      # note eliminate latest peak if equal
      d1 <- bcp[i - 1,1]
      d2 <- bcp[i,1]
      if (x[d1] >= x[d2]) s[i,1] <- TRUE else s[i - 1,1] <- TRUE
      bcp <- bcp[!s,,drop = F]
      nv <- 1
      if (nrow(bcp) == 0 | nrow(bct) == 0) {
        return(list(peak = bcp, trough = bct))
      }
      i <- i - 1
    }
    i <- i + 1
  }
  if (nv != 0) {
    bc <- alter2(bcp, bct, x)
    bcp <- bc$peak
    bct <- bc$trough
    rm(bc)
  }

  nv <- 0
  # Check Trough to Trough
  i <- 2
  while (i <= nrow(bct)) {
    if ((bct[i,1] - bct[i - 1,1]) < 15) {
      s <- as.matrix(rep(FALSE, nrow(bct)))
      d1 <- bct[i - 1,1]
      d2 <- bct[i,1]
      # Note: eliminate first trough if equal
      if (x[d1] < x[d2]) s[i,1] <- TRUE else s[i - 1,1] <- TRUE
      bct <- bct[!s,,drop = F]
      nv <- 1
      if (nrow(bcp) == 0 | nrow(bct) == 0) {
        return(list(peak = bcp, trough = bct))
      }
      i <- i - 1
    }
    i <- i + 1
  }
  if (nv != 0) {
    bc <- alter2(bcp, bct, x)
    bcp <- bc$peak
    bct <- bc$trough
    rm(bc)
  }
  return(list(peak = bcp, trough = bct))
}


# Function enfve ----------------------------------------------------------

enfve <- function(bcp, bct, x){
  # This follows the BB code, subroutine TP3 (lines 0250-580)
  # If a violation is found the two turning points are eliminated
  # If the violation occurs at the last turning point, then only
  # the last turn is eliminated

  # Skip this procedure if either bcp or bct contain missing values
  if (nrow(bcp) == 0 | nrow(bct) == 0) {
    warning("Missing values in the peak or trough")
    return(list(peak = bcp, trough = bct))
  }

  bcp <- cbind(bcp,1)
  bct <- cbind(bct,0)
  anmat <- rbind(bcp, bct)
  anmat <- anmat[order(anmat[,1]),]

  nv <- 0

  # Check all phases except last
  j <- 1
  while (j <= (nrow(anmat) - 2)) {
    # bcp <- cbind(bcp,1)
    # bct <- cbind(bct,0)
    # anmat <- rbind(bcp, bct)
    # anmat <- anmat[order(anmat[,1]),]
    # CODIGO ACIMA PARECIA ERRADO NO ORIGINAL...
    if ((j + 1) > (nrow(anmat) - 1)) break # This is included because dates
    # and anmat is modified in loop
    cind <- anmat[j,1]
    nind <- anmat[j + 1,1]
    if ((nind - cind) < 5) {
      nv <- 1
      s <- as.matrix(rep(FALSE, nrow(anmat)))
      s[j,1] <- TRUE
      s[j + 1,1] <- TRUE
      anmat <- anmat[!s,,drop = F]
      bcp <- anmat[anmat[,2] == 1,1,drop = F]
      bct <- anmat[anmat[,2] == 0,1,drop = F]
      if ((nrow(bcp) == 0) | (nrow(bct) == 0)) {
        return(list(peak = bcp, trough = bct))
      }
    }
    j <- j + 1
  }

  # Check final phase
  cind <- anmat[nrow(anmat) - 1,1]
  nind <- anmat[nrow(anmat),1]
  if ((nind - cind) < 5) {
    nv <- 1
    anmat <- anmat[1:(nrow(anmat) - 1),, drop = F]
    bcp <- anmat[anmat[,2] == 1,1,drop = F]
    bct <- anmat[anmat[,2] == 0,1,drop = F]
    if ((nrow(bcp) == 0) | (nrow(bct) == 0)) {
      return(list(peak = bcp, trough = bct))
    }
  }
  return(list(peak = bcp[,1], trough = bct[,1]))
}


# Function MA -------------------------------------------------------------

MA <- function(x, n) {
  if (n == 1) return(x)
  xsp <- forecast::ma(x, n)
  # First nm observations (As in BB for MA12 -- modify for mcd)
  n2 <- n / 2
  n2t <- trunc(n2)
  nm <-
    np <- n2t # no need for the above code as MA is always centered
  xsp[1:nm] <-
    (cumsum(x[1:(nm + np)]) / 1:(nm + np))[(np + 1):(nm + np)]

  # Last np observations (As in BB for MA12 -- modify for mcd)
  nd <- length(x)
  xsp[(nd - np + 1):nd] <-
    rev((cumsum(x[rev((nd - (np + nm - 1)):nd)]) / 1:(nm + np))[1:nm])
  return(xsp)
}


# Function mcd1, mcd2 (combined in mcd) -----------------------------------

mcd <- function(series, h, ilog) {
  # This procedure calculates the months of cyclical dominance of a series
  # With logs differences are used

  c  <- spencer(series)
  nd <- length(series)
  if (ilog == 1) {
    i <- series - c
  } else {
    i <- series / c
  }

  mcdv <- vector('numeric', length = h)

  if (ilog == 1) {
    for (j in 1:h) {
      num <- sum(abs(diff(i,j)))
      den <- sum(abs(diff(c,j)))
      mcdv[j] <- num/den
    }
  } else {
    for (j in 1:h) {
      num <- sum(abs(i[(1 + j):nd] / i[1:(nd - j)] - 1))
      den <- sum(abs(c[(1 + j):nd] / c[1:(nd - j)] - 1))
      mcdv[j] <- num/den
    }

  }


  if (is.na(which(mcdv < 1)[1])) {
    stop("MCD beyond given h-limit. Increase h !!!")
  } else {
    mcd <- which(mcdv < 1)[1]
  }
  return(mcd)
}


# Function refine ---------------------------------------------------------

refine <- function(bcp, bct, x, n) {
  # Refine turning point dates
  # The original dates are bcp (peaks) and bct (troughs)
  # these are refined using the series x. The new dates are
  # those highest (for peaks) or lowest (for troughs) which are
  # within n periods of the old dates
  #
  # Notes, peaks and troughs must alternate, and this is checked and
  # the dates modified if necessary

  # for each peak date find highest n months
  bcpn <- bcp
  for (i in 1:nrow(bcp)) {
    d  <- bcp[i,1]
    i1 <- max(c(n / 2, d - n))
    i1 <- ceiling(i1)
    i2 <- min(c(length(x) - n / 2, d + n))
    i2 <- floor(i2)
    tt <- i1:i2
    b  <- cbind(x[i1:i2], tt)
    b  <- b[order(b[, 1]), ]
    # check for ties - chose earliest of the ties (earliest peak)
    b <- b[b[, 1] == b[nrow(b), 1], , drop = F]
    bcpn[i,1] <- min(b[, 2])
  }

  # for each trough date find lowest in spencer curve n 5 months
  bctn <- bct
  for (i in 1:nrow(bct)) {
    d  <- bct[i,1]
    #tt <- -n:n  ?? bug in the code?
    i1 <- max(c(n / 2, d - n))
    i1 <- ceiling(i1)
    i2 <- min(c(length(x) - n / 2, d + n))
    i2 <- floor(i2)
    tt <- i1:i2
    b  <- cbind(x[i1:i2], tt)
    b  <- b[order(b[, 1]), ]
    # check for ties -- choose latest of the ties (latest trough)
    b <- b[b[, 1] == b[1, 1], , drop = F] ## DIFERENTE DO loop anterior...
    bctn[i,1] <- max(b[, 2])
  }

  # make sure dates alternate
  return(alter2(bcpn, bctn, x))
}



### END END END

#
# myfun <- function(x) list(a = x, b = 2*x)
# lresult <- structure(NA,class = "result")
# "[<-.result" <- function(x,...,value) {
#   args <- as.list(match.call())
#   args <- args[-c(1:2,length(args))]
#   length(value) <- length(args)
#   for (i in seq(along = args)) {
#     a <- args[[i]]
#     if (!missing(a)) eval.parent(substitute(a <- v,list(a = a,v = value[[i]])))
#   }
#   x
# }
