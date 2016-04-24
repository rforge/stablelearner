### -- user interface ----------------------------------------------------------

stabletree <- function(x, data = NULL, sampler = bootstrap,
  applyfun = NULL, cores = NULL, ...)
{
  ## process sampler
  if (is.character(sampler)) 
    sampler <- get(sampler, mode = "function", envir = parent.frame())
  if (is.function(sampler)) {
    samp <- try(sampler(...), silent = TRUE)
    if (!inherits(samp, "try-error") && is.list(samp) && c("method", "sampler") %in% 
      names(samp)) {
      sampler <- samp
    } else {
      sampler <- list(method = "User-defined", sampler = sampler)
    }
  }
  
  ## data extraction
  if (is.null(data)) 
    data <- eval(getCall(x)$data)  #FIXME# more elegant default?
  n <- nrow(data)
  
  ## draw bootstrap samples
  bix <- sampler$sampler(n)
  B <- ncol(bix)
  
  ## facilitate parallelization
  if (is.null(applyfun)) {
    applyfun <- if (is.null(cores)) {
      lapply
    } else {
      function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
    }
  }
  
  ## bootstrap trees (and omit data copy from tree)
  xx <- applyfun(1L:B, function(i) {
    datai <- data[na.omit(bix[, i]), , drop = FALSE]
    xi <- update(x, data = datai)
    if (!inherits(xi, "party")) 
      xi <- partykit::as.party(xi)
    xi$data <- xi$data[0L, , drop = FALSE]
    return(xi)
  })
  
  ## extract names of all variables and omit response (FIXME: currently assuming a
  ## single response)
  nm <- all.vars(terms(x))
  cl <- sapply(data[, nm], function(x) class(x)[1L])
  yi <- attr(terms(x), "response")
  nm <- nm[-yi]
  
  ## convert original tree to party (if necessary)
  if (!inherits(x, "party")) 
    x <- as.party(x)
  
  ## function for extracting splits from trees
  extract_split <- function(x) {
    ids <- nodeids(x)
    ids <- ids[-nodeids(x, terminal = TRUE)]
    nodeapply(x, ids = ids, FUN = split_node)
  }
  
  ## function for extracting variable id from trees
  extract_varid <- function(x) {
    nm <- attr(terms(x), "term.labels")
    sp <- extract_split(x)
    if (!is.null(sp)) {
      vi <- sapply(sp, "[[", "varid") - 1L
      vi <- sort(unique(vi))
    } else vi <- NULL
    vi <- as.numeric(nm %in% nm[vi])
    names(vi) <- nm
    return(vi)
  }
  
  ## function for split information from trees
  extract_breaks <- function(x) {
    sp <- extract_split(x)
    if (!is.null(sp)) {
      vi <- sapply(sp, "[[", "varid") - 1L
    } else vi <- NULL
    br <- lapply(sp, "[[", "breaks")
    id <- lapply(sp, "[[", "index")
    names(id) <- names(br) <- nm[vi]
    br <- lapply(nm, function(n) {
      brs <- br[names(br) == n]
      ids <- id[names(id) == n]
      if (length(brs) > 0L || length(ids) > 0L) {
        if (is.null(brs[[n]])) {
          ans <- do.call("rbind", ids)
          if (!is.null(ans)) {
          rownames(ans) <- NULL
          colnames(ans) <- levels(data[, n])
          }
        } else {
          ans <- unlist(brs)
          names(ans) <- NULL
        }
      } else ans <- NULL
      return(ans)
    })
    names(br) <- nm
    return(br)
  }
  
  ## add levels to list with breakpoints
  add_levels <- function(x) {
    nm <- names(x)
    ans <- lapply(nm, function(n) {
      br <- x[[n]]
      if (!is.null(br)) {
        if (cl[n] == "ordered") 
          br <- ordered(br, levels = 1L:nlevels(data[, n]), labels = levels(data[, 
          n]))
        br
      } else NULL
    })
    names(ans) <- nm
    return(ans)
  }
  
  node_level <- function(x, terminal = FALSE) {
    nl <- function(x, level = 1) {
      if (!is.terminal(x)) {
        l <- lapply(kids_node(x), function(x) nl(x, level = level + 1))
        level <- c(level, l)
      }
      level
    }
    ans <- unlist(nl(node_party(x)))
    ids <- nodeids(x)
    names(ans) <- ids
    if (!terminal) {
      term <- unlist(nodeapply(x, ids = ids, is.terminal))
      ans[!term]
    } else ans
  }
  
  extract_splitinfo <- function(x) {
    sp <- extract_split(x)
    if (!is.null(sp)) {
      br <- add_levels(extract_breaks(x))
      vi <- sapply(sp, "[[", "varid") - 1L
      lv <- node_level(x)
      ids <- nodeids(x)
      ids <- ids[-nodeids(x, terminal = TRUE)]
      names(lv) <- nm[vi]
      names(ids) <- nm[vi]
      ninfo <- lapply(names(br), function(n) {
        if (!is.null(br[[n]])) {
          list(nodeids = unname(ids[names(ids) == n]), levels = unname(lv[names(lv) == 
          n]), breaks = br[[n]])
        } else NULL
      })
      names(ninfo) <- names(br)
      ninfo
    } else NULL
  }
  
  ## selection proportions
  vi <- applyfun(xx, FUN = extract_varid)
  vi_mat <- do.call("rbind", vi)
  
  ## breakpoints
  br <- applyfun(xx, FUN = extract_breaks)
  br <- lapply(nm, function(n) {
    if (cl[n] == "factor") {
      do.call("rbind", lapply(br, "[[", n))
    } else {
      unlist(lapply(br, "[[", n))
    }
  })
  names(br) <- nm
  
  ## collect observed and bootstrapped results
  rval <- list(
    call = getCall(x),
    B = B,
    sampler = sampler,
    vs0 = extract_varid(x), 
    br0 = extract_splitinfo(x),
    vs = vi_mat,
    br = add_levels(br),
    classes = cl[-yi]
  )
  class(rval) <- "stabletree"
  return(rval)
}

### -- simple standard methods -------------------------------------------------

print.stabletree <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nSampler:\n")
  cat("B =", x$B, "\n")
  cat("Method =", x$sampler$method)
  cat("\n")
}

summary.stabletree <- function(object, show.breaks = FALSE, digits = 3, ...)
{
  ans <- list()
  ans$B <- object$B
  ans$call <- object$call
  ans$method <- object$sampler$method

  varselect <- function(x, digits = 3) {
    nm <- names(x$classes)
    avsc <- round(sapply(nm, function(n) {
      bri <- x$br[[n]]
      if (is.matrix(bri)) nrow(bri)/x$B else length(bri)/x$B
    }), digits)
    avsc.star <- round(sapply(nm, function(n) {
      bri <- x$br0[[n]]$br
      if (is.matrix(bri)) nrow(bri) else length(bri)
    }), digits)
    ans <- data.frame(round(colMeans(x$vs), digits), x$vs0, avsc, avsc.star)
    colnames(ans) <- c("freq", "*", "mean", "*")
    ans[order(ans[, "freq"], decreasing = TRUE), ]
  }
  ans$vstab <- varselect(object, digits = digits)

  if (show.breaks) {
    breaks <- function(x, digits) {
      nm <- names(x$br)
      ans <- lapply(nm, function(n) {
        br <- x$br[[n]]
        if (is.numeric(br)) 
          br <- round(br, digits)
        sort(table(br, dnn = NULL), decreasing = TRUE)
      })
      names(ans) <- nm
      ans
    }
    ans$br <- breaks(object, digits = digits)
  }

  class(ans) <- "summary.stabletree"
  ans
}

print.summary.stabletree <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nSampler:\n")
  cat("B =", x$B, "\n")
  cat("Method =", x$method)
  cat("\n\nVariable selection overview:\n\n")
  print(x$vstab)
  cat("(* = complete data tree)\n")
  if (!is.null(x$br)) {
    cat("\n\n")
    print(x$br)
  }
}

### -- graphical methods -------------------------------------------------------

barplot.stabletree <- function(height, main = "Variable selection frequencies",
  xlab = "", ylab = "Relative frequency (in %)", horiz = FALSE, col = gray.colors(2),
  names.arg = NULL, names.uline = TRUE, names.diag = TRUE,
  cex.names = 0.9, 
  ylim = if (horiz) NULL else c(0, 100), xlim = if (horiz) c(0, 100) else NULL, ...)
{  
  vsp <- colMeans(height$vs)
  ord <- order(vsp, decreasing = TRUE)
  vsp <- vsp[ord] * 100
  vs0 <- height$vs0[ord]
  
  if (is.null(names.arg)) 
    names.arg <- names(vsp)
  
  labs <- sapply(seq_along(names.arg), function(i) format_labels(names.arg[i], 
    uline = if (names.uline) 
      vs0[i] else 0))
  
  par0 <- par()
  
  if (horiz) {
    mai <- par0$mai
    mai[2] <- max(strwidth(labs, "inches")) + 0.3
    par(mai = mai)
  }
  
  b <- barplot(vsp, col = col[1L + (vs0 < 1)], names.arg = if (names.diag & !horiz) 
    NA else labs, horiz = horiz, ylim = ylim, xlim = xlim, cex.names = cex.names, main = main, 
    xlab = xlab, ylab = ylab, axes = FALSE, las = ifelse(horiz, 2, 1), ...)
  
  if (horiz) {
    axis(1)
  } else {
    axis(2, las = 2)
    if (names.diag) {
      # draw_labels <- function(x, labels, line = 1, ...) {
      #   par0 <- par()
      #   dy_us <- diff(par0$usr[3:4])
      #   dy_in <- par0$pin[2]
      #   line <- line + (1 - par0$ylbias)
      #   line_inch <- line * par0$cin[2]
      #   yshift <- line_inch/(dy_in/dy_us)
      #   text(x, par0$usr[3] - yshift, labels, xpd = TRUE, ...)
      # }
      # draw_labels(b, labels = labs, srt = 45, adj = c(1,1), cex = cex.names)
      text(b, par0$usr[3] - 2, labels = labs, srt = 45, adj = c(1, 1), xpd = TRUE, 
        cex = cex.names)
    }
  }
  
  par(mai = par0$mai)
  
}

image.stabletree <- function(x, main = "Variable selections",
  ylab = "Repetitions", xlab = "", col = gray.colors(2),
  names.arg = NULL, names.uline = TRUE, names.diag = TRUE, 
  cex.names = 0.9, xaxs = "i", yaxs = "i",
  col.tree = 2, lty.tree = 2,
  xlim = c(0, length(x$vs0)), ylim = c(0, x$B), ...)
{
  ord <- ordermat(x$vs)
  z <- 1L - ord$x
  
  nr <- nrow(z)
  nc <- ncol(z)
  
  plot(xlim, ylim, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xaxs = xaxs, 
    yaxs = yaxs, xlab = xlab, ylab = ylab, main = main, ...)
  
  sapply(1:nc, function(j) {
    r <- rle(z[, j])
    y <- c(0, cumsum(r$lengths), x$B)
    sapply(2:length(y), function(k) rect(j - 1, y[k - 1], j, y[k], col = col[r$values[k - 
      1] + 1], border = NA))
  })
  
  axis(2)
  box()
  grid(nx = nc, ny = NA, col = "darkgray", lty = "solid")
  
  if (is.null(names.arg)) 
    names.arg <- colnames(z)
  
  labs <- sapply(seq_along(names.arg), function(i) format_labels(names.arg[i], 
    uline = if (names.uline) 
      x$vs0[ord$colind][i] else 0))
  
  if (names.diag) {
    text(seq(nc) - 0.5, par("usr")[3] - 0.02 * x$B, labels = labs, srt = 45, 
      adj = c(1, 1), xpd = TRUE, cex = cex.names)
  } else {
    axis(1, at = seq(nc) - 0.5, labels = labs, lty = 0, cex.axis = cex.names, 
      ...)
  }
  
  if (!is.null(x$vs0)) {
    vs <- ord$x
    rownames(vs) <- apply(vs, 1, paste, collapse = "-")
    vs0 <- paste(x$vs0[ord$colind], collapse = "-")
    eq <- rownames(vs) %in% vs0 + 0L
    yy <- which(abs(diff(c(0, eq, 0))) > 0) - 1
    abline(h = yy, col = col.tree, lty = lty.tree)
    axis(4, at = yy, labels = NA, col = col.tree, lwd = 1, line = 0.2)
  }

}

plot.stabletree <- function(x, select = order(colMeans(x$vs), decreasing = TRUE), 
  type.breaks = "levels", col.breaks = "red", lty.breaks = "dashed", cex.breaks = 0.7, 
  col.main = c("black", "gray50"), main.uline = TRUE, args.numeric = NULL, args.factor = NULL, 
  args.ordered = NULL, main = NULL, ...)
{
  br <- x$br
  cl <- x$classes
  if (is.character(select)) 
    select <- sapply(select, function(n) which(names(br) == n))
  nplt <- length(select)
  if (nplt < 1L) {
    message("Nothing to plot!")
    return(invisible(NULL))
  } else {
    par(mfrow = n2mfrow(nplt))
  }
  for (i in select) {
    if (sum(x$vs[, i]) > 0L) {
      bri <- br[[i]]
      br0 <- x$br0[[i]]$br
      tx0 <- x$br0[[i]][[type.breaks]]
      args <- list(bri = bri, br0 = br0, tx0 = tx0, B = x$B, col.breaks = col.breaks, 
        lty.breaks = lty.breaks, cex.breaks = cex.breaks)
      switch(cl[i],
        "numeric" = do.call("breaks_hist", c(args, args.numeric)), 
        "integer" = do.call("breaks_hist", c(args, args.numeric)),
      	"factor" = do.call("breaks_image", c(args, args.factor)),
	      "ordered" = do.call("breaks_barplot", c(args, args.ordered)),
      	NULL
      )
      abline(h = x$B, col = "black", lty = "dotted")
    } else {
      ## empty plot
      plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")
      text(0, 0, labels = "nothing to plot")
    }
    t_i <- if(is.null(main)) names(br)[i] else main[i]
    if(nchar(t_i)>0) {
      title(main = format_labels(t_i,
        uline = if (main.uline) x$vs0[i] else 0, bold = TRUE), 
        col.main = col.main[2L - x$vs0[i]], ...)
    }
  }
  par(mfrow = c(1, 1))
}

### -- graphical auxiliary functions -------------------------------------------

format_labels <- function(n, uline = FALSE, bold = FALSE)
{
  n <- paste0(n, "[]")
  if (bold) 
    n <- paste0("bold(", n, ")")
  if (uline) 
    n <- paste0("underline(", n, ")")
  return(parse(text = n))
}

ordermat <- function(x, order.rows = TRUE, order.cols = TRUE)
{
  if (order.cols) {
    colind <- order(colMeans(x, na.rm = TRUE), decreasing = TRUE)
  } else colind <- 1:ncol(x)
  if (order.rows) {
    rowind <- order(apply(x[, colind], 1, paste0, collapse = ""), decreasing = TRUE)
  } else rowind <- 1:nrow(x)
  return(list(x = x[rowind, colind], rowind = rowind, colind = colind))
}

breaks_barplot <- function(bri, br0 = NULL, tx0 = NULL, B = NULL,
  ylab = "Counts", xlab = "", col.breaks = "red", lty.breaks = "dashed",
  cex.breaks = 0.7, col = "#E6E6E6", ...)
{
  
  tb <- table(bri)
  
  xlim <- c(0.1, (length(tb) - 1) * 1.2 + 0.1)
  
  b <- barplot.default(tb[-length(tb)], ylim = c(0, max(B, max(tb))), col = col, 
    names.arg = NA, main = "", ylab = ylab, xlab = xlab, xlim = xlim, ...)
  at <- 0.1 + seq(0, length(tb) - 1) * 1.2
  axis(1, at = at, labels = names(tb), lwd = 0, lwd.ticks = 1)
  
  if (length(br0) > 0) {
    br0. <- as.numeric(br0) * 1.2 - 0.5
    abline(v = unique(br0.), col = col.breaks, lty = lty.breaks)
    if (!is.null(tx0)) {
      tx <- tapply(tx0, br0., paste0, collapse = "\n")
      mtext(tx, side = 3, at = names(tx), col = col.breaks, cex = cex.breaks)
    }
  }
  
}

breaks_hist <- function(bri, br0 = NULL, tx0 = NULL, B = NULL, breaks = "Sturges", 
  col = "#E6E6E6", ylab = "Counts", xlab = "", col.breaks = "red", lty.breaks = "dashed", 
  cex.breaks = 0.7, ...)
{
  
  if (length(bri) < 1L) 
    bri <- 0
  
  h <- hist.default(bri, breaks = breaks, plot = FALSE)  
  plot(h, main = "", ylab = ylab, xlab = xlab, col = col,
    ylim = c(0, max(B, max(h$counts))), ...)  
  rug(bri)
  
  if (length(br0) > 0) {
    abline(v = unique(br0), col = col.breaks, lty = lty.breaks)
    if (!is.null(tx0)) {
      tx <- tapply(tx0, br0, paste0, collapse = "\n")
      mtext(tx, side = 3, at = names(tx), col = col.breaks, cex = cex.breaks)
    }
  }
  
}

breaks_image <- function(bri, br0 = NULL, tx0 = NULL, B = NULL, ylab = "Repetitions", 
  xlab = "", col = c("#97BDE1", "#ECD1A5"), col.na = "#E6E6E6", col.breaks = "red", 
  lty.breaks = "dashed", cex.breaks = 0.7, xaxs = "i", yaxs = "i", 
  xlim = c(0, ncol(bri)), ylim = c(0, max(B, nrow(bri))), ...)
{
  
  z <- ordermat(bri, order.rows = TRUE, order.cols = FALSE)$x
  rownames(z) <- apply(z, 1, paste, collapse = "-")
  z[is.na(z)] <- 0L
  
  nr <- nrow(bri)
  nc <- ncol(bri)
  
  col <- c(col.na, col)

  plot(xlim, ylim, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, 
       xaxs = xaxs, yaxs = yaxs, xlab = xlab, ylab = ylab, ...)
  
  sapply(1:nc, function(j) {
    r <- rle(z[, j])
    y <- c(0, cumsum(r$lengths), B)
    sapply(2:length(y), function(k) rect(j - 1, y[k - 1], j, y[k], col = col[r$values[k - 1] + 1], border = NA))
  })
  
#   image(x = seq(nc), y = seq(0, nr), z = t(z), axes = FALSE, ylim = c(0, max(B, 
#     nr)), col = col, ylab = ylab, xlab = xlab)
  
  grid(nx = nc, ny = NA, col = "#4D4D4D", lty = "solid")
  axis(1, at = seq(nc) - 0.5, labels = colnames(bri), lwd = 0, lwd.ticks = 1)
  axis(2)
  
  if (!is.null(br0)) {
    rownames(br0) <- apply(br0, 1, paste, collapse = "-")
    if (is.matrix(tx0)) 
      tx0 <- rownames(br0)
    tx0 <- tapply(tx0, rownames(br0), paste0, collapse = "\n")
    br0 <- unique(br0)
    sapply(seq(nrow(br0)), function(i) {
      eq <- rownames(z) %in% rownames(br0)[i] + 0L
      yy <- which(abs(diff(c(0, eq, 0))) > 0) - 1
      abline(h = yy, col = col.breaks, lty = lty.breaks)
      axis(4, at = yy, labels = NA, col = col.breaks, lwd = 1, line = 0.2, 
        cex = cex.breaks)
      if (!is.null(tx0)) {
        mtext(tx0[i], side = 4, at = mean(yy), line = 0.5, col = col.breaks, 
          cex = cex.breaks, las = 2)
      }
    })
  }
  
}

### ----------------------------------------------------------------------------
