Abs_diff_multi2 = function (prop.real, prop.est, method.name = NULL, title = NULL, 
  eval = TRUE, ...) 
{
  # QWL: added the following 3 statements to avoid warnings: no visible binding for global variable
  CellType = NULL
  Sub = NULL
  Abs.Diff = NULL
  
  ct.real = colnames(prop.real)
  sub.real = rownames(prop.real)
  if (is.list(prop.est)) {
    L = length(prop.est)
    if (is.null(method.name)) {
      method.name = paste0("Est.Method", 1:L)
    }
    else {
      if (length(method.name) < L) {
        method.name = c(method.name, paste0("Est.Method", 
          1:(L - length(method.name))))
      }
      else {
        method.name = method.name[1:L]
      }
    }
    l.sub.est = lapply(prop.est, rownames)
    l.ct.est = lapply(prop.est, colnames)
    sub.est = Reduce(intersect, l.sub.est)
    ct.est = Reduce(intersect, l.ct.est)
    celltype = intersect(ct.real, ct.est)
    sub = intersect(sub.real, sub.est)
    N = length(sub)
    K = length(celltype)
    if (N < 1) {
      stop("No common Subjects! Check rowname!")
    }
    if (K <= 1) {
      stop("Not enough cell types!")
    }
    cm.prop.real = prop.real[match(sub, sub.real), match(celltype, 
      ct.real)]
    m.prop.real = data.frame(Prop = c(cm.prop.real), CellType = factor(rep(celltype, 
      each = N), levels = celltype), Sub = factor(rep(sub, 
      K), levels = sub), Method = rep("Real", N * K))
    abs.diff = NULL
    ann = NULL
    for (l in 1:L) {
      abs.diff.temp = m.prop.real
      colnames(abs.diff.temp)[1] = "Abs.Diff"
      cm.prop.est = prop.est[[l]][match(sub, l.sub.est[[l]]), 
        match(celltype, l.ct.est[[l]])]
      abs.diff.temp$Abs.Diff = abs(c(cm.prop.est) - m.prop.real$Prop)
      abs.diff.temp$Method = factor(rep(method.name[l], 
        K * N), levels = method.name)
      ann = c(ann, paste0("RMSD = ", round(sqrt(mean((cm.prop.real - 
        cm.prop.est)^2)), digits = 5), "\n mAD = ", 
        round(mean(abs.diff.temp$Abs.Diff), digits = 5), "(", method.name[l], ")\n"))
      #abs.diff = rbind(abs.diff, abs.diff.temp)
      # QWL: added "as.data.frame"
      abs.diff = as.data.frame(rbind(abs.diff, abs.diff.temp))
      
    }
  }
  else {
    ct.est = colnames(prop.est)
    sub.est = rownames(prop.est)
    celltype = intersect(ct.real, ct.est)
    sub = intersect(sub.real, sub.est)
    N = length(sub)
    K = length(celltype)
    if (N < 1) {
      stop("No common Subjects! Check rowname!")
    }
    if (K <= 1) {
      stop("Not enough cell types!")
    }
    cm.prop.real = prop.real[match(sub, sub.real), match(celltype, 
      ct.real)]
    cm.prop.est = prop.est[match(sub, sub.real), match(celltype, 
      ct.real)]
    m.prop.real = data.frame(Prop = c(cm.prop.real), CellType = factor(rep(celltype, 
      each = N), levels = celltype), Sub = factor(rep(sub, 
      K), levels = sub), Method = rep("Real", N * K))
    abs.diff = m.prop.real
    colnames(abs.diff)[1] = "Abs.Diff"
    abs.diff$Abs.Diff = abs(c(cm.prop.est - cm.prop.real))
    abs.diff$Method = rep(method.name, K * N)
    ann = paste0("RMSD = ", round(sqrt(mean((cm.prop.real - 
      cm.prop.est)^2)), digits = 5), 
		 "\n mAD = ", round(mean(abs.diff$Abs.Diff), 
      digits = 5))
  }
  if (is.null(title)) {
    title = "Heatmap of Absolute Difference |p - Est.p|"
  }
  ann = paste(ann, collapse = "")
  title = paste(title, "\n(", ann, ")", sep="")
  if (eval) {
    ggplot(abs.diff, aes(CellType, Sub)) + geom_tile(aes(fill = Abs.Diff), 
      colour = "white") + scale_fill_gradient(low = "white", 
      high = "steelblue", name = "Abs.Diff\n") + facet_wrap(~Method, 
      nrow = 1) + theme(axis.text.x = element_text(angle = -90, 
      size = 10, vjust = 0)) + ggtitle(title) #+ annotate("text", 
      #label = ann, x = round(4 * K/5), y = N - 0.5, size = 2.5, 
      #colour = "black")
  }
  else {
    ggplot(abs.diff, aes(CellType, Sub)) + geom_tile(aes(fill = Abs.Diff), 
      colour = "white") + scale_fill_gradient(low = "white", 
      high = "steelblue", name = "Abs.Diff\n") + facet_wrap(~Method, 
      nrow = 1) + theme(axis.text.x = element_text(angle = -90, 
      size = 10, vjust = 0)) + ggtitle(title)
  }
}

