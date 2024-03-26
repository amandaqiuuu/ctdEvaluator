Prop_comp_multi2=function (prop.real, prop.est, method.name = NULL, title = NULL, 
  eval = TRUE, ...) 
{
  # QWL: added the following 3 statements to avoid warnings: no visible binding for global variable
  CellType = NULL
  Sub = NULL
  Prop = NULL
  
  
  ct.real = colnames(prop.real)
  sub.real = rownames(prop.real)
  if (!is.list(prop.est)) {
    prop.est = list(prop.est)
  }
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
  m.prop.est = NULL
  ann = ""
  for (l in 1:L) {
    m.prop.temp = m.prop.real
    cm.prop.est = prop.est[[l]][match(sub, l.sub.est[[l]]), 
      match(celltype, l.ct.est[[l]])]
    m.prop.temp$Prop = c(cm.prop.est)
    m.prop.temp$Method = factor(rep(method.name[l], K * 
      N), levels = method.name)
    ann = c(ann, paste0("R = ", round(cor(c(cm.prop.real), 
      c(cm.prop.est)), digits = 4), "(", method.name[l], ") "))
    m.prop.est = rbind(m.prop.est, m.prop.temp)
  }
  #m.prop = rbind(m.prop.real, m.prop.est)
  # QWL: added 'as.data.frame'
  m.prop = as.data.frame(rbind(m.prop.real, m.prop.est))
  if (is.null(title)) {
    title = "Heatmap of Estimated and Real Proportion"
  }
  ann = paste(ann, collapse = "")
  title = paste(title, "\n(", ann, ")", sep="")
  if (eval) {
    ggplot(m.prop, aes(CellType, Sub)) + geom_tile(aes(fill = Prop), 
      colour = "white") + scale_fill_gradient2(low = "steelblue", 
      high = "red", mid = "white", midpoint = 0.5, limit = c(0, 
        1), name = "Est Prop\n") + facet_wrap(~Method, 
      nrow = 1) + theme(axis.text.x = element_text(angle = -90, 
      size = 10, vjust = 0)) + ggtitle(title) #+ annotate("text", 
      #label = ann, x = round(4 * K/5), y = N, size = 2.5, 
      #colour = "black")
  }
  else {
    ggplot(m.prop, aes(CellType, Sub)) + geom_tile(aes(fill = Prop), 
      colour = "white") + scale_fill_gradient2(low = "steelblue", 
      high = "red", mid = "white", midpoint = 0.5, limit = c(0, 
        1), name = "Est Prop\n") + facet_wrap(~Method, 
      nrow = 1) + theme(axis.text.x = element_text(angle = -90, 
      size = 10, vjust = 0)) + ggtitle(title)
  }
}

