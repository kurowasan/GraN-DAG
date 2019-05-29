# Code copyright by Prof. Peter Spirtes
# Used in this package with permission.

require("combinat")
require("ggm")

DEBUG = FALSE

NO_EDGE = 0
EMP = 1
ARR = 2
CIR = 3

setClass("dataset", representation(size = "numeric", variable_names =
                                   "character"))
setClass("continuousdata", contains = "dataset", representation =
         representation(data = "data.frame"))
setClass("covariancedata", contains = "dataset", representation =
         representation(cov_matrix = "data.frame"))

setGeneric("isIndependent", function(d, x, y, given = NULL, prgt = .05) {
#  if (DEBUG) cat("x:", x, "y:", y, "given:", given, "\n")
  standardGeneric("isIndependent")
})


setMethod("isIndependent", "continuousdata",
          function(d, x, y, given = NULL, prgt = .05) {
            model = lm(d@data[,y]~., data =  data.frame(d@data[,c(x,given)]))
            summary(model)$coefficients[2, "Pr(>|t|)"] > prgt
          })

setMethod("isIndependent", "covariancedata",
          function(d, x, y, given = NULL, prgt = .05) {
            p = pcor(c(x,y,given), d@cov_matrix)
            pcor.test(p, length(given), d@size)$pvalue > prgt
          })


make_continuousdata = function(d) {
  new("continuousdata", size = nrow(d), variable_names = dimnames(d)[[2]],
      data = d)
}

make_covariancedata = function(d, size) {
  new("covariancedata", size = size, variable_names = dimnames(d)[[2]],
      cov_matrix = d)
}

import_tetrad_file = function(filename) {
  f = file(filename, open = "r")
  firstline = readLines(f, 1)
  if (length(grep("/continuousdata", firstline)) == 1) {
    d = read.table(f, header = T)
    #return(make_continuousdata(d))
    return(make_covariancedata(data.frame(cov(d)),nrow(d)))
  } else if (length(grep("/Covariance", firstline)) == 1) {
    size = as.numeric(readLines(f, 1))
    d = read.table(f, header = T, fill = T)
    for (j in 2:ncol(d)) {
      for (i in 1:(j-1)) {
        d[i,j] = d[j,i]
      }}
    return(make_covariancedata(d, size))
  } else {
    stop("First line of file was ", firstline,
         "This format of file is unknown.")
  }
}

read_tetrad_output = function(filename,rout) {
  f = file(filename, open = "r")
  repeat {
      firstline = readLines(f, 1)
      if (length(grep("Graph Edges:", firstline)) == 1) break }
  tab = read.table(f,stringsAsFactors=FALSE)
  nm = dimnames(rout)[[1]]
  adj = matrix(nrow=length(nm),ncol=length(nm),data=0)
  attr(adj,'headnames') <- nm
  for (i in 1:dim(tab1)[1]){
         first = which(nm == tab[i,2])
         second = which(nm == tab[i,4])
         if (substr(tab[i,3],3,3) == "-") {end2 = CIR} else {if (substr(tab[i,3],3,3) == ">") {end2 = ARR}}
         if (end2 == ARR) end1 = EMP else end1 = CIR
         adj[first,second] = end2
         adj[second,first] = end1
 }
 return(adj)
}

printlist = function(ambiguous,res){
  n = dim(res)[[2]]
  nn = 1:n
  nm = dimnames(res)[[1]]
  if (length(ambiguous) != 0){
    cat("Ambiguous triples","\n")
    for (i in 1:dim(ambiguous)[[2]]){
        cat(nm[ambiguous[1,i]]," ",nm[ambiguous[2,i]]," ",nm[ambiguous[3,i]],"\n")}}}

printadj = function(res){
  n = dim(res)[[2]]
  nn = 1:n
  nm = dimnames(res)[[1]]
  for (i in 1:(n-1)){
     for (j in nn[res[i,] != 0]){
       if (j > i){
        if ((res[j,i] == ARR) && (res[i,j] != ARR)){
          first = j
          second = i} else {
            first = i
            second = j}
        cat(nm[first])
        if (res[second,first] == CIR) cat("o-") else
        if (res[second,first] == ARR) cat ("<-") else
        if (res[second,first] == EMP) cat(" -")
        if (res[first,second] == CIR) cat("o") else
        if (res[first,second] == ARR) cat(">") else
        if (res[first,second] == EMP) cat("-")
        cat(nm[second],"\n")}}}}


inambiguous = function(x,ambiguous){
  foundit = FALSE
  if (length(ambiguous) != 0){
    for (i in 1:dim(ambiguous)[[2]]){
      if (((x[1] == ambiguous[1,i]) && (x[2] == ambiguous[2,i]) && (x[3] == ambiguous[3,i])) ||
          ((x[1] == ambiguous[3,i]) && (x[2] == ambiguous[2,i]) && (x[3] == ambiguous[1,i]))){
        foundit = TRUE
      break}}}
  return(foundit)}

#added pc_ to skeleton by Antti
pc_skeleton = function(data, prgt = .05) {
  n = length(data@variable_names)
  skel = matrix(1, nrow=n, ncol=n)
  dimnames(skel) <- list(data@variable_names, data@variable_names)
  diag(skel) <- 0
  sepset = vector(length = n, mode = "list")
  for (i in 1:n) {
    sepset[[i]] = vector(length = n, mode = "list")
  }

  nn = 1:n

  for (x in 1:(n-1)) {
    for (y in (x+1):n) {
      if (isIndependent(data, x, y, prgt = prgt)) {
          skel[x,y] = skel[y,x] = 0
          sepset[[x]][[y]] = sepset[[y]][[x]] = numeric(0)
  }}}

  i=1
  while (max(skel %*% rep(1,n)) >= i) {
    for (x in 1:n) {
      if (sum(skel[x,]) <= i) next
      for (y in nn[skel[x,] == 1]) {
        if (sum(skel[x,]) <= i) break
        ## The combn function is obnoxiously inconsistent.  It behaves
        ## differently when the first argument has length one and
        ## produces a non-matrix output when there is only one combination
        ## in the output.
        tf = setdiff(nn[skel[x,] == 1], y)
        sepsets = if (length(tf) == i) matrix(tf) else combn(tf, i)

        for (idx in 1:ncol(sepsets)) {
          if (isIndependent(data, x, y, sepsets[,idx], prgt = prgt)) {
              skel[x,y] = skel[y,x] = 0
              sepset[[x]][[y]] = sepset[[y]][[x]] = sepsets[,idx]
              break
          }}}}

    i = i + 1
  }

  list(skel, sepset)
}

 reach = function(a,b,c,adjacency){
        #reachable      set of vertices;
        #edgeslist      array[1..maxvertex] of list of edges
        #more           Boolean
        #reachable      list of vertices
        #numvertex      integer
        #labeled        array (by depth) of list of edges that have been labeled

 makeedge = function(x,y)(list(list(x,y)))

 legal <- function(t,...) {UseMethod("legal")}

 legal.possibledsep = function(t,u,v,r,s) {
   if (((adjacency[r[[1]],r[[2]]] == ARR) &&
        (adjacency[s,r[[2]]] == ARR) &&
        (r[[1]] != s)) ||
       ((adjacency[r[[1]],s] != 0) &&
        (r[[1]] != s))){
           edgeslist[[r[[2]]]]  <<- setdiff(edgeslist[[r[[2]]]],s)
           makeedge(r[[2]],s)}}

  legal.discriminating = function(t,u,v,r,s) {
     if ((length(intersect(s,c(t,u,v))) == 0) &&   # s not in the triangle t,u,v
         (adjacency[s,r[[2]]] == ARR) &&            # s collides with r edge at r[[2]]
         (r[[1]] != s) &&                           # s is not on the r edge
         (adjacency[u,r[[2]]] == EMP) &&
         (adjacency[r[[2]],u] == ARR)){
           edgeslist[[r[[2]]]]  <<- setdiff(edgeslist[[r[[2]]]],s)
           makeedge(r[[2]],s)}}

 initialize <- function(x,...) UseMethod("initialize")

 initialize.possibledsep <- function(x,y,z) {mapply(makeedge,x=a,y=edgeslist[[a]])}

 initialize.discriminating <- function(x,y,z) {mapply(makeedge,x=a,y=setdiff(which(adjacency[,a] == ARR),c(b,c)))}

 labeled = list()
 numvertex = dim(adjacency)[1]
 edgeslist = list()
 for (i in 1:numvertex) edgeslist = c(edgeslist,list(which(adjacency[,i] != 0)))
 labeled[[1]] = initialize(a,b,c)
 edgeslist[[a]] = list()
 depth = 2
 repeat
       {more = FALSE
        labeled[[depth]] = list()
        for (i in 1:length(labeled[[depth-1]])) {
           edgestemp = edgeslist[[labeled[[depth-1]][[i]][[2]]]]
           if (length(edgestemp) == 0) break
           for (j in 1:length(edgestemp))
            {labeled[[depth]] =
              union(legal(a,b,c,labeled[[depth-1]][[i]],edgestemp[[j]]),labeled[[depth]])}}
        if (length(labeled[[depth]]) != 0){
              more = TRUE
              depth = depth  + 1} else break}
 return(unique(unlist(labeled)))
}

cpc_colliders = function(data,res,prgt,ambiguous) {
  n = nrow(res)
  nn = 1:n
  nm = dimnames(res)[[1]]
  all = FALSE
  none = FALSE
  for (x in 1:n) {
    for (y in nn[res[x,] != 0]) {
      for (z in nn[res[y,] != 0]) {
        if (x < z && res[x,z] == 0) {
           if (isIndependent(data, x, z, prgt = prgt)) {
             all = FALSE
             none = TRUE} else {
               all = TRUE
               none = TRUE}
           tf1 = nn[res[x,] != 0]                                                # adjacent to x
           tf2 = nn[res[z,] != 0]                                                # adjacent to z
           maxlength = max(length(tf1),length(tf2))
           for (i in 1:maxlength){
             condsets1 = list()                                                  # combinations of adjacent to x of length i
             if (length(tf1) >= i)
               condsets1 = if (i == length(tf1)) {matrix(tf1)} else {combn(tf1,i)}
             if (length(condsets1) == 0) xlength = 0 else xlength = dim(condsets1)[[2]]
             condsets2 = list()                                                  # combinations of adjacent to z of length i
             if (length(tf2) >= i)
               condsets2 = if (i == length(tf2)) {matrix(tf2)} else {combn(tf2,i)}
             if (length(condsets1) != 0){                                         # combinations of adjacent to x or
               if (length(condsets2 != 0)) condsets1 = cbind(condsets1,condsets2)} else
             if (length(condsets2) != 0) condsets1 = condsets2
             for (idx in 1:dim(condsets1)[[2]]) {
                                                                                 # skip if subset already looked at
               if ((idx > xlength) && (length(setdiff(condsets1[,idx],tf1)) == 0)) next
               if (isIndependent(data, x, z, condsets1[,idx], prgt = prgt)) {
                 if (length(which(condsets1[,idx] == y)) == 0){                  # y is not in independence set
                     all = FALSE
                     if (!(all || none)) {
                       ambiguous = cbind(ambiguous,c(x,y,z))
                       if (DEBUG) {cat(nm[x], " ",nm[y]," ", nm[z],"ambiguous", "\n")}
                       break}
                   } else {                                                      # y is in independence set
                          none = FALSE
                          if (!(all || none)) {
                          ambiguous = cbind(ambiguous,c(x,y,z))
                          if (DEBUG) {cat(nm[x], " ",nm[y]," ", nm[z],"ambiguous", "\n")}
                          break}}}
               if (!(all || none)) break}
             if (!(all || none)) break}
           if (all && none) {
              ambiguous = cbind(ambiguous,c(x,y,z))
              if (DEBUG) {cat(nm[x], " ",nm[y]," ", nm[z],"ambiguous", "\n")}}  # there is no edge, but not independent conditional on adjacent
           if (none && (!all)){
             if ((res[x,y] != EMP) && (res[z,y] != EMP)){
               res[x,y] = res[z,y] = ARR
               res[y,z] = res[y,x] = EMP
               if (DEBUG) {cat(nm[x], "*->", nm[y], "<-*", nm[z],"by Colliders","\n")}} else {
                 cat(nm[x], "*->", nm[y], "<-*", nm[z],"already oriented","\n")}

           }}}}}
return(list(res,ambiguous))}

orient_colliders = function(res, sepset) {
  n = nrow(res)
  nn = 1:n
  nm = dimnames(res)[[1]]
  for (x in 1:n) {
    for (y in nn[res[x,] != 0]) {
      for (z in nn[res[y,] != 0]) {
        if (x != z && res[x,z] == 0) {
          if (!is.element(y, sepset[[x]][[z]])) {
            if ((res[x,y] != EMP) && (res[z,y] != EMP)){
              res[z,y] = res[x,y] = ARR
              res[y,x] = res[y,z] = EMP
              if (DEBUG) cat(nm[x], "*->", nm[y], "<-*", nm[z],"by Colliders with sepset = {",nm[sepset[[x]][[z]]], "}\n")} else {
                if (DEBUG) cat(nm[x]," ", nm[y], " ", nm[z], "already oriented with sepset = {", nm[sepset[[x]][[z]]], "}\n")}}
          }}}}
  res
}

pcinternal = function(data, prgt = .05, noambiguous) {
  #added pc_skeleton by antti
  l = pc_skeleton(data, prgt = prgt)
  adjacencies = l[[1]]
  nm = dimnames(adjacencies)[[1]]
  sepset = l[[2]]
  n = length(data@variable_names)
  nn = 1:n

  res = adjacencies
## res[res == 0] = NO_EDGE
  res[res == 1] = CIR
  ambiguous = NULL

  if (noambiguous) res = orient_colliders(res,sepset) else {
    results = cpc_colliders(data=data,res=res,prgt=prgt,ambiguous)
    res = results[[1]]
    ambiguous = results[[2]]}
  res[t(res) == ARR] = EMP
  if (DEBUG) printadj(res)

  MAKING_PROGRESS = T
  while(MAKING_PROGRESS) {
    MAKING_PROGRESS = F
    for (a in 1:n) {
      for (b in nn[res[a,] == ARR]) {
        for (c in nn[res[,b] == CIR]) {
          if (a != c && adjacencies[a,c] == 0 && !(inambiguous(c(a,b,c),ambiguous))) {
              res[b,c] = ARR
              res[c,b] = EMP
              MAKING_PROGRESS = T
              if (DEBUG) cat(nm[b], "->", nm[c], "by Away From Collider since ",
                           nm[a], "->", nm[b], "o-o", nm[c], "\n")
            }}}}}

    for (a in 1:n) {
      for (c in nn[res[a,] == ARR]) {
        for (b in nn[res[c,] == ARR]) {
          if (res[a,b] == CIR) {
            res[a,b] = ARR
            res[b,a] = EMP
            MAKING_PROGRESS = T
            if (DEBUG) cat(nm[a], "->", nm[b], "by Away From Cycle since ", nm[a], "->", nm[c], "->", nm[b], "\n")
          }}}}

    for (a in 1:n) {
      for (b in nn[res[a,] == ARR]) {
        for (c in nn[res[b,] == EMP]) {
          for (d in nn[res[,b] == CIR]) {
            if ((a != c) && (res[c,d] != 0) && (res[a,d] != 0)) {
              if ((adjacencies[a,c] == 0) && !(inambiguous(c(a,d,c),ambiguous)) || (res[a,d] == EMP) || (res[c,d] == EMP)){
                res[d,b] = ARR
                res[b,d] = EMP
                MAKING_PROGRESS = T
                if (DEBUG) cat(nm[d], "->", nm[b], "by Double Triangle ",
                           nm[a], "->", nm[b], "<-", nm[c], "  ",nm[a], "o-o",nm[c], "*-*",nm[d], "o-o", nm[b], "\n")}}
            }}}}
i = 1
#  while (i <= length(ambiguous)){
#    if ((res[ambiguous[[i]][[1]],ambiguous[[i]][[2]]] == EMP) ||
#        (res[ambiguous[[i]][[3]],ambiguous[[i]][[2]]] == EMP) ||
#        ((res[ambiguous[[i]][[1]],ambiguous[[i]][[2]]] == ARR) && (res[ambiguous[[i]][[3]],ambiguous[[i]][[2]]] == ARR))){
#          browser()
#          ambiguous[[i]] = NULL} else i = i + 1}
  if (DEBUG) {
    printadj(res)
    printlist(ambiguous,data,res)}
  if (noambiguous) return(res) else return(list(res,ambiguous))
}


pc = function(data, prgt = 0.05) {pcinternal(data,prgt,noambiguous = TRUE)}



