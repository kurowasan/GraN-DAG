# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 


## To source() a bunch of files in a folder
sourceDir = function(path, trace = FALSE)
{
    for (nm in list.files(path, pattern = "\\.[R]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm))
        if(trace) cat("\n")
    }
}


library(clue)
sourceDir("./bayeslingam/rlingam/", trace = FALSE)
sourceDir("./bayeslingam/common/")
#sourceDir("./bayeslingam/ETUDAG/")

source("../code/inferDAG/lingamWrap.R")
