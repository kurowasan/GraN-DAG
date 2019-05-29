source("../startups/startupSID.R")
p <- 4; G1 <- randomDAG(p,0.5); G2 <- randomDAG(p,1)
a <- structIntervDist2old(G1,dag2cpdagAdj(G2),FALSE);show(c(a$sid,a$sidLowerBound,a$sidUpperBound))
a <- structIntervDist(G1,dag2cpdagAdj(G2),FALSE);show(c(a$sid,a$sidLowerBound,a$sidUpperBound))
