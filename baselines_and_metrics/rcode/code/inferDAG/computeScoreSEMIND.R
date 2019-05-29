computeScoreSEMIND <- function(State, pars, penFactor)
    # Copyright (c) 2013 Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    SigmaHat <- pars$SigmaHat
    n <- State$n
    p <- State$p
    
    
    score <- sum(State$SingleScores) + 0.7*sum(State$Adj)
    scoreStat <- sum(State$SingleScoresStat)
    # score <- - indtestMutual(State$Res, pars$indtest.method, pars$indtest.pars)$p.value
    # score <- indtestMutualHsic(State$Res, pars$indtest.method, pars$indtest.pars)$statistic

    return(list(score=score, scoreStat=scoreStat))
}
