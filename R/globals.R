## The rgeneric engine reads these from its definition environment (populated by
## inla.rgeneric.define); they are not true globals. Declared to quiet R CMD check.
utils::globalVariables(c("Dfunc", "Vinv", "rankdef", "graph.pattern",
                         "theta.initial", "logprior", "lognc"))
