# Methods Update: Class Weighting and Data-Driven Panel Reduction

The primary binary modeling cohort contained 10 fast progressors and 22 slow progressors. To reduce bias toward the majority slow-progressor class, primary logistic-regression and random-forest models used class weighting. Unweighted versions of each algorithm were evaluated as sensitivity analyses using the same cross-validation splits and feature sets.

The revised model set includes balanced and unweighted variants for both logistic regression and random forest. Logistic-regression models use L2 regularization, median imputation, and standard scaling inside the modeling pipeline. Random-forest models use median imputation and constrained tree complexity.

Exploratory data-driven protein models were reduced to a single top-10 panel. For each outer training fold, proteins were ranked using only training-fold data based on a combined rank of Mann-Whitney p-value and absolute Cohen's d, after excluding proteins marked as age- or sex-associated in healthy controls. The top 10 ranked proteins from that training fold were then evaluated on the held-out fold. This data-driven panel is reported as an exploratory comparator rather than as the preferred manuscript panel.

No SMOTE or synthetic oversampling was used.
