# survivalAnalysis 0.2.0

* add the methods pluck_survival_analysis and pluck_multivariate_analysis
  providing for the first clean API to access the "black-box" result objects
* add the multivariate_as_data_frame method to access the multivariate result object
* fix bug when cleaning column names in analyse_ functions

# survivalAnalysis 0.1.3

* fix correct extraction of p value from survdiff in case degree-of-freedom is >1 (thanks to Nolan A. Wages)

# survivalAnalysis 0.1.2

* fix crash with upcoming dplyr 1.0.0, and some rlang deprecations
* fix x scaling and breaking in KM plot

# survivalAnalysis 0.1.1

* export identity_order, which is used by forest plots
* work around upstream name collision (flatten_raw in purrr and rlang)

# survivalAnalysis 0.1.0

* initial release
