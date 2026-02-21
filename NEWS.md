# bhmbasket 1.0.0

### New and Altered Features

* Added structured unit tests for core functionality

* Integrated code coverage

* Replaced validation logic with checkmate assertions and removed unused helper functions

* Replaced dependency on R2jags to rjags

* Added Github Actions workflows for automated R CMD check and CI validation

* Minor changes in code

# bhmbasket 0.9.6

### Fixed Bugs

* Fixed a bug in getEstimates(), which used to calculate variances instead of the standard deviations of response rates estimated with the methods pooled and stratified

### New & Altered Features

* Added print methods for scenario_list, analysis_list, and decision_list objects

* Added function getAverageNSubjects()

* Switched from R2jags to rjags for performance increase

* Performance increases in performAnalyses(), getEstimates(), and getGoDecisions()

* getGoDecisions() saves the decisions rules used to derive the decisions

* In continueRecruitment(), the argument method_name can be omitted if only one method has been used in performAnalyses()

* Caution message in performAnalyses() if no parallel backend is detected disabled for methods stratified and pooled

* Changed message from performAnalyses()

* Disabled loading package messages when calling functions from other packages

* Updated Imports in DESCRIPTION

* Minor changes in code

# bhmbasket 0.9.5

### Fixed Bugs

* Fixed a bug that occured in performAnalyses() using R-devel due to a recent change in stats::aggregate()

### New & Altered Features

* Introduced nested parallelization for better usage of HPC resources

* Introduced chunking of tasks for better performance in parallel environments

* Updated vignette on HPC environment

* Update documentation of performAnalyses()

* Recommended doFuture and future over doParallel and parallel when no parallel backend is detected

* Updated Suggests and SystemRequirements in DESCRIPTION

* Added a WORDLIST

* Minor changes in code

# bhmbasket 0.9.4

### Fixed Bugs

* Fixed a bug in continueRecruitment() that could result in additional subjects to be recruited when the overall decision for a trial realization is NoGo but some cohorts of that trial realization have Go decisions.

* Fixed a bug in performAnalyses() that would occur if all trial realizations of a scenario had a previous overall NoGo decision and would result in performAnalyses() to return an empty list for that scenario's posterior quantiles.

* Specified R2jags package version requirement in DESCRIPTION to prevent 'unused argument' bug in performAnalyses()

* Fixed warning message not showing when specifying deprecated arguments 'seed' and 'n_cores' in performAnalyses() 

### New & Altered Features

* Usage of doRNG package for fully reproducible results in parallel execution

* Usage of hash tables for mapping unique trial realizations to scenario trial realizations resulting in performance improvement

* Added a vignette that provides a short example on how to use bhmbasket in a high performance computing environment

* JAGS model files stored in package instead of writing to temporary files

* Updated Imports in DESCRIPTION

* Updated CITATION

* Minor changes in code

* Updated documentation

* Removed superfluous files

# bhmbasket 0.9.3

### Fixed Bugs

* Fixed a bug in continueRecruitment() that could result in an error message and prevent the function from running although all conditions were met.

* Fixed a bug leading to an error messages in getEstimates(). It would occur if in a previous call of performAnalyses() differences between cohorts were calculated, but not in the current call of performAnalyses().

### New & Altered Features

* Registration of parallel backend is now the responsibility of the user to allow for flexibility. A respective message is displayed in performAnalyses() if no parallel backend is registered.

* Functions simulateScenarios() and performAnalyses() look up for their arguments 'n_trials' and 'n_mcmc_iterations', respectively, in the global environment if not provided by user.

* Arguments 'n_subjects_list' and 'n_subjects_add_list' of the functions simulateScenarios() and continueRecruitment(), respectively, can be provided as a single vector for the case when all scenarios recruit the same number of subjects.

* Arguments 'seed' and 'n_cores' were deprecated in performAnalyses(). The argument 'seed' was not used in the call of R2jags::jags(), see also the the documentation of 'jags.seed' in ?R2jags::jags. The argument 'n_cores' is no longer needed as the registration of the parallel backend is now the responsibility of the user.

* Added the argument 'overall_min_nogos' to the function negateGoDecisions().

* Added an input check in function getEstimates() regarding the argument 'add_parameters'.

* Changed date format in message from performAnalyses()

* Updated documentation

* Updated Imports and Suggests in DESCRIPTION

* Minor changes in code

# bhmbasket 0.9.2

* Fixed a rare bug that could result in wrong overall decision in getGoDecisions() in the presence of previous go decisions

* Updated documentation

* Updated Description in DESCRIPTION

* Minor changes in code

# bhmbasket 0.9.1

* Measures to fulfill CRAN policies and recommendations

# bhmbasket 0.9.0

* Initial release
