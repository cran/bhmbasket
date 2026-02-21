## 1. Helper scenario_list and analysis_list for decision tests
## ------------------------------------------------------------------

helper_n_subjects <- c(10, 10)
helper_rr         <- c(0.5, 0.5)

helper_scenarios <- simulateScenarios(
  n_subjects_list     = list(helper_n_subjects),
  response_rates_list = list(helper_rr),
  scenario_numbers    = 1L,
  n_trials            = 5L
)

helper_analyses <- performAnalyses(
  scenario_list     = helper_scenarios,
  target_rates      = c(0.5, 0.5),
  method_names      = "berry",
  n_mcmc_iterations = 100,
  verbose           = FALSE
)

decisions <- getGoDecisions(
  analyses_list   = helper_analyses,
  cohort_names    = c("p_1", "p_2"),
  evidence_levels = c(0.5, 0.5),
  boundary_rules  = quote(c(TRUE, TRUE)),
  overall_min_gos = 1L
)

decisions_list <- decisions


# Tests for simulateScenarios --------------------------------------------------

# ------------------------------------------------------------------
# Test: basic structure and handling of historical cohorts
# Input:
#   - n_subjects = c(10,20,30), multiple response rate vectors (rr_1..rr_4),
#     n_trials = 25 for 4 scenarios.
# Behaviour:
#   - simulateScenarios() returns a scenario_list with expected matrices,
#     names, dimensions, and fixed responders for historical cohorts
#     (rr <= 0 or >= 1).
# Expectations:
#   - Class "scenario_list" with names scenario_1..scenario_4.
#   - n_subjects, response_rates, n_responders are matrices.
#   - Column names are n_*, r_*, rr_*.
#   - For rr_3 historical columns, n_responders are constant across trials.
#   - rr_5 outside [0,1] triggers an error.
# Why:
#   - Verifies core construction and the special handling of historical arms.
# ------------------------------------------------------------------
test_that("simulateScenarios has correct structure", {
  n_subjects <- c(10, 20, 30)
  rr_1 <- c(0.1, 0.1, 0.1)
  rr_2 <- c(0.9, 0.9, 0.9)
  rr_3 <- c(0.9, 0.9, 1)
  rr_4 <- c(0.9, 0.9, 0)
  rr_5 <- c(0.9, 0.9, 1.1)
  
  scenarios <- simulateScenarios(
    n_subjects_list     = list(n_subjects, n_subjects, n_subjects, n_subjects),
    response_rates_list = list(rr_1, rr_2, rr_3, rr_4),
    scenario_numbers    = c(1, 2, 3, 4),
    n_trials            = 25
  )
  
  expect_s3_class(scenarios, "scenario_list")
  expect_equal(names(scenarios), paste0("scenario_", 1:4))
  
  s1 <- scenarios[[1]]
  s3 <- scenarios[[3]]
  
  expect_true(is.matrix(s1$n_subjects))
  expect_true(is.matrix(s1$response_rates))
  expect_equal(s1$n_trials, 25)
  
  hist_cols <- which(rr_3 <= 0 | rr_3 >= 1)
  
  expect_true(
    all(apply(
      s3$n_responders[, hist_cols, drop = FALSE],
      2,
      function(x) length(unique(x)) == 1
    ))
  )
  
  expect_equal(ncol(s1$n_subjects), length(n_subjects))
  expect_equal(ncol(s1$response_rates), length(n_subjects))
  
  expect_equal(
    colnames(s1$n_subjects),
    paste0("n_", seq_len(length(n_subjects)))
  )
  
  expect_equal(
    colnames(s1$n_responders),
    paste0("r_", seq_len(length(n_subjects)))
  )
  
  expect_equal(
    colnames(s1$response_rates),
    paste0("rr_", seq_len(length(n_subjects)))
  )
  
  expect_error(
    simulateScenarios(
      n_subjects_list     = list(n_subjects),
      response_rates_list = list(rr_5),
      scenario_numbers    = 5,
      n_trials            = 5
    )
  )
})


# ------------------------------------------------------------------
# Test: argument validation, list consistency, and default n_trials
# Input:
#   - Various invalid n_trials, mismatched list lengths, cohort lengths,
#     missing arguments, and a global n_trials in .GlobalEnv.
# Behaviour:
#   - simulateScenarios() should:
#     • error on missing required args,
#     • error on invalid n_trials values,
#     • error on length mismatches / cohort mismatches,
#     • accept vector n_subjects_list (recycled to list),
#     • use .GlobalEnv$n_trials when argument is omitted.
# Expectations:
#   - Correct errors mentioning relevant argument names.
#   - Positive cases return scenario_list.
#   - When n_trials omitted and global n_trials = 7L, output n_trials = 7L.
# Why:
#   - Ensures robust front-end validation and backward-compatible defaults.
# ------------------------------------------------------------------
test_that("simulateScenarios has no mismatch and gives appropriate error messages", {
  n_subjects <- list(c(10, 20, 30))
  rr_1       <- c(0.1, 0.1, 0.1)
  rr_ok      <- c(0.2, 0.3, 0.4)
  
  expect_true(checkmate::test_int(5L, lower = 1))
  
  expect_error(
    simulateScenarios(n_subjects_list = n_subjects),
    "response_rates_list"
  )
  
  expect_error(
    simulateScenarios(response_rates_list = list(rr_1)),
    "n_subjects_list"
  )
  
  expect_s3_class(
    simulateScenarios(
      n_subjects_list     = c(10, 20, 30),
      response_rates_list = list(rr_1, rr_1),
      n_trials            = 2
    ),
    "scenario_list"
  )
  
  invalid_trials <- list(
    zero     = 0L,
    negative = -1L,
    nonint   = 2.5,
    na       = NA_integer_,
    nan      = NaN,
    inf      = Inf,
    ninf     = -Inf,
    length2  = c(2L, 3L)
  )
  
  for (case in invalid_trials) {
    expect_error(
      simulateScenarios(
        n_subjects_list     = list(c(10, 20, 30)),
        response_rates_list = list(rr_1),
        n_trials            = case
      ),
      "n_trials",
      info = paste("invalid n_trials case:", paste(case, collapse = ","))
    )
  }
  
  expect_error(
    simulateScenarios(
      n_subjects_list     = list(c(10, 20, 30)),
      response_rates_list = list(rr_1, rr_1)
    ),
    "same length"
  )
  
  expect_error(
    simulateScenarios(
      n_subjects_list     = list(c(10, 20, 30), c(10, 20, 30)),
      response_rates_list = list(c(0.1, 0.2, 0.3), c(0.2, 0.3, 0.4)),
      scenario_numbers    = c(1)
    ),
    "scenario_numbers"
  )
  
  expect_error(
    simulateScenarios(
      n_subjects_list     = list(10),
      response_rates_list = list(0.1),
      n_trials            = 2
    ),
    "at least 2 cohorts"
  )
  
  expect_error(
    simulateScenarios(
      n_subjects_list     = list(c(10, 20), c(10, 20)),
      response_rates_list = list(c(0.1, 0.2), 0.3),
      n_trials            = 2
    ),
    "same number of cohorts"
  )
  
  expect_error(
    simulateScenarios(
      n_subjects_list     = list(c(10, 20), 10),
      response_rates_list = list(c(0.1, 0.2), c(0.1, 0.3)),
      n_trials            = 2
    ),
    "same number of cohorts"
  )
  
  old <- if (exists("n_trials", .GlobalEnv)) get("n_trials", .GlobalEnv) else NULL
  
  on.exit({
    if (is.null(old)) {
      rm(list = "n_trials", envir = .GlobalEnv)
    } else {
      assign("n_trials", old, .GlobalEnv)
    }
  }, add = TRUE)
  
  assign("n_trials", 7L, envir = .GlobalEnv)
  
  out <- simulateScenarios(
    n_subjects_list     = list(c(10, 20, 30)),
    response_rates_list = list(rr_ok)
  )
  
  expect_equal(out[[1]]$n_trials, 7L)
})

# Tests for saveScenarios ------------------------------------------------------

# ------------------------------------------------------------------
# Test: saveScenarios writes RDS files and returns metadata
# Input:
#   - scenario_list with 2 simple mock scenarios (scenario_number 1 and 2),
#     save_path = temp directory.
# Behaviour:
#   - For each scenario, saveScenarios() writes a scenario_data_<num>.rds file
#     and returns a list with the scenario_numbers and the path used.
# Expectations:
#   - Files scenario_data_1.rds and scenario_data_2.rds exist.
#   - Returned scenario_numbers == c(1,2) and path == temp_dir.
# Why:
#   - Checks that serialization of scenarios and meta-return are correct.
# ------------------------------------------------------------------
test_that("saveScenarios saves files correctly", {
  mock_scenario <- function(num) {
    structure(
      list(scenario_number = num, data = paste("Scenario", num)),
      class = "scenario"
    )
  }
  scenario_list <- list(mock_scenario(1), mock_scenario(2))
  class(scenario_list) <- "scenario_list"
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  result <- saveScenarios(scenario_list, save_path = temp_dir)
  
  expect_true(file.exists(file.path(temp_dir, "scenario_data_1.rds")))
  expect_true(file.exists(file.path(temp_dir, "scenario_data_2.rds")))
  
  expect_equal(result$scenario_numbers, c(1, 2))
  expect_equal(result$path, temp_dir)
  
  unlink(temp_dir, recursive = TRUE)
})


# ------------------------------------------------------------------
# Test: saveScenarios creates directory if it does not exist
# Input:
#   - scenario_list with one scenario, save_path = non-existent temp path.
# Behaviour:
#   - Function must create the directory, save the file and return path.
# Expectations:
#   - dir.exists(save_path) becomes TRUE.
#   - Returned path equals the requested save_path.
# Why:
#   - Ensures convenience behaviour when user passes a new directory.
# ------------------------------------------------------------------
test_that("saveScenarios creates directory when save_path does not exist", {
  temp_dir <- tempfile()
  expect_false(dir.exists(temp_dir))
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
  
  sl <- structure(list(list(scenario_number = 1)), class = "scenario_list")
  
  res <- saveScenarios(sl, save_path = temp_dir)
  
  expect_true(dir.exists(temp_dir))
  expect_equal(res$path, temp_dir)
})

# Tests for print.scenario_list ------------------------------------------------

# ------------------------------------------------------------------
# Test: print.scenario_list header shows scenario and cohort count
# Input:
#   - scenario_list with 2 scenarios, each with 2 cohorts.
# Behaviour:
#   - print.scenario_list() prints a header with counts plus scenario names.
# Expectations:
#   - Output contains “scenario_list of 2 scenarios with 2 cohorts”.
#   - Output mentions scenario_1 and scenario_2.
# Why:
#   - Verifies the informative summary header and scenario labels.
# ------------------------------------------------------------------
test_that("print.scenario_list header is correct and shows cohort count", {
  n_subjects <- c(5, 6)
  rr1 <- c(0.1, 0.2)
  rr2 <- c(0.3, 0.4)
  
  scenarios <- simulateScenarios(
    n_subjects_list     = list(n_subjects, n_subjects),
    response_rates_list = list(rr1, rr2),
    scenario_numbers    = c(1, 2),
    n_trials            = 3
  )
  
  expect_output(print(scenarios), "scenario_list of 2 scenarios with 2 cohorts")
  expect_output(print(scenarios), "  - scenario_1", fixed = FALSE)
  expect_output(print(scenarios), "  - scenario_2", fixed = FALSE)
})


# ------------------------------------------------------------------
# Test: print.scenario_list prints table section labels
# Input:
#   - scenario_list with 1 scenario and 3 cohorts.
# Behaviour:
#   - print.scenario_list() prints labeled rows for “true response rates”
#     and “average number of subjects”.
# Expectations:
#   - Output contains those labels.
# Why:
#   - Confirms that the summary table describes what is being shown.
# ------------------------------------------------------------------
test_that("print.scenario_list prints table row labels", {
  n_subjects <- c(7, 8, 9)
  rr <- c(0.2, 0.3, 0.4)
  
  scenarios <- simulateScenarios(
    n_subjects_list     = list(n_subjects),
    response_rates_list = list(rr),
    scenario_numbers    = 1,
    n_trials            = 2
  )
  
  expect_output(print(scenarios), "true response rates:")
  expect_output(print(scenarios), "average number of subjects:")
})


# ------------------------------------------------------------------
# Test: print.scenario_list uses cohort column names c_1, c_2, ...
# Input:
#   - scenario_list with 1 scenario and 3 cohorts.
# Behaviour:
#   - Columns in the printed tables are labeled c_1, c_2, c_3.
# Expectations:
#   - Output contains literal c_1, c_2, c_3 as word boundaries.
# Why:
#   - Ensures cohort naming scheme is user-friendly and consistent.
# ------------------------------------------------------------------
test_that("print.scenario_list uses cohort column names c_1, c_2, ...", {
  n_subjects <- c(7, 8, 9)
  rr <- c(0.2, 0.3, 0.4)
  
  scenarios <- simulateScenarios(
    n_subjects_list     = list(n_subjects),
    response_rates_list = list(rr),
    scenario_numbers    = 1,
    n_trials            = 2
  )
  
  expect_output(print(scenarios), "\\bc_1\\b")
  expect_output(print(scenarios), "\\bc_2\\b")
  expect_output(print(scenarios), "\\bc_3\\b")
})


# ------------------------------------------------------------------
# Test: print.scenario_list footer shows realizations info
# Input:
#   - scenario_list with n_trials = 5.
# Behaviour:
#   - Footer prints number of realizations per scenario and overall.
# Expectations:
#   - Output mentions “5 trial realizations per scenario”.
#   - Output mentions “unique trial realizations overall”.
# Why:
#   - Provides context for how many simulations underpin the summaries.
# ------------------------------------------------------------------
test_that("print.scenario_list prints footer lines with realizations info", {
  n_subjects <- c(10, 10)
  rr <- c(0.1, 0.2)
  
  scenarios <- simulateScenarios(
    n_subjects_list     = list(n_subjects),
    response_rates_list = list(rr),
    scenario_numbers    = 1,
    n_trials            = 5
  )
  
  expect_output(print(scenarios), "5 trial realizations per scenario")
  expect_output(print(scenarios), "unique trial realizations overall")
})

# Tests for loadScenarios ------------------------------------------------------

# ------------------------------------------------------------------
# Test: loadScenarios round-trip with saveScenarios
# Input:
#   - Two scenarios generated via simulateScenarios(), saved to temp_dir.
# Behaviour:
#   - loadScenarios() reads back the RDS files into a scenario_list.
# Expectations:
#   - Returned object has class "scenario_list".
#   - Loaded object equals the original scenario_list.
# Why:
#   - Validates persistence and recovery of scenario definitions.
# ------------------------------------------------------------------
test_that("loadScenarios loads saved scenarios and returns proper class/names", {
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  scenario_list_to_save <- simulateScenarios(
    n_subjects_list     = list(c(10, 20, 30), c(10, 20, 30)),
    response_rates_list = list(c(0.1, 0.2, 0.3), c(0.5, 0.6, 0.7)),
    scenario_numbers    = c(1, 2),
    n_trials            = 5
  )
  
  saveScenarios(scenario_list_to_save, save_path = temp_dir)
  
  loaded <- loadScenarios(c(1, 2), load_path = temp_dir)
  
  expect_s3_class(loaded, "scenario_list")
  expect_equal(loaded, scenario_list_to_save)
  
  unlink(temp_dir, recursive = TRUE, force = TRUE)
})


# ------------------------------------------------------------------
# Test: loadScenarios validates load_path argument
# Input:
#   - Incorrect load_path values: vector, numeric, NULL, NA_character_.
# Behaviour:
#   - loadScenarios() should error when load_path is not a single valid string.
# Expectations:
#   - Each invalid call errors with message mentioning 'load_path'.
# Why:
#   - Ensures the file location argument is properly validated.
# ------------------------------------------------------------------
test_that("loadScenarios points out invalid paths", {
  expect_error(loadScenarios(1L, load_path = c(tempdir(), tempdir())), "load_path")
  expect_error(loadScenarios(1L, load_path = 123), "load_path")
  expect_error(loadScenarios(1L, load_path = NULL), "load_path")
  expect_error(loadScenarios(1L, load_path = NA_character_), "load_path")
  # Optional: empty string behaviour may be handled elsewhere
  # expect_error(loadScenarios(1L, load_path = ""), "load_path")
})


# ------------------------------------------------------------------
# Test: loadScenarios errors when requested files are missing
# Input:
#   - Valid directory without any scenario_data_*.rds files,
#     scenario_numbers = c(999).
# Behaviour:
#   - loadScenarios() attempts to read missing files and should error.
# Expectations:
#   - Error is thrown.
# Why:
#   - Confirms that missing files are not silently ignored.
# ------------------------------------------------------------------
test_that("loadScenarios errors when files do not exist", {
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  expect_error(
    suppressWarnings(loadScenarios(c(999), load_path = temp_dir))
  )
  
  unlink(temp_dir, recursive = TRUE, force = TRUE)
})

# Tests for is.scenario_list ---------------------------------------------------

# ------------------------------------------------------------------
# Test: is.scenario_list TRUE for valid scenario_list and FALSE otherwise
# Input:
#   - Valid scenario_list (proper class and scenario elements),
#     plus several invalid objects.
# Behaviour:
#   - is.scenario_list() returns TRUE only for the valid scenario_list.
# Expectations:
#   - TRUE for valid_list.
#   - FALSE for list(), character, NULL.
# Why:
#   - Ensures the helper detects correct class/structure.
# ------------------------------------------------------------------
test_that("is.scenario_list returns TRUE for valid input and FALSE otherwise", {
  mock_scenario <- function(num) {
    structure(list(scenario_number = num), class = "scenario")
  }
  valid_list <- list(mock_scenario(1), mock_scenario(2))
  class(valid_list) <- "scenario_list"
  
  expect_true(is.scenario_list(valid_list))
  
  expect_false(is.scenario_list(list()))
  expect_false(is.scenario_list("not a scenario list"))
  expect_false(is.scenario_list(NULL))
})


# ------------------------------------------------------------------
# Test: is.scenario_list errors when x is missing
# Input:
#   - Call is.scenario_list() with no argument.
# Behaviour:
#   - Function should throw an informative error about missing x.
# Expectations:
#   - Error message mentions "argument 'x'" (or similar custom message).
# Why:
#   - Enforces explicit use and avoids accidental calls without input.
# ------------------------------------------------------------------
test_that("is.scenario_list errors when argument x is missing", {
  expect_error(
    is.scenario_list(),
    "Please provide an object for the argument 'x'"
  )
})

# Tests for createTrial --------------------------------------------------------

# ------------------------------------------------------------------
# Test: createTrial returns a trial object for valid input
# Input:
#   - n_subjects = c(10,20), n_responders = c(5,15), with simulateScenarios
#     mocked to avoid heavy computation.
# Behaviour:
#   - createTrial() should call simulateScenarios() internally and
#     return a non-NULL trial object.
# Expectations:
#   - Result is not NULL (shape-specific checks may be done elsewhere).
# Why:
#   - Confirms the wrapper-path from trial inputs to a scenario object.
# ------------------------------------------------------------------
test_that("createTrial returns trial object with valid input", {
  with_mocked_bindings({
    result <- createTrial(n_subjects = c(10, 20), n_responders = c(5, 15))
    expect_true(!is.null(result))
  }, simulateScenarios = function(n_subjects_list, response_rates_list, n_trials) {
    list(trial = "mocked_trial")
  })
})


# ------------------------------------------------------------------
# Test: createTrial errors on invalid n_subjects or n_responders
# Input:
#   - n_subjects containing NA, or n_responders being non-integer.
# Behaviour:
#   - Argument checks should fail and throw informative errors.
# Expectations:
#   - Error mentioning 'n_subjects' for NA in n_subjects.
#   - Error mentioning 'n_responders' for non-integer responders.
# Why:
#   - Ensures basic input validation for trial setup.
# ------------------------------------------------------------------
test_that("createTrial errors with invalid input", {
  expect_error(
    createTrial(n_subjects = c(10, NA), n_responders = c(5, 15)),
    "n_subjects"
  )
  
  expect_error(
    createTrial(n_subjects = c(10, 20), n_responders = c(5.5, 15)),
    "n_responders"
  )
})

# Tests for continueRecruitment ------------------------------------------------

# ------------------------------------------------------------------
# Test: continueRecruitment errors if n_subjects_add_list is missing
# Input:
#   - decisions_list provided, n_subjects_add_list omitted.
# Behaviour:
#   - Function must complain about missing n_subjects_add_list.
# Expectations:
#   - Error message mentions 'n_subjects_add_list'.
# Why:
#   - Additional recruitment per scenario must be specified explicitly.
# ------------------------------------------------------------------
test_that("Error is thrown if n_subjects_add_list is missing", {
  expect_error(
    continueRecruitment(
      decisions_list = list(scenario_1 = list(decisions_list = list()))
    ),
    "n_subjects_add_list"
  )
})


# ------------------------------------------------------------------
# Test: continueRecruitment errors if decisions_list is missing
# Input:
#   - n_subjects_add_list provided, decisions_list omitted.
# Behaviour:
#   - Function must complain about missing decisions_list.
# Expectations:
#   - Error message mentions 'decisions_list'.
# Why:
#   - Decisions drive which cohorts continue recruiting.
# ------------------------------------------------------------------
test_that("Error is thrown if decisions_list is missing", {
  expect_error(
    continueRecruitment(
      n_subjects_add_list = list(10)
    ),
    "decisions_list"
  )
})


# ------------------------------------------------------------------
# Test: continueRecruitment errors when decisions_list has wrong class
# Input:
#   - decisions_list is a plain list, not class 'decision_list'.
# Behaviour:
#   - Class assertion should fail.
# Expectations:
#   - Error mentioning 'decisions_list'.
# Why:
#   - Ensures only properly constructed decision_list objects are used.
# ------------------------------------------------------------------
test_that("Error is thrown if decisions_list is not of class 'decision_list'", {
  expect_error(
    continueRecruitment(
      n_subjects_add_list = list(10),
      decisions_list      = list()
    ),
    "decisions_list"
  )
})


# ------------------------------------------------------------------
# Test: continueRecruitment errors for invalid method_name
# Input:
#   - method_name = "invalid_method" not present in decisions_list.
# Behaviour:
#   - Function should reject unknown method names.
# Expectations:
#   - Error is thrown.
# Why:
#   - Prevents referencing non-existent analysis methods.
# ------------------------------------------------------------------
test_that("Error for invalid method_name", {
  expect_error(
    continueRecruitment(
      n_subjects_add_list = list(c(10, 10)),
      decisions_list      = decisions,
      method_name         = "invalid_method"
    )
  )
})


# ------------------------------------------------------------------
# Test: continueRecruitment errors when method_name is NULL and multiple methods
# Input:
#   - decisions modified to contain two methods ("berry", "exnex"),
#     method_name not supplied.
# Behaviour:
#   - Function should require method_name when >1 method is available.
# Expectations:
#   - Error is thrown.
# Why:
#   - Disambiguates which method’s decisions to follow.
# ------------------------------------------------------------------
test_that("Error when method_name is NULL and multiple methods are present", {
  decisions$scenario_1$analysis_data$analysis_parameters$method_names <- c("berry", "exnex")
  decisions$scenario_1$decisions_list <- list(berry = data.frame(), exnex = data.frame())
  
  expect_error(
    continueRecruitment(
      n_subjects_add_list = list(c(10, 10)),
      decisions_list      = decisions
    )
  )
})


# ------------------------------------------------------------------
# Test: continueRecruitment errors if list lengths do not match
# Input:
#   - n_subjects_add_list length 2 vs decisions_list length 1.
# Behaviour:
#   - Length consistency check should fail.
# Expectations:
#   - Error is thrown.
# Why:
#   - Each scenario’s decisions need corresponding n_subjects_add vector.
# ------------------------------------------------------------------
test_that("Error is thrown if n_subjects_add_list and decisions_list lengths do not match", {
  expect_error(
    continueRecruitment(
      n_subjects_add_list = list(10, 20),
      decisions_list      = decisions_list,
      method_name         = "berry"
    )
  )
})


# ------------------------------------------------------------------
# Test: continueRecruitment errors when selected method not analyzed
# Input:
#   - decisions scenario_1 only analysed with 'berry',
#     requested method_name = 'exnex'.
# Behaviour:
#   - Function should detect that exnex is not in method_names.
# Expectations:
#   - Error containing "Selected method_name not analyzed".
# Why:
#   - Prevents accidentally requesting an unsupported method.
# ------------------------------------------------------------------
test_that("Error when selected method_name not analyzed in scenario", {
  decisions$scenario_1$analysis_data$analysis_parameters$method_names <- "berry"
  decisions$scenario_1$decisions_list <- list(berry = data.frame())
  
  expect_error(
    continueRecruitment(
      n_subjects_add_list = list(c(10, 10)),
      decisions_list      = decisions,
      method_name         = "exnex"
    ),
    "Selected method_name not analyzed"
  )
})


# ------------------------------------------------------------------
# Test: continueRecruitment errors when all cohorts are historical
# Input:
#   - scenario_1 response_rates <= 0 or >= 1, making all cohorts historical.
# Behaviour:
#   - continueRecruitment() must error with "Only historical cohorts in scenario".
# Expectations:
#   - Error message contains that text.
# Why:
#   - There is nothing to recruit if all cohorts are purely historical.
# ------------------------------------------------------------------
test_that("Error when response rate <=0 or >=1", {
  decisions$scenario_1$analysis_data$analysis_parameters$method_names <- "berry"
  decisions$scenario_1$decisions_list <- list()
  decisions$scenario_1$scenario_data$response_rates <-
    matrix(c(0, 1), nrow = 1, dimnames = list(NULL, c("rr_hist1", "rr_hist2")))
  decisions$scenario_1$scenario_data$n_trials <- 1
  decisions$scenario_1$scenario_data$n_responders <-
    matrix(c(0, 1), nrow = 1, dimnames = list(NULL, c("rr_hist1", "rr_hist2")))
  decisions$scenario_1$scenario_data$n_subjects <-
    matrix(c(0, 1), nrow = 1, dimnames = list(NULL, c("rr_hist1", "rr_hist2")))
  decisions$scenario_1$scenario_data$scenario_number <- 1
  
  expect_error(
    continueRecruitment(
      n_subjects_add_list = list(c(10, 10)),
      decisions_list      = decisions,
      method_name         = "berry"
    ),
    "Only historical cohorts in scenario 1"
  )
})


# ------------------------------------------------------------------
# Test: continueRecruitment errors when n_subjects_add length mismatches
# Input:
#   - Two recruiting cohorts and one historical, but n_subjects_add_list[[1]]
#     has length 1 instead of 2.
# Behaviour:
#   - Function should detect mismatch between recruiting cohorts and
#     n_subjects_add length.
# Expectations:
#   - Error is thrown.
# Why:
#   - Ensures user supplies one recruitment size per recruiting cohort.
# ------------------------------------------------------------------
test_that("Error when n_subjects_add length does not match recruiting cohorts", {
  decisions$scenario_1$analysis_data$analysis_parameters$method_names <- "berry"
  decisions$scenario_1$decisions_list <- list(berry = data.frame())
  decisions$scenario_1$scenario_data$response_rates <-
    matrix(c(0.3, 0.4, 1), nrow = 1,
           dimnames = list(NULL, c("rr_A", "rr_B", "rr_hist")))
  
  expect_error(
    continueRecruitment(
      n_subjects_add_list = list(c(10)),
      decisions_list      = decisions,
      method_name         = "berry"
    )
  )
})


# ------------------------------------------------------------------
# Test: continueRecruitment updates only rows where overall == TRUE
# Input:
#   - One recruiting cohort rr_A (0.5) and one historical rr_hist (1),
#     baseline n_subjects/n_responders, decisions overall=TRUE,
#     mocked getScenario that returns updated counts for rr_A only.
# Behaviour:
#   - continueRecruitment() should:
#     • add new subjects/responders for rr_A,
#     • leave rr_hist unchanged.
# Expectations:
#   - New n_subjects["rr_A"] == 10, n_responders["rr_A"] == 5.
#   - n_subjects["rr_hist"] and n_responders["rr_hist"] unchanged at 1.
# Why:
#   - Verifies that recruitment only happens in Go rows and only for
#     non-historical cohorts.
# ------------------------------------------------------------------
test_that("Only overall==TRUE rows are updated", {
  decisions$scenario_1$scenario_data$response_rates <-
    matrix(c(0.5, 1), nrow = 1,
           dimnames = list(NULL, c("rr_A", "rr_hist")))
  
  decisions$scenario_1$scenario_data$n_subjects <-
    matrix(c(0, 1), nrow = 1,
           dimnames = list(NULL, c("rr_A", "rr_hist")))
  
  decisions$scenario_1$scenario_data$n_responders <-
    matrix(c(0, 1), nrow = 1,
           dimnames = list(NULL, c("rr_A", "rr_hist")))
  
  decisions$scenario_1$decisions_list <-
    list(berry = data.frame(decision_1 = TRUE, decision_2 = TRUE, overall = TRUE))
  decisions$scenario_1$analysis_data$analysis_parameters$method_names <- "berry"
  
  testthat::with_mocked_bindings(
    getScenario = function(n_subjects, response_rates, cohort_names, n_trials) {
      list(
        n_subjects   = matrix(10, nrow = 1, dimnames = list(NULL, cohort_names)),
        n_responders = matrix(5, nrow = 1, dimnames = list(NULL, cohort_names))
      )
    }, {
      out <- continueRecruitment(
        n_subjects_add_list = list(c(10)),
        decisions_list      = decisions,
        method_name         = "berry"
      )
      
      expect_equal(unname(out$scenario_1$n_subjects[, "rr_A"]), 10)
      expect_equal(unname(out$scenario_1$n_responders[, "rr_A", drop = TRUE]), 5)
      expect_equal(unname(out$scenario_1$n_subjects[, "rr_hist"]), 1)
      expect_equal(unname(out$scenario_1$n_responders[, "rr_hist", drop = TRUE]), 1)
    }
  )
})