# Tests for chunkVector --------------------------------------------------------

# ------------------------------------------------------------------
# Test: chunkVector returns a single chunk when n_chunks <= 1
# Input:
#   - x = 1:10
#   - n_chunks = 1
# Behaviour:
#   - returns list(x) (one chunk containing all elements)
# Expectations:
#   - output is a list of length 1
#   - the only element equals x exactly
# Why:
#   - ensures the "no chunking" path is stable and predictable.
# ------------------------------------------------------------------
test_that("chunkVector returns single chunk when n_chunks <= 1", {
  
  x <- 1:10
  
  out <- chunkVector(x, n_chunks = 1)
  
  expect_type(out, "list")
  
  expect_length(out, 1)
  
  expect_identical(out[[1]], x)
  
})

# ------------------------------------------------------------------
# Test: chunkVector splits into multiple chunks and preserves all elements
# Input:
#   - x = 1:10
#   - n_chunks = 2, 3
# Behaviour:
#   - splits x into n_chunks chunks
# Expectations:
#   - output is a list with requested length
#   - unlisted output contains exactly the original elements
# Why:
#   - verifies chunking doesn't drop or duplicate elements.
# ------------------------------------------------------------------
test_that("chunkVector splits into multiple chunks and preserves elements", {
  
  x <- 1:10
  
  out2 <- chunkVector(x, n_chunks = 2)
  
  expect_type(out2, "list")
  expect_length(out2, 2)
  expect_identical(sort(unlist(out2)), x)
  
  out3 <- chunkVector(x, n_chunks = 3)
  
  expect_type(out3, "list")
  expect_length(out3, 3)
  expect_identical(sort(unlist(out3)), x)
  
})


# Tests for convertVector2Matrix -----------------------------------------------

# ------------------------------------------------------------------
# Test: convertVector2Matrix converts a vector to a 1-row matrix
# Input:
#   - vector = c(1,2,3,4)
# Behaviour:
#   - wraps into a matrix and transposes to 1 x length(vector)
# Expectations:
#   - result is a matrix with dim (1,4)
#   - content equals the original vector
# Why:
#   - many downstream functions assume matrix input.
# ------------------------------------------------------------------
test_that("convertVector2Matrix converts vector to 1-row matrix", {
  
  v <- 1:4
  
  out <- convertVector2Matrix(v)
  
  expect_true(is.matrix(out))
  
  expect_equal(dim(out), c(1, 4))
  
  expect_identical(as.integer(out[1, ]), v)
  
})

# ------------------------------------------------------------------
# Test: convertVector2Matrix returns matrix unchanged
# Input:
#   - matrix input
# Behaviour:
#   - should not modify an existing matrix
# Expectations:
#   - identical to the original matrix
# Why:
#   - prevents accidental reshaping when caller already provides a matrix.
# ------------------------------------------------------------------
test_that("convertVector2Matrix leaves matrices unchanged", {
  
  m <- matrix(1:6, nrow = 2)
  
  out <- convertVector2Matrix(m)
  
  expect_true(is.matrix(out))
  
  expect_identical(out, m)
  
})


# Tests for createHashTable + getHashValues ------------------------------------

# ------------------------------------------------------------------
# Test: createHashTable stores key-value pairs retrievable via getHashValues
# Input:
#   - keys = c("a","b","c")
#   - values = list(10,20,30)
# Behaviour:
#   - creates hash and assigns values to keys
# Expectations:
#   - returned object is an environment
#   - getHashValues returns values in the requested key order
# Why:
#   - mapUniqueTrials relies on retrieval from hashed keys.
# ------------------------------------------------------------------
test_that("createHashTable and getHashValues store and retrieve values", {
  
  keys <- c("a", "b", "c")
  vals <- list(10, 20, 30)
  
  h <- createHashTable(keys, vals)
  
  expect_true(is.environment(h))
  
  out <- getHashValues(c("b", "a"), h)
  
  expect_type(out, "list")
  
  expect_equal(out, list(20, 10))
  
})


# Tests for getHashKeys ---------------------------------------------------------

# ------------------------------------------------------------------
# Test: getHashKeys uses default x1..xp names if colnames are missing
# Input:
#   - 2x2 matrix without colnames
# Behaviour:
#   - returns character keys "x1:<val>|x2:<val>" per row
# Expectations:
#   - correct key strings for each row
# Why:
#   - hashing must be stable even when matrices lack colnames.
# ------------------------------------------------------------------
test_that("getHashKeys uses default x* names without colnames", {
  
  x <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  
  out <- getHashKeys(x)
  
  expect_type(out, "character")
  expect_length(out, 2)
  
  expect_identical(out[1], "x1:1|x2:2")
  expect_identical(out[2], "x1:3|x2:4")
  
})

# ------------------------------------------------------------------
# Test: getHashKeys uses existing colnames when present
# Input:
#   - 2x2 matrix with colnames c("A","B")
# Behaviour:
#   - returns "A:<val>|B:<val>" per row
# Expectations:
#   - correct key strings for each row
# Why:
#   - ensures keys remain consistent when names exist.
# ------------------------------------------------------------------
test_that("getHashKeys uses colnames when present", {
  
  x <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  colnames(x) <- c("A", "B")
  
  out <- getHashKeys(x)
  
  expect_identical(out[1], "A:1|B:2")
  expect_identical(out[2], "A:3|B:4")
  
})


# Tests for cummulativeMovingAverage -------------------------------------------

# ------------------------------------------------------------------
# Test: cummulativeMovingAverage returns running mean
# Input:
#   - x = c(2,4,6,8)
# Behaviour:
#   - returns cumsum(x) / seq_along(x)
# Expectations:
#   - equals c(2,3,4,5)
# Why:
#   - check for helper used in summaries.
# ------------------------------------------------------------------
test_that("cummulativeMovingAverage returns running mean", {
  
  x <- c(2, 4, 6, 8)
  
  out <- cummulativeMovingAverage(x)
  
  expect_equal(out, c(2, 3, 4, 5))
  
})


# Tests for firstUpper ----------------------------------------------------------

# ------------------------------------------------------------------
# Test: firstUpper capitalizes only the first character
# Input:
#   - "hello", "Hello", "h"
# Behaviour:
#   - replaces first character with uppercase
# Expectations:
#   - correct transformed strings
# Why:
#   - ensures string formatting is deterministic.
# ------------------------------------------------------------------
test_that("firstUpper capitalizes first character", {
  
  expect_identical(firstUpper("hello"), "Hello")
  expect_identical(firstUpper("Hello"), "Hello")
  expect_identical(firstUpper("h"), "H")
  
})


# Tests for getBlankString ------------------------------------------------------

# ------------------------------------------------------------------
# Test: getBlankString returns a string of spaces of requested length
# Input:
#   - length = 0, 3, 10
# Behaviour:
#   - concatenates that many " " characters
# Expectations:
#   - correct string and nchar
# Why:
#   - used for formatting and alignment in text outputs.
# ------------------------------------------------------------------
test_that("getBlankString returns requested number of spaces", {
  out0 <- getBlankString(0)
  expect_type(out0, "character")
  expect_length(out0, 0)
  
  out5 <- getBlankString(5)
  expect_identical(out5, "     ")
  expect_type(out5, "character")
  expect_length(out5, 1)
})


# Tests for getListLevel --------------------------------------------------------

# ------------------------------------------------------------------
# Test: getListLevel returns nesting depth
# Input:
#   - non-list and nested lists
# Behaviour:
#   - returns 0 for non-list; otherwise 1 + depth of first element
# Expectations:
#   - correct depths for 0..3 levels
# Why:
#   - scaleRoundList relies on correct depth to traverse lists.
# ------------------------------------------------------------------
test_that("getListLevel returns correct nesting depth", {
  
  expect_identical(getListLevel(1), 0)
  expect_identical(getListLevel(list(1)), 1)
  expect_identical(getListLevel(list(list(1))), 2)
  expect_identical(getListLevel(list(list(list(1)))), 3)
  
})


# Tests for getRowIndexOfVectorInMatrix ----------------------------------------

# ------------------------------------------------------------------
# Test: getRowIndexOfVectorInMatrix finds exact matching row
# Input:
#   - matrix with known rows
#   - vectors equal to those rows
# Behaviour:
#   - returns the indices of rows where all columns match
# Expectations:
#   - correct row index returned; empty integer for no match
# Why:
#   - used to map trial patterns to unique rows.
# ------------------------------------------------------------------
test_that("getRowIndexOfVectorInMatrix finds matching rows", {
  
  m <- matrix(
    c(1, 2,
      3, 4,
      1, 9),
    ncol  = 2,
    byrow = TRUE
  )
  
  expect_identical(getRowIndexOfVectorInMatrix(c(1, 2), m), 1L)
  expect_identical(getRowIndexOfVectorInMatrix(c(3, 4), m), 2L)
  expect_identical(getRowIndexOfVectorInMatrix(c(1, 9), m), 3L)
  
  expect_length(getRowIndexOfVectorInMatrix(c(7, 7), m), 0)
  
})

# ------------------------------------------------------------------
# Test: getRowIndexOfVectorInMatrix errors if vector length != ncol(matrix)
# Input:
#   - vector of wrong length
# Behaviour:
#   - stops with informative error
# Expectations:
#   - error is thrown
# Why:
#   - prevents silent partial matching / wrong indexing.
# ------------------------------------------------------------------
test_that("getRowIndexOfVectorInMatrix errors on length mismatch", {
  
  m <- matrix(1:4, ncol = 2)
  
  expect_error(
    getRowIndexOfVectorInMatrix(c(1, 2, 3), m),
    "length of the vector must be equal"
  )
  
})


# Tests for invLogit ------------------------------------------------------------

# ------------------------------------------------------------------
# Test: invLogit matches stats::binomial()$linkinv and inverts logit
# Input:
#   - theta = 0 and theta = logit(p)
# Behaviour:
#   - returns inverse logit
# Expectations:
#   - invLogit(0) = 0.5
#   - invLogit(logit(p)) = p
# Why:
#   - core link utilities must be numerically correct.
# ------------------------------------------------------------------
test_that("invLogit returns correct inverse-logit values", {
  
  expect_equal(invLogit(0), 0.5)
  
  p <- c(0.1, 0.25, 0.8)
  
  expect_equal(invLogit(logit(p)), p, tolerance = 1e-12)
  
})

# ------------------------------------------------------------------
# Test: invLogit validates inputs
# Input:
#   - missing theta, non-numeric theta
# Behaviour:
#   - stops with the documented error
# Expectations:
#   - error thrown
# Why:
#   - input validation prevents silent misuse.
# ------------------------------------------------------------------
test_that("invLogit errors for missing or non-numeric input", {
  
  expect_error(invLogit(), "Please provide a numeric")
  expect_error(invLogit("x"), "Please provide a numeric")
  
})


# Tests for logit ---------------------------------------------------------------

# ------------------------------------------------------------------
# Test: logit matches stats::binomial()$linkfun and inverts invLogit
# Input:
#   - p = 0.5 and random p values in (0,1)
# Behaviour:
#   - returns log(p/(1-p))
# Expectations:
#   - logit(0.5) = 0
#   - invLogit(logit(p)) = p
# Why:
#   - correctness of link function is fundamental for Berry/ExNex models.
# ------------------------------------------------------------------
test_that("logit returns correct log-odds and inverts invLogit", {
  
  expect_equal(logit(0.5), 0)
  
  p <- c(0.1, 0.25, 0.8)
  
  expect_equal(invLogit(logit(p)), p, tolerance = 1e-12)
  
})

# ------------------------------------------------------------------
# Test: logit validates inputs
# Input:
#   - missing p, non-numeric p, p outside [0,1]
# Behaviour:
#   - stops with the documented error
# Expectations:
#   - error thrown
# Why:
#   - prevents invalid probability inputs leaking into models.
# ------------------------------------------------------------------
test_that("logit errors for invalid input", {
  
  expect_error(logit(), "Please provide a numeric")
  expect_error(logit("x"), "Please provide a numeric")
  expect_error(logit(c(-0.1, 0.2)), "Please provide a numeric")
  expect_error(logit(c(0.2, 1.2)), "Please provide a numeric")
  
})


# Tests for listPerMethod ------------------------------------------------------

# ------------------------------------------------------------------
# Test: listPerMethod reshapes scenario-first lists into method-first lists
# Input:
#   - analyses_list created via simulateScenarios() + performAnalyses()
#   - estimates_per_method from getEstimates() (method-first)
#   - scenario-first object reconstructed from estimates_per_method
# Behaviour:
#   - listPerMethod() should convert scenario-first -> method-first
# Expectations:
#   - Output equals estimates_per_method (same methods, same scenarios, same matrices)
# Why:
#   - Verifies listPerMethod performs pure reshaping on real package outputs.
# ------------------------------------------------------------------
test_that("listPerMethod reshapes scenario-first lists into method-first lists", {
  
  skip_if_not_installed("foreach")
  foreach::registerDoSEQ()
  set.seed(123)
  
  method_names <- c("pooled", "stratified")
  target_rates <- c(0.5, 0.5)
  
  scen <- simulateScenarios(
    n_subjects_list     = list(c(10, 10), c(20, 20)),
    response_rates_list = list(c(0.4, 0.4), c(0.6, 0.6)),
    n_trials = 2
  )
  
  analyses <- performAnalyses(
    scenario_list      = scen,
    evidence_levels    = c(0.05, 0.5, 0.95),
    target_rates       = target_rates,
    method_names       = method_names,
    n_mcmc_iterations  = 200,
    verbose            = FALSE
  )
  
  estimates_per_method <- getEstimates(
    analyses_list   = analyses,
    point_estimator = "median",
    alpha_level     = 0.10
  )
  
  method_names_from_obj <- names(estimates_per_method)
  scenario_names        <- names(estimates_per_method[[1]])
  
  expect_true(all(grepl("^scenario_[0-9]+$", scenario_names)))
  
  scenario_first <- vector("list", length = length(scenario_names))
  names(scenario_first) <- scenario_names
  
  for (s in scenario_names) {
    scenario_first[[s]] <- lapply(
      method_names_from_obj,
      function(m) estimates_per_method[[m]][[s]]
    )
    names(scenario_first[[s]]) <- method_names_from_obj
  }
  
  reshaped <- listPerMethod(scenario_first)
  
  expect_setequal(names(reshaped), method_names_from_obj)
  expect_setequal(names(reshaped[[method_names_from_obj[1]]]), scenario_names)
  
  for (m in method_names_from_obj) {
    for (s in scenario_names) {
      expect_equal(reshaped[[m]][[s]], estimates_per_method[[m]][[s]])
    }
  }
})


# Tests for roundList -----------------------------------------------------------

# ------------------------------------------------------------------
# Test: roundList rounds numerics for nesting levels 1..3
# Input:
#   - list structures with depths 1, 2, 3
# Behaviour:
#   - rounds each numeric entry to round_digits
# Expectations:
#   - correct rounded outputs
# Why:
#   - ensures list traversal is correct for supported nesting depths.
# ------------------------------------------------------------------
test_that("roundList rounds list entries for levels 1..3", {
  
  l1 <- list(1.234, 2.345)
  expect_equal(roundList(l1, round_digits = 1, list_levels = 1), list(1.2, 2.3))
  
  l2 <- list(a = list(1.234), b = list(2.345))
  expect_equal(roundList(l2, round_digits = 2, list_levels = 2), list(a = list(1.23), b = list(2.35)))
  
  l3 <- list(a = list(b = list(1.234)))
  expect_equal(roundList(l3, round_digits = 0, list_levels = 3), list(a = list(b = list(1))))
  
})

# ------------------------------------------------------------------
# Test: roundList errors for non-numeric entries and unsupported depth
# Input:
#   - list containing non-numeric
#   - list_levels > 3
# Behaviour:
#   - stops with informative error
# Expectations:
#   - error thrown
# Why:
#   - avoids silently returning wrong types.
# ------------------------------------------------------------------
test_that("roundList errors for non-numerics and unsupported depth", {
  
  expect_error(roundList(list("x"), round_digits = 1, list_levels = 1), "must contain numerics")
  
  expect_error(roundList(list(1), round_digits = 1, list_levels = 4), "greater than 3")
  
})


# Tests for scaleList -----------------------------------------------------------

# ------------------------------------------------------------------
# Test: scaleList scales numerics for nesting levels 1..3
# Input:
#   - list structures with depths 1, 2, 3
# Behaviour:
#   - multiplies each numeric by scale_param
# Expectations:
#   - correct scaled outputs
# Why:
#   - used by scaleRoundList; must behave deterministically.
# ------------------------------------------------------------------
test_that("scaleList scales list entries for levels 1..3", {
  
  l1 <- list(1, 2)
  expect_equal(scaleList(l1, scale_param = 10, list_levels = 1), list(10, 20))
  
  l2 <- list(a = list(1), b = list(2))
  expect_equal(scaleList(l2, scale_param = 0.5, list_levels = 2), list(a = list(0.5), b = list(1)))
  
  l3 <- list(a = list(b = list(2)))
  expect_equal(scaleList(l3, scale_param = 3, list_levels = 3), list(a = list(b = list(6))))
  
})

# ------------------------------------------------------------------
# Test: scaleList errors for non-numeric entries and unsupported depth
# Input:
#   - list containing non-numeric
#   - list_levels > 3
# Behaviour:
#   - stops with informative error
# Expectations:
#   - error thrown
# Why:
#   - protects against silent type coercion.
# ------------------------------------------------------------------
test_that("scaleList errors for non-numerics and unsupported depth", {
  
  expect_error(scaleList(list("x"), scale_param = 2, list_levels = 1), "must contain numerics")
  
  expect_error(scaleList(list(1), scale_param = 1, list_levels = 4), "greater than 3")
  
})


# Tests for scaleRoundList ------------------------------------------------------

# ------------------------------------------------------------------
# Test: scaleRoundList applies scaling then rounding
# Input:
#   - list of numerics
# Behaviour:
#   - scales by scale_param then rounds to round_digits
# Expectations:
#   - correct scaled+rounded values
# Why:
#   - user-facing helper; expected to behave consistently.
# ------------------------------------------------------------------
test_that("scaleRoundList scales then rounds list entries", {
  
  skip_if_not_installed("checkmate")
  
  l <- list(1.234, 2.345)
  
  out <- scaleRoundList(list = l, scale_param = 10, round_digits = 1)
  
  expect_equal(out, list(12.3, 23.5))
  
})

# ------------------------------------------------------------------
# Test: scaleRoundList returns unchanged list if round_digits is NULL
# Input:
#   - list of numerics, round_digits = NULL
# Behaviour:
#   - scales by scale_param, does no rounding when round_digits is NULL
# Expectations:
#   - equals scaled values (or original if scale_param = 1)
# Why:
#   - matches documented behaviour and avoids surprising rounding.
# ------------------------------------------------------------------
test_that("scaleRoundList skips rounding when round_digits is NULL", {
  
  skip_if_not_installed("checkmate")
  
  l <- list(1.234, 2.345)
  
  out <- scaleRoundList(list = l, scale_param = 1, round_digits = NULL)
  
  expect_equal(out, l)
  
})

# ------------------------------------------------------------------
# Test: scaleRoundList validates inputs
# Input:
#   - missing list / non-list
#   - invalid scale_param
#   - invalid round_digits
# Behaviour:
#   - stops with informative errors
# Expectations:
#   - error thrown
# Why:
#   - ensures correct user-facing input checks.
# ------------------------------------------------------------------
test_that("scaleRoundList errors for invalid inputs", {
  
  skip_if_not_installed("checkmate")
  
  expect_error(scaleRoundList(), "Please provide a list")
  expect_error(scaleRoundList(list = 1), "Please provide a list")
  
  l <- list(1.234)
  
  expect_error(scaleRoundList(list = l, scale_param = -1))
  expect_error(scaleRoundList(list = l, round_digits = -1))
  
})


# Tests for substrRight ---------------------------------------------------------

# ------------------------------------------------------------------
# Test: substrRight returns rightmost n characters
# Input:
#   - x = "abcdef", n = 1,2,6
# Behaviour:
#   - returns substring from the right
# Expectations:
#   - correct suffixes
# Why:
#   - used for parsing parameter/cohort suffix indices.
# ------------------------------------------------------------------
test_that("substrRight returns rightmost n characters", {
  
  expect_identical(substrRight("abcdef", 1), "f")
  expect_identical(substrRight("abcdef", 2), "ef")
  expect_identical(substrRight("abcdef", 6), "abcdef")
  
})