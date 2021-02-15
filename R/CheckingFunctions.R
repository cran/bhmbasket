
### numeric ####

is.positive.numeric <- function (x) {

  return (is.numeric(x) && all(x > 0))

}

is.numeric.in.zero.one <- function (x) {

  return (is.numeric(x) && all(x > 0) && all(x < 1))

}

is.single.numeric <- function (x) {

  return (is.numeric(x) && length(x) == 1)

}

is.single.positive.numeric <- function (x) {

  return (is.single.numeric(x) && x > 0)

}

is.single.numeric.in.zero.one <- function (x) {

  return (is.single.numeric(x) && is.numeric.in.zero.one(x))

}

### whole number ####

is.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  if (is.numeric(x)) {

    return (abs(x - round(x)) < tol)

  } else {

    return (FALSE)

  }

}

is.positive.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  return (all(is.wholenumber(x, tol = tol)) && all(x > 0))

}

is.non.negative.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  return (all(is.wholenumber(x, tol = tol)) && all(x >= 0))

}

is.single.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  return (all(is.wholenumber(x, tol = tol)) && length(x) == 1)

}

is.single.non.negative.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  return (is.single.wholenumber(x, tol = tol) && all(is.non.negative.wholenumber(x, tol = tol)))

}

is.single.positive.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  return (is.single.wholenumber(x, tol = tol) && all(is.positive.wholenumber(x, tol = tol)))

}
