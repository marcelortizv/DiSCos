#-----------------------------------------------------------------------------
# permutation determinism (perm_seed) and quantile sub-range (perm_q_range)
#-----------------------------------------------------------------------------

test_that("perm_seed makes the permutation p-value independent of num.cores", {
  Ts <- 3
  t0 <- 3
  set.seed(123)
  df <- ex_gmm(Ts=Ts, num.con=6)

  a <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, permutation=TRUE,
             perm_seed=1L, num.cores=1, graph=FALSE)
  b <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, permutation=TRUE,
             perm_seed=1L, num.cores=2, graph=FALSE)

  # perm_seed pins the RNG of each placebo iteration, so the placebo distances
  # (distp) — and hence the permutation distribution — are identical at any core
  # count. (The main fit's own RNG, controlled by `seed`, still uses per-core
  # L'Ecuyer streams, so distt is not covered by this guarantee.)
  expect_equal(a$perm$distp, b$perm$distp)
  expect_equal(a$perm$p_overall, b$perm$p_overall)
})


test_that("perm_q_range spanning the full distribution reproduces the default", {
  Ts <- 3
  t0 <- 3
  set.seed(123)
  df <- ex_gmm(Ts=Ts, num.con=6)

  full <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, permutation=TRUE,
                perm_seed=1L, num.cores=1, graph=FALSE)
  spanned <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, permutation=TRUE,
                   perm_seed=1L, perm_q_range=c(0, 1), num.cores=1, graph=FALSE)

  expect_equal(full$perm$p_overall, spanned$perm$p_overall)
  expect_equal(full$perm$distt, spanned$perm$distt)
})


test_that("perm_q_range sub-range returns a valid p-value without changing the fit", {
  Ts <- 3
  t0 <- 3
  set.seed(123)
  df <- ex_gmm(Ts=Ts, num.con=6)

  full <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, permutation=TRUE,
                perm_seed=1L, num.cores=1, graph=FALSE)
  upper <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, permutation=TRUE,
                 perm_seed=1L, perm_q_range=c(0.5, 1), num.cores=1, graph=FALSE)

  expect_true(upper$perm$p_overall >= 0 && upper$perm$p_overall <= 1)
  # the synthetic-control fit is not affected by the sub-range restriction
  expect_equal(full$weights, upper$weights)
})
