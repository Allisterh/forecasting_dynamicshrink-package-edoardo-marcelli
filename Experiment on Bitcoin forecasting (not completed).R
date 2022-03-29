# BITCOIN ESTIMATION WITH EXCHANGE-RATES, VIX, GOOGLE TRENDS,

library(BatchGetSymbols)
library(hciR)

# DATA-SET BUILDING

BTC_USD = BatchGetSymbols(
  tickers="BTC-USD",
  first.date = "2015-12-31",
  last.date = "2020-12-31",
  thresh.bad.data = 0.75,
  bench.ticker = "BTC-USD",
  type.return = "arit",
  freq.data = "weekly",
  how.to.aggregate = "last",
  do.complete.data = FALSE,
  do.fill.missing.prices = TRUE,
  do.cache = TRUE,
  cache.folder = file.path(tempdir(), "BGS_Cache"),
  do.parallel = FALSE,
  be.quiet = FALSE
)

USD_EUR = BatchGetSymbols(
  tickers="USDEUR=X",
  first.date = "2015-12-31",
  last.date = "2020-12-31",
  thresh.bad.data = 0.75,
  bench.ticker = "USDEUR=X",
  type.return = "arit",
  freq.data = "weekly",
  how.to.aggregate = "last",
  do.complete.data = FALSE,
  do.fill.missing.prices = TRUE,
  do.cache = TRUE,
  cache.folder = file.path(tempdir(), "BGS_Cache"),
  do.parallel = FALSE,
  be.quiet = FALSE
)

VIX = BatchGetSymbols(
  tickers="^VIX",
  first.date = "2015-12-31",
  last.date = "2020-12-31",
  thresh.bad.data = 0.75,
  bench.ticker = "^VIX",
  type.return = "arit",
  freq.data = "weekly",
  how.to.aggregate = "last",
  do.complete.data = FALSE,
  do.fill.missing.prices = TRUE,
  do.cache = TRUE,
  cache.folder = file.path(tempdir(), "BGS_Cache"),
  do.parallel = FALSE,
  be.quiet = FALSE
)

CNY_USD = BatchGetSymbols(
  tickers="CNYUSD=X",
  first.date = "2015-12-31",
  last.date = "2020-12-31",
  thresh.bad.data = 0.75,
  bench.ticker = "CNYUSD=X",
  type.return = "arit",
  freq.data = "weekly",
  how.to.aggregate = "last",
  do.complete.data = FALSE,
  do.fill.missing.prices = TRUE,
  do.cache = TRUE,
  cache.folder = file.path(tempdir(), "BGS_Cache"),
  do.parallel = FALSE,
  be.quiet = FALSE
)


Gold = BatchGetSymbols(
  tickers="GC=F",
  first.date = "2015-12-31",
  last.date = "2020-12-31",
  thresh.bad.data = 0.75,
  bench.ticker = "GC=F",
  type.return = "arit",
  freq.data = "weekly",
  how.to.aggregate = "last",
  do.complete.data = FALSE,
  do.fill.missing.prices = TRUE,
  do.cache = TRUE,
  cache.folder = file.path(tempdir(), "BGS_Cache"),
  do.parallel = FALSE,
  be.quiet = FALSE
)



Extractcolumn <- function(x,column="price.close"){
  findcolumn = as_matrix(x$df.tickers[,column])
  mat = matrix(as.numeric(rownames(findcolumn)))
  colnames(mat) = deparse(substitute(x))
  mat
}