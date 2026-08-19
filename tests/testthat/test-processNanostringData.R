test_that("RUVg handles technical replicate averaging", {
  
  example_data <- system.file(
    "extdata",
    "GSE117751_RAW",
    package = "NanoTube"
  )
  
  sample_data <- system.file(
    "extdata",
    "GSE117751_sample_data_replicates.csv",
    package = "NanoTube"
  )
  
  warnings <- character()
  
  dat <- withCallingHandlers(
    processNanostringData(
      nsFiles = example_data,
      sampleTab = sample_data,
      groupCol = "Sample_Diagnosis",
      replicateCol = "Replicate_ID",
      normalization = "RUVg",
      bgType = "t.test",
      bgPVal = 0.01,
      output.format = "ExpressionSet"
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  # The RUVg "does not contain counts" warning should disappear.
  expect_false(any(grepl(
    "expression matrix does not contain counts",
    warnings,
    ignore.case = TRUE
  )))
  
  expect_s4_class(dat, "ExpressionSet")
  
  expr <- Biobase::exprs(dat)
  
  # RUVg should return count-scale, integer-valued normalized counts.
  expect_true(all(expr >= 0))
  expect_true(all(expr == round(expr)))
})