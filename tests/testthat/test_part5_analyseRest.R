library(GGIR)
context("g.part5.analyseRest")
tz = "Europe/Amsterdam"
ds_names = rep("", 40)
dsummary = matrix("", 1, 40)
startday = as.POSIXlt(x = "2022-06-02 08:00:00", tz = tz)
time = seq(startday, startday + (16 * 3600), by = 60)



test_that("Overlap 1 nap and 1 sib", {
  fi = 1
  di = 1
  sibreport = data.frame(ID = rep("test123", 2), type = c("nap", "sib"),
                         start = c("2022-06-02 14:00:00", "2022-06-02 14:05:00"),
                         end = c("2022-06-02 14:20:00", "2022-06-02 14:20:00"))
  sibreport$duration = as.numeric(difftime(time1 = sibreport$end, time2 = sibreport$start, units = "mins", tz = tz))
  restAnalyses = g.part5.analyseRest(sibreport = sibreport, dsummary = dsummary,
                                     ds_names = ds_names, fi = fi, di = di,
                                     time = time,
                                     tz = tz)
  fi = restAnalyses$fi
  di = restAnalyses$di
  dsummary = restAnalyses$dsummary
  ds_names = restAnalyses$ds_names
  dsummary = as.numeric(dsummary[, which(ds_names != "")])
  ds_names = ds_names[ds_names != ""]
  names(dsummary) = ds_names
  dsummary = as.data.frame(t(dsummary))
  
  expect_equal(dsummary$nbouts_day_sib, 1)
  expect_equal(dsummary$nbouts_day_srnap, 1)
  expect_equal(dsummary$frag_mean_dur_sib_day, 15)
  expect_equal(dsummary$frag_mean_dur_srnap_day, 20)
  expect_equal(dsummary$perc_sib_overl_srnap, 100)
  expect_equal(dsummary$perc_srnap_overl_sib, 75)
  expect_equal(sum(dsummary), 323)
})

test_that("Overlap 1 nonwear and 1 sib", {
  fi = 1
  di = 1
  sibreport = data.frame(ID = rep("test123", 2), type = c("nonwear", "sib"),
                         start = c("2022-06-02 14:00:00", "2022-06-02 14:05:00"),
                         end = c("2022-06-02 14:20:00", "2022-06-02 14:20:00"))
  sibreport$duration = as.numeric(difftime(time1 = sibreport$end, time2 = sibreport$start, units = "mins", tz = tz))
  restAnalyses = g.part5.analyseRest(sibreport = sibreport, dsummary = dsummary,
                                     ds_names = ds_names, fi = fi, di = di,
                                     time = time,
                                     tz = tz)
  fi = restAnalyses$fi
  di = restAnalyses$di
  dsummary = restAnalyses$dsummary
  ds_names = restAnalyses$ds_names
  dsummary = as.numeric(dsummary[, which(ds_names != "")])
  ds_names = ds_names[ds_names != ""]
  names(dsummary) = ds_names
  dsummary = as.data.frame(t(dsummary))
  
  expect_equal(dsummary$nbouts_day_sib, 1)
  expect_equal(dsummary$nbouts_day_srnonw, 1)
  expect_equal(dsummary$frag_mean_dur_sib_day, 15)
  expect_equal(dsummary$frag_mean_dur_srnonw_day, 20)
  expect_equal(dsummary$perc_sib_overl_srnonw, 100)
  expect_equal(dsummary$perc_srnonw_overl_sib, 75)
  expect_equal(sum(dsummary), 323)
})


test_that("No overlap 1 nonwear, 1 nap, and 1 sib", {
  fi = 1
  di = 1
  sibreport = data.frame(ID = rep("test123", 3), type = c("nonwear", "nap", "sib"),
                         start = c("2022-06-02 12:00:00", "2022-06-02 13:00:00", "2022-06-02 15:00:00"),
                         end = c("2022-06-02 12:20:00", "2022-06-02 13:20:00", "2022-06-02 15:20:00"))
  sibreport$duration = as.numeric(difftime(time1 = sibreport$end, time2 = sibreport$start, units = "mins", tz = tz))
  restAnalyses = g.part5.analyseRest(sibreport = sibreport, dsummary = dsummary,
                                     ds_names = ds_names, fi = fi, di = di,
                                     time = time,
                                     tz = tz)
  fi = restAnalyses$fi
  di = restAnalyses$di
  dsummary = restAnalyses$dsummary
  ds_names = restAnalyses$ds_names
  dsummary = as.numeric(dsummary[, which(ds_names != "")])
  ds_names = ds_names[ds_names != ""]
  names(dsummary) = ds_names
  dsummary = as.data.frame(t(dsummary))
  
  expect_equal(dsummary$nbouts_day_sib, 1)
  expect_equal(dsummary$nbouts_day_srnap, 1)
  expect_equal(dsummary$nbouts_day_srnonw, 1)
  expect_equal(dsummary$frag_mean_dur_sib_day, 20)
  expect_equal(dsummary$frag_mean_dur_srnap_day, 20)
  expect_equal(dsummary$frag_mean_dur_srnonw_day, 20)
  expect_equal(dsummary$perc_sib_overl_srnap, 0)
  expect_equal(dsummary$perc_sib_overl_srnonw, 0)
  expect_equal(dsummary$perc_srnonw_overl_sib, 0)
  expect_equal(sum(dsummary), 129)
})
