library(GGIR)
context("gloadlog")
test_that("gloadlog is able to load different log formats", {
  tempdir = "mytestdir"
  if (dir.exists(tempdir) == FALSE) dir.create(tempdir)
  ID = "123A"

  # test that a sleep log with any separator from the set of [,\t |;] can be successfully parsed
  for (sep in c(",", "\t", " ", "|", ";")) {
    create_test_sleeplog_csv(Nnights = 7, advanced = TRUE, sep = sep)
    fn = "testsleeplogfile.csv"
    expect_true(file.exists(fn))

    # rec_starttime = "2016-06-25T20:20:20+0200"
    rec_starttime = format(as.POSIXct("2016-06-25 20:20:20"), "%Y-%m-%dT%H:%M:%S%z")
    save(ID, rec_starttime, file = "mytestdir/dummyms3.RData")
    cat(paste0("\nLocal timezone: ", Sys.timezone()))
    logs1 = g.loadlog(loglocation = fn, coln1 = 2, colid = 1, #nnights = 7,
                      meta.sleep.folder = tempdir , desiredtz = "")
    # after deprecating nnights, the whole sleeplog is read (it contains 8 nights)
    # ammend tests as needed
    expect_equal(nrow(logs1$sleeplog), 8) 
    expect_equal(ncol(logs1$sleeplog), 5)
    expect_equal(logs1$sleeplog$night, as.character(1:8))
    expect_equal(logs1$sleeplog$sleepwake, c("7:0:2", "7:0:3", "7:0:4", "7:0:5", "7:0:6", "7:0:7", "7:0:8", "7:0:9"))
    expect_equal(nrow(logs1$nonwearlog), 9)
    expect_equal(ncol(logs1$nonwearlog), 4)
    expect_equal(nrow(logs1$naplog), 9)
    expect_equal(ncol(logs1$naplog), 6)
    
    # start of accelerometer recording 6 days earlier
    rec_starttime = format(as.POSIXct("2016-06-20 20:20:20"), "%Y-%m-%dT%H:%M:%S%z")
    save(ID, rec_starttime, file = "mytestdir/dummyms3.RData")
    logs2 = g.loadlog(loglocation = fn, coln1 = 2, colid = 1, #nnights = 7,
                     meta.sleep.folder = tempdir , desiredtz = "")
    expect_equal(nrow(logs2$sleeplog), 8)
    expect_equal(ncol(logs2$sleeplog), 5)
    expect_equal(logs2$sleeplog$night, as.character(6:13)) # these are nights in the acc recording
    expect_equal(logs2$sleeplog$sleepwake, c("7:0:2", "7:0:3", "7:0:4", "7:0:5", "7:0:6", "7:0:7", "7:0:8", "7:0:9"))
    expect_equal(nrow(logs2$nonwearlog), 9)
    expect_equal(ncol(logs2$nonwearlog), 4)
    expect_equal(nrow(logs2$naplog), 9)
    expect_equal(ncol(logs2$naplog), 6)
    
    # start of accelerometer recording 3 days later
    rec_starttime = format(as.POSIXct("2016-06-29 20:20:20"), "%Y-%m-%dT%H:%M:%S%z")
    save(ID, rec_starttime, file = "mytestdir/dummyms3.RData")
    logs3 = g.loadlog(loglocation = fn, coln1 = 2, colid = 1, #nnights = 7,
                      meta.sleep.folder = tempdir , desiredtz = "")
    expect_equal(nrow(logs3$sleeplog), 4)
    expect_equal(ncol(logs3$sleeplog), 5)
    expect_equal(logs3$sleeplog$night, as.character(1:4)) # 1:3 because these are the first three nights of the acc recording
    expect_equal(logs3$sleeplog$sleeponset, c("23:0:5", "23:0:6", "23:0:7", "23:0:8"))
    expect_equal(nrow(logs3$nonwearlog), 5)
    expect_equal(ncol(logs3$nonwearlog), 4)
    expect_equal(nrow(logs3$naplog), 5)
    expect_equal(ncol(logs3$naplog), 6)
  }

  if (dir.exists(tempdir)) unlink(tempdir, recursive = TRUE)
})
