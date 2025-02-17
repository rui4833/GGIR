convertEpochData = function(datadir = c(), metadatadir = c(),
                            params_general = c(), verbose = TRUE) {
  # Function to convert externally derived epoch data
  # to a format that allows GGIR to process them as if they were part 1 output.
  
  # Capture current local settings
  Syslocale = Sys.getenv("LC_TIME")
  Sys.setlocale("LC_TIME", "C") # set language to English
  tz = params_general[["desiredtz"]]
  epSizeShort = params_general[["windowsizes"]][1]
  epSizeLong = params_general[["windowsizes"]][2]
  # Identify input data file extensions
  if (dir.exists(datadir) == FALSE) {
    stop("\nWhen working with external data, argument datadir is expected to be a directory")
  }
  fnames_csv = dir(datadir,full.names = TRUE,pattern = "[.]csv")
  fnames_awd = dir(datadir,full.names = TRUE,pattern = "[.]awd|[.]AWD")
  if (length(fnames_csv) > 0 & length(fnames_awd) > 0) {
    stop("Do not mix csv and awd files in the same data directory")
  } else {
    if (length(fnames_awd) > 0) {
      if (params_general[["dataFormat"]] != "actiwatch_awd") {
        stop("Specified dataFormat does not match the data")
      }
      fnames = fnames_awd
    } else {
      if (params_general[["dataFormat"]] == "actiwatch_awd") {
        stop("Specified dataFormat does not match the data")
      }
      fnames = fnames_csv
    }
  }
  #-------------
  # Create output folder, normally with raw data g.part1 would do this:
  if (!file.exists(metadatadir)) {
    dir.create(metadatadir)
    dir.create(file.path(metadatadir, "meta"))
    dir.create(file.path(metadatadir, "meta", "basic"))
    dir.create(file.path(metadatadir, "results"))
    dir.create(file.path(metadatadir, "results", "QC"))
  }
  
  #============
  # Based on knowledge about data format
  # we can already assign monitor names and codes
  # dataformat names and codes
  # and sample rate
  actiwatchData = length(grep(pattern = "actiwatch", x = params_general[["dataFormat"]], ignore.case = TRUE)) > 0
  if (params_general[["dataFormat"]] == "ukbiobank_csv") {
    deviceName = "Axivity"
    monn = "axivity"
    monc = 4
    dformc = 4
    dformn = "cwa"
    sf = 100
    epSizeShort = 5
  } else if (params_general[["dataFormat"]] == "actigraph_csv") {
    deviceName = "ActiGraph"
    monn = "actigraph"
    monc = 3
    dformc = 2
    dformn = "csv"
    sf = 0 # extract from file header
    epSizeShort = 0 # extract from file header
  } else if (actiwatchData == TRUE) {
    deviceName = "Actiwatch"
    monn = "actiwatch"
    monc = 99
    dformc = 99
    dformn = "epochdata"
    sf = 100 # <= EXTRACT FROM FILE?
  }
  
  # Before we look inside the epoch files we can already create templates
  # with dummy data
  C = list(cal.error.end = 0, cal.error.start = 0)
  C$scale = c(1, 1, 1)
  C$offset = c(0, 0, 0)
  C$tempoffset =  c(0, 0, 0)
  C$QCmessage = "Autocalibration not done"
  C$npoints = 0
  C$nhoursused = 0
  C$use.temp = TRUE
  M = list(filecorrupt = FALSE, filetooshort = FALSE,
           metalong = data.frame(
             timestamp = rep(0, 3),
             nonwearscore = rep(0, 3),
             clippingscore = rep(0, 3),
             lightmean = rep(0, 3),
             lightpeak = rep(0, 3),
             temperaturemean = rep(0, 3),
             EN = rep(0, 3)
           ),
           metashort = data.frame(timestamp = rep(0, 3),
                                  accmetric = rep(0, 3)), 
           wday = 0,
           wdayname = "Sunday",
           windowsizes = params_general[["windowsizes"]],
           bsc_qc = data.frame(A = 1:3,
                               B = 1:3))
  dummyheader = data.frame(uniqueSerialCode = 0, 
                           frequency = 100,
                           start = "", # if relevant, this will be filled later on in code
                           device = deviceName,
                           firmwareVersion = "unknown",
                           block = 0)
  dummyheader = t(dummyheader)
  colnames(dummyheader) = "value"
  dummyheader = as.data.frame(dummyheader)
  I = list(
    header = dummyheader,
    monc = monc,
    monn = monn,
    dformc = dformc,
    dformn = dformn,
    sf = sf,
    filename = "unknown",
    deviceSerialNumber = "unknown")
  
  detectQuote = function(fn, index) {
    # From experience we know that on some machine the quotes in the files are
    # poorly recognised, to catch this first try to check whether this is the case:
    quote = "\""
    Dtest = NULL
    try(expr = {Dtest = data.table::fread(input = fn,
                                          header = FALSE, sep = ",", skip = index,
                                          nrows = 20, quote = quote)}, silent = TRUE)
    if (length(Dtest) == 0) {
      quote = ""
    } else {
      if (nrow(Dtest) <= 1) {
        quote = "" 
      }
    }
    return(quote)
  }
  I_bu = I # backup version of I (for now only used in actigraph_csv)
  for (i in 1:length(fnames)) { # loop over all epoch files
    # include verbose info
    if (verbose == TRUE) {
      if (i  == 1) {
        cat(paste0("\nLoading file: ", i))
      } else {
        cat(paste0(" ", i))
      }
    }
    # filename
    fname = basename(fnames[i])
    
    outputFileName = paste0(metadatadir, "/meta/basic/meta_", fname, ".RData")
    
    skip = TRUE
    if (params_general[["overwrite"]] == TRUE) {
      skip = FALSE
    } else {
      if (!file.exists(outputFileName)) {
        skip = FALSE
      }
    }
    if (skip == FALSE) {
      I = I_bu
      I$filename = filename_dir = fname
      if (params_general[["dataFormat"]] == "ukbiobank_csv") {
        # read data
        D = data.table::fread(input = fnames[i],
                              header = TRUE,
                              data.table = FALSE,
                              sep = ",")
        header = as.character(data.table::fread(input = fnames[i], header = FALSE, nrows = 1, data.table = FALSE, sep = ",")[1, 1])
        # extract date/timestamp from fileheader
        timestamp = unlist(strsplit(header," - "))[2]
        timestamp_POSIX = as.POSIXlt(timestamp, tz = tz)
      } else if (params_general[["dataFormat"]] == "actigraph_csv") {
        # check if file was exported with header:
        header_test = FALSE
        header = data.table::fread(input = fnames[i], header = FALSE, nrows = 10, data.table = FALSE, sep = ",")
        splitHeader = function(x) {
          tmp = unlist(strsplit(x, " "))
          variable = gsub(pattern = ":| ", replacement = "", x = paste0(tmp[1:(length(tmp) - 1)], collapse = ""))
          df = data.frame(variable = tolower(variable), value = tmp[length(tmp)])
          return(df)
        }
        
        AGh = NULL
        for (hh in header[2:9,1]) {
          AGh = rbind(AGh, splitHeader(hh))
        }
        
        if (any(grepl("serialnumber", AGh$variable))) header_test = TRUE
        
        # rows to skip:
        if (header_test == TRUE) {
          skip = 10
        } else {
          tmp = data.table::fread(input = fnames[i],
                                  header = FALSE,
                                  data.table = FALSE,
                                  skip = 0,
                                  nrows = 1)
          if (any(grepl("data|scoring", tmp[1,]))) {
            skip = 1
          } else {
            skip = 0
          }
        }
        
        # check if file was exported with column names:
        colnames = FALSE
        colnames_test = data.table::fread(input = fnames[i],
                                        header = FALSE,
                                        data.table = FALSE,
                                        skip = skip,
                                        nrows = 1)
        if (any(grepl("Axis|vector magnitude|vm", colnames_test[1,], ignore.case = TRUE))) {
          colnames = TRUE
        }
        
        # read data
        D = data.table::fread(input = fnames[i],
                              header = colnames,
                              data.table = FALSE,
                              skip = skip)
        
        # Identify time and acceleration columns 
        acccol = vmcol = NA
        if (colnames == TRUE) {
          acccol = grep("axis", colnames(D), ignore.case = TRUE)
          vmcol = grep("vector magnitude|vm", colnames(D), ignore.case = TRUE)
        } else {
          # Then assume 
          # first 3 columns are axis1, axis2, axis3 if ncol(D) >= 3
          # first column is VM if ncol(D) < 3
          # NOTE: this could not be true, from actilife, the user can select
          # the columns to export (e.g, it could be "Axis1", "Vector Magnitude", "Steps")
          if (ncol(D) >= 3) {
            acccol = 1:3
          } else {
            vmcol = 1
          }
        }
        # D colnames and formatting
        colnames(D)[acccol] = c("NeishabouriCount_y", "NeishabouriCount_x", "NeishabouriCount_z")
        colnames(D)[vmcol] = c("NeishabouriCount_vm")
        keep = c(acccol, vmcol)[!is.na(c(acccol, vmcol))]
        D = D[, keep]
        if (ncol(D) == 3 & is.na(vmcol)) {
          D$NeishabouriCount_vm = sqrt(D[,1]^2 + D[,2]^2 + D[,3]^2)
        }
        
        # extract information from header
        
        if (header_test == TRUE) {
          I$deviceSerialNumber = AGh$value[grep(pattern = "serialnumber", x = AGh$variable)]
          # Add serial number to header object because g.extractheadersvars for Actigraph uses this
          I$header[nrow(I$header) + 1, 1] = NA
          I$header$value[nrow(I$header)] = as.character(I$deviceSerialNumber)
          row.names(I$header)[nrow(I$header)] = "Serial Number:"
          
          epochSize = AGh$value[grep(pattern = "epochperiod", x = AGh$variable)]
          epSizeShort = sum(as.numeric(unlist(strsplit(epochSize, ":"))) * c(3600, 60, 1))
          
          # extract date/timestamp from fileheader
          starttime = AGh$value[grep(pattern = "starttime", x = AGh$variable)]
          startdate = AGh$value[grep(pattern = "startdate", x = AGh$variable)]
          
          timestamp = paste0(startdate, " ", starttime)
          I$header[which(rownames(I$header) == "start"), 1] = timestamp
          timestamp_POSIX = as.POSIXlt(timestamp, tz = tz,
                                       format = params_general[["extEpochData_timeformat"]])
        } else if (header_test == FALSE) {
          # extract date/timestamp if included in the columns
          tmp = data.table::fread(input = fnames[i],
                                  header = colnames,
                                  data.table = FALSE,
                                  skip = skip,
                                  nrows = 2)
          if (colnames == TRUE) {
            datecol = grep("date", colnames(tmp), ignore.case = TRUE)
            timecol = grep("time|epoch", colnames(tmp), ignore.case = TRUE)
            time = tmp[, timecol]
            date = tmp[, datecol]
            
            starttime = time[1]
            starttime = date[1]
            timestamp = paste0(date, " ", time)
            I$header[which(rownames(I$header) == "start"), 1] = timestamp[1]
            format = params_general[["extEpochData_timeformat"]]
            timestamp_POSIX = as.POSIXlt(timestamp, tz = tz, format = format)
            if (all(is.na(timestamp_POSIX))) {
              stop(paste0("\nTimestamps are not available in hte file, neither it has
                          a header to extract the timestamps from. Therefore, the file
                          cannot be processed.\n"))
            }
            
            epochSize = difftime(timestamp_POSIX[2], timestamp_POSIX[1], 
                                 units = "secs")
            epSizeShort = as.numeric(epochSize)
            timestamp_POSIX = timestamp_POSIX[1]
          }
        }
        # CHECKS
        if (all(is.na(timestamp_POSIX))) {
          stop(paste0("\nTime format in data ", timestamp, " does not match with time format ",
                      params_general[["extEpochData_timeformat"]],
                      " as specified by argument extEpochData_timeformat, please correct.\n"))
        }
        if (epSizeShort != params_general[["windowsizes"]][1]) {
          stop(paste0("\nThe short epoch size as specified by the user as the first value of argument windowsizes (",
                      params_general[["windowsizes"]][1],
                      " seconds) does NOT match the short epoch size we see in the data (", epSizeShort),
               " seconds). Please correct.", call. = FALSE)
        }
        # mode = AGh$value[grep(pattern = "mode", x = AGh$variable)]
      } else if (length(grep(pattern = "actiwatch", x = params_general[["dataFormat"]], ignore.case = TRUE)) > 0) {
        if (params_general[["dataFormat"]] == "actiwatch_csv") {
          # ! Assumptions that timeseries start before line 1000
          index = 1000
          while (index > 0) {
            quote = detectQuote(fn = fnames[i], index = index)
            testraw = data.table::fread(input = fnames[i],
                                        header = FALSE, sep = ",", skip = index,
                                        nrows = 2, data.table = FALSE, quote = quote)
            if (length(testraw) > 0) {
              if (nrow(testraw) == 2) {
                if (testraw$V1[2] == testraw$V1[1] + 1) {
                  break
                }
              }
            }
            index = index - 100
          }
          # ! Assumption that first column are the epoch numbers
          delta = 1 - testraw$V1[1]
          index = index + delta
          startFound = FALSE
          while (startFound == FALSE) {
            Dtest = data.table::fread(input = fnames[i], sep = ",", skip = index, quote = quote, nrows = 1)  
            if (Dtest$V1[1] == 1) {
              startFound = TRUE
            } else {
              # This happens when file is has an empty row between each measurement point is stored
              index = index - 1
              if (index < 1) stop("Could not find start of recording", call. = FALSE)
            }
          }
          
          D = data.table::fread(input = fnames[i], sep = ",", skip = index, quote = quote)
          # ! Assumption that column names are present 2 lines prior to timeseries
          colnames = data.table::fread(input = fnames[i],
                                       header = FALSE, sep = ",",
                                       skip = index - 2, nrows = 1, quote = quote)
          if (all(is.na(colnames))) {
            colnames = data.table::fread(input = fnames[i],
                                         header = FALSE, sep = ",",
                                         skip = index - 4, nrows = 1, quote = quote)
          }
          colnames(D) = as.character(colnames)[1:ncol(D)]
          
          # ! Assumptions about columns names
          colnames(D) = gsub(pattern = "datum|date", replacement = "date", x = colnames(D), ignore.case = TRUE)
          colnames(D) = gsub(pattern = "tijd|time", replacement = "time", x = colnames(D), ignore.case = TRUE)
          colnames(D) = gsub(pattern = "activiteit|activity", replacement = "ZCY", x = colnames(D), ignore.case = TRUE)
          I$header[which(rownames(I$header) == "start"), 1] = paste(D$date[1], D$time[1], sep = " ")
          
          timestamp_POSIX = as.POSIXct(x = paste(D$date[1:4], D$time[1:4], sep = " "),
                                       format = params_general[["extEpochData_timeformat"]],
                                       tz = tz)
          if (all(is.na(timestamp_POSIX))) {
            stop(paste0("\nTime format in data ", D$date[1], " does not match with time format ",
                        params_general[["extEpochData_timeformat"]],
                        " as specified by argument extEpochData_timeformat, please correct.\n"))
          }
          epSizeShort = mean(diff(as.numeric(timestamp_POSIX)))
          if (epSizeShort != params_general[["windowsizes"]][1]) {
            stop(paste0("\nThe short epoch size as specified by the user as the first value of argument windowsizes (",
                        params_general[["windowsizes"]][1],
                        " seconds) does NOT match the short epoch size we see in the data (", epSizeShort),
                 " seconds). Please correct.", call. = FALSE)
          }
          timestamp_POSIX = timestamp_POSIX[1]
          D = D[, "ZCY"]
          if (quote == "") D$ZCY = as.numeric(D$ZCY)
        } else if (params_general[["dataFormat"]] == "actiwatch_awd") {
          # ! Assumption that first data row equals the first row with 3 columns
          index = 0
          
          quote = detectQuote(fn = fnames[i], index = 50)
          NC = 1
          while (NC >= 3) {
            
            testraw = data.table::fread(input = fnames[i],
                                        header = FALSE, sep = ",", skip = index,
                                        nrows = 1, data.table = TRUE, quote = quote)
            NC = ncol(testraw)
            if (NC >= 3) {
              break()
            } else {
              index = index + 1
            }
          }
          D = data.table::fread(input = fnames[i],
                                header = FALSE, sep = ",", skip = index, quote = quote)
          D = D[,1]
          colnames(D)[1] = "ZCY"
          header = data.table::fread(input = fnames[i],
                                     header = FALSE, sep = ",", nrows =  7, quote = quote)
          # Get epoch size
          optionalEpochs = data.frame(code = c("1", "2", "4", "8", "20", "81", "C1", "C2"),
                                      size = c(15, 30, 60, 120, 300, 2, 5, 10))
          epSizeShort = optionalEpochs$size[which(optionalEpochs$code == as.character(header[4]))]
          # Get starttime 
          timestampFormat = paste0(unlist(strsplit(params_general[["extEpochData_timeformat"]], " "))[1], " %H:%M")
          I$header[which(rownames(I$header) == "start"), 1] =  paste(header[2], header[3], sep = " ")
          timestamp_POSIX = as.POSIXct(x = paste(header[2], header[3], sep = " "),
                                       format = timestampFormat, tz = tz)
          if (is.na(timestamp_POSIX)) {
            stop(paste0("\nTime format in data ", header[2], " does not match with time format ",
                        params_general[["extEpochData_timeformat"]],
                        " as specified by argument extEpochData_timeformat, please correct.\n"))
          }
          if (epSizeShort != params_general[["windowsizes"]][1]) {
            stop(paste0("\nThe short epoch size as specified by the user as the first value of argument windowsizes (",
                        params_general[["windowsizes"]][1],
                        " seconds) does NOT match the short epoch size we see in the data (", epSizeShort),
                 " seconds). Please correct.", call. = FALSE)
          }
          if (quote == "") D$ZCY = as.numeric(D$ZCY)
        }
      }
      
      quartlystart = (ceiling(as.numeric(timestamp_POSIX) / epSizeLong)) * epSizeLong
      
      ts_longEp_POSIX = as.POSIXlt(quartlystart, tz = tz, origin = "1970-01-01") # POSIX time
      M$wdayname = weekdays(x = ts_longEp_POSIX, abbreviate = FALSE) # extract weekday as day name
      M$wday = ts_longEp_POSIX$wday # extract weekday as a daynumber
      
      if (length(ts_longEp_POSIX) > 2) {
        # calculate time difference relative to original data
        # to know how much epoch to ignore at the beginning of the file
        deltastart = as.numeric(difftime(ts_longEp_POSIX, timestamp_POSIX,units = "sec") / 5)
        # shorten the data (object D) accordingly
        D = D[(deltastart + 1):nrow(D),]
      }
      # Define compatible metashort and metalong dimensions
      LMS = nrow(D)
      LML = LMS * epSizeShort / epSizeLong
      expected_LML = floor(LML)
      if (expected_LML != LML) {
        # this means LMS and LML are not compatible
        # keep only until last epSizeLong available
        # redefine metashort length
        LMS = expected_LML * epSizeLong / epSizeShort
        D = D[1:LMS,]
        LML = expected_LML
      } 
      
      #create timeseries for metashort
      time_shortEp_num = seq(quartlystart, quartlystart + ((nrow(D) - 1) * epSizeShort), by = epSizeShort)
      time_longEp_num = seq(quartlystart, quartlystart + ((nrow(D) - 1) * epSizeShort), by = epSizeLong)
      
      if (LML * 15 > LMS) {
        time_longEp_num = time_longEp_num[1:floor(length(time_shortEp_num) / (epSizeLong/epSizeShort))]
        time_shortEp_num = time_shortEp_num[1:(length(time_longEp_num) * (epSizeLong/epSizeShort))]
      }
      starttime = as.POSIXlt(time_longEp_num[1], origin = "1970-1-1", tz = tz)
      time_shortEp_8601 = POSIXtime2iso8601(x = as.POSIXlt(time_shortEp_num, tz = tz,
                                                           origin = "1970-01-01"),
                                            tz = tz)
      time_longEp_8601 = POSIXtime2iso8601(x = as.POSIXlt(time_longEp_num, tz = tz, origin = "1970-01-01"),
                                           tz = tz)
      if (params_general[["dataFormat"]] != "actigraph_csv") {
        M$metashort = data.frame(timestamp = time_shortEp_8601,
                                 accmetric = D[1:length(time_shortEp_8601),1],stringsAsFactors = FALSE)
      } else if (params_general[["dataFormat"]] == "actigraph_csv") {
        M$metashort = as.data.frame(cbind(time_shortEp_8601,
                                          D[1:length(time_shortEp_8601), ]))
        colnames(M$metashort) = c("timestamp", colnames(D))
      }
      
      if (length(which(is.na(M$metashort$ZCY) == TRUE)) > 0) {
        # impute missing data by zero
        # if it is a lot then this will be detected as non-wear
        M$metashort$ZCY[is.na(M$metashort$ZCY)] = 0 
      }
      
      LML = length(time_longEp_8601)
      if (params_general[["dataFormat"]] == "ukbiobank_csv") {
        names(M$metashort)[2] = "LFENMO"
        #Collapse second column of input data to use as non-wear score
        imp = D[, 2]
        M$metashort$LFENMO = M$metashort$LFENMO/1000
      } else if (length(grep(pattern = "actiwatch", x = params_general[["dataFormat"]], ignore.case = TRUE)) > 0) {
        imp = unlist(D[, 1])
      } else if (params_general[["dataFormat"]] == "actigraph_csv") {
        imp = unlist(D[, 1])
      }
      navalues = which(is.na(imp) == TRUE)
      if (length(navalues) > 0) imp[navalues] = 1
      
      
      if (params_general[["dataFormat"]] == "ukbiobank_csv") {
        # Take long epoch mean of UK Biobank based invalid data indicater per 5 seconds
        imp2 = cumsum(imp)
        imp3 = diff(imp2[seq(1, length(imp2),
                             by = ((60/epSizeShort) * (epSizeLong/60)))]) / ((60/epSizeShort) * (epSizeLong/60)) # rolling mean
        imp4 = round(imp3 * 3) # create three level nonwear score from it, not really necessary for GGIR, but helps to retain some of the information
      } else if (length(grep(pattern = "actiwatch", x = params_general[["dataFormat"]], ignore.case = TRUE)) > 0 |
                 params_general[["dataFormat"]] == "actigraph_csv") {
        # Using rolling 60 minute sum to indicate whether it is nonwear
        imp2 = zoo::rollapply(imp, width = (1*3600) / epSizeShort, FUN = sum, fill = 0)
        imp4 = imp2
        imp4[which(imp2 > 0)] = 0
        imp4[which(imp2 == 0)] = 3 # If rolling average is zero then consider it nonwear
        imp5 = cumsum(imp4)
        step = (60/epSizeShort) * (epSizeLong/60)
        imp4 = diff(imp5[seq(1, length(imp5) + step, by = step)]) / step # rolling mean
        if (length(imp4) > length(time_longEp_8601)) imp4 = imp4[1:length(time_longEp_8601)]
        imp4 = round(imp4)
        if (any(is.na(imp4))) imp4[which(is.na(imp4))] = 3
      }
      
      if (length(imp4) < LML) {
        imp4 = c(imp4, rep(0, LML - length(imp4)))
      }
      # create data.frame for metalong, note that light and temperature are just set at zero
      M$metalong = data.frame(timestamp = time_longEp_8601,nonwearscore = imp4, #rep(0,LML)
                              clippingscore = rep(0,LML), lightmean = rep(0,LML),
                              lightpeak = rep(0,LML), temperaturemean = rep(0,LML), EN = rep(0,LML))
      
      # update weekday name and number based on actual data
      M$wdayname = weekdays(x = starttime, abbreviate = FALSE)
      M$wday = as.POSIXlt(starttime)$wday + 1
      
      # Save these files as new meta-file
      filefoldername = NA
      save(M, C, I, filename_dir, filefoldername,
           file = outputFileName)
      Sys.setlocale("LC_TIME", Syslocale) # set language to English
    }
  }
}
