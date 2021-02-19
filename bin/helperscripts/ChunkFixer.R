args <- commandArgs(trailingOnly = TRUE)

DivideSNPsPerChunks <- function(InputPath, NeededChunkSize) {
  library(data.table)
  library(stringr)

  # data.table uses 4 threads
  setDTthreads(4)

  # String of files:
  Files <- list.files(InputPath)

  # Sort naturally
  Files <- str_sort(Files, numeric = TRUE)

  NeededChunkSize <- as.numeric(NeededChunkSize)
  FileIndex <- 0
  NrOfFiles <- length(list.files(InputPath))

  print(Files)
  print(NrOfFiles)
  print(NeededChunkSize)

  for (f in Files[-NrOfFiles]) {
    FileIndex <- FileIndex + 1
    CheckRows <- as.numeric(str_trim(system(paste0("cat ", InputPath, "/", f, " | wc -l"), intern = TRUE)))

    if (!file.exists(paste0(InputPath, "/", f))) {
      message("Program finishes...")
      break
    }
    if (CheckRows == NeededChunkSize) {
      message(paste0("Input file ", f, " already has specified size."))
      next
    }
    if (CheckRows < NeededChunkSize) {
      message(paste0("Fixing ", f, ", missing ", NeededChunkSize - CheckRows, " rows."))
      MissingRows <- NeededChunkSize - CheckRows
      orig_file <- fread(paste0(InputPath, "/", f), header = FALSE)

      # Pre-allocate empty data.table for adding missing rows
      MoveChunk <- orig_file[-c(1:nrow(orig_file))]

      NextFileIndex <- FileIndex

      while (MissingRows > 0 & NextFileIndex <= (NrOfFiles - 1)) {
        NextFileIndex <- NextFileIndex + 1

        # If next file from current has not been removed, continue with program
        NextFile <- fread(paste0(InputPath, "/", Files[NextFileIndex]), header = FALSE)
        NrowFile <- nrow(NextFile)

        if (NrowFile >= MissingRows) {
          message(paste0("Taking ", MissingRows, " rows from: ", Files[NextFileIndex]))
          MoveChunk <- rbindlist(list(MoveChunk, NextFile[1:MissingRows]))
          file.remove(paste0(InputPath, "/", Files[NextFileIndex]))
          fwrite(NextFile[-c(1:MissingRows)], paste0(InputPath, "/", Files[NextFileIndex]), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
          MissingRows <- MissingRows - nrow(NextFile)
          rm(NextFile)
          gc()
        }

        if (NrowFile < MissingRows) {
          if (nrow(NextFile) == 0) {
            message(paste0("There are ", nrow(NextFile), " rows in: ", Files[NextFileIndex], ", skipping."))
            next
          } else {
            message(paste0("Taking ", nrow(NextFile), " rows from: ", Files[NextFileIndex]))
            MoveChunk <- rbindlist(list(MoveChunk, NextFile[1:nrow(NextFile), ]))
            # Remove next file
            file.remove(paste0(InputPath, "/", Files[NextFileIndex]))
            # Write empty file as a placeholder
            fwrite(NextFile[-c(1:nrow(NextFile))], paste0(InputPath, "/", Files[NextFileIndex]), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
            MissingRows <- MissingRows - nrow(NextFile)
            rm(NextFile)
            gc()
          }
        }
      }
      message(paste("Chunk size to move:", nrow(MoveChunk)))
      if (nrow(MoveChunk) < MissingRows) {
        message("No SNPs in remaining files. Removing remaining empty files...")
        for (empty_file in Files[FileIndex:NrOfFiles]) {
          file.remove(paste0(InputPath, "/", empty_file))
        }
      }

      orig_file <- rbindlist(list(orig_file, MoveChunk))

      if (nrow(orig_file) > 0) {
        file.remove(paste0(InputPath, "/", f))
        fwrite(orig_file, paste0(InputPath, "/", f), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        rm(orig_file)
        gc()
      }
    }
    if (CheckRows > NeededChunkSize) {
      stop("Something is wrong, there should not be any chunks larger than specified chunk size! Aborting...")
    }
    }
}

suppressWarnings(DivideSNPsPerChunks(args[1], args[2]))
