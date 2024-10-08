rename_fastqs <- function(samplesheet, fastq_dir, out_dir) {
    sheet <- read.delim(samplesheet, sep="\t")

    ifelse(!dir.exists(out_dir), dir.create(file.path(out_dir), recursive = TRUE), FALSE)

    for (samp in sheet$sample) {
        print(paste("Renaming PE reads for:", samp))
        
        # Rename R1 reads
        file_1 <- paste0(fastq_dir, "/", unlist(strsplit(subset(sheet, sample == samp)$fq1, split = ",")))
        print(paste("Found", length(file_1), "R1 files for", samp))
        
        for (i in 1:length(file_1)) {
            if (file.exists(file_1[i])) {
                new_name <- paste0(out_dir, "/", samp, "-", i, "-R1.fastq.gz")
                system(paste0("ln -sr ", file_1[i], " ", new_name))
                print(paste(file_1[i], "successfully symlinked to", new_name))
            } else {
                stop(paste("Error:", file_1[i], "not found in", fastq_dir))
            }
        }
        
        # Rename R2 reads
        if (is.na(subset(sheet, sample == samp)$fq2) | subset(sheet, sample == samp)$fq2 == "") {
            print(paste("No R2 reads found for", samp))
            next
        }
        
        file_2 <- paste0(fastq_dir, "/", unlist(strsplit(subset(sheet, sample == samp)$fq2, split = ",")))
        print(paste("Found", length(file_2), "R2 files for", samp))
        
        for (i in 1:length(file_2)) {
            if (file.exists(file_2[i])) {
                new_name <- paste0(out_dir, "/", samp, "-", i, "-R2.fastq.gz")
                system(paste0("ln -sr ", file_2[i], " ", new_name))
                print(paste(file_2[i], "successfully symlinked to", new_name))
            } else {
                stop(paste("Error:", file_2[i], "not found in", fastq_dir))
            }
        }
    }
}

rename_fastqs(snakemake@params[["samplesheet"]], snakemake@params[["fastq_dir"]], snakemake@output[["symlink_dir"]])
