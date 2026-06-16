# ======== 1. 加载依赖包 ========
suppressPackageStartupMessages({
  library(data.table)
  library(DNAcopy)
})

# ======== 2. 参数设置 ========
bdg_dir    <- "C:/Users/10260/Desktop/AA_Osat_jap_bdg_files"   # bdg目录
output_dir <- "AA_Osat_jap_segmentation"                        # 输出目录
undo_sd    <- 1

if (!dir.exists(output_dir)) dir.create(output_dir)

# 从文件名自动获取所有染色体列表（取第一个bdg文件的chr列）
sample_files <- list.files(bdg_dir, pattern = "\\.bdg$", full.names = TRUE)
if (length(sample_files) == 0) stop("❌ 未找到任何 .bdg 文件，请检查 bdg_dir 路径")

cat(">>> 读取染色体列表...\n")
chroms <- unique(fread(sample_files[1], col.names = c("chr", "start", "end", "log2ratio"))$chr)
cat(">>> 共检测到", length(chroms), "条染色体：", paste(chroms, collapse = ", "), "\n")

# ======== 3. 循环每条染色体输出 TSV ========
for (cur_chr in chroms) {
  cat("\n>>> 正在处理：", cur_chr, "\n")

  for (hap in c("hap1", "hap2")) {
    for (rep in 1:3) {
      fname <- sprintf("AA_Osat_jap_%s.sample%d.CENH3.bdg", hap, rep)
      fpath <- file.path(bdg_dir, fname)

      if (!file.exists(fpath)) {
        warning("❌ 缺失文件：", fpath)
        next
      }

      # 读取并筛选目标染色体
      dat <- fread(fpath, col.names = c("chr", "start", "end", "log2ratio"))
      dat <- dat[chr == cur_chr]
      if (nrow(dat) == 0) {
        warning("⚠️ 文件中未找到目标染色体 ", cur_chr, " 数据：", fpath)
        next
      }
      dat[, pos := floor((start + end) / 2)]
      cat(fname, " 数据点数：", nrow(dat), "\n")

      # CBS 分段
      cna    <- CNA(genomdat  = dat$log2ratio,
                    chrom     = rep(1, nrow(dat)),
                    maploc    = dat$pos,
                    data.type = "logratio",
                    sampleid  = fname)
      cna_sm <- smooth.CNA(cna)
      seg    <- segment(cna_sm, undo.splits = "sdundo", undo.SD = undo_sd, verbose = 0)
      seg_df <- as.data.frame(seg$output)
      names(seg_df) <- gsub("\\.", "_", names(seg_df))

      # 将 chrom 列还原为真实染色体名
      if ("chrom" %in% names(seg_df)) {
        seg_df$chrom <- cur_chr
      } else if ("chromosome" %in% names(seg_df)) {
        seg_df$chromosome <- cur_chr
      }

      # 输出 segmentation TSV
      seg_file <- file.path(output_dir,
                            sprintf("AA_Osat_jap_%s.%s.rep%d.segmentation.tsv", cur_chr, hap, rep))
      fwrite(seg_df, seg_file, sep = "\t")
      cat("    ▶ 已写出 TSV:", basename(seg_file), "\n")
    }
  }
}

cat("\n✅ 所有染色体 CBS 分段 TSV 输出完成！\n")
