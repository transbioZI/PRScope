# Calculate polygenic scores using ldpred2
# this script is an adaptation of the demo script available at the bigsnpr homepage
library(bigsnpr, quietly = T)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(tools)
library(argparser, quietly=T)
library(stringr)
library(data.table)
library(doParallel)

### Maybe there's some environment variable availble to determine the location of the script instead
coms <- commandArgs()
coms <- coms[substr(coms, 1, 7) == '--file=']
dirScript <- dirname(substr(coms, 8, nchar(coms)))
source(paste0(dirScript, '/fun.R'))

par <- arg_parser('Calculate polygenic scores using ldpred2')
# Mandatory arguments (files)
par <- add_argument(par, "--geno-file-rds", help="Input .rds (bigSNPR) file with genotypes. A similarly named .bk (backing) file must exist in the same location")
par <- add_argument(par, "--sumstats", help="Input file with GWAS summary statistics")
par <- add_argument(par, "--out", help="Output file with calculated PGS")

# Optional files
par <- add_argument(par, "--out-merge", flag=T, help="Merge output with existing file.")
par <- add_argument(par, "--out-merge-ids", nargs=2, help='Pass ID columns in this order: [family ID] [individual ID]')
par <- add_argument(par, "--file-keep-snps", help="File with RSIDs of SNPs to keep")
par <- add_argument(par, "--ld-file", default="/ldpred2_ref/ldref_hm3_plus/LD_with_blocks_chr@.rds", help="LD reference files, split per chromosome; chr label should be indicated by '@' symbol")
par <- add_argument(par, "--ld-meta-file", default="/ldpred2_ref/map_hm3_plus.rds", help="list of variants in --ld-file")

# Genotype
par <- add_argument(par, "--geno-impute-zero", help="Set missing genotypes to zero.", flag=T)
# Sumstats file.
par <- add_argument(par, "--chr2use", help="list of chromosomes to use (by default it uses chromosomes 1 to 22)", nargs=Inf)
par <- add_argument(par, "--col-chr", help="CHR number column", default="CHR", nargs=1)
par <- add_argument(par, "--col-snp-id", help="SNP ID (RSID) column", default="SNP", nargs=1)
par <- add_argument(par, "--col-A1", help="Effective allele column", default="A1", nargs=1)
par <- add_argument(par, "--col-A2", help="Noneffective allele column", default="A2", nargs=1)
par <- add_argument(par, "--col-bp", help="SNP position column. Will be ignored if --merge-by-rsid flag is used", default="BP", nargs=1)
par <- add_argument(par, "--col-stat", help="Effect estimate column", default="BETA", nargs=1)
par <- add_argument(par, "--col-stat-se", help="Effect estimate standard error column", default="SE", nargs=1)
par <- add_argument(par, "--col-pvalue", help="P-value column", default="P", nargs=1)
par <- add_argument(par, "--col-n", help="Effective sample size. Override with --effective-sample-size", default="N", nargs=1)
par <- add_argument(par, "--stat-type", help="Effect estimate type (BETA for linear, OR for odds-ratio)", default="BETA", nargs=1)
par <- add_argument(par, "--effective-sample-size", help="Effective sample size, if unavailable in sumstats (--col-n)", nargs=1)
par <- add_argument(par, "--n-cases", help="Nr cases when phenotype is binary. For calculating effective sample size.", nargs=1)
par <- add_argument(par, "--n-controls", help="Nr controls when phenotype is binary. For calculating effective sample size.", nargs=1)
# Polygenic score
par <- add_argument(par, "--name-score", help="Set column name for the created score", nargs=1, default='score')
# Parameters to LDpred
par <- add_argument(par, "--hyper-p-length", help="Length of hyperparameter p sequence to use for --ldpred-mode auto", default=50)
par <- add_argument(par, "--hyper-p-max", help="Maximum (<1) of hyperparameter p sequence to use for --ldpred-mode auto", default=0.2)
# Others
par <- add_argument(par, "--ldpred-mode", help='Ether "auto" or "inf" (infinitesimal)', default="inf")
par <- add_argument(par, "--cores", help="Number of CPU cores to use, otherwise use the available number of cores minus 1", default=nb_cores())
par <- add_argument(par, '--set-seed', help="Set a seed for reproducibility", nargs=1)
par <- add_argument(par, "--merge-by-rsid", help="Merge using rsid (the default is to merge by chr:bp:a1:a2 codes).", flag=TRUE)
par <- add_argument(par, "--genomic-build", help="Genomic build to use. Either hg19, hg18 or hg38", default="hg19", nargs=1)
par <- add_argument(par, "--tmp-dir", help="Directory to store temporary files. Default is output of base::tempdir()", default=tempdir())

parsed <- parse_args(par)

### Mandatory
fileGeno <- parsed$geno_file_rds
fileSumstats <- parsed$sumstats
fileOutput <- parsed$out
fileOutputPlot <- fileOutput # diagnostic plot
fileOutputHer <- fileOutput # diagnostic plot
fileOutputLDReg <- fileOutput
fileLD <- parsed$ld_file
fileMetaLD <- parsed$ld_meta_file

### Optional
fileKeepSNPs <- parsed$file_keep_snps
fileOutputMerge <- parsed$out_merge
fileOutputMergeIDs <- parsed$out_merge_ids
# Vectors as defaults causes a warning for argparse, so the default is set this way instead
if (fileOutputMerge & isVarNA(fileOutputMergeIDs)) fileOutputMergeIDs <- COLNAMES_ID_PLINK
### Genotype
genoImputeZero <- parsed$geno_impute_zero

# Sumstats file
chr2use <- parsed$chr2use
if (any(is.na(chr2use))) chr2use <- 1:22

colChr <- parsed$col_chr
colSNPID <- parsed$col_snp_id
colA1 <- parsed$col_A1
colA2 <- parsed$col_A2
colBP <- parsed$col_bp
colStat <- parsed$col_stat
colStatSE <- parsed$col_stat_se
colPValue <- parsed$col_pvalue
colN <- parsed$col_n
mergeByRsid <- parsed$merge_by_rsid

# unset colBP in case of merging by rsid
if (mergeByRsid) {
  cat(paste('The --merge-by-rsid flag is used; --col-bp', colBP, 'will be ignored\n'))
  colBP <- NULL
}

# Parameters to LDpred
parHyperPLength <- parsed$hyper_p_length
parHyperPMax <- parsed$hyper_p_max
# Others
argEffectiveSampleSize <- parsed$effective_sample_size
# These two have to be accessed differently due to some argparser bug
argNCases <- parsed$`n-cases`
argNControls <- parsed$`n_controls`

argLdpredMode <- parsed$ldpred_mode
validModes <- c('inf', 'auto')
if (!argLdpredMode %in% validModes) stop("--ldpred-mode should be one of: ", paste0(validModes, collapse=', '))
argStatType <- parsed$stat_type
if (!argStatType %in% c('BETA', 'OR')) stop('--stat-type should be one of "BETA" or "OR"')
setSeed <- parsed$set_seed

# These vectors are used to convert headers in the sumstat files to those
# used by bigsnpr
if (mergeByRsid) {
  colSumstatsOld <- c(colChr, colSNPID, colA1, colA2, colStat, colStatSE)
  colSumstatToGeno <- c("chr",  "rsid",  "a1",  "a0",  "beta",  "beta_se")
} else {
  colSumstatsOld <- c(colChr, colSNPID, colBP, colA1, colA2, colStat, colStatSE)
  colSumstatToGeno <- c("chr",  "rsid",  "pos",  "a1",  "a0",  "beta",  "beta_se")
}

cat('Loading backingfile:', fileGeno ,'\n')
obj.bigSNP <- snp_attach(fileGeno)

# Store some key variables
G <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
NCORES <- parsed$cores

# Check genotype data for missingness
if (genoImputeZero) {
  cat('### Imputing missing genotypes with zero\n')
  G <- zeroMissingGenotypes(G)
}

reference_LD_list = list()
for (chr in chr2use) {
  fileLD_chr <- str_replace(fileLD, "@", toString(chr))
  reference_LD_list[[chr]] <- readRDS(fileLD_chr)
}

studies = fread("",data.table = F)
tmp_path = "" # todo change
fileOutput = "" # todo change
res_result = "" # todo change
number_of_threads = 20

map_ldref_hg19 <- readRDS(fileMetaLD)
map_ldref_hg38 <- readRDS(fileMetaLD)
map_ldref_hg38$pos <- map_ldref_hg38$pos_hg38
map_ldref_hg38$pos_hg38 <- NULL

cl <- makeCluster(number_of_threads)
registerDoParallel(cl)

res = foreach(st=1:nrow(studies)) %dopar% {
    genomic_build = studies$genome_build[st]
    if(genomic_build == "hg37") {
        map_ldref = map_ldref_hg19
    } else {
        map_ldref = map_ldref_hg38
    }
    results_fam = obj.bigSNP$fam
    nameScore = studies$study_id[st]
    NEFF = studies$neff[st]
    dirPlot <- dirname(fileOutput)
    fileOutputPlot <- paste0(dirPlot, '/', nameScore, '.png')
    fileOutputResults <- paste0(dirPlot, '/', nameScore, '.tsv')
    fileSumstats = paste0("path_to_gwas/",nameScore,"/ldpred/",nameScore,".qced.h.tsv.gz")

    cat('\n### Reading summary statistics', fileSumstats,'\n')
    sumstats <- bigreadr::fread2(fileSumstats)
    cat('Loaded', nrow(sumstats), 'SNPs\n')
    # Rename columns in bigSNP object
    colMap <- c('chr', 'rsid', 'pos', 'a1', 'a0')
    map <- setNames(obj.bigSNP$map[-3], colMap)
    tmpdir = paste0(tmp_path,"/",nameScore)

    if (!file.exists(tmpdir)){
        dir.create(tmpdir)
    }

    h2_est <- studies$observed_h2[st]
    tmp_file <- tempfile(tmpdir=tmpdir)
    # Rename headers in sumstats file if necessary
    colSumstats <- colnames(sumstats)
    # Lowercase them
    colSumstats <- tolower(colSumstats)
    colSumstatsOld <- tolower(colSumstatsOld)
    # Check that all necessary columns are present
    colReplacements <- match(colSumstatsOld, colSumstats)
    if (sum(is.na(colReplacements)) > 0) {
        cat('The following necessary columns could not be found in', basename(fileSumstats), '\n')
        for (i in 1:length(colReplacements)) {
            if (is.na(colReplacements[i])) {
            cat('\t', colSumstatsOld[i],'\n')
            }
        }
        cat('Columns in ', basename(fileSumstats), ': ', paste0(colSumstats, collapse=' | '), '\n', sep='')
        stop('Necessary columns not found')
    }

    # Replace columns in sumstat data so the match bigSNP
    colSumstats[colReplacements] <- colSumstatToGeno
    colnames(sumstats) <- colSumstats
    sumstats$a0 <- toupper(sumstats$a0)
    sumstats$a1 <- toupper(sumstats$a1)

    # Check that there are no characters in chromosome column (causes snp_match to fail)
    if (!isOnlyNumeric(sumstats$chr)) {
        cat('Removing rows with non-integers in column', colChr, '\n')
        numeric <- getNumericIndices(sumstats$chr)
        cat('Removing', nrow(sumstats) - length(numeric), 'SNPs...\n')
        sumstats <- sumstats[numeric,]
        cat('Retained', nrow(sumstats), 'SNPs\n')
        sumstats$chr <- as.numeric(sumstats$chr)
    }
    sumstats$n_eff = NEFF

    if (!is.na(fileKeepSNPs)) {
        cat('Filtering SNPs using', fileKeepSNPs, '\n')
        keepSNPs <- read.table(fileKeepSNPs)
        nSNPsBefore <- nrow(sumstats)
        sumstats <- sumstats[sumstats$rsid %in% keepSNPs[,1],]
        nSNPsAfter <- nrow(sumstats)
        cat('Retained', nSNPsAfter, 'out of', nSNPsBefore,'\n')
    }

    cat('Filtering SNPs based on --chr2use\n')
    nSNPsBefore <- nrow(sumstats)
    sumstats <- sumstats[sumstats$chr %in% chr2use,]
    nSNPsAfter <- nrow(sumstats)
    cat('Retained', nSNPsAfter, 'out of', nSNPsBefore,'\n')
    # match sumstats to genotypes
    # (df_beta is quite an ugly name for GWAS sumstats, but let it be so
    # for consistency with original LDpred2 tutorial)
    cat('Matching sumstats to genotypes\n')
    df_beta <- snp_match(sumstats, map, join_by_pos=!mergeByRsid, match.min.prop=0)
    drops <- c("_NUM_ID_.ss", "_NUM_ID_", "rsid.ss")
    df_beta <- df_beta[ , !(names(df_beta) %in% drops)]

    cat('Matching sumstats to LD reference\n')
    df_beta <- snp_match(df_beta, map_ldref, join_by_pos=!mergeByRsid, match.min.prop=0)  # this adds af_UKBB and ld columns
    drops <- c("_NUM_ID_.ss", "rsid.ss", 'block_id', 'pos_hg18', 'pos_hg38')
    df_beta <- df_beta[ , !(names(df_beta) %in% drops)]

    cat('\n### Loading LD reference from ', fileLD, '\n')

    ld_size <- 0; corr <- NULL
    for (chr in chr2use) {
    ## indices in 'df_beta' corresponding to a particular 'chr'
        ind.chr <- which(df_beta$chr == chr)
        if (length(ind.chr) == 0) next
        ## indices in 'map_ldref'
        ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
        ## indices in 'corr_chr'
        ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

        num_ldref_snps <- sum(map_ldref$chr == chr)
        ld_size <- ld_size + num_ldref_snps

        corr_chr <- reference_LD_list[[chr]][ind.chr3, ind.chr3]

        if (is.null(corr)) {
            corr <- as_SFBM(corr_chr, tmp_file, compact = TRUE)
        } else {
            corr$add_columns(corr_chr, nrow(corr))
        }
    }

    ldsc <- with(df_beta, snp_ldsc(ld, ld_size, chi2=(beta/beta_se)^2, sample_size=NEFF, blocks = NULL, ncores=NCORES))
    h2_est <- ldsc[["h2"]]
    if(h2_est < 0) {
        cat('\n### h2_init < 0, calculation finished with an empty result\n')
        results_fam[,nameScore] = rep(NA,dim(results_fam)[1])
    } else {
        cat('Running LDPRED2 auto model\n')
        if (!is.na(setSeed)) set.seed(setSeed)

        sh_cor = sort(runif(dim(map_ldref)[1], min = 0.6, max = 0.99))[dim(df_beta)[1]]

        multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init=h2_est, vec_p_init=seq_log(1e-4, parHyperPMax, length.out=parHyperPLength),
                                    allow_jump_sign=F, shrink_corr=sh_cor, ncores=NCORES, burn_in = 500, num_iter = 250)

        #multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init=h2_est, vec_p_init=seq_log(1e-4, parHyperPMax, length.out=parHyperPLength),
        #                           allow_jump_sign=F, shrink_corr=0.95, ncores=NCORES, burn_in = 800, num_iter = 400)

        cat('Plotting diagnostics: ', fileOutputPlot, '\n', sep='')
        library(ggplot2)
        auto <- multi_auto[[1]]
        dta <- data.frame(path_p_est=auto$path_p_est, path_h2_est=auto$path_h2_est, x=1:length(auto$path_p_est))
        plt <- plot_grid(
        ggplot(dta, aes(y=path_p_est, x=x)) + geom_point() + theme_bigstatsr() +
            geom_hline(aes(yintercept=auto$p_est), col="blue") +
            scale_y_log10() + labs(y="p"),
        ggplot(dta, aes(y=path_h2_est, x=x)) + geom_point() + theme_bigstatsr() +
            geom_hline(aes(yintercept=auto$h2_est), col="blue") + labs(y="h2"),
        ncol=1, align="hv"
        )
        ggsave(plt, file=fileOutputPlot)
        cat('Filtering chains\n')
        range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
        # Keep chains that pass the filtering below
        keep <- (range > (0.95 * quantile(range, 0.95, na.rm = T)))
        nas <- sum(is.na(keep))
        if (nas > 0) cat('Omitting', nas, 'chains out of', length(keep), ' due to missing values in correlation range\n')
        keep[is.na(keep)] <- F
        beta <- rowMeans(as.data.frame(sapply(multi_auto[keep], function (auto) auto$beta_est)))

        h2_est_med <- sapply(multi_auto[keep], function(auto) auto$h2_est)[keep]
        p_est <-sapply(multi_auto[keep], function(auto) auto$p_est)[keep]

        cat('Scoring all individuals...\n')
        # find which SNPs to use, and whether we need to flip their sign

        map_pgs <- df_beta[1:4]; map_pgs$beta <- 1
        map_pgs2 <- snp_match(map_pgs, map, join_by_pos=!mergeByRsid, match.min.prop=0)

        tryCatch(
            pred_all <- big_prodVec(G, beta * map_pgs2$beta, ind.col=map_pgs2[['_NUM_ID_']], ncores=NCORES),
            error=function(er) {
                cat('bigstatsr::big_prodVec threw an error:\n')
                message(er)
                cat('\n\nErrors regarding "missingness in X" may be solved by imputing genotype data or passing --geno-impute-zero\n')
                er
                quit(status=1, save='no')
        },
            warning=function(er) {
                cat('bigstatsr::big_prodVec threw a warning:\n')
                message(er)
            }
        )

        results_fam[,nameScore] = pred_all
    }

    write.table(results_fam, fileOutputResults, row.names = F, col.names = F, quote = F, sep = "\t")
    file.remove(paste0(tmp_file, '.sbk'))
}

writeLines(as.character(res), res_result)
stopCluster(cl)
