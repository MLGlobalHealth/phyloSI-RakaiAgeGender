require(data.table)
require(here)

#########
# PATHS #
#########

usr <- Sys.info()[["user"]]
indir <- here::here()

if (usr == "andrea") {
    indir.deepsequencedata <- "/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live"
    indir.deepanalyses <- "/home/andrea/HPC/project/ratmann_deepseq_analyses/live"
    indir.deepanalyses_xiaoyue <- "/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI"
}

# outdir is actually in the git repository
outdir.data <- file.path(indir, "data")

file.path.tsiestimates <- file.path(
    indir.deepanalyses,
    "PANGEA2_RCCS_MRC_UVRI_TSI",
    "2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_001_rla_T_zla_T",
    "aggregated_TSI_with_estimated_dates.csv"
)

cat("\n---- Extract non-confidential columns ----\n")

cols <- c( "AID", "RENAME_ID",
    "RF_pred_sqrt", "RF_std", "RF_cc025", "RF_cc975", "RF_pred_MAE",
    "RF_pred_linear", "RF_pred_min_linear", "RF_pred_max_linear")
dtsi <- fread(file.path.tsiestimates, select = cols)


cat("\n---- Save in gitdir.data ----\n")
filename <- file.path(outdir.data, "TSI_estimates.csv")
fwrite(dtsi, file=filename)
