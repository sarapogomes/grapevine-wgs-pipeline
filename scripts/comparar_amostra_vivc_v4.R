# =============================================================
# COMPARAÇÃO DE PERFIL SNP COM BASE DE DADOS VIVC
# =============================================================

rm(list = ls())  # limpa variáveis de sessões anteriores
library(readxl)
library(writexl)

# =============================================================
# CONFIGURAÇÃO — editar apenas esta secção
# =============================================================

FICHEIRO_AMOSTRA <- "C:/Users/User/Desktop/apresentação/TN.xlsx"
NOME_AMOSTRA     <- "TN"   # igual ao nome do separador no Excel
FICHEIRO_DB      <- "C:/Users/User/Desktop/apresentação/vivc_snp_db.txt"
PASTA_OUTPUT     <- "C:/Users/User/Desktop/apresentação/"

# =============================================================
# PASSO 1 — Ler os dados
# =============================================================

cat("\n[1/4] A ler ficheiros...\n")

db <- read.delim(FICHEIRO_DB, sep = "\t", header = TRUE,
                 stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("      Base de dados: %d variedades, %d colunas\n", nrow(db), ncol(db)))

amostra <- read_excel(FICHEIRO_AMOSTRA, sheet = NOME_AMOSTRA)
cat(sprintf("      Amostra %s: %d SNPs\n", NOME_AMOSTRA, nrow(amostra)))

# =============================================================
# PASSO 2 — Preparar os dados
# =============================================================

cat("\n[2/4] A preparar dados...\n")

meta_cols   <- c("snp_id", "akzessionsname", "kenn_nr", "acc_number",
                 "gespeichert_am", "stempel")
snp_cols_db <- setdiff(colnames(db), meta_cols)

gt_amostra  <- setNames(amostra$GT_CORRECT_2, amostra$SNP_ID)

snps_comuns <- intersect(snp_cols_db, names(gt_amostra))
n_snps      <- length(snps_comuns)
cat(sprintf("      SNPs em comum: %d\n", n_snps))

# =============================================================
# PASSO 3 — Função de normalização
# =============================================================

normalizar_gt <- function(gt) {
  gt <- trimws(as.character(gt))
  if (is.na(gt) || gt %in% c("", "--", "nd", "NA")) return("")
  paste(sort(strsplit(gt, "")[[1]]), collapse = "")
}

gt_amostra_norm <- sapply(gt_amostra[snps_comuns], normalizar_gt)

# =============================================================
# PASSO 4 — Comparar com cada variedade
# =============================================================

cat(sprintf("\n[3/4] A comparar com %d variedades...\n", nrow(db)))

match_matrix <- matrix(FALSE, nrow = nrow(db), ncol = n_snps,
                       dimnames = list(NULL, snps_comuns))
score_vec    <- integer(nrow(db))
md_vec       <- integer(nrow(db))

for (i in seq_len(nrow(db))) {
  score    <- 0
  md_count <- 0

  for (j in seq_along(snps_comuns)) {
    snp   <- snps_comuns[j]
    gt_db <- normalizar_gt(db[i, snp])
    gt_am <- gt_amostra_norm[snp]

    if (gt_db == "" || gt_am == "") {
      md_count <- md_count + 1
    } else if (gt_db == gt_am) {
      score <- score + 1
      match_matrix[i, j] <- TRUE
    }
  }

  score_vec[i] <- score
  md_vec[i]    <- md_count
}

pct_md_vec <- md_vec / n_snps

# =============================================================
# PASSO 5 — Montar e guardar o resultado
# =============================================================

resultado <- data.frame(
  score      = score_vec,
  md_count   = md_vec,
  pct_md     = pct_md_vec,
  accsession = db[["akzessionsname"]],
  match_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

resultado <- resultado[order(-resultado$score), ]
rownames(resultado) <- NULL

cat("\n[4/4] A guardar ficheiros...\n")

ficheiro_output <- paste0(PASTA_OUTPUT, NOME_AMOSTRA, "_results.txt")

write.table(resultado,
            file      = ficheiro_output,
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)

cat(sprintf("      .txt guardado: %s\n", ficheiro_output))

# =============================================================
# PASSO 6 — Exportar para Excel
# =============================================================

ficheiro_excel <- paste0(PASTA_OUTPUT, NOME_AMOSTRA, "_results.xlsx")

# Separar em duas sheets: resumo (top 50) e resultado completo
resumo <- resultado[1:min(50, nrow(resultado)),
                    c("score", "md_count", "pct_md", "accsession")]

write_xlsx(
  list(
    "Top50"    = resumo,
    "Completo" = resultado
  ),
  path = ficheiro_excel
)

cat(sprintf("      .xlsx guardado: %s\n", ficheiro_excel))
cat(sprintf("\n✓ Concluído! Ficheiros guardados em: %s\n", PASTA_OUTPUT))

cat("\nTop 10 resultados:\n")
print(resultado[1:10, c("score", "md_count", "pct_md", "accsession")],
      row.names = FALSE)
