# VDZ_ADA_Vascular_SPP1_lineplots.R
# Objective: plot Olink NPX expression and ssGSEA signature trajectories for
# Adalimumab/Vedolizumab responders vs non-responders (raw, change, percent change);

# paths are defined up front for clarity.

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
project_dir <- getwd()
fig2_dir <- file.path(project_dir, "Fig2")
olink_data_dir <- file.path(fig2_dir, "Olink_data")
dir.dat <- olink_data_dir
dir.code <- file.path(fig2_dir, "code")
dir.res <- file.path(fig2_dir, "results")
dir.figures <- file.path(dir.res, "figures")

dir.create(dir.code, recursive = TRUE, showWarnings = FALSE)
dir.create(dir.res, recursive = TRUE, showWarnings = FALSE)
dir.create(dir.figures, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# Libraries and helper functions
# ------------------------------------------------------------------
library(DT)
print_datatable <- function(data, caption, col_names = colnames(data), digits = 4, rownames = FALSE) {
  require(DT)
  for (i in 1:ncol(data)) {
    if (is.numeric(data[, i])) {
      if (!is.integer(abs(data[, i]))) {
        data[, i][abs(data[, i]) >= 0.001 & !is.na(abs(data[, i]) >= 0.001)] <- round(data[, i][abs(data[, i]) >= 0.001 & !is.na(abs(data[, i]) >= 0.001)], 3)
        data[, i][abs(data[, i]) < 0.001 & !is.na(abs(data[, i]) < 0.001)] <- formatC(data[, i][abs(data[, i]) < 0.001 & !is.na(abs(data[, i]) < 0.001)], format = "e", digits = 1)
      }
    }
  }
  if (knitr::is_html_output()) {
    datatable(data, caption = caption, rownames = rownames,
      extensions = c("Buttons"),
      options = list(scrollX = TRUE, scrollcollapse = TRUE, dom = "Bfrtip", buttons = c("copy", "csv", "excel"))
    )
  } else {
    if (nrow(data) > 10) {
      sel <- 1:10
    } else {
      sel <- 1:nrow(data)
    }
    knitr::kable(data[sel, ], format = "pandoc", digits = digits, caption = caption, col.names = col_names, row.names = rownames)
  }
}

knitr::opts_chunk$set(echo = TRUE, eval = TRUE)
DT::datatable(NULL)

tidy.dataframe <- function(df) {
  for (i in 1:ncol(df)) {
    if (is.numeric(df[, i])) {
      if (!is.integer(abs(df[, i]))) {
        df[, i][(!is.na((df[, i])) & (abs(df[, i]) >= 0.001))] <- round(df[, i][(!is.na((df[, i])) & (abs(df[, i]) >= 0.001))], 3)
        df[, i][!is.na((df[, i])) & (abs(df[, i]) < 0.001)] <- formatC(df[, i][!is.na((df[, i])) & (abs(df[, i]) < 0.001)], format = "e", digits = 1)
      }
    }
  }
  return(df)
}

rowCallback <- c(
  "function(row, data){",
  "  for(var i=0; i<data.length; i++){",
  "    if(data[i] === null){",
  "      $('td:eq('+i+')', row).html('-')",
  "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
  "    }",
  "  }",
  "}"
)

render <- JS(
  "function(data, type, row) {",
  "  if(type === 'sort' && data === null) {",
  "    return 999999;",
  "  }",
  "  return data;",
  "}"
)
knitr::opts_chunk$set(echo = TRUE, eval = TRUE)

library(here)
library(tidyverse)
library(openxlsx)
library(edgeR)
library(ggfortify)
library(ggpubr)
library(reshape2)
library(glue)
library(grid)
library(GSVA)

# ------------------------------------------------------------------
# Genes and signatures
# ------------------------------------------------------------------
prots_interest <- c("SPP1")
signatures <- list(
  Vascular_inflammation = c(
    "SPP1", "STC1", "TIMP1", "CXCL6", "OLR1", "THBD",
    "PGLYRP1", "PRG2", "FSTL1", "KCNJ8", "FGFR1", "PF4",
    "CCND2", "MSX1"
  )
)

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
md <- readxl::read_xlsx(file.path(dir.dat, "olink samples.xlsx"), sheet = 1)
response <- read.table(file.path(dir.dat, "outcome_response.csv"), header = TRUE, sep = ",")

md_clinical_all <- read.table(file.path(dir.dat, "PX_md_clinical.txt"), sep = "\t", header = TRUE)
md_clinical_all$SampleID <- paste0("0", md_clinical_all$SampleID)
row.names(md_clinical_all) <- md_clinical_all$SampleID
md_clinical_all$Treatment <- gsub(" .*", "", md_clinical_all$TRT01P.Olink)
md_clinical_all <- md_clinical_all %>%
  dplyr::mutate(PlateID = factor(case_when(
    PlateID == "Q-13795__flowcellposition_PA_SS240201_SP240234" ~ "Plate1",
    PlateID == "Q-13795__flowcellposition_PB_SS240201_SP240235" ~ "Plate2",
    PlateID == "Q-13795__flowcellposition_PA_SS240202_SP240236" ~ "Plate3",
    PlateID == "Q-13795__flowcellposition_PB_SS240202_SP240237" ~ "Plate4",
    PlateID == "Q-13795__flowcellposition_PA_SS240203_SP240238" ~ "Plate5",
    PlateID == "Q-13795__flowcellposition_PB_SS240203_SP240239" ~ "Plate6",
    PlateID == "Q-13795__flowcellposition_PA_SS240204_SP240240" ~ "Plate7",
    PlateID == "Q-13795__flowcellposition_PB_SS240204_SP240241" ~ "Plate8",
    PlateID == "Q-13795__flowcellposition_PA_SS240205_SP240242" ~ "Plate9",
    PlateID == "Q-13795__flowcellposition_PB_SS240205_SP240243" ~ "Plate10"
  )))
md_clinical_all$PlateID <- factor(md_clinical_all$PlateID, levels = c("Plate1", "Plate2", "Plate3", "Plate4", "Plate5", "Plate6", "Plate7", "Plate8", "Plate9", "Plate10"))

match <- read.table(file.path(dir.dat, "matching_IDs.txt"), header = TRUE, sep = "\t")
response <- response %>% dplyr::left_join(match, by = c("Originating_ID" = "Originating.ID"))

npx_mat <- read.table(file.path(dir.dat, "NPX_expression_all.txt"), sep = "\t", header = TRUE, row.names = 1)
samples_in <- colnames(npx_mat)
npx_mat <- npx_mat[, samples_in]

nn <- colnames(npx_mat)
nn <- gsub("X", "", nn)
colnames(npx_mat) <- nn

md_clinical_all <- md_clinical_all %>% dplyr::left_join(response, by = "Subject.ID")
md_clinical_all <- md_clinical_all[md_clinical_all$SampleID %in% colnames(npx_mat), ]
md_clinical_all$Group <- paste0(md_clinical_all$Treatment, "_", md_clinical_all$Week)

keep.vars <- c("SampleID", "Subject.ID", "PlateID", "AGE", "SEX", "RACE", "COUNTRY", "TNFGR1", "CORTGR1", "CRPN", "BMI", "Group", "Treatment", "Week", "Response_W14", "Remission_W52", "Sustained_Rvs_NR", "Histo_Remission_W52", "Endo_W52let1")

md_clinical_all <- md_clinical_all[, keep.vars]
md_clinical_all <- md_clinical_all %>% dplyr::distinct()

md_clinical_all$R_W14 <- gsub("Yes", "Resp", md_clinical_all$Response_W14)
md_clinical_all$R_W14 <- gsub("No", "Non_Resp", md_clinical_all$R_W14)
md_clinical_all$R_W52 <- gsub("Yes", "Resp", md_clinical_all$Remission_W52)
md_clinical_all$R_W52 <- gsub("No", "Non_Resp", md_clinical_all$R_W52)
md_clinical_all$Histo_RW52 <- gsub("Yes", "Resp", md_clinical_all$Histo_Remission_W52)
md_clinical_all$Histo_RW52 <- gsub("No", "Non_Resp", md_clinical_all$Histo_RW52)
md_clinical_all$SustR <- gsub("C_R", "Resp", md_clinical_all$Sustained_Rvs_NR)
md_clinical_all$SustR <- gsub("C_NR", "Non_Resp", md_clinical_all$SustR)
md_clinical_all$Endo_RW52 <- gsub("Yes", "Resp", md_clinical_all$Endo_W52let1)
md_clinical_all$Endo_RW52 <- gsub("No", "Non_Resp", md_clinical_all$Endo_RW52)

INVN <- function(x) {
  return(qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
}

npx_mat.invn <- data.frame(Assay = rownames(npx_mat), sapply(npx_mat, INVN)) %>% column_to_rownames(var = "Assay")
colnames(npx_mat.invn) <- colnames(npx_mat)

responses <- c("R_W52")
covariates_cols <- c("PlateID", "Subject.ID", "AGE", "SEX", "RACE", "COUNTRY", "TNFGR1", "CORTGR1", "CRPN", "BMI")

md_clinical_list <- list()
for (r in c("R_W14", "R_W14_CORT", "R_W52", "R_W52_CORT", "Histo_RW52", "Histo_RW52_CORT", "SustR", "SustR_CORT", "Endo_RW52", "Endo_RW52_CORT")) {
  t_name <- gsub("_CORT", "", r)

  md_clinical <- md_clinical_all
  md_clinical$Group_Resp <- paste0(md_clinical$Treatment, "_", md_clinical$Week, "_", md_clinical[, t_name])
  md_clinical <- md_clinical[md_clinical$SampleID %in% colnames(npx_mat.invn), ]
  md_clinical_resp <- subset(md_clinical, grepl("Resp|Non_Resp", Group_Resp))

  if (grepl("_CORT", r) == TRUE) {
    md_clinical_resp <- md_clinical_resp[md_clinical_resp$CORTGR1 == "No", ]
    md_clinical_resp <- md_clinical_resp[!is.na(md_clinical_resp$SampleID), ]
  }

  md_clinical_resp <- md_clinical_resp %>% dplyr::select(all_of(covariates_cols), Group_Resp, SampleID)
  md_clinical_resp <- md_clinical_resp[complete.cases(md_clinical_resp), ]
  row.names(md_clinical_resp) <- md_clinical_resp$SampleID

  md_clinical_list[[r]] <- md_clinical_resp
}

df_titles <- data.frame(
  "t_names" = c("R_W14", "R_W14_CORT", "R_W52", "R_W52_CORT", "Histo_RW52", "Histo_RW52_CORT", "SustR", "SustR_CORT", "Endo_RW52", "Endo_RW52_CORT"),
  "title" = c("Clinical Response W14", "Clinical Response W14 (corticosteroid free)", "Clinical Remission W52", "Clinical Remission W52 (corticosteroid free)", "Histological Remission W52", "Histological Remission W52 (corticosteroid free)", "Sustained Clinical Remission", "Sustained Clinical Remission (corticosteroid free)", "Endoscopic Remission W52", "Endoscopic Remission W52 (corticosteroid free)"),
  "cap" = c("Clinical", "Clinical", "Clinical", "Clinical", "Histological", "Histological", "Clinical Sustained", "Clinical Sustained", "Endoscopic", "Endoscopic"),
  "week" = c("W14", "W14", "W52", "W52", "W52", "W52", "W14 and W52", "W14 and W52", "W52", "W52")
)

responses <- c("R_W52")

# ------------------------------------------------------------------
# Line plots of proteins of interest
# ------------------------------------------------------------------
std.error <- function(x) sd(x) / sqrt(length(x))

for (p in intersect(prots_interest, rownames(npx_mat.invn))) {
  cat("####", p, " {.tabset} \n")

  for (r in responses) {
    dft <- df_titles %>% dplyr::filter(t_names == r)
    cat("#####", dft$title, "\n")

    meta <- md_clinical_list[[r]]
    p_exp <- list()
    p_baseline <- list()
    p_pch <- list()

    for (drug in c("Adalimumab", "Vedolizumab")) {
      samples_names <- md_clinical_list[[r]] %>% dplyr::filter(grepl(drug, Group_Resp)) %>% pull(SampleID)

      if (length(which(rownames(npx_mat.invn) == p)) > 0) {
        expr <- npx_mat.invn[p, samples_names] %>% as.data.frame()

        df <- as.data.frame(t(expr))
        df$Sample <- rownames(df)

        df_melted <- melt(expr, value.name = "Expression")
        colnames(df_melted) <- c("SampleID", "Expression")
        df_melted$SampleID <- as.character(df_melted$SampleID)

        df_melted <- df_melted %>% left_join(meta[, c("SampleID", "Subject.ID", "Group_Resp")], by = "SampleID")
        df_melted$Resp <- gsub(".*W[0-9]+_", "", df_melted$Group_Resp)
        df_melted$Week <- gsub(".*W([0-9]+)_.*", "Week \\1", df_melted$Group_Resp)

        paired_subj <- df_melted %>% dplyr::select(Subject.ID, Week, Resp) %>% dplyr::distinct() %>% dplyr::group_by(Subject.ID) %>%
          dplyr::summarize(Num = n(), .groups = "drop") %>% dplyr::filter(Num > 1) %>% dplyr::pull(Subject.ID)

        df_melted <- df_melted[df_melted$Subject.ID %in% paired_subj, ]

        df_sum <- df_melted %>% group_by(Week, Resp) %>%
          summarize(
            mean = mean(Expression),
            se = std.error(Expression),
            .groups = "drop"
          ) %>%
          dplyr::mutate(
            Resp = gsub("Resp", "Responders", Resp),
            Resp = gsub("_", "-", Resp)
          )

        resp <- df_melted[df_melted$Resp == "Resp", ]
        nonresp <- df_melted[df_melted$Resp == "Non_Resp", ]
        res.resp <- wilcox.test(resp[resp$Week == "Week 0", "Expression"], resp[resp$Week == "Week 14", "Expression"], conf.int = TRUE)$p.value
        res.nonresp <- wilcox.test(nonresp[nonresp$Week == "Week 0", "Expression"], nonresp[nonresp$Week == "Week 14", "Expression"], conf.int = TRUE)$p.value

        baseline <- compare_means(Expression ~ Resp, df_melted[df_melted$Week == "Week 0", ], method = "t.test", p.adjust.method = "fdr")
        wk14 <- compare_means(Expression ~ Resp, df_melted[df_melted$Week == "Week 14", ], method = "t.test", p.adjust.method = "fdr")

        n_resp <- length(unique(resp$SampleID))
        n_nonresp <- length(unique(nonresp$SampleID))
        n_uni_resp <- length(unique(resp$Subject.ID))
        n_uni_nonresp <- length(unique(nonresp$Subject.ID))

        plot1 <- ggplot(df_sum, aes(x = Week, y = mean, group = Resp, colour = Resp)) +
          geom_line() +
          geom_point(fill = "white") +
          geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
          theme_minimal() +
          scale_colour_manual(values = c("red", "darkblue")) +
          xlab("") +
          ylab("Expression") +
          labs(
            title = paste0(p, " in ", drug), colour = "",
            caption = paste0(
              "Non-responder subjects = ", n_uni_nonresp, "\nResponder subjects = ", n_uni_resp,
              "\nBaseline p-val = ", round(baseline$p, 2), "\nWeek 14 p-val = ", round(wk14$p, 2)
            )
          ) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 15),
            plot.caption = element_text(size = 15),
            plot.title = element_text(size = 16),
            legend.text = element_text(size = 14)
          ) +
          annotation_custom(grob = textGrob(paste0("p-val= ", round(res.resp, 2)), x = unit(0.5, "npc"), y = unit(.85, "npc"), gp = gpar(col = "darkblue"))) +
          annotation_custom(grob = textGrob(paste0("p-val= ", round(res.nonresp, 2)), x = unit(0.5, "npc"), y = unit(.95, "npc"), gp = gpar(col = "red")))

        p_exp[[drug]] <- plot1

        a <- df_melted %>%
          dplyr::select(-SampleID, -Group_Resp) %>%
          pivot_wider(names_from = c("Week", "Resp"), values_from = "Expression")

        nonresponders <- a %>% dplyr::select(Subject.ID, contains("_Non_Resp")); nonresponders <- nonresponders[complete.cases(nonresponders), ]
        responders <- a %>% dplyr::select(Subject.ID, setdiff(names(a), names(nonresponders))); responders <- responders[complete.cases(responders), ]

        nonresponders_change <- data.frame(
          SubjectID = nonresponders$Subject.ID,
          "Week 0" = 0,
          "Week 14" = nonresponders$`Week 14_Non_Resp` - nonresponders$`Week 0_Non_Resp`
        )

        responders_change <- data.frame(
          SubjectID = responders$Subject.ID,
          "Week 0" = 0,
          "Week 14" = responders$`Week 14_Resp` - responders$`Week 0_Resp`
        )

        new.df <- data.frame(
          rbind(
            cbind(reshape2::melt(responders_change), "Resp" = "Responders"),
            cbind(reshape2::melt(nonresponders_change), "Resp" = "Non-Responders")
          )
        )
        new.df$variable <- gsub("\\.", " ", new.df$variable)
        new.df.backup <- new.df
        new.df <- new.df %>% group_by(variable, Resp) %>%
          summarize(
            mean = mean(value),
            se = std.error(value),
            .groups = "drop"
          )

        res.resp2 <- wilcox.test(responders_change$Week.0, responders_change$Week.14, paired = TRUE, conf.int = TRUE)$p.value
        res.nonresp2 <- wilcox.test(nonresponders_change$Week.0, nonresponders_change$Week.14, paired = TRUE, conf.int = TRUE)$p.value
        wk14 <- compare_means(value ~ Resp, data = new.df.backup[new.df.backup$variable == "Week 14", ])

        plot2 <- ggplot(new.df, aes(x = variable, y = mean, group = Resp, colour = Resp)) +
          geom_line() +
          geom_point(fill = "white") +
          geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
          theme_minimal() +
          scale_colour_manual(values = c("red", "darkblue")) +
          xlab("") +
          ylab("Expression change from baseline") +
          labs(
            title = paste0(p, " in ", drug), colour = "",
            caption = paste0(
              "Non-responder subjects = ", n_uni_nonresp, "\nResponder subjects = ", n_uni_resp,
              "\nWeek 14 p-val = ", round(wk14$p, 2)
            )
          ) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 15),
            plot.caption = element_text(size = 15),
            plot.title = element_text(size = 16),
            legend.text = element_text(size = 14)
          ) +
          annotation_custom(grob = textGrob(paste0("p-val= ", round(res.resp2, 2)), x = unit(0.5, "npc"), y = unit(.85, "npc"), gp = gpar(col = "darkblue"))) +
          annotation_custom(grob = textGrob(paste0("p-val= ", round(res.nonresp2, 2)), x = unit(0.5, "npc"), y = unit(.95, "npc"), gp = gpar(col = "red")))

        p_baseline[[drug]] <- plot2

        nonresponders_pch <- data.frame(
          SubjectID = nonresponders$Subject.ID,
          "Week 0" = 0,
          "Week 14" = (exp(nonresponders$`Week 14_Non_Resp` - nonresponders$`Week 0_Non_Resp`) - 1) * 100
        )

        responders_pch <- data.frame(
          SubjectID = responders$Subject.ID,
          "Week 0" = 0,
          "Week 14" = (exp(responders$`Week 14_Resp` - responders$`Week 0_Resp`) - 1) * 100
        )

        new.df.pch <- data.frame(
          rbind(
            cbind(reshape2::melt(responders_pch), "Resp" = "Responders"),
            cbind(reshape2::melt(nonresponders_pch), "Resp" = "Non-Responders")
          )
        )
        new.df.pch$variable <- gsub("\\.", " ", new.df.pch$variable)

        new.df.pch.backup <- new.df.pch

        new.df.pch <- new.df.pch %>% group_by(variable, Resp) %>%
          summarize(
            mean = mean(value),
            se = std.error(value),
            .groups = "drop"
          )

        res.resp2.pch <- wilcox.test(responders_pch$Week.0, responders_pch$Week.14, paired = TRUE, conf.int = TRUE)$p.value
        res.nonresp2.pch <- wilcox.test(nonresponders_pch$Week.0, nonresponders_pch$Week.14, paired = TRUE, conf.int = TRUE)$p.value
        wk14 <- compare_means(value ~ Resp, data = new.df.pch.backup[new.df.pch.backup$variable == "Week 14", ])

        plot3 <- ggplot(new.df.pch, aes(x = variable, y = mean, group = Resp, colour = Resp)) +
          geom_line() +
          geom_point(fill = "white") +
          geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
          theme_minimal() +
          scale_colour_manual(values = c("red", "darkblue")) +
          xlab("") +
          ylab("Expression change from baseline") +
          labs(
            title = paste0(p, " in ", drug), colour = "",
            caption = paste0(
              "Non-responder subjects = ", n_uni_nonresp, "\nResponder subjects = ", n_uni_resp,
              "\nWeek 14 p-val = ", round(wk14$p, 2)
            )
          ) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 15),
            plot.caption = element_text(size = 15),
            plot.title = element_text(size = 16),
            legend.text = element_text(size = 14)
          ) +
          annotation_custom(grob = textGrob(paste0("p-val= ", round(res.resp2.pch, 2)), x = unit(0.5, "npc"), y = unit(.85, "npc"), gp = gpar(col = "darkblue"))) +
          annotation_custom(grob = textGrob(paste0("p-val= ", round(res.nonresp2.pch, 2)), x = unit(0.5, "npc"), y = unit(.95, "npc"), gp = gpar(col = "red")))

        p_pch[[drug]] <- plot3
      }
    }

    if (length(p_exp) > 1) {
      cat("###### Protein expression \n\n")
      limits <- c(ggplot_build(p_exp[[1]])$layout$panel_params[[1]]$y.range, ggplot_build(p_exp[[2]])$layout$panel_params[[1]]$y.range)
      plot1 <- ggarrange(p_exp[[1]] + ylim(range(limits)), p_exp[[2]] + ylim(range(limits)), common.legend = TRUE, legend = "bottom") +
        theme(legend.text = element_text(size = 14))
      print(plot1)
      cat("\n\n")

      cat("###### Protein expression change from baseline \n\n")
      limits <- c(ggplot_build(p_baseline[[1]])$layout$panel_params[[1]]$y.range, ggplot_build(p_baseline[[2]])$layout$panel_params[[1]]$y.range)
      plot2 <- ggarrange(p_baseline[[1]] + ylim(range(limits)), p_baseline[[2]] + ylim(range(limits)), common.legend = TRUE, legend = "bottom") +
        theme(legend.text = element_text(size = 14))
      print(plot2)
      cat("\n\n")

      cat("###### Protein percent change from baseline \n\n")
      limits <- c(ggplot_build(p_pch[[1]])$layout$panel_params[[1]]$y.range, ggplot_build(p_pch[[2]])$layout$panel_params[[1]]$y.range)
      plot3 <- ggarrange(p_pch[[1]] + ylim(range(limits)), p_pch[[2]] + ylim(range(limits)), common.legend = TRUE, legend = "bottom") +
        theme(legend.text = element_text(size = 14))
      print(plot3)
      cat("\n\n")
    }
    cat("\n\n")
  }
}

# ------------------------------------------------------------------
# Gene signatures (ssGSEA)
# ------------------------------------------------------------------
npx.ssgsea <- gsva(ssgseaParam(exprData = as.matrix(npx_mat.invn), geneSets = signatures), verbose = TRUE)

signature.gsva <- as.data.frame(t(npx.ssgsea)); names(npx.ssgsea) <- names(signatures)
signature.gsva <- signature.gsva %>% rownames_to_column("Sample")

for (sig in names(signature.gsva)[2:ncol(signature.gsva)]) {
  cat("####", sig, " {.tabset} \n")
  for (r in responses) {
    dft <- df_titles %>% dplyr::filter(t_names == r)
    cat("#####", dft$title, "\n")

    meta <- md_clinical_list[[r]]
    p_exp <- list()
    p_baseline <- list()
    p_pch <- list()

    for (drug in c("Adalimumab", "Vedolizumab")) {
      samples_names <- md_clinical_list[[r]] %>% dplyr::filter(grepl(drug, Group_Resp)) %>% pull(SampleID)

      sig.df <- signature.gsva %>% dplyr::select(Sample, sig)
      sig.df <- sig.df[sig.df$Sample %in% samples_names, ]
      colnames(sig.df) <- c("SampleID", "signature")
      sig.df$SampleID <- as.character(sig.df$SampleID)

      sig.df <- sig.df %>% left_join(meta[, c("SampleID", "Subject.ID", "Group_Resp")], by = "SampleID")
      sig.df$Resp <- gsub(".*W[0-9]+_", "", sig.df$Group_Resp)
      sig.df$Week <- gsub(".*W([0-9]+)_.*", "Week \\1", sig.df$Group_Resp)

      paired_subj <- sig.df %>% dplyr::select(Subject.ID, Week, Resp) %>% dplyr::distinct() %>% dplyr::group_by(Subject.ID) %>%
        dplyr::summarize(Num = n(), .groups = "drop") %>% dplyr::filter(Num > 1) %>% dplyr::pull(Subject.ID)

      sig.df <- sig.df[sig.df$Subject.ID %in% paired_subj, ]

      df_sum <- sig.df %>% group_by(Week, Resp) %>%
        summarize(
          mean = mean(signature),
          se = std.error(signature),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          Resp = gsub("Resp", "Responders", Resp),
          Resp = gsub("_", "-", Resp)
        )

      resp <- sig.df[sig.df$Resp == "Resp", ]
      nonresp <- sig.df[sig.df$Resp == "Non_Resp", ]
      res.resp <- wilcox.test(resp[resp$Week == "Week 0", "signature"], resp[resp$Week == "Week 14", "signature"], conf.int = TRUE)$p.value
      res.nonresp <- wilcox.test(nonresp[nonresp$Week == "Week 0", "signature"], nonresp[nonresp$Week == "Week 14", "signature"], conf.int = TRUE)$p.value
      baseline <- compare_means(signature ~ Resp, sig.df[sig.df$Week == "Week 0", ], method = "t.test", p.adjust.method = "fdr")
      wk14 <- compare_means(signature ~ Resp, sig.df[sig.df$Week == "Week 14", ], method = "t.test", p.adjust.method = "fdr")

      n_resp <- length(unique(resp$SampleID))
      n_nonresp <- length(unique(nonresp$SampleID))
      n_uni_resp <- length(unique(resp$Subject.ID))
      n_uni_nonresp <- length(unique(nonresp$Subject.ID))

      plot1 <- ggplot(df_sum, aes(x = Week, y = mean, group = Resp, colour = Resp)) +
        geom_line() +
        geom_point(fill = "white") +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
        theme_minimal() +
        scale_colour_manual(values = c("red", "darkblue")) +
        xlab("") +
        ylab(paste0(sig, " ssGSEA")) +
        labs(
          title = paste0(sig, " in ", drug), colour = "",
          caption = paste0(
            "Non-responder subjects = ", n_uni_nonresp, "\nResponder subjects = ", n_uni_resp,
            "\nBaseline p-val = ", round(baseline$p, 2), "\nWeek 14 p-val = ", round(wk14$p, 2)
          )
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 15),
          plot.caption = element_text(size = 15),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 14)
        ) +
        annotation_custom(grob = textGrob(paste0("p-val= ", round(res.resp, 2)), x = unit(0.5, "npc"), y = unit(.85, "npc"), gp = gpar(col = "darkblue"))) +
        annotation_custom(grob = textGrob(paste0("p-val= ", round(res.nonresp, 2)), x = unit(0.5, "npc"), y = unit(.95, "npc"), gp = gpar(col = "red")))
      p_exp[[drug]] <- plot1

      a <- sig.df %>%
        dplyr::select(-SampleID, -Group_Resp) %>%
        pivot_wider(names_from = c("Week", "Resp"), values_from = "signature")

      nonresponders <- a %>% dplyr::select(Subject.ID, contains("_Non_Resp")); nonresponders <- nonresponders[complete.cases(nonresponders), ]
      responders <- a %>% dplyr::select(Subject.ID, setdiff(names(a), names(nonresponders))); responders <- responders[complete.cases(responders), ]

      nonresponders_change <- data.frame(
        SubjectID = nonresponders$Subject.ID,
        "Week 0" = 0,
        "Week 14" = nonresponders$`Week 14_Non_Resp` - nonresponders$`Week 0_Non_Resp`
      )

      responders_change <- data.frame(
        SubjectID = responders$Subject.ID,
        "Week 0" = 0,
        "Week 14" = responders$`Week 14_Resp` - responders$`Week 0_Resp`
      )

      new.df <- data.frame(rbind(
        cbind(reshape2::melt(responders_change), "Resp" = "Responders"),
        cbind(reshape2::melt(nonresponders_change), "Resp" = "Non-Responders")
      ))
      new.df$variable <- gsub("\\.", " ", new.df$variable)
      new.df.backup <- new.df
      new.df <- new.df %>% group_by(variable, Resp) %>%
        summarize(
          mean = mean(value),
          se = std.error(value),
          .groups = "drop"
        )

      res.resp2 <- wilcox.test(responders_change$Week.0, responders_change$Week.14, paired = TRUE, conf.int = TRUE)$p.value
      res.nonresp2 <- wilcox.test(nonresponders_change$Week.0, nonresponders_change$Week.14, paired = TRUE, conf.int = TRUE)$p.value
      wk14 <- compare_means(value ~ Resp, data = new.df.backup[new.df.backup$variable == "Week 14", ])

      plot2 <- ggplot(new.df, aes(x = variable, y = mean, group = Resp, colour = Resp)) +
        geom_line() +
        geom_point(fill = "white") +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
        theme_minimal() +
        scale_colour_manual(values = c("red", "darkblue")) +
        xlab("") +
        ylab(paste0(sig, " ssGSEA change from baseline")) +
        labs(
          title = paste0(sig, " in ", drug), colour = "",
          caption = paste0(
            "Non-responder subjects = ", n_uni_nonresp, "\nResponder subjects = ", n_uni_resp,
            "\nWeek 14 p-val = ", round(wk14$p, 2)
          )
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 15),
          plot.caption = element_text(size = 15),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 14)
        ) +
        annotation_custom(grob = textGrob(paste0("p-val= ", round(res.resp2, 2)), x = unit(0.5, "npc"), y = unit(.85, "npc"), gp = gpar(col = "darkblue"))) +
        annotation_custom(grob = textGrob(paste0("p-val= ", round(res.nonresp2, 2)), x = unit(0.5, "npc"), y = unit(.95, "npc"), gp = gpar(col = "red")))

      p_baseline[[drug]] <- plot2

      nonresponders_pch <- data.frame(
        SubjectID = nonresponders$Subject.ID,
        "Week 0" = 0,
        "Week 14" = (exp(nonresponders$`Week 14_Non_Resp` - nonresponders$`Week 0_Non_Resp`) - 1) * 100
      )

      responders_pch <- data.frame(
        SubjectID = responders$Subject.ID,
        "Week 0" = 0,
        "Week 14" = (exp(responders$`Week 14_Resp` - responders$`Week 0_Resp`) - 1) * 100
      )

      new.df.pch <- data.frame(
        rbind(
          cbind(reshape2::melt(responders_pch), "Resp" = "Responders"),
          cbind(reshape2::melt(nonresponders_pch), "Resp" = "Non-Responders")
        )
      )
      new.df.pch$variable <- gsub("\\.", " ", new.df.pch$variable)
      new.df.pch.backup <- new.df.pch
      new.df.pch <- new.df.pch %>% group_by(variable, Resp) %>%
        summarize(
          mean = mean(value),
          se = std.error(value),
          .groups = "drop"
        )

      res.resp2.pch <- wilcox.test(responders_pch$Week.0, responders_pch$Week.14, paired = TRUE, conf.int = TRUE)$p.value
      res.nonresp2.pch <- wilcox.test(nonresponders_pch$Week.0, nonresponders_pch$Week.14, paired = TRUE, conf.int = TRUE)$p.value
      wk14 <- compare_means(value ~ Resp, data = new.df.pch.backup[new.df.pch.backup$variable == "Week 14", ])

      plot3 <- ggplot(new.df.pch, aes(x = variable, y = mean, group = Resp, colour = Resp)) +
        geom_line() +
        geom_point(fill = "white") +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
        theme_minimal() +
        scale_colour_manual(values = c("red", "darkblue")) +
        xlab("") +
        ylab(paste0(sig, " ssGSEA PCH from baseline")) +
        labs(
          title = paste0(sig, " in ", drug), colour = "",
          caption = paste0(
            "Non-responder subjects = ", n_uni_nonresp, "\nResponder subjects = ", n_uni_resp,
            "\nWeek 14 p-val = ", round(wk14$p, 2)
          )
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 15),
          plot.caption = element_text(size = 15),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 14)
        ) +
        annotation_custom(grob = textGrob(paste0("p-val= ", round(res.resp2.pch, 2)), x = unit(0.5, "npc"), y = unit(.85, "npc"), gp = gpar(col = "darkblue"))) +
        annotation_custom(grob = textGrob(paste0("p-val= ", round(res.nonresp2.pch, 2)), x = unit(0.5, "npc"), y = unit(.95, "npc"), gp = gpar(col = "red")))

      p_pch[[drug]] <- plot3
    }

    if (length(p_exp) > 1) {
      cat("###### Protein expression \n\n")
      limits <- c(ggplot_build(p_exp[[1]])$layout$panel_params[[1]]$y.range, ggplot_build(p_exp[[2]])$layout$panel_params[[1]]$y.range)
      plot1 <- ggarrange(p_exp[[1]] + ylim(range(limits)), p_exp[[2]] + ylim(range(limits)), common.legend = TRUE, legend = "bottom") +
        theme(legend.text = element_text(size = 14))
      print(plot1)
      cat("\n\n")

      cat("###### Protein expression change from baseline \n\n")
      limits <- c(ggplot_build(p_baseline[[1]])$layout$panel_params[[1]]$y.range, ggplot_build(p_baseline[[2]])$layout$panel_params[[1]]$y.range)
      plot2 <- ggarrange(p_baseline[[1]] + ylim(range(limits)), p_baseline[[2]] + ylim(range(limits)), common.legend = TRUE, legend = "bottom") +
        theme(legend.text = element_text(size = 14))
      print(plot2)
      cat("\n\n")

      cat("###### Protein percent change from baseline \n\n")
      limits <- c(ggplot_build(p_pch[[1]])$layout$panel_params[[1]]$y.range, ggplot_build(p_pch[[2]])$layout$panel_params[[1]]$y.range)
      plot3 <- ggarrange(p_pch[[1]] + ylim(range(limits)), p_pch[[2]] + ylim(range(limits)), common.legend = TRUE, legend = "bottom") +
        theme(legend.text = element_text(size = 14))
      print(plot3)
      cat("\n\n")
    }
    cat("\n\n")
  }
}
