# Objective: load GEMINI survival data, fit Kaplan-Meier curves for hospital/procedure outcomes by group,
# and export risk-table survival plots as PDFs.
library(ggsurvfit)
library(survminer) # for survfit2 function (optional, but recommended)

#df <- read.csv("C:\\Users\\ibi9245\\OneDrive - Takeda\\Projects\\Entyvio\\C13006\\gemini survival.csv")
df <- read.csv("C:\\Users\\ibi9245\\OneDrive - Takeda\\Projects\\Entyvio\\C13006\\gemini survival2.csv")
#df <- read.csv("C:\\Users\\ibi9245\\OneDrive - Takeda\\Projects\\Entyvio\\C13006\\wk6mh survival.csv")
# Sample data
#df <- data.frame(
#  time = sample(0:10, 100, replace = TRUE),
#  status = sample(0:1, 100, replace = TRUE),
#  sex = sample(0:1, 100, replace = TRUE)
#)

# Create survival curves
fit1 <- survfit2(Surv(time_hosp_or_proc, event_hosp_or_proc) ~ group, data = df)
fit2 <- survfit2(Surv(time_hosp, event_hosp) ~ group, data = df)
fit3 <- survfit2(Surv(time_proc, event_proc) ~ group, data = df)

# Create Kaplan-Meier plot with ggsurvfit
#ggsurvfit(fit) + 
#  scale_y_continuous(limits = c(0, 1)) + # Set y-axis limits to 0-1
#  ggtitle("Kaplan-Meier Plot") # Add a title

# Add risk table (optional)
p1<-ggsurvfit(fit1 ) + 
 add_risktable() + # Add risk table with default settings
  scale_y_continuous(limits = c(0, 1)) + 
  ggtitle("Kaplan-Meier Plot with Risk Table") + scale_ggsurvfit(x_scales = list(breaks = seq(0, 3000, by = 365)))+add_confidence_interval()

ggsave("C:\\Users\\ibi9245\\OneDrive - Takeda\\Projects\\Entyvio\\C13006\\km_hosp_or_proc_v1v2_v3v4.pdf", plot = p1, width = 6, height = 6, units = "in")

p2<-ggsurvfit(fit2 ) + 
  add_risktable() + # Add risk table with default settings
  scale_y_continuous(limits = c(0, 1)) + 
  ggtitle("Kaplan-Meier Plot with Risk Table") + scale_ggsurvfit(x_scales = list(breaks = seq(0, 3000, by = 365)))+add_confidence_interval()

ggsave("C:\\Users\\ibi9245\\OneDrive - Takeda\\Projects\\Entyvio\\C13006\\km_hosp_wk6endo.pdf", plot = p2, width = 6, height = 6, units = "in")

p3<-ggsurvfit(fit3 ) + 
 add_risktable() + # Add risk table with default settings
  scale_y_continuous(limits = c(0, 1)) + 
  ggtitle("Kaplan-Meier Plot with Risk Table") + scale_ggsurvfit(x_scales = list(breaks = seq(0, 3000, by = 365)))+add_confidence_interval()

ggsave("C:\\Users\\ibi9245\\OneDrive - Takeda\\Projects\\Entyvio\\C13006\\km_proc_wk6endo.pdf", plot = p3, width = 6, height = 6, units = "in")

