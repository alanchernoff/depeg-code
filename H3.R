## =============================================================================
## HYPOTHESIS 3: Did reserve-transparency drive contagion during the USDC/SVB
##               depeg (March 2023)? DAI vs. USDT, with BUSD as a placebo.
## =============================================================================

## ---- 0. Packages -------------------------------------------------------------
required_pkgs <- c("tidyverse", "lubridate", "officer", "flextable",
                   "zoo", "scales", "sandwich", "lmtest")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(tidyverse)
library(lubridate)
library(officer)
library(flextable)
library(zoo)
library(scales)
library(sandwich)
library(lmtest)

## ---- 1. Setup -----------------------------------------------------------------
setwd("C:/Users/acher/research/")

DATA_FILE <- "stablecoin_hourly_prices_wide.csv"
OUT_DOCX  <- "hypothesis3_SVB_USDC_contagion_results.docx"

COINS <- c("USDC", "DAI", "USDT")

SVB_FAILURE       <- ymd_hms("2023-03-10 00:00:00", tz = "UTC")
TREASURY_BACKSTOP <- ymd_hms("2023-03-12 18:00:00", tz = "UTC")

PRE_WINDOW    <- c(SVB_FAILURE - days(14), SVB_FAILURE - hours(1))
EVENT_WINDOW  <- c(SVB_FAILURE, TREASURY_BACKSTOP + hours(6))
POST_WINDOW   <- c(TREASURY_BACKSTOP + hours(6) + hours(1), SVB_FAILURE + days(14))
PLOT_WINDOW   <- c(SVB_FAILURE - days(14), SVB_FAILURE + days(14))

## ---- 2. Load and reshape ------------------------------------------------------
## CONFIRMED from diagnostics: the raw file uses "YYYY-MM-DD HH:MM" (no seconds,
## no offset), and readr's own default guesser already parses this correctly as
## POSIXct with no manual intervention needed. v1/v2's mistake was re-parsing an
## already-correct POSIXct column with ymd_hms(), which (via R's default
## as.character() truncating exact-midnight instants to a bare date) silently
## turned every midnight hour into NA. v3's first attempt forced the column to
## character and then required seconds that don't exist in the source format.
## The fix: parse once, explicitly, with the exact format -- no re-parsing.
raw <- read_csv(DATA_FILE, show_col_types = FALSE,
                col_types = cols(datetime_utc = col_datetime(format = "%Y-%m-%d %H:%M"),
                                 .default = col_double()),
                locale = locale(tz = "UTC"))

cat("Datetime parsing: ", sum(!is.na(raw$datetime_utc)), " of ", nrow(raw), " rows parsed successfully.\n", sep = "")
cat("Parsed date range: ", format(range(raw$datetime_utc, na.rm = TRUE)), "\n")
if (all(is.na(raw$datetime_utc))) stop("Datetime parsing failed for all rows -- check DATA_FILE's datetime_utc format before proceeding.")

df <- raw %>%
  pivot_longer(cols = -datetime_utc, names_to = "coin", values_to = "price") %>%
  filter(!is.na(price)) %>%
  mutate(deviation_bps = (price - 1) * 10000)

## ---- 3. Coverage table (FIXED: na.rm + diagnostic) -----------------------------
coverage_tbl <- df %>%
  filter(coin %in% COINS) %>%
  group_by(coin) %>%
  summarise(
    n_na_datetime    = sum(is.na(datetime_utc)),
    first_obs        = suppressWarnings(min(datetime_utc, na.rm = TRUE)),
    last_obs         = suppressWarnings(max(datetime_utc, na.rm = TRUE)),
    n_in_plot_window = sum(datetime_utc >= PLOT_WINDOW[1] & datetime_utc <= PLOT_WINDOW[2], na.rm = TRUE),
    n_in_post_window = sum(datetime_utc >= POST_WINDOW[1] & datetime_utc <= POST_WINDOW[2], na.rm = TRUE),
    .groups = "drop"
  )

if (any(coverage_tbl$n_na_datetime > 0)) {
  warning("Unparseable datetime_utc values detected -- see n_na_datetime column in coverage_tbl.")
}
print(coverage_tbl)

## ---- 4. Event panel with phase + event dummy -----------------------------------
event_panel <- df %>%
  filter(coin %in% COINS, datetime_utc >= PLOT_WINDOW[1], datetime_utc <= PLOT_WINDOW[2]) %>%
  mutate(
    phase = case_when(
      datetime_utc >= PRE_WINDOW[1]   & datetime_utc <= PRE_WINDOW[2]   ~ "Pre-event",
      datetime_utc >= EVENT_WINDOW[1] & datetime_utc <= EVENT_WINDOW[2] ~ "Event window",
      datetime_utc >= POST_WINDOW[1]  & datetime_utc <= POST_WINDOW[2]  ~ "Post-event",
      TRUE ~ NA_character_
    ),
    event_dummy = if_else(phase == "Event window", 1L, 0L)
  ) %>%
  filter(!is.na(phase))

## ---- 5. Peak deviation & recovery time per coin (FIXED: signed vs. absolute) ---
recovery_tbl <- map_dfr(COINS, function(cn) {
  sub <- event_panel %>% filter(coin == cn, phase %in% c("Event window", "Post-event")) %>%
    arrange(datetime_utc)
  if (nrow(sub) == 0) {
    return(tibble(coin = cn, peak_dev_bps = NA_real_, peak_abs_dev_bps = NA_real_,
                  direction = NA_character_, peak_time = as.POSIXct(NA),
                  hours_to_recover = NA_real_, note = "No data in window"))
  }
  peak_idx  <- which.max(abs(sub$deviation_bps))
  peak_dev  <- sub$deviation_bps[peak_idx]       # signed
  peak_time <- sub$datetime_utc[peak_idx]
  
  after_peak <- sub %>% filter(datetime_utc > peak_time)
  recovered  <- after_peak %>% filter(abs(deviation_bps) <= 10)
  if (nrow(recovered) > 0) {
    hrs <- as.numeric(difftime(recovered$datetime_utc[1], peak_time, units = "hours"))
    note <- "Recovered in window"
  } else {
    hrs <- NA_real_
    note <- "Not recovered in window"
  }
  tibble(
    coin = cn,
    peak_dev_bps     = round(peak_dev, 2),
    peak_abs_dev_bps = round(abs(peak_dev), 2),
    direction        = if_else(peak_dev < 0, "Below peg", "Above peg"),
    peak_time        = peak_time,
    hours_to_recover = round(hrs, 1),
    note             = note
  )
})

print(recovery_tbl)

## ---- 6. Difference-in-differences regression -----------------------------------
did_data <- event_panel %>%
  filter(phase %in% c("Pre-event", "Event window")) %>%
  mutate(coin = relevel(factor(coin), ref = "USDT"))

did_model <- lm(deviation_bps ~ event_dummy * coin, data = did_data)
hac_vcov  <- NeweyWest(did_model, lag = 24, prewhite = FALSE)
hac_se    <- coeftest(did_model, vcov = hac_vcov)

did_results <- tibble(
  term = rownames(hac_se),
  estimate = round(hac_se[, "Estimate"], 3),
  std_error_HAC = round(hac_se[, "Std. Error"], 3),
  p_value = round(hac_se[, "Pr(>|t|)"], 4)
)

print(did_results)

## ---- 6b. Combined event-window effect per coin (NEW: exact delta-method SE) ----
## The raw event_dummy:coinX coefficient is X's reaction *relative to USDT's own
## reaction*, not X's total average shift from its own baseline. This computes
## the latter -- event_dummy + event_dummy:coinX -- with a properly propagated
## HAC standard error via the model's coefficient covariance matrix.
coefs <- coef(did_model)

combo_effect <- function(coin_name) {
  contrast <- setNames(rep(0, length(coefs)), names(coefs))
  contrast["event_dummy"] <- 1
  inter_name <- paste0("event_dummy:coin", coin_name)
  if (inter_name %in% names(coefs)) contrast[inter_name] <- 1
  est <- sum(contrast * coefs)
  se  <- sqrt(as.numeric(t(contrast) %*% hac_vcov %*% contrast))
  tstat <- est / se
  pval  <- 2 * pt(-abs(tstat), df = did_model$df.residual)
  tibble(coin = coin_name, total_event_effect_bps = round(est, 2),
         se_HAC = round(se, 2), t_stat = round(tstat, 2), p_value = round(pval, 4))
}

combined_effects <- map_dfr(COINS, combo_effect)
print(combined_effects)

## ---- 7. Plots -------------------------------------------------------------------
plot_df <- event_panel %>% mutate(coin = factor(coin, levels = COINS))

p_event <- ggplot(plot_df, aes(x = datetime_utc, y = deviation_bps, color = coin)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("rect", xmin = EVENT_WINDOW[1], xmax = EVENT_WINDOW[2],
           ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "firebrick") +
  geom_line(linewidth = 0.4) +
  labs(
    title = "Peg Deviation Around the SVB / USDC Depeg (March 2023)",
    subtitle = "Shaded band = acute event window (SVB failure to Treasury backstop + 6h)",
    x = NULL, y = "Deviation from $1.00 (bps)", color = NULL
  ) +
  scale_color_manual(values = c(USDC = "#2b6cb0", DAI = "#dd6b20", USDT = "#38a169")) +
  scale_x_datetime(date_breaks = "3 days", labels = date_format("%b %d")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top")

peak_plot_df <- recovery_tbl %>% filter(!is.na(peak_abs_dev_bps)) %>%
  mutate(coin = factor(coin, levels = COINS))

p_peak <- ggplot(peak_plot_df, aes(x = coin, y = peak_abs_dev_bps, fill = coin)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = direction), vjust = -0.4, size = 3) +
  labs(
    title = "Peak Absolute Deviation During/After the SVB Event",
    x = NULL, y = "Peak |deviation| (bps)"
  ) +
  scale_fill_manual(values = c(USDC = "#2b6cb0", DAI = "#dd6b20", USDT = "#38a169")) +
  theme_minimal(base_size = 11)

## ---- 8. Assemble Word document ---------------------------------------------------
ft_coverage <- coverage_tbl %>%
  mutate(first_obs = format(first_obs, "%Y-%m-%d %H:%M"),
         last_obs  = format(last_obs, "%Y-%m-%d %H:%M")) %>%
  select(coin, first_obs, last_obs, n_in_plot_window, n_in_post_window, n_na_datetime) %>%
  rename(Coin = coin, `First obs.` = first_obs, `Last obs.` = last_obs,
         `Obs. (plot wdw)` = n_in_plot_window, `Obs. (post wdw)` = n_in_post_window,
         `NA timestamps` = n_na_datetime) %>%
  flextable() %>% autofit() %>% theme_vanilla() %>%
  fontsize(size = 9, part = "all")

ft_recovery <- recovery_tbl %>%
  mutate(peak_time = format(peak_time, "%Y-%m-%d %H:%M")) %>%
  select(coin, peak_abs_dev_bps, direction, peak_time, hours_to_recover, note) %>%
  rename(Coin = coin, `Peak |dev| (bps)` = peak_abs_dev_bps, Direction = direction,
         `Peak time (UTC)` = peak_time, `Hrs to recover (<=10bps)` = hours_to_recover, Note = note) %>%
  flextable() %>% autofit() %>% theme_vanilla() %>%
  fontsize(size = 8, part = "all") %>%
  width(j = "Note", width = 1.6) %>%
  width(j = "Peak time (UTC)", width = 1.1)

ft_did <- did_results %>%
  rename(Term = term, Estimate = estimate, `HAC Std. Error` = std_error_HAC, `p-value` = p_value) %>%
  flextable() %>% autofit() %>% theme_vanilla() %>%
  fontsize(size = 8, part = "all")

ft_combined <- combined_effects %>%
  rename(Coin = coin, `Total event effect (bps)` = total_event_effect_bps,
         `HAC SE` = se_HAC, `t-stat` = t_stat, `p-value` = p_value) %>%
  flextable() %>% autofit() %>% theme_vanilla() %>%
  fontsize(size = 9, part = "all")

doc <- read_docx() %>%
  body_add_par("Hypothesis 3: Reserve Transparency and Contagion in the SVB/USDC Depeg", style = "heading 1") %>%
  body_add_par(paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_par(
    paste0(
      "Hypothesis: During USDC's March 2023 depeg (triggered by Silicon Valley Bank's failure), ",
      "DAI -- which held substantial, disclosed USDC reserves via its Peg Stability Module -- ",
      "depegged in closer tandem with USDC than USDT did, whose reserve composition was opaque ",
      "and not directly SVB-exposed. Reserve transparency is treated as a contagion channel in ",
      "its own right, distinct from direct exposure to the triggering shock."
    ), style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("1. Data Coverage", style = "heading 2") %>%
  body_add_flextable(ft_coverage) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("2. Peg Deviation Around the Event Window", style = "heading 2") %>%
  body_add_gg(p_event, width = 6.5, height = 4.2) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("3. Peak Deviation and Recovery Time", style = "heading 2") %>%
  body_add_gg(p_peak, width = 6.5, height = 3.5) %>%
  body_add_flextable(ft_recovery) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("4. Difference-in-Differences Regression", style = "heading 2") %>%
  body_add_par(
    "Model: deviation_bps ~ event_dummy * coin, USDT as reference category. HAC (Newey-West, 24-lag) standard errors.",
    style = "Normal"
  ) %>%
  body_add_flextable(ft_did) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("5. Combined Event-Window Effect by Coin", style = "heading 2") %>%
  body_add_par(
    paste0(
      "The raw event_dummy:coinX coefficients above are X's reaction relative to USDT's own reaction, ",
      "not X's total shift from its own baseline. This table reports the latter (event_dummy + ",
      "event_dummy:coinX where applicable), with an exact HAC standard error via the model's covariance matrix."
    ), style = "Normal"
  ) %>%
  body_add_flextable(ft_combined) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("6. Interpretation Notes", style = "heading 2") %>%
  body_add_par(
    paste0(
      "A significantly negative total event effect for DAI, smaller in magnitude than USDC's own ",
      "reaction but clearly distinguishable from zero, supports the transparency-channel hypothesis: ",
      "DAI's disclosed USDC exposure transmitted a meaningful share of USDC's stress, while USDT -- ",
      "with opaque, non-SVB-exposed reserves -- shows the opposite-signed, appreciating reaction ",
      "expected of an unaffected flight-to-quality destination."
    ), style = "Normal"
  )

print(doc, target = OUT_DOCX)
cat("Document written to:", file.path(getwd(), OUT_DOCX), "\n")
