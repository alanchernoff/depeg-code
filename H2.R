## ---- 0. Packages -----------------------------------------------------------
required_pkgs <- c("tidyverse", "lubridate", "officer", "flextable",
                   "zoo", "scales")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(tidyverse)
library(lubridate)
library(officer)
library(flextable)
library(zoo)
library(scales)

## ---- 1. Setup ---------------------------------------------------------------
setwd("C:/Users/acher/research/")
set.seed(1)

DATA_FILE   <- "stablecoin_hourly_prices_wide.csv"
OUT_DOCX    <- "hypothesis2_UST_contagion_results.docx"

CONTAGION_COINS <- c("USDC", "USDT", "DAI", "FRAX", "BUSD")

EVENT_DATE     <- ymd_hms("2022-05-09 00:00:00", tz = "UTC")
PRE_WINDOW     <- c(EVENT_DATE - days(21), EVENT_DATE - hours(1))
ACUTE_WINDOW   <- c(EVENT_DATE, EVENT_DATE + days(7))
POST_WINDOW    <- c(EVENT_DATE + days(7) + hours(1), EVENT_DATE + days(35))
PLOT_WINDOW    <- c(EVENT_DATE - days(21), EVENT_DATE + days(35))

## ---- 2. Load and reshape data -----------------------------------------------
## CONFIRMED from diagnostics: the raw file uses "YYYY-MM-DD HH:MM" (no seconds,
## no offset), and readr's own default guesser already parses this correctly as
## POSIXct with no manual intervention needed. v1/v2's mistake was re-parsing an
## already-correct POSIXct column with ymd_hms(), which (via R's default
## as.character() truncating exact-midnight instants to a bare date) silently
## turned every midnight hour into NA. v3's mistake was forcing the column to
## character and then requiring seconds that don't exist in the source format.
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
  filter(!is.na(price))

df <- df %>% mutate(deviation_bps = (price - 1) * 10000)

## ---- 3. Coverage check for the event window (FIXED: na.rm + diagnostic) ----
coverage_tbl <- df %>%
  filter(coin %in% c(CONTAGION_COINS, "UST")) %>%
  group_by(coin) %>%
  summarise(
    n_na_datetime    = sum(is.na(datetime_utc)),
    first_obs        = suppressWarnings(min(datetime_utc, na.rm = TRUE)),
    last_obs         = suppressWarnings(max(datetime_utc, na.rm = TRUE)),
    n_in_plot_window = sum(datetime_utc >= PLOT_WINDOW[1] & datetime_utc <= PLOT_WINDOW[2], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(coin)

if (any(coverage_tbl$n_na_datetime > 0)) {
  warning("Unparseable datetime_utc values detected -- see n_na_datetime column in coverage_tbl.")
}
print(coverage_tbl)

## ---- 4. Build event-window panel --------------------------------------------
event_panel <- df %>%
  filter(coin %in% CONTAGION_COINS,
         datetime_utc >= PLOT_WINDOW[1], datetime_utc <= PLOT_WINDOW[2]) %>%
  mutate(
    phase = case_when(
      datetime_utc >= PRE_WINDOW[1]   & datetime_utc <= PRE_WINDOW[2]   ~ "Pre-event",
      datetime_utc >= ACUTE_WINDOW[1] & datetime_utc <= ACUTE_WINDOW[2] ~ "Acute collapse",
      datetime_utc >= POST_WINDOW[1]  & datetime_utc <= POST_WINDOW[2]  ~ "Post-event",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(phase))

## ---- 4b. Clipping / data-quality diagnostic (NEW) ---------------------------
clipping_diag <- event_panel %>%
  group_by(coin, phase) %>%
  summarise(
    max_abs_dev  = max(abs(deviation_bps)),
    n_at_max     = sum(abs(deviation_bps) >= max_abs_dev - 1e-9),
    n_obs        = n(),
    share_at_max = round(n_at_max / n_obs, 3),
    .groups = "drop"
  )

flagged <- clipping_diag %>% filter(share_at_max > 0.05, max_abs_dev < 50)
if (nrow(flagged) > 0) {
  message("Possible clipping/quantization detected (>5% of hours sitting exactly at the phase max, ceiling under 50 bps):")
  print(flagged)
}

## ---- 5. Pre vs. post summary stats per coin (long format) -------------------
summary_stats_long <- event_panel %>%
  group_by(coin, phase) %>%
  summarise(
    n_obs           = n(),
    mean_dev_bps    = round(mean(deviation_bps), 2),
    sd_dev_bps      = round(sd(deviation_bps), 2),
    max_abs_dev_bps = round(max(abs(deviation_bps)), 2),
    .groups = "drop"
  ) %>%
  mutate(phase = factor(phase, levels = c("Pre-event", "Acute collapse", "Post-event"))) %>%
  arrange(coin, phase)

print(summary_stats_long)

## ---- 6. Statistical tests: pre vs. acute, and pre vs. post (long format) ----
run_tests <- function(coin_name) {
  pre_vals   <- event_panel %>% filter(coin == coin_name, phase == "Pre-event") %>% pull(deviation_bps)
  acute_vals <- event_panel %>% filter(coin == coin_name, phase == "Acute collapse") %>% pull(deviation_bps)
  post_vals  <- event_panel %>% filter(coin == coin_name, phase == "Post-event") %>% pull(deviation_bps)
  
  safe_test <- function(x, y, window_label) {
    if (length(x) < 2 || length(y) < 2) {
      return(tibble(window = window_label, mean_shift_bps = NA_real_, p_ttest = NA_real_,
                    p_wilcox = NA_real_, var_ratio = NA_real_, p_vartest = NA_real_))
    }
    tt <- t.test(y, x)
    vt <- var.test(y, x)
    wp <- suppressWarnings(wilcox.test(y, x)$p.value)
    tibble(
      window = window_label,
      mean_shift_bps = round(unname(diff(tt$estimate)), 4),
      p_ttest        = round(tt$p.value, 4),
      p_wilcox       = round(wp, 4),
      var_ratio      = round(unname(vt$estimate), 4),
      p_vartest      = round(vt$p.value, 4)
    )
  }
  
  bind_rows(
    safe_test(pre_vals, acute_vals, "Acute collapse"),
    safe_test(pre_vals, post_vals,  "Post-event")
  ) %>% mutate(coin = coin_name, .before = 1)
}

test_results_long <- map_dfr(CONTAGION_COINS, run_tests) %>%
  mutate(window = factor(window, levels = c("Acute collapse", "Post-event"))) %>%
  arrange(coin, window)

print(test_results_long)

## ---- 7. Plot: peg deviation around the event window -------------------------
plot_df <- event_panel %>% mutate(coin = factor(coin, levels = CONTAGION_COINS))

p_event <- ggplot(plot_df, aes(x = datetime_utc, y = deviation_bps)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = as.numeric(EVENT_DATE), linetype = "dotted", color = "firebrick") +
  geom_line(color = "#2b6cb0", linewidth = 0.3) +
  facet_wrap(~coin, scales = "free_y", ncol = 1) +
  labs(
    title = "Peg Deviation Around the UST Collapse (May 2022)",
    subtitle = "Dotted red line = event date (2022-05-09). Dashed line = perfect peg.",
    x = NULL, y = "Deviation from $1.00 (bps)"
  ) +
  scale_x_datetime(date_breaks = "1 week", labels = date_format("%b %d")) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

## ---- 8. Bar chart: mean absolute deviation by phase --------------------------
bar_df <- event_panel %>%
  group_by(coin, phase) %>%
  summarise(mean_abs_dev = mean(abs(deviation_bps)), .groups = "drop") %>%
  mutate(
    coin  = factor(coin, levels = CONTAGION_COINS),
    phase = factor(phase, levels = c("Pre-event", "Acute collapse", "Post-event"))
  )

p_bar <- ggplot(bar_df, aes(x = coin, y = mean_abs_dev, fill = phase)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  labs(
    title = "Mean Absolute Peg Deviation by Event Phase",
    x = NULL, y = "Mean |deviation| (bps)", fill = NULL
  ) +
  scale_fill_manual(values = c("Pre-event" = "#90cdf4", "Acute collapse" = "#e53e3e", "Post-event" = "#68d391")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

## ---- 9. Assemble Word document ----------------------------------------------
ft_coverage <- coverage_tbl %>%
  mutate(first_obs = format(first_obs, "%Y-%m-%d %H:%M"),
         last_obs  = format(last_obs, "%Y-%m-%d %H:%M")) %>%
  select(coin, first_obs, last_obs, n_in_plot_window, n_na_datetime) %>%
  rename(Coin = coin, `First obs.` = first_obs, `Last obs.` = last_obs,
         `Obs. in plot window` = n_in_plot_window, `NA timestamps (full history)` = n_na_datetime) %>%
  flextable() %>% autofit() %>% theme_vanilla()

ft_summary <- summary_stats_long %>%
  rename(Coin = coin, Phase = phase, `N` = n_obs, `Mean dev. (bps)` = mean_dev_bps,
         `SD dev. (bps)` = sd_dev_bps, `Max |dev.| (bps)` = max_abs_dev_bps) %>%
  flextable() %>% autofit() %>% theme_vanilla() %>%
  fontsize(size = 9, part = "all")

ft_tests <- test_results_long %>%
  rename(Coin = coin, Window = window, `Mean shift (bps)` = mean_shift_bps,
         `p (t-test)` = p_ttest, `p (Wilcoxon)` = p_wilcox,
         `Variance ratio` = var_ratio, `p (F-test)` = p_vartest) %>%
  flextable() %>% autofit() %>% theme_vanilla() %>%
  fontsize(size = 9, part = "all")

ft_clipping <- clipping_diag %>%
  rename(Coin = coin, Phase = phase, `Max |dev.| (bps)` = max_abs_dev, `N at max` = n_at_max,
         N = n_obs, `Share at max` = share_at_max) %>%
  flextable() %>% autofit() %>% theme_vanilla() %>%
  fontsize(size = 9, part = "all")

doc <- read_docx() %>%
  body_add_par("Hypothesis 2: Contagion Effects from the UST Collapse (May 2022)", style = "heading 1") %>%
  body_add_par(paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_par(
    paste0(
      "Hypothesis: Following UST's de-peg beginning 2022-05-09, other stablecoins ",
      "(USDT, USDC, DAI, FRAX, BUSD) exhibited abnormally elevated peg deviation ",
      "and/or volatility relative to their pre-event baseline, consistent with ",
      "market-wide contagion rather than UST-idiosyncratic risk."
    ), style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("1. Data Coverage", style = "heading 2") %>%
  body_add_flextable(ft_coverage) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("2. Peg Deviation Around the Event Window", style = "heading 2") %>%
  body_add_gg(p_event, width = 6.5, height = 8) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("3. Mean Absolute Deviation by Phase", style = "heading 2") %>%
  body_add_gg(p_bar, width = 6.5, height = 3.8) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("4. Descriptive Statistics by Phase", style = "heading 2") %>%
  body_add_flextable(ft_summary) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("5. Statistical Tests (Pre-event baseline vs. Acute / Post windows)", style = "heading 2") %>%
  body_add_par(
    "Welch t-test and Wilcoxon rank-sum test on deviation levels; F-test (var.test) on variance ratio.",
    style = "Normal"
  ) %>%
  body_add_flextable(ft_tests) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("6. Data-Quality Diagnostic: Clipping Check", style = "heading 2") %>%
  body_add_par(
    paste0(
      "For each coin and phase, the share of hourly observations sitting exactly at that phase's ",
      "maximum |deviation| value. A high share with a low absolute ceiling (e.g., repeatedly hitting ",
      "exactly 10 bps) suggests the underlying price feed may be rounded, quantized, or clipped rather ",
      "than reflecting genuine market prices, and should be verified against the raw data before drawing ",
      "conclusions about that coin's resilience."
    ), style = "Normal"
  ) %>%
  body_add_flextable(ft_clipping) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("7. Interpretation Notes", style = "heading 2") %>%
  body_add_par(
    paste0(
      "A significant positive mean shift with p < 0.05 in the Acute collapse or Post-event window ",
      "indicates the coin's peg deviation moved abnormally away from baseline during/after the UST ",
      "collapse. A variance ratio > 1 with a significant p (F-test) indicates elevated volatility. ",
      "Interpret BUSD cautiously outside its own data window (BUSD data begins 2022-04-27, shortly ",
      "before the event). Interpret any coin flagged in Section 6 cautiously until the clipping ",
      "concern is resolved against the raw source data."
    ), style = "Normal"
  )

print(doc, target = OUT_DOCX)
cat("Document written to:", file.path(getwd(), OUT_DOCX), "\n")
