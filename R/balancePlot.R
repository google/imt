# Copyright 2024 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Create a Baseline Balance Plot.
#'
#' @param data tibble produced with `checkBaseline`.
#'
#' @return ggplot2 baseline balance plot.
#' @export
#'
#' @examples
#' \donttest{
#' library(imt)
#' set.seed(123)  # for reproducibility
#' N <- 1000
#' fake_data <- tibble::tibble(x1 = rnorm(N), x2 = rnorm(N), t = rbinom(N, 1, 0.5))
#' baseline <- checkBaseline(data = fake_data, variables = c("x1", "x2"), treatment = "t")
#' balancePlot(data = baseline)
#' }
balancePlot <-
  function(data) {
    # The following values come from the standard definition for baseline
    # equivalence. The goal is to define the range [min_x, max_x] for plot.
    # More information can be found at
    # https://ies.ed.gov/ncee/wwc/Docs/referenceresources/wwc_brief_baseline_080715.pdf
    # Compute the x-axis interval [min_x, max_x] so it is symmetric around 0.
    a <- -0.5
    b <- 0.5
    min_x <- round(min(data$std_diff, na.rm = TRUE), 1) - 0.1
    max_x <- round(max(data$std_diff, na.rm = TRUE), 1) + 0.1
    if (min_x > a) {
      min_x <- a
    }
    if (max_x < b) {
      max_x <- b
    }
    absolute_max <- max(abs(c(min_x, max_x)))
    min_x <- -absolute_max
    max_x <- absolute_max

    ncovs <- nrow(data)
    regions <- tibble::tibble(
      xmin = c(min_x, -0.25, -0.05, 0.05, 0.25),
      xmax = c(-0.25, -0.05, 0.05, 0.25, max_x),
      balance = factor(
        c(
          "Very Concerned",
          "Concerned",
          "Not concerned",
          "Concerned",
          "Very Concerned"
        ),
        levels = c("Not concerned", "Concerned", "Very Concerned")
      )
    )
    std_diff_plot <- ggplot2::ggplot() +
      # regions
      ggplot2::geom_rect(
        data = regions,
        ggplot2::aes(
          xmin = xmin,
          xmax = xmax,
          ymin = -Inf,
          ymax = Inf,
          fill = balance
        ),
        alpha = 0.9
      ) +
      ggplot2::scale_fill_manual(values = c("#f0f0f0", "#bdbdbd", "#636363")) +
      ggplot2::scale_alpha(guide = "none") +
      # scatters
      ggplot2::geom_point(
        data = data,
        ggplot2::aes(
          x = std_diff,
          y = variables
        ),
        size = 4
      ) +
      ggplot2::xlab("Standardized Difference") +
      ggplot2::scale_x_continuous(
        breaks = c(-0.9, -0.7, -0.5, -0.25, -0.05, 0.05, 0.25, 0.5, 0.7, 0.9),
        limits = c(-absolute_max, absolute_max)
      ) +
      ggplot2::ylab("Variable") +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "longdash",
        colour = "gray"
      ) +
      ggplot2::coord_cartesian(xlim = c(min_x, max_x)) +
      # Base
      ggplot2::theme_bw(base_size = 14, base_family = "sans") +
      # Set the entire chart region to a light gray color
      ggplot2::theme(
        panel.background = ggplot2::element_rect(
          fill = "#FFFFFF",
          color = "#FFFFFF"
        ),
        plot.background = ggplot2::element_rect(
          fill = "#FFFFFF",
          color = "#FFFFFF"
        ),
        panel.border = ggplot2::element_rect(color = "#F0F0F0"),
        strip.background = ggplot2::element_rect(
          fill = "#F0F0F0",
          color = "#F0F0F0"
        )
      ) +
      # Grid
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(
          color = "#D9D9D9",
          size = 0.25
        ),
        panel.grid.minor = ggplot2::element_line(
          color = "#D9D9D9",
          size = 0.15,
          linetype = "dashed"
        ),
        axis.ticks = ggplot2::element_blank()
      ) +
      # Legend
      ggplot2::theme(
        legend.background = ggplot2::element_rect(
          fill = "#FFFFFF",
          color = "#FFFFFF"
        ),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.key = ggplot2::element_rect(color = "#FFFFFF", fill = "#FFFFFF")
      ) +
      # Title, and axis labels
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          color = "#000000",
          hjust = 0,
          size = ggplot2::rel(1.5),
          face = "bold"
        ),
        axis.text.x = ggplot2::element_text(color = "#737373"),
        axis.text.y = ggplot2::element_text(color = "#737373"),
        axis.title.x = ggplot2::element_text(color = "#525252"),
        axis.title.y = ggplot2::element_text(color = "#525252", angle = 0)
      ) +
      # Margins
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(1, 1, 1, 1), "lines")
      )

    # Facet if "group" column exists (NEW)
    if ("group" %in% names(data)) {
      std_diff_plot <- std_diff_plot +
        ggplot2::facet_wrap(~group, ncol = 2)
    }

    return(std_diff_plot)
  }
