# ns_report.R: the look of the charts and tables of som_next_steps.Rmd, and the drawing of maps
# (the same palette and chart style as som_variants_summary.Rmd)

PAL  <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948", "#7a5c3a", "#00a3c4",
          "#b0b000", "#c04ad6", "#5a8a00", "#d66a00", "#2a2a9a", "#9a2a2a", "#2a9a9a")
INK  <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
FONT <- "system-ui, -apple-system, 'Segoe UI', sans-serif"

ax <- function(title = NULL, ...) list(title = list(text = title, font = list(color = INK2)), gridcolor = GRID,
                                       zerolinecolor = AXIS, linecolor = AXIS, tickfont = list(color = INK2), ...)
chart <- function(p, ...) {
  lay <- modifyList(list(font = list(family = FONT, size = 13, color = INK2), paper_bgcolor = "#ffffff", plot_bgcolor = "#ffffff",
                         hoverlabel = list(bgcolor = "#ffffff", bordercolor = AXIS, font = list(family = FONT, size = 13, color = INK)),
                         legend = list(orientation = "h", x = 0, y = -0.18, font = list(color = INK2))), list(...))
  do.call(plotly::layout, c(list(p), lay)) |>
    plotly::config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d", "toggleSpikelines"))
}
tab <- function(d, search = nrow(d) > 12, digits = NULL, page = FALSE, ...) {
  dt <- DT::datatable(d, rownames = FALSE, class = "compact hover stripe", escape = FALSE,
                      options = list(dom = if (page) "ftip" else if (search) "ft" else "t", paging = page, pageLength = 15,
                                     autoWidth = FALSE, ...))
  if (!is.null(digits)) for (cn in names(digits)) dt <- DT::formatRound(dt, cn, digits[[cn]])
  dt
}
# a table with a row of dropdown filters above it, one per chosen column, wired to the DataTables API
# (no selectize or slider library needed, so it renders the same everywhere); the filter columns must
# hold plain text
filtered_table <- function(d, filter_cols, id, signif_cols = NULL, round_cols = NULL, page_length = 15, dom = "tip",
                           escape = TRUE, column_defs = NULL, ...) {
  sel <- lapply(filter_cols, function(cn) {
    vals <- unique(as.character(d[[cn]])); if (is.factor(d[[cn]])) vals <- levels(droplevels(d[[cn]])) else vals <- sort(vals)
    htmltools::tags$label(style = "margin: 0 16px 6px 0; font-size: 13px; color: #52514e; display: inline-block;", paste0(cn, ": "),
      htmltools::tags$select(`data-col` = match(cn, names(d)) - 1, style = "font-size: 13px; padding: 2px 4px; max-width: 260px;",
                             htmltools::tags$option(value = "", "all"), lapply(vals, function(v) htmltools::tags$option(value = v, v))))
  })
  bar <- htmltools::tags$div(id = paste0(id, "-filters"), style = "margin: 6px 0 4px;", sel)
  cb <- DT::JS(sprintf(paste0(
    "$('#%s-filters select').on('change', function() {",
    "  var v = $(this).val(); var col = $(this).data('col');",
    "  table.column(col).search(v ? '^' + $.fn.dataTable.util.escapeRegex(v) + '$' : '', true, false).draw();",
    "}); return table;"), id))
  dt <- DT::datatable(d, rownames = FALSE, elementId = id, class = "compact hover stripe", escape = escape, callback = cb,
                      options = c(list(dom = dom, pageLength = page_length, autoWidth = FALSE, scrollX = TRUE, columnDefs = column_defs), list(...)))
  if (!is.null(signif_cols)) dt <- DT::formatSignif(dt, signif_cols, 3)
  if (!is.null(round_cols)) dt <- DT::formatRound(dt, round_cols, 3)
  htmltools::tagList(bar, dt)
}
f2 <- function(x) sprintf("%.2f", x); f3 <- function(x) sprintf("%.3f", x); pct <- function(x) sprintf("%.0f%%", 100 * x)
menu <- function(buttons, x = 0, y = 1.03, type = "buttons", active = 0)
  list(type = type, direction = "right", x = x, y = y, xanchor = "left", yanchor = "bottom", showactive = TRUE, active = active,
       buttons = buttons, pad = list(r = 4, t = 2), bgcolor = "#ffffff", bordercolor = AXIS, font = list(size = 12, color = INK))

# ------------------------------------------------------------------------------------------- map drawing
# the outline of every cell as one x/y path (NA lifts the pen), hexagons on a hex lattice, squares otherwise
cell_paths <- function(grid) {
  pos <- grid$pos
  if (grid$topo == "hex") { ang <- (0:6 * 60 + 30) * pi / 180; shape <- cbind(cos(ang), sin(ang)) / sqrt(3) }
  else { shape <- cbind(c(-1, 1, 1, -1, -1), c(-1, -1, 1, 1, -1)) * 0.5 }
  x <- y <- c()
  for (k in seq_len(nrow(pos))) { x <- c(x, pos[k, 1] + shape[, 1], NA); y <- c(y, pos[k, 2] + shape[, 2], NA) }
  data.frame(x = x, y = y)
}
# the prototypes as small curves in their cells, one channel at a time; at most n_show points per curve
# and three decimals, to keep the page small
proto_curves <- function(M, grid, L, nch, channel = 1, y_range = range(M), width = 0.74, height = 0.62, n_show = 40) {
  keep <- unique(round(seq(1, L, length.out = min(L, n_show))))
  tt <- seq(-width / 2, width / 2, length.out = L)[keep]
  Y <- M[, (channel - 1) * L + keep, drop = FALSE]
  x <- y <- c()
  for (k in seq_len(nrow(M))) {
    x <- c(x, grid$pos[k, 1] + tt, NA)
    y <- c(y, grid$pos[k, 2] + ((Y[k, ] - y_range[1]) / diff(y_range) - 0.5) * height, NA)
  }
  data.frame(x = round(x, 3), y = round(y, 3))
}
# colour of every cell: the class of most of the items it wins, white when it wins none
cell_fill <- function(bmu, cls, K, class_col, alpha = 0.55) {
  main <- tapply(cls, bmu, function(v) names(which.max(table(v))))
  fill <- rep("#ffffff", K)
  fill[as.integer(names(main))] <- grDevices::adjustcolor(class_col[main], alpha)
  fill
}
class_colours <- function(cls) setNames(PAL[seq_len(nlevels(cls))], levels(cls))
