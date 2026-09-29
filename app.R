# Lotka-Volterra explorer: two-species competition and mutualism
# BES 550 Advanced Ecology, Species interactions
#
# To run: open this file in RStudio and click "Run App",
# or run shiny::runApp("path/to/BES550-LotkaVolterra_App_v1") in the R console.
# The only package needed is shiny: install.packages("shiny")
#
# Model (same sign convention as the lecture slides):
#   dN1/dt = r1 N1 (K1 - N1 - alpha N2) / K1
#   dN2/dt = r2 N2 (K2 - N2 - beta  N1) / K2
# alpha = per-capita effect of species 2 on species 1, beta = effect of 1 on 2.
# alpha, beta > 0: competition.  alpha, beta < 0: mutualism.

library(shiny)

teal <- "#169C9A"; orange <- "#D9622B"; navy <- "#0E2841"; grey <- "#8A8A8A"

# ---------------------------------------------------------------------------
# Scenarios (illustrative values for teaching, not data)
# ---------------------------------------------------------------------------
scenarios <- list(
  exclusion   = list(r1 = 0.8, r2 = 0.6, K1 = 100, K2 = 80,  a = 0.6,  b = 1.4,  N1 = 10, N2 = 10),
  coexist     = list(r1 = 0.5, r2 = 0.9, K1 = 100, K2 = 90,  a = 0.5,  b = 0.6,  N1 = 10, N2 = 10),
  priority    = list(r1 = 0.7, r2 = 0.7, K1 = 100, K2 = 100, a = 1.3,  b = 1.4,  N1 = 30, N2 = 20),
  mutualism   = list(r1 = 0.5, r2 = 0.5, K1 = 100, K2 = 80,  a = -0.3, b = -0.4, N1 = 10, N2 = 10),
  runaway     = list(r1 = 0.5, r2 = 0.5, K1 = 100, K2 = 80,  a = -0.9, b = -1.2, N1 = 10, N2 = 10)
)

# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
lv_rates <- function(N, p) {
  c(p$r1 * N[1] * (p$K1 - N[1] - p$a * N[2]) / p$K1,
    p$r2 * N[2] * (p$K2 - N[2] - p$b * N[1]) / p$K2)
}

# Fourth-order Runge-Kutta. Stops early if either population passes `cap`
# (20 times the larger K by default), which only happens in runaway mutualism.
# Stopping there also avoids the numerical blow-up of the explicit solver.
lv_solve <- function(p, N0, tmax = 100, dt = 0.05, cap = 20 * max(p$K1, p$K2)) {
  steps <- ceiling(tmax / dt)
  out <- matrix(NA_real_, nrow = steps + 1, ncol = 3)
  N <- pmax(N0, 0)
  out[1, ] <- c(0, N)
  runaway <- FALSE
  for (i in seq_len(steps)) {
    k1 <- lv_rates(N, p)
    k2 <- lv_rates(N + dt / 2 * k1, p)
    k3 <- lv_rates(N + dt / 2 * k2, p)
    k4 <- lv_rates(N + dt * k3, p)
    Nnew <- N + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    if (any(!is.finite(Nnew)) || any(Nnew > cap)) { runaway <- TRUE; break }
    N <- pmax(Nnew, 0)
    out[i + 1, ] <- c(i * dt, N)
  }
  out <- out[stats::complete.cases(out), , drop = FALSE]
  list(traj = data.frame(t = out[, 1], N1 = out[, 2], N2 = out[, 3]), runaway = runaway)
}

# Interior equilibrium and outcome, from the invasion criteria.
lv_outcome <- function(p) {
  inv1 <- p$K1 > p$a * p$K2     # species 1 can grow when rare, with species 2 at K2
  inv2 <- p$K2 > p$b * p$K1     # species 2 can grow when rare, with species 1 at K1
  den  <- 1 - p$a * p$b
  eq <- if (abs(den) > 1e-9) c((p$K1 - p$a * p$K2) / den, (p$K2 - p$b * p$K1) / den) else c(NA, NA)

  if (p$a < 0 && p$b < 0 && den <= 0) {
    return(list(code = "runaway", eq = c(NA, NA),
                text = "Mutualism with α × β ≥ 1: each species boosts the other faster than self-limitation can hold it back, so both populations grow without bound."))
  }
  if (inv1 && inv2 && den > 0 && all(eq > 0)) {
    txt <- "Stable coexistence. Each species limits itself more than it limits the other, so each can invade when rare."
    if (p$a < 0 || p$b < 0)
      txt <- "Stable coexistence. At least one species benefits from the other, and both settle at an equilibrium."
    return(list(code = "coexist", eq = eq, text = txt))
  }
  if (!inv1 && !inv2) {
    return(list(code = "priority", eq = eq,
                text = "Priority effect. Neither species can invade the other, so the winner depends on starting densities."))
  }
  if (inv1 && !inv2) return(list(code = "sp1", eq = c(p$K1, 0), text = "Species 1 excludes species 2."))
  if (inv2 && !inv1) return(list(code = "sp2", eq = c(0, p$K2), text = "Species 2 excludes species 1."))
  list(code = "other", eq = eq, text = "Outcome is at a boundary between cases. Try nudging a parameter.")
}

# Interaction type from the signs of the effects, as on Pringle's (2016) interaction compass.
# The effect of species 2 on species 1 is -alpha; the effect of species 1 on species 2 is -beta.
interaction_type <- function(a, b, tol = 1e-9) {
  s <- function(x) if (abs(x) < tol) "0" else if (-x > 0) "+" else "-"
  key <- paste0(s(a), "/", s(b))
  lab <- switch(key,
    "+/+" = "Mutualism: each species benefits the other.",
    "-/-" = "Competition: each species harms the other.",
    "+/-" = "Exploitation (like predation or parasitism): species 1 benefits, species 2 is harmed.",
    "-/+" = "Exploitation (like predation or parasitism): species 2 benefits, species 1 is harmed.",
    "+/0" = "Commensalism: species 1 benefits, species 2 is unaffected.",
    "0/+" = "Commensalism: species 2 benefits, species 1 is unaffected.",
    "-/0" = "Amensalism: species 1 is harmed, species 2 is unaffected.",
    "0/-" = "Amensalism: species 2 is harmed, species 1 is unaffected.",
    "0/0" = "Neutral: neither species affects the other.")
  list(key = gsub("-", "\u2212", key, fixed = TRUE), label = lab)
}

fmt <- function(x, d = 1) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
num <- function(id, label, value, step) numericInput(id, label, value = value, step = step)

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-size: 15px; }
    .stat-note { color: #4B5563; font-size: 14px; }
    .scen .btn { margin: 0 6px 6px 0; }
    h4 { margin-top: 16px; }
    .outcome { font-size: 17px; font-weight: 600; color: #0E2841; }
  "))),
  titlePanel("Lotka-Volterra explorer: competition and mutualism"),

  sidebarLayout(
    sidebarPanel(width = 4,
      h4("Scenarios"),
      div(class = "scen",
        actionButton("sc_excl", "Exclusion"),
        actionButton("sc_coex", "Coexistence"),
        actionButton("sc_prio", "Priority effect"),
        actionButton("sc_mut", "Mutualism"),
        actionButton("sc_run", "Runaway mutualism")
      ),
      p(class = "stat-note", "Scenarios fill in the values below. The numbers are illustrative."),

      h4("Species 1"),
      fluidRow(column(4, num("r1", "r₁", 0.8, 0.05)),
               column(4, num("K1", "K₁", 100, 5)),
               column(4, num("N10", "N₁(0)", 10, 1))),
      h4("Species 2"),
      fluidRow(column(4, num("r2", "r₂", 0.6, 0.05)),
               column(4, num("K2", "K₂", 80, 5)),
               column(4, num("N20", "N₂(0)", 10, 1))),
      h4("Interaction coefficients"),
      fluidRow(column(6, num("a", "α (effect of 2 on 1)", 0.6, 0.05)),
               column(6, num("b", "β (effect of 1 on 2)", 1.4, 0.05))),
      p(class = "stat-note", "Positive values = competition. Negative values = mutualism."),
      div(class = "stat-note", style = "background:#F2F5F7; padding:8px 10px; border-radius:6px; margin-bottom:10px;",
        HTML(paste0(
          "<b>r</b>: intrinsic rate of increase, the per-capita growth rate when the population is small and alone. It sets how fast a population changes, not who wins.<br>",
          "<b>K</b>: carrying capacity, the population size a species reaches on its own.<br>",
          "<b>α</b>: per-capita effect of species 2 on species 1, in species 1 equivalents (α = 0.5: one individual of species 2 limits species 1 as much as half an individual of species 1).<br>",
          "<b>β</b>: per-capita effect of species 1 on species 2, in species 2 equivalents.<br>",
          "Full list in the Definitions tab."))),

      h4("Organismal modifier: parasite cost"),
      sliderInput("par1", "Infection lowers K₁ by (%)", min = 0, max = 90, value = 0, step = 1),
      sliderInput("par2", "Infection lowers K₂ by (%)", min = 0, max = 90, value = 0, step = 1),
      numericInput("tmax", "Years to simulate", value = 100, min = 10, max = 500, step = 10)
    ),

    mainPanel(width = 8,
      fluidRow(
        column(6, h4("Outcome"), uiOutput("outcome")),
        column(6, h4("Invasion checks"), tableOutput("checks"))
      ),
      tabsetPanel(
        tabPanel("Phase plane", plotOutput("phase", height = "460px"),
                 p(class = "stat-note", "Isoclines show where each species stops growing. Grey arrows show the direction of change; the orange path is the simulated trajectory from your starting densities.")),
        tabPanel("Interaction compass", plotOutput("compass", height = "460px"),
                 p(class = "stat-note", "After Pringle (2016). Each axis is the per-capita effect one species has on the other. Because Lotka-Volterra subtracts \u03b1 and \u03b2, the effect of species 2 on species 1 is \u2212\u03b1, so positive \u03b1 (competition) plots on the negative side. Points farther from the center are stronger interactions.")),
        tabPanel("Over time", plotOutput("timeplot", height = "420px")),
        tabPanel("Definitions", br(), tableOutput("defs"),
                 p(class = "stat-note", "Subscript 1 or 2 refers to species 1 or species 2.")),
        tabPanel("About", uiOutput("about"))
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  set_sc <- function(s) {
    updateNumericInput(session, "r1", value = s$r1); updateNumericInput(session, "r2", value = s$r2)
    updateNumericInput(session, "K1", value = s$K1); updateNumericInput(session, "K2", value = s$K2)
    updateNumericInput(session, "a",  value = s$a);  updateNumericInput(session, "b",  value = s$b)
    updateNumericInput(session, "N10", value = s$N1); updateNumericInput(session, "N20", value = s$N2)
    updateSliderInput(session, "par1", value = 0);   updateSliderInput(session, "par2", value = 0)
  }
  observeEvent(input$sc_excl, set_sc(scenarios$exclusion))
  observeEvent(input$sc_coex, set_sc(scenarios$coexist))
  observeEvent(input$sc_prio, set_sc(scenarios$priority))
  observeEvent(input$sc_mut,  set_sc(scenarios$mutualism))
  observeEvent(input$sc_run,  set_sc(scenarios$runaway))

  val <- function(x, default) if (is.null(x) || is.na(x)) default else x

  pars <- reactive({
    K1 <- max(val(input$K1, 100), 1) * (1 - val(input$par1, 0) / 100)
    K2 <- max(val(input$K2, 80), 1)  * (1 - val(input$par2, 0) / 100)
    list(r1 = max(val(input$r1, 0.8), 0.001), r2 = max(val(input$r2, 0.6), 0.001),
         K1 = max(K1, 0.1), K2 = max(K2, 0.1),
         a = val(input$a, 0.6), b = val(input$b, 1.4))
  })
  N0  <- reactive(c(max(val(input$N10, 10), 0), max(val(input$N20, 10), 0)))
  sim <- reactive(lv_solve(pars(), N0(), tmax = min(max(val(input$tmax, 100), 10), 500)))
  res <- reactive(lv_outcome(pars()))

  # --- Outcome ---
  output$outcome <- renderUI({
    o <- res(); s <- sim(); p <- pars()
    last <- tail(s$traj, 1)
    it <- interaction_type(p$a, p$b)
    tagList(
      p(HTML(paste0("<b>Interaction (effect on species 1 / species 2): ", it$key, "</b>. ", it$label))),
      p(class = "outcome", o$text),
      if (o$code %in% c("coexist", "priority") && all(!is.na(o$eq)) && all(o$eq > 0))
        p(HTML(paste0("Interior equilibrium: N₁* = <b>", fmt(o$eq[1]), "</b>, N₂* = <b>", fmt(o$eq[2]), "</b>",
                      if (o$code == "priority") " (unstable)" else ""))),
      p(class = "stat-note",
        if (s$runaway) paste0("The simulation stopped at year ", fmt(last$t, 1), " because a population passed 20 times the larger carrying capacity.")
        else paste0("After ", fmt(last$t, 0), " years: N₁ = ", fmt(last$N1), ", N₂ = ", fmt(last$N2), ".")),
      if (input$par1 > 0 || input$par2 > 0)
        p(class = "stat-note", paste0("With the parasite cost, K₁ = ", fmt(p$K1), " and K₂ = ", fmt(p$K2), "."))
    )
  })

  output$checks <- renderTable({
    p <- pars()
    data.frame(
      Check = c("Species 1 invades? K₁ > αK₂", "Species 2 invades? K₂ > βK₁", "α × β"),
      Values = c(paste0(fmt(p$K1), " vs ", fmt(p$a * p$K2)),
                 paste0(fmt(p$K2), " vs ", fmt(p$b * p$K1)),
                 fmt(p$a * p$b, 2)),
      Result = c(ifelse(p$K1 > p$a * p$K2, "Yes", "No"),
                 ifelse(p$K2 > p$b * p$K1, "Yes", "No"),
                 ifelse(p$a * p$b < 1, "< 1", "≥ 1")),
      check.names = FALSE
    )
  }, striped = TRUE, spacing = "s")

  # --- Phase plane ---
  output$phase <- renderPlot({
    p <- pars(); s <- sim(); o <- res(); tr <- s$traj
    xmax <- max(p$K1, p$K2 / max(p$b, 1e-6) * (p$b > 0), max(tr$N1, na.rm = TRUE), N0()[1], 10) * 1.1
    ymax <- max(p$K2, p$K1 / max(p$a, 1e-6) * (p$a > 0), max(tr$N2, na.rm = TRUE), N0()[2], 10) * 1.1
    xmax <- min(xmax, 5 * max(p$K1, p$K2)); ymax <- min(ymax, 5 * max(p$K1, p$K2))
    par(mar = c(4.5, 4.5, 1.5, 1), cex = 1.05)
    plot(NA, xlim = c(0, xmax), ylim = c(0, ymax), xaxs = "i", yaxs = "i",
         xlab = expression(N[1]), ylab = expression(N[2]))
    # direction field
    gx <- seq(xmax * 0.05, xmax * 0.95, length.out = 12)
    gy <- seq(ymax * 0.05, ymax * 0.95, length.out = 12)
    for (x in gx) for (y in gy) {
      d <- lv_rates(c(x, y), p)
      u <- d[1] / xmax; v <- d[2] / ymax; m <- sqrt(u^2 + v^2)   # direction in plot units
      if (is.finite(m) && m > 0)
        arrows(x, y, x + 0.035 * xmax * u / m, y + 0.035 * ymax * v / m,
               length = 0.05, col = "grey75")
    }
    # isoclines: N1 = K1 - a N2 ; N2 = K2 - b N1
    yy <- seq(0, ymax, length.out = 200); xx <- seq(0, xmax, length.out = 200)
    lines(p$K1 - p$a * yy, yy, col = navy, lwd = 3)
    lines(xx, p$K2 - p$b * xx, col = teal, lwd = 3)
    lines(tr$N1, tr$N2, col = orange, lwd = 2.5)
    points(N0()[1], N0()[2], pch = 21, bg = orange, col = "white", cex = 1.8)
    if (all(!is.na(o$eq)) && all(o$eq >= 0) && o$code != "runaway")
      points(o$eq[1], o$eq[2], pch = if (o$code == "priority") 1 else 19, cex = 2, lwd = 2, col = navy)
    legend("topright", c("Species 1 isocline", "Species 2 isocline", "Trajectory"),
           col = c(navy, teal, orange), lwd = c(3, 3, 2.5), bty = "n", bg = "white")
  })

  # --- Interaction compass ---
  output$compass <- renderPlot({
    p <- pars()
    x <- -p$a; y <- -p$b
    lim <- max(1.5, abs(x), abs(y)) * 1.15
    par(mar = c(4.5, 4.5, 1.5, 1), cex = 1.05, pty = "s")
    plot(NA, xlim = c(-lim, lim), ylim = c(-lim, lim), asp = 1,
         xlab = expression("Effect of species 2 on species 1  (" * -alpha * ")"),
         ylab = expression("Effect of species 1 on species 2  (" * -beta * ")"))
    rect(0, 0, lim * 2, lim * 2, col = adjustcolor(teal, 0.10), border = NA)
    rect(-lim * 2, -lim * 2, 0, 0, col = adjustcolor(navy, 0.08), border = NA)
    rect(0, -lim * 2, lim * 2, 0, col = adjustcolor(orange, 0.08), border = NA)
    rect(-lim * 2, 0, 0, lim * 2, col = adjustcolor(orange, 0.08), border = NA)
    for (rr in seq(0.5, lim, by = 0.5))
      symbols(0, 0, circles = rr, inches = FALSE, add = TRUE, fg = "grey85")
    abline(h = 0, v = 0, col = "grey40", lwd = 1.5)
    k <- lim * 0.62
    text( k,  k, "Mutualism\n(+/+)", col = teal, font = 2)
    text(-k, -k, "Competition\n(\u2212/\u2212)", col = navy, font = 2)
    text( k, -k, "Exploitation\n(+/\u2212)\nsp. 1 gains", col = orange, font = 2)
    text(-k,  k, "Exploitation\n(\u2212/+)\nsp. 2 gains", col = orange, font = 2)
    text(lim * 0.97, 0, "Commensalism /\namensalism on axes", adj = c(1, -0.4), cex = 0.8, col = "grey30")
    arrows(0, 0, x, y, length = 0.12, lwd = 2, col = "grey30")
    points(x, y, pch = 21, bg = orange, col = "white", cex = 2.6, lwd = 2)
    text(x, y, paste0("(", fmt(x, 2), ", ", fmt(y, 2), ")"), pos = 4, cex = 0.95)
  })

  # --- Time series ---
  output$timeplot <- renderPlot({
    tr <- sim()$traj
    par(mar = c(4.5, 4.5, 1.5, 1), cex = 1.05)
    yl <- c(0, max(c(tr$N1, tr$N2, 1), na.rm = TRUE) * 1.1)
    plot(tr$t, tr$N1, type = "l", lwd = 3, col = navy, ylim = yl, xlab = "Years", ylab = "Population size")
    lines(tr$t, tr$N2, lwd = 3, col = teal)
    p <- pars()
    abline(h = p$K1, col = navy, lty = 3); abline(h = p$K2, col = teal, lty = 3)
    legend("right", c("Species 1", "Species 2", "K₁", "K₂"), col = c(navy, teal, navy, teal),
           lty = c(1, 1, 3, 3), lwd = c(3, 3, 1, 1), bty = "n")
  })

  # --- Definitions ---
  output$defs <- renderTable(data.frame(
    Symbol = c("N₁, N₂", "t", "dN/dt", "r₁, r₂", "K₁, K₂",
               "α", "β", "N₁(0), N₂(0)", "N₁*, N₂*",
               "α × β", "Isocline", "Invasion check", "Parasite cost", "Interaction compass"),
    Meaning = c(
      "Population size (number of individuals) of each species.",
      "Time, in years.",
      "Rate of change in population size per year.",
      "Intrinsic rate of increase: per-capita growth rate when the population is small and has no competitors. Sets how fast a population changes, not who wins.",
      "Carrying capacity: the population size each species reaches alone, set by its own resources and self-limitation.",
      "Per-capita effect of species 2 on species 1, measured in species 1 equivalents. α = 0.5 means one individual of species 2 limits species 1 as much as half an individual of species 1. Positive = competition, negative = mutualism.",
      "Per-capita effect of species 1 on species 2, measured in species 2 equivalents. Same sign convention as α.",
      "Starting population sizes.",
      "Equilibrium population sizes, where neither population is changing.",
      "Product of the two interaction coefficients. For competition, αβ < 1 is needed for stable coexistence. For mutualism, αβ ≥ 1 means runaway growth.",
      "Line in the phase plane where one species' growth is zero (dN/dt = 0). Species 1: N₁ = K₁ − αN₂. Species 2: N₂ = K₂ − βN₁.",
      "Whether a species can grow when rare while the other sits at its carrying capacity: K₁ > αK₂ for species 1, K₂ > βK₁ for species 2.",
      "Percentage reduction in a species' carrying capacity, a simple stand-in for the physiological cost of infection.",
      "Plot of the two per-capita effects, \u2212\u03b1 and \u2212\u03b2 (Pringle 2016). The sign of each effect sets the type of interaction, and the distance from the center sets its strength."),
    check.names = FALSE
  ), striped = TRUE, spacing = "s")

  # --- About ---
  output$about <- renderUI({
    tagList(
      h4("The model"),
      p("dN₁/dt = r₁N₁(K₁ − N₁ − αN₂)/K₁ and dN₂/dt = r₂N₂(K₂ − N₂ − βN₁)/K₂. α is the per-capita effect of species 2 on species 1, measured in units of species 1 individuals; β is the reverse."),
      h4("Reading the outcome"),
      p("A species can invade when rare if its K is larger than the competitive pressure from the other species at its K (K₁ > αK₂ for species 1). Both invade: stable coexistence. Neither: a priority effect, where the starting densities decide. One: that species wins. With mutualism (negative α and β), populations settle above their K values if αβ < 1 and grow without bound if αβ ≥ 1."),
      h4("The parasite cost slider"),
      p("Infection is represented, very simply, as a reduction in a species' carrying capacity. It is a stand-in for the physiological costs of infection, such as the anemia in Schall's Anolis."),
      h4("Questions to try"),
      tags$ol(
        tags$li("Click Exclusion. Double r₂. Does the outcome change? What changes?"),
        tags$li("Stay on Exclusion and raise the parasite cost on species 1. At what cost does species 2 persist? At what cost does species 1 disappear?"),
        tags$li("Click Priority effect, then swap the starting densities. Why does the winner change?"),
        tags$li("Click Mutualism, then Runaway mutualism. What changed in α × β? Why don't real mutualists behave like the runaway case?")
      ),
      h4("Reference"),
      p("Pringle, E. G. 2016. Orienting the interaction compass: resource availability as a major driver of context dependence. PLoS Biology 14:e2000891."),
      p("Holland, J. N. & DeAngelis, D. L. 2010. A consumer-resource approach to the density-dependent population dynamics of mutualism. Ecology 91:1286-1295.")
    )
  })
}

shinyApp(ui, server)
