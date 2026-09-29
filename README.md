# Lotka–Volterra explorer

This is an interactive teaching app for two-species competition and mutualism, built for BES 550 Advanced Ecology at the University of Rhode Island. It accompanies the lecture on species interactions (competition, mutualism, and parasitism).

**Open the app:** https://YOUR-USERNAME.github.io/lotka-volterra/

The first load takes 10–20 seconds because the app runs entirely in your browser.

## What it does

The model uses the same sign convention as the lecture:

dN₁/dt = r₁N₁(K₁ − N₁ − αN₂)/K₁  
dN₂/dt = r₂N₂(K₂ − N₂ − βN₁)/K₂

α is the per-capita effect of species 2 on species 1, and β is the reverse. Positive values mean competition and negative values mean mutualism.

The app shows:

- **Outcome.** Exclusion, stable coexistence, a priority effect, or runaway mutualism, worked out from the invasion checks (K₁ > αK₂ and K₂ > βK₁).
- **Phase plane.** Both isoclines, a direction field, and the simulated trajectory.
- **Over time.** Population trajectories for both species.
- **Parasite cost.** A slider that lowers either species' carrying capacity, a simple stand-in for the physiological cost of infection.

## How to use it

1. Click a scenario to load an example, or type in your own values.
2. Adjust r, K, α, β and the starting densities.
3. Use the parasite-cost sliders to see how infection changes the outcome.

The *About* tab has questions to work through. Scenario values are illustrative and are not data from real systems.

## Notes

- The simulation uses a fourth-order Runge–Kutta solver with a time step of 0.05 years.
- In runaway mutualism (α and β both negative, with αβ ≥ 1), the simulation stops once a population passes 20 times the larger carrying capacity.

## Running it locally

```r
install.packages("shiny")
shiny::runApp("BES550-LotkaVolterra_App_v1")
```

The web version is built with Shinylive:
`shinylive::export("BES550-LotkaVolterra_App_v1", "lotka-volterra-site")`

## Reference

Holland, J. N. & DeAngelis, D. L. (2010). A consumer–resource approach to the density-dependent population dynamics of mutualism. *Ecology* 91:1286–1295.

## Contact

Johanna Harvey, University of Rhode Island ([Avian Disease Lab](https://aviandiseaselab.com)).
