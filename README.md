# Lotka–Volterra explorer

This is an interactive teaching app for two-species competition and mutualism, built for BES 550 Advanced Ecology at the University of Rhode Island. It accompanies the lecture on species interactions (competition, mutualism, and parasitism).

**Open the app:** https://jaharvey.github.io/lotka-volterra/

The first load takes 10–20 seconds because the app runs entirely in your browser.

## What it does

The model uses the same sign convention as the lecture:

dN₁/dt = r₁N₁(K₁ − N₁ − αN₂)/K₁  
dN₂/dt = r₂N₂(K₂ − N₂ − βN₁)/K₂

α is the per-capita effect of species 2 on species 1, and β is the reverse. Positive values mean competition and negative values mean mutualism.

## Definitions

Subscript 1 or 2 refers to species 1 or species 2.

| Symbol | Meaning |
|---|---|
| N₁, N₂ | Population size (number of individuals) of each species |
| t | Time, in years |
| dN/dt | Rate of change in population size per year |
| r₁, r₂ | Intrinsic rate of increase: the per-capita growth rate when the population is small and has no competitors. It sets how fast a population changes, not who wins. |
| K₁, K₂ | Carrying capacity: the population size each species reaches on its own, set by its own resources and self-limitation |
| α | Per-capita effect of species 2 on species 1, in species 1 equivalents. α = 0.5 means one individual of species 2 limits species 1 as much as half an individual of species 1. |
| β | Per-capita effect of species 1 on species 2, in species 2 equivalents. Same sign convention as α. |
| N₁(0), N₂(0) | Starting population sizes |
| N₁\*, N₂\* | Equilibrium population sizes, where neither population is changing |
| α × β | Product of the two interaction coefficients. For competition, αβ < 1 is needed for stable coexistence. For mutualism, αβ ≥ 1 means runaway growth. |
| Isocline | Line in the phase plane where one species' growth is zero (dN/dt = 0). Species 1: N₁ = K₁ − αN₂. Species 2: N₂ = K₂ − βN₁. |
| Invasion check | Whether a species can grow when rare while the other sits at its carrying capacity: K₁ > αK₂ for species 1, K₂ > βK₁ for species 2 |
| Parasite cost | Percentage reduction in a species' carrying capacity, a simple stand-in for the physiological cost of infection |

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
