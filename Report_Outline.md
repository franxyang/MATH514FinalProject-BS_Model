# Report Outline (Methods & Results)

## Phenomenon & Model Choice
- Problem: Price European call options under Black–Scholes with continuous dividend yield; interpret price and Greeks as hedging signals.
- Why MOL (Chapter 12 context): Spatial discretization converts the PDE into a stiff ODE system, letting us apply implicit solvers (Backward Euler, BDF2) studied in class.

## Discretization & Solvers
- Space: Central differences on \(S \in [0, S_{\max}]\) with boundary conditions \(V(0,\tau)=0\), \(V(S_{\max},\tau)=S_{\max}e^{-q\tau}-Ke^{-r\tau}\).
- Time-stepping: Backward Euler (1st order, A-stable) for startup; BDF2 (2nd order, A-stable multistep) for production. Linear solves via sparse tridiagonal systems.
- Parameter baseline (synthetic): \(S_{\max}=400\), \(K=100\), \(T=1\), \(r=0.03\), \(\sigma=0.20\), \(q=0.01\), grids \(M=400\), \(N=200\).

## Stability & Error Expectations
- Stiffness: Large \(M\) yields tight eigenvalues; implicit schemes avoid step-size restrictions of explicit FDM.
- Consistency: First-order startup, then second-order BDF2; expect temporal EOC \(\approx 2\) once startup transient decays.
- Boundary sensitivity: Error grows if \(S_{\max}\) too small; note effect on deep ITM pricing and Gamma near \(S_{\max}\).

## Convergence Evidence (include in paper)
- Temporal EOC table (fixed \(M=400\), BDF2): report errors for \(N = 50, 100, 200, 400\); compute observed order \(\log(e_i/e_{i+1})/\log(\Delta\tau_i/\Delta\tau_{i+1})\).
- Spatial sensitivity: show one plot/table varying \(M = 100, 200, 400, 800\) at fixed \(N\), note diminishing returns vs. cost.
- Figure: log–log plot of error vs. \(\Delta\tau\) with order-1 and order-2 reference slopes.

## Results Interpretation
- Prices: Compare numerical \(V(S,\tau=T)\) to analytic Black–Scholes; comment on max error location and the role of \(q\).
- Greeks: Discuss Delta monotonicity and Gamma peak near \(S \approx K\); note grid spacing influence on Gamma smoothness.
- Market case study: Report one cached snapshot (AAPL/MSFT/GOOGL/NVDA) with \(S_0, K, T, \sigma, q\) and market mid; show error with/without dividend correction and highlight which scenario reduces error.

## Reproducibility Notes
- Use pinned parameters above for synthetic runs; rely on offline snapshot in `02_RealWorld_CaseStudy.ipynb` when live data is unavailable.
- Cite: Süli & Mayers (2003), Chapter 12 for stiff solvers; note use of AI assistance where applicable.
