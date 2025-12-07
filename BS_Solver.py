"""
-------------------------------------------------------------------------------
File Name: BS_Solver.py
Author: Yifan Yang
Date Revised: December 7, 2025
Course: MATH 514 - Numerical Analysis (Fall 2025)
Institution: University of Wisconsin–Madison

Description
-----------
Core numerical engine for pricing European Call Options using the
Method of Lines (MOL). The generalized Black–Scholes PDE with
continuous dividend yield q

    ∂V/∂τ = 1/2 σ² S² ∂²V/∂S² + (r − q) S ∂V/∂S − r V

is semi–discretized in space (central differences on [0, S_max]) to
obtain a stiff ODE system

    dV/dτ = A V + b(τ),

which is integrated in (time–to–expiry) τ using either

    • Backward Euler  (order 1, A–stable), or
    • BDF2            (order 2, A–stable, with BE start-up).

Contents
--------
1. BSParams         : Data class for model parameters.
2. exact_solution   : Closed-form Black–Scholes formula (with dividends).
3. construct_system : Spatial discretization (central differences).
4. get_boundary_vector : Right-hand side vector b(τ) from boundary data.
5. solve_ode_system : Time-stepping solver (Backward Euler / BDF2).

Typical Usage
-------------
    from BS_Solver import BSParams, exact_solution, solve_ode_system

    params = BSParams(S_max=500.0, K=100.0, T=1.0, r=0.03, sigma=0.30, q=0.0)
    S_inner, V_num, dS = solve_ode_system("BackwardEuler", N=400, M=400, p=params)
-------------------------------------------------------------------------------
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.stats import norm


# --------------------------------------------------------------------------- #
# 1. Parameter container
# --------------------------------------------------------------------------- #

@dataclass
class BSParams:
    """
    Container for Black–Scholes model and numerical domain parameters.

    Attributes
    ----------
    S_max : float
        Maximum underlying price used as the right boundary S_max.
    K : float
        Strike price.
    T : float
        Time to expiry (years). Present time is τ = T, maturity is τ = 0.
    r : float
        Risk-free interest rate.
    sigma : float
        Volatility of the underlying.
    q : float
        Continuous dividend yield.
    """
    S_max: float = 500.0
    K: float = 100.0
    T: float = 1.0
    r: float = 0.03
    sigma: float = 0.30
    q: float = 0.0


# --------------------------------------------------------------------------- #
# 2. Closed-form solution (for verification)
# --------------------------------------------------------------------------- #

def exact_solution(S: np.ndarray | float, tau: float, p: BSParams) -> np.ndarray | float:
    """
    Generalized Black–Scholes call option closed-form solution with dividends.

    Parameters
    ----------
    S : array_like or float
        Underlying asset price(s).
    tau : float
        Time to expiry. For pricing at the current time, tau = p.T.
        Note: tau >= 0.  The limit tau → 0 is handled explicitly.
    p : BSParams
        Parameter object.

    Returns
    -------
    V : ndarray or float
        Option value(s) at time-to-expiry tau.
    """
    # Handle maturity (tau = 0) explicitly to avoid division-by-zero issues.
    if tau <= 0.0:
        if isinstance(S, np.ndarray):
            return np.maximum(S - p.K, 0.0)
        else:
            return max(float(S) - p.K, 0.0)

    # Avoid log(0) issues near S = 0
    S_arr = np.asarray(S, dtype=float)
    S_safe = np.where(S_arr > 1e-12, S_arr, 1e-12)

    sqrt_tau = np.sqrt(tau)
    d1 = (np.log(S_safe / p.K) + (p.r - p.q + 0.5 * p.sigma**2) * tau) / (p.sigma * sqrt_tau)
    d2 = d1 - p.sigma * sqrt_tau

    # Call price with continuous dividends:
    #     V = S e^{-q τ} N(d1) − K e^{-r τ} N(d2)
    V_arr = S_safe * np.exp(-p.q * tau) * norm.cdf(d1) - p.K * np.exp(-p.r * tau) * norm.cdf(d2)

    # Enforce V(0, τ) = 0 exactly
    V_arr[S_arr <= 1e-12] = 0.0

    if np.isscalar(S):
        # Convert single value back to Python float
        return float(V_arr.item())
    return V_arr


# --------------------------------------------------------------------------- #
# 3. Spatial discretization: construct MOL system
# --------------------------------------------------------------------------- #

def construct_system(
    M: int,
    p: BSParams
) -> Tuple[diags, float, np.ndarray, np.ndarray, np.ndarray]:
    """
    Discretize the spatial domain using central differences (Method of Lines).

    The spatial grid is [0, S_max] with M uniform intervals of size ΔS.
    We work on the interior nodes S_inner = {ΔS, 2ΔS, ..., S_max − ΔS}.

    Parameters
    ----------
    M : int
        Number of spatial intervals. There are M+1 grid points and M−1
        interior unknowns.
    p : BSParams
        Model parameters.

    Returns
    -------
    A : scipy.sparse.csc_matrix
        Tridiagonal matrix for the interior ODE system.
    dS : float
        Spatial step size ΔS.
    S_inner : ndarray of shape (M-1,)
        Interior grid points.
    alpha : ndarray
        Lower diagonal coefficients alpha_j (including the left boundary
        row, which is later separated into the boundary vector).
    gamma : ndarray
        Upper diagonal coefficients gamma_j (including the right boundary
        row, which is later separated into the boundary vector).
    """
    dS = p.S_max / M
    S_inner = np.linspace(dS, p.S_max - dS, M - 1)

    drift = p.r - p.q

    # Coefficients for central differences.
    # alpha_j: coefficient of V_{j-1}
    # beta_j : coefficient of V_j
    # gamma_j: coefficient of V_{j+1}
    ratio = S_inner / dS
    alpha = 0.5 * p.sigma**2 * ratio**2 - 0.5 * drift * ratio
    beta = -p.sigma**2 * ratio**2 - p.r
    gamma = 0.5 * p.sigma**2 * ratio**2 + 0.5 * drift * ratio

    # Build sparse tridiagonal matrix A (size (M-1) x (M-1)).
    # NOTE:
    #   • alpha[0] multiplies the left boundary V_0 (S=0) which is known and
    #     identically zero for a call, so it is not included in A.
    #   • gamma[-1] multiplies the right boundary V_M (S = S_max) and is
    #     handled via the b(τ) vector.
    A = diags(
        diagonals=[alpha[1:], beta, gamma[:-1]],
        offsets=[-1, 0, 1],
        format="csc",
    )

    return A, dS, S_inner, alpha, gamma


# --------------------------------------------------------------------------- #
# 4. Boundary vector b(τ)
# --------------------------------------------------------------------------- #

def get_boundary_vector(
    tau: float,
    p: BSParams,
    M: int,
    gamma_full: np.ndarray
) -> np.ndarray:
    """
    Construct the boundary condition vector b(τ) for the ODE system.

    For a European call:
        • Left boundary  (S = 0):        V(0, τ) = 0       ⇒ no contribution.
        • Right boundary (S = S_max):   V(S_max, τ) ≈
                                        S_max e^{-q τ} − K e^{-r τ}.

    Only the last interior node (adjacent to S_max) receives a contribution
    from the right boundary.

    Parameters
    ----------
    tau : float
        Time to expiry at which the boundary is evaluated.
    p : BSParams
        Model parameters.
    M : int
        Number of spatial intervals.
    gamma_full : ndarray
        The full gamma_j coefficients from construct_system (size M-1).

    Returns
    -------
    b : ndarray of shape (M-1,)
        Column vector for the boundary contribution in dV/dτ = A V + b(τ).
    """
    b = np.zeros(M - 1, dtype=float)

    # Right Dirichlet boundary at S_max:
    V_M = p.S_max * np.exp(-p.q * tau) - p.K * np.exp(-p.r * tau)

    # Contribution enters only the last interior node through gamma_{M-1} V_M.
    b[-1] += gamma_full[-1] * V_M

    return b


# --------------------------------------------------------------------------- #
# 5. Time-stepping solver
# --------------------------------------------------------------------------- #

def solve_ode_system(
    method: str,
    N: int,
    M: int,
    p: BSParams
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Solve the ODE system dV/dτ = A V + b(τ) using the chosen time-stepping method.

    Parameters
    ----------
    method : str
        Time-stepping scheme. Accepted values (case-insensitive):
            • "BackwardEuler", "backward_euler", "BE"
            • "BDF2"
    N : int
        Number of time steps in τ ∈ [0, T]. The time step is Δτ = T / N.
    M : int
        Number of spatial intervals. There are M−1 interior unknowns.
    p : BSParams
        Model parameters.

    Returns
    -------
    S_inner : ndarray of shape (M-1,)
        Interior spatial grid points in [0, S_max].
    V : ndarray of shape (M-1,)
        Numerical solution at τ = T on the interior grid.
    dS : float
        Spatial step size ΔS.

    Notes
    -----
    • The initial condition at τ = 0 is the option payoff:
          V(S, 0) = max(S − K, 0).
    • Backward Euler:
          (I − Δτ A) V^{n+1} = V^n + Δτ b^{n+1}
    • BDF2 (with Backward Euler start-up):
          (I − (2/3) Δτ A) V^{n+1}
              = (4/3) V^n − (1/3) V^{n−1} + (2/3) Δτ b^{n+1}
    """
    # --------------------- basic input validation --------------------- #
    if N <= 0:
        raise ValueError(f"N must be a positive integer, got N = {N}.")
    if M < 3:
        raise ValueError(
            f"M must be at least 3 so that there are interior points; got M = {M}."
        )

    method_key = method.strip().lower()
    if method_key in {"backwardeuler", "backward_euler", "be"}:
        scheme = "BackwardEuler"
    elif method_key == "bdf2":
        scheme = "BDF2"
    else:
        raise ValueError(
            f"Unknown method '{method}'. "
            "Use 'BackwardEuler' (or 'BE') or 'BDF2'."
        )

    # ---------------------- assemble MOL system ----------------------- #
    dtau = p.T / N
    A, dS, S_inner, alpha, gamma = construct_system(M, p)

    # Identity matrix on the interior space
    I = diags(
        diagonals=[np.ones(M - 1, dtype=float)],
        offsets=[0],
        format="csc",
    )

    # Initial condition at τ = 0: payoff max(S − K, 0)
    V = np.maximum(S_inner - p.K, 0.0)

    # --------------------------- Backward Euler ------------------------ #
    if scheme == "BackwardEuler":
        LHS = I - dtau * A

        for n in range(N):
            tau_next = (n + 1) * dtau
            b_next = get_boundary_vector(tau_next, p, M, gamma)
            V = spsolve(LHS, V + dtau * b_next)

    # ----------------------------- BDF2 -------------------------------- #
    elif scheme == "BDF2":
        # Step 1: Backward Euler start-up
        LHS_BE = I - dtau * A
        tau_1 = dtau
        b_1 = get_boundary_vector(tau_1, p, M, gamma)

        V_prev = V.copy()
        V = spsolve(LHS_BE, V_prev + dtau * b_1)

        # Steps 2..N: BDF2 with constant time step Δτ
        LHS_BDF = I - (2.0 / 3.0) * dtau * A

        for n in range(1, N):
            tau_next = (n + 1) * dtau
            b_next = get_boundary_vector(tau_next, p, M, gamma)

            RHS = (4.0 / 3.0) * V - (1.0 / 3.0) * V_prev + (2.0 / 3.0) * dtau * b_next
            V_prev = V.copy()
            V = spsolve(LHS_BDF, RHS)

    return S_inner, V, dS
