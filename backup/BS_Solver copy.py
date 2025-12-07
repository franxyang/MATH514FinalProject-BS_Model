"""
-------------------------------------------------------------------------------
File Name:      BS_Solver.py
Author:         Yifan Yang
Date Revised:   December 2, 2025
Course:         MATH 514 - Numerical Analysis (Fall 2025)
Institution:    University of Wisconsin-Madison

Description:
    This module implements the core numerical engine for pricing European Call
    Options using the Method of Lines (MOL). It discretizes the generalized 
    Black-Scholes PDE (with dividend yield) into a stiff system of ODEs.

    It includes:
    1.  BSParams: A data class for model parameters.
    2.  construct_system: Spatial discretization (Central Differences).
    3.  solve_ode_system: Time-stepping solvers (Backward Euler & BDF2).
    4.  exact_solution: Closed-form BS solution for validation.

Usage:
    from BS_Solver import BSParams, solve_ode_system
-------------------------------------------------------------------------------
"""

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.stats import norm

class BSParams:
    """
    Data class to hold Black-Scholes parameters.
    """
    def __init__(self, S_max=500.0, K=100.0, T=1.0, r=0.03, sigma=0.30, q=0.0):
        self.S_max = S_max  # Maximum asset price (domain boundary)
        self.K = K          # Strike price
        self.T = T          # Time to expiry (years)
        self.r = r          # Risk-free interest rate
        self.sigma = sigma  # Volatility
        self.q = q          # Dividend yield

def exact_solution(S, tau, p):
    """
    Generalized Black-Scholes Call Option Closed-form Solution (with dividends).
    
    Args:
        S (array or float): Underlying asset price.
        tau (float): Time to expiry.
        p (BSParams): Parameter object.
    """
    # Avoid log(0)
    S_safe = np.where(S > 1e-8, S, 1e-8)
    
    d1 = (np.log(S_safe / p.K) + (p.r - p.q + 0.5 * p.sigma**2) * tau) / (p.sigma * np.sqrt(tau))
    d2 = d1 - p.sigma * np.sqrt(tau)
    
    # Call Price formula with dividends: S * e^(-q*tau) * N(d1) - K * e^(-r*tau) * N(d2)
    V = S * np.exp(-p.q * tau) * norm.cdf(d1) - p.K * np.exp(-p.r * tau) * norm.cdf(d2)
    
    # Correction for S=0
    if isinstance(S, np.ndarray):
        V[S <= 1e-8] = 0.0
    else:
        if S <= 1e-8: V = 0.0
    return V

def construct_system(M, p):
    """
    Discretize the spatial domain using Central Differences (Method of Lines).
    Returns Matrix A, step size dS, grid points S_inner, and diag vectors.
    """
    dS = p.S_max / M
    S_inner = np.linspace(dS, p.S_max - dS, M - 1)
    
    drift = p.r - p.q
    
    # Coefficients for Central Difference
    # alpha (Lower): Coeff of V_{i-1}
    alpha = 0.5 * p.sigma**2 * (S_inner/dS)**2 - 0.5 * drift * (S_inner/dS)
    # beta (Main):   Coeff of V_{i}
    beta  = -p.sigma**2 * (S_inner/dS)**2 - p.r
    # gamma (Upper): Coeff of V_{i+1}
    gamma = 0.5 * p.sigma**2 * (S_inner/dS)**2 + 0.5 * drift * (S_inner/dS)
    
    # Construct Sparse Matrix
    A = diags([alpha[1:], beta, gamma[:-1]], [-1, 0, 1], format='csc')
    
    return A, dS, S_inner, alpha, gamma

def get_boundary_vector(tau, p, M, gamma_full):
    """
    Construct the boundary condition vector b(tau).
    """
    b = np.zeros(M - 1)
    # Right Boundary: V(S_max) = S_max * e^(-q*tau) - K * e^(-r*tau)
    V_M = p.S_max * np.exp(-p.q * tau) - p.K * np.exp(-p.r * tau)
    b[-1] += gamma_full[-1] * V_M
    return b

def solve_ode_system(method, N, M, p):
    """
    Solve the ODE system dV/dtau = AV + b using specified time-stepping method.
    
    Args:
        method (str): 'BackwardEuler' or 'BDF2'.
        N (int): Number of time steps.
        M (int): Number of spatial steps.
        p (BSParams): Parameters.
    """
    dtau = p.T / N
    A, dS, S_inner, alpha, gamma = construct_system(M, p)
    I = diags([1] * (M - 1), 0, format='csc')
    
    # Initial Condition: Payoff at tau=0
    V = np.maximum(S_inner - p.K, 0)
    
    if method == 'BackwardEuler':
        LHS = I - dtau * A
        for n in range(N):
            tau_next = (n + 1) * dtau
            b_next = get_boundary_vector(tau_next, p, M, gamma) # corrected arg name
            V = spsolve(LHS, V + dtau * b_next)
            
    elif method == 'BDF2':
        # Step 1: Backward Euler (Startup)
        LHS_BE = I - dtau * A
        tau_1 = 1 * dtau
        b_1 = get_boundary_vector(tau_1, p, M, gamma)
        V_prev = V.copy()
        V = spsolve(LHS_BE, V_prev + dtau * b_1)
        
        # Step 2 to N: BDF2
        LHS_BDF = I - (2.0/3.0) * dtau * A
        for n in range(1, N):
            tau_next = (n + 1) * dtau
            b_next = get_boundary_vector(tau_next, p, M, gamma)
            RHS = (4.0/3.0)*V - (1.0/3.0)*V_prev + (2.0/3.0)*dtau*b_next
            V_prev = V.copy()
            V = spsolve(LHS_BDF, RHS)
            
    return S_inner, V, dS