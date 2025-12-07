# Repository Guidelines

## Project Structure & Module Organization
- Core solver: `BS_Solver.py` holds the Method of Lines implementation, parameter container, discretization, and solvers (Backward Euler, BDF2, exact solution).
- Notebooks: `01_Verification_Analysis.ipynb` (convergence, Greeks) and `02_RealWorld_CaseStudy.ipynb` (market data validation). `draft/` keeps scratch notebooks; `backup/` contains older copies.
- Assets: `images/` stores generated plots; keep rendered outputs small and replace rather than append.
- Dependencies: pinned in `requirements.txt`; notebooks expect Python 3.8+ with NumPy/SciPy/Matplotlib.

## Build, Test, and Development Commands
- Create env and install deps:
  - `python -m venv .venv && source .venv/bin/activate`
  - `pip install -r requirements.txt`
- Run notebooks locally: `jupyter notebook 01_Verification_Analysis.ipynb` (or the case-study notebook). Restart kernels when changing solver code.
- Quick solver check in a REPL:
  ```bash
  python - <<'PY'
  from BS_Solver import BSParams, solve_ode_system
  p = BSParams(T=1.0, K=100, r=0.03, sigma=0.2)
  S, V, dS = solve_ode_system("BDF2", N=200, M=400, p=p)
  print(float(V[len(V)//2]))
  PY
  ```

## Coding Style & Naming Conventions
- Follow PEP 8: 4-space indentation, snake_case for functions/variables, PascalCase for classes. Use descriptive financial names (`S_max`, `sigma`, `tau`).
- Keep functions side-effect free; prefer NumPy vectorization over loops. Document numeric choices (grid sizes, boundary handling) in docstrings or short comments.
- Notebooks: title cells with numbered headings; keep outputs collapsed before committing.

## Testing Guidelines
- No formal automated suite yet. Validate changes by rerunning `01_Verification_Analysis.ipynb` (EOC, Greeks) and `02_RealWorld_CaseStudy.ipynb` (market comparison) and confirm plots/metrics remain reasonable.
- For new features, add lightweight checks (e.g., compare `solve_ode_system` against `exact_solution` for coarse grids) and prefer deterministic parameters so results are reproducible.
- If adding automated tests, place them under `tests/` using `pytest`; avoid network calls and large artifacts.

## Commit & Pull Request Guidelines
- Commit messages: short imperative summaries (e.g., `Refine BDF2 startup step`, `Add convergence notebook run`). Group related file changes together.
- PRs: include a brief problem statement, what changed, and how it was validated (notebooks run, plots regenerated). Attach before/after figures if visuals changed and note any parameter defaults updated.
- Link to relevant coursework prompt or issue ID when applicable, and mention any dependencies added/removed.***
