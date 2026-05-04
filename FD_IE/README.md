# FD_IE

Finite-difference Cahn-Hilliard test with explicit Euler and fully implicit Euler time stepping.

## Run

In MATLAB:

```matlab
restoredefaultpath
cd('/path/to/FD_IE')
addpath(pwd)
which ch_fd_solve -all
run_fd_implicit_reviewer_tests
```

`which ch_fd_solve -all` should print one file, from this folder.

## Expected pass message

```text
Reviewer implicit-FD tests passed.
```

## Output files

The script writes:

```text
output/implicit_fd_reviewer/
  implicit_fd_reviewer_summary.csv
  implicit_fd_reviewer_results.mat
  N64_full_gate_fields.png
  N64_full_gate_diagnostics.png
  N128_early_fields.png
  N128_early_diagnostics.png
```

The field PNGs show:

```text
Initial condition | Explicit reference | Implicit Euler FD
```

The explicit reference uses a 100x smaller time step at the same physical time. The implicit and reference panels should look nearly identical.

The diagnostics PNGs show normalized energy and mass error. Newton convergence is checked by the script and recorded in the CSV/MAT output, but it is not plotted.

## Expected numbers

Typical MATLAB R2024a output:

```text
N64_full_gate:
  explicit large-step failure step: 5
  implicit mass error: ~1e-17
  energy ratio: ~0.535
  RMSE(implicit, reference): ~4.14e-4
  correlation: ~0.999999

N128_early:
  explicit large-step failure step: 5
  implicit mass error: ~1e-17
  energy ratio: ~0.563
  RMSE(implicit, reference): ~1.80e-3
  correlation: ~0.999952
```

The script also checks Newton convergence. The maximum Newton residual should be below `1e-8`.

MATLAB may print a `CodeCache is full` warning on some machines. If the final pass message appears, the test passed.
