# Codes for Paper *Local Second-Order Limit Dynamics of the Alternating Direction Method of Multipliers for Semidefinite Programming*

All codes are ready to run. No need to install anything except MATLAB.

## 1. $\sigma$ updating rules

The $\sigma$ updating rules for the first $20,000$ iterations can be found in `admm-single-light/update_sig_1.m`.

## 2. File execution order

### Toy examples

Directly run `toy-examples/toy1.m`, `toy-examples/toy2.m`, and `toy-examples/toy3.m` separately. All results are already contained in the `toy-examples/figs` folder.

### `Mittelmann` dataset

For any instances out of the selected $25$ `Mittelmann` problems (single block, matrix size no greater than $3,000$):

1. Run `sig-wave/generate_start_points.m` to generate starting points for both the three-step ADMM formal test and $\sigma$ tuning test. 
2. After (1), run `acceleration/main_baseline.m` to generate the three-step ADMM long trajectories $|| \Delta Z^{(k)} ||_F$, $ r_{\max}^{(k)}$, and $ \angle (\Delta Z^{(k)}, \Delta Z^{(k+1)})$. The maximum iteration number is set to $1,000,000$.
3. After (1), we can also directly run `sig-wave/study_sig_wave.m` to generate $\sigma$ tuning trajectories $r_p^{(k)}$, $r_{d}^{(k)}$, $|| \Delta X^{(k)} ||_F$, and $|| \Delta S^{(k)} ||_F$.

Please change the `type` to the data instance you want to test in each folder. After running all $25$ instances, you can then use `sig-wave/plot_sig_wave.m` and `acceleration/plot_all.m` to generate figures in the paper. The results are also already contained in the `acceleration/figs` folder.
