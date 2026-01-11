# Magnetic Pose Estimation Workspace Instructions

## Architecture Overview
This project specializes in magnetic localization using magnetometer sensor arrays and rotation invariant methods.
- **Entry Points**: `exp_random_point_test.m`, `exp_point_test.m`, and `exp2_initial_robustness.m`.
- **Batch Processing**: `utils/run_batch_experiments.m` orchestrates multiple test points and trials.
- **Trial Execution**: `utils/run_single_experiment.m` manages individual solver calls and error calculations.
- **Core Algorithms**: Located in `Functions/`, following the naming pattern `estimate_pose_[method].m`.

## Developer Workflows
- **Running Experiments**: Run the top-level scripts (e.g., `exp_random_point_test.m`). These scripts setup paths via `addpath` and configure parameters.
- **Data Initialization**: Parameters are centrally managed via `utils/get_experiment_params.m`. Use this struct to propagate settings.

## Coding Conventions
- **Variable Naming**:
    - `p` (3x1): Position vector.
    - `R` (3x3): Rotation matrix.
    - `theta` (3x1): Exponential coordinates for rotation (axis-angle).
- **Sensor Configurations**: `d_list` is a `3xN` matrix representing the relative positions of $N$ sensors in the array.
- **Error Metrics**:
    - Position error: `norm(p_est - p_true)` (meters).
    - Rotation error: Derived from the axis-angle representation of the mismatch matrix.

## Key Files & Directories
- `Functions/`: Contains the implementation of various pose estimation algorithms (LM, ELM, ours, Fischer).
- `utils/`: Core utilities for data generation, batch execution, and parameter management.
- `tools/`: Visualization and post-processing tools like `plot_batch_results.m`.
- `./exp/`: Contains experimental datasets and optimized parameters in `.mat` and `.json` formats.
