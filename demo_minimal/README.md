# Magnetic Sensor Array Pose Estimation - Minimal Demo

A minimal reproducible example for magnetic sensor array pose estimation using rotation invariants.

## Overview

This demo implements four pose estimation algorithms for a magnetic sensor array:
- **Ours**: The proposed two-stage method (eigenvalue-based position + iterative rotation estimation)
- **LM**: Levenberg-Marquardt optimization on magnetic field residuals
- **ELM**: Optimization on magnetic field difference residuals
- **Fischer**: Grid search with eigenvalue invariants

## Directory Structure

```
demo_minimal/
├── README.md
├── run_batch_demo.m          # Main entry script
├── get_params.m             # Parameter configuration
├── dipole_b_and_gradb.m     # Dipole magnetic field calculation
├── calcFieldAndGradient.m    # Field and gradient summation
├── lc_grad_tensor_estimator.m  # Local gradient tensor estimation
├── generate_test_points.m   # Random test point generation
├── generate_magnetic_data.m # Magnetic data generation
├── estimateR.m              # Rotation estimation (eigen decomposition)
├── estimateR_iter.m         # Rotation estimation (iterative PPI)
├── estimate_pose_ours.m     # Proposed algorithm
├── estimate_pose_lm.m       # LM algorithm
├── estimate_pose_elm.m      # ELM algorithm
├── estimate_pose_fischer.m   # Fischer algorithm
└── helper/
    ├── MatrixExp3.m         # SO(3) exponential map
    └── VecToso3.m           # Vector to skew-symmetric matrix
```

## Dependencies

- MATLAB R2018b or later
- Optimization Toolbox (for `lsqnonlin`, `fmincon`, `fminunc`)

No Robotics Toolbox required - `MatrixExp3` and `VecToso3` are implemented locally.

## Usage

1. Add the `demo_minimal` directory to your MATLAB path:

```matlab
addpath demo_minimal
cd demo_minimal
```

2. Run the batch demo:

```matlab
run_batch_demo
```

## Algorithm Description

### Stage 1: Position Estimation
Uses eigenvalue matching between the measured gradient tensor and the model gradient tensor. The position is found by minimizing the eigenvalue differences through grid search followed by gradient refinement.

### Stage 2: Rotation Estimation
1. Initialize rotation using eigen decomposition
2. Refine using iterative Procrustes-like optimization on the manifold

## Test Configuration

- **Magnets**: 3 magnets with specific positions and orientations
- **Sensors**: 7 sensors in a cross configuration
- **Workspace**: Spherical region with radius 0.15m
- **Test points**: Randomly generated within the workspace (upper hemisphere)

## Output

The script prints per-test-point errors and summary statistics:
- Position error: Euclidean distance in millimeters
- Rotation error: Frobenius norm of rotation difference in radians

## Citation

If you use this code, please cite our paper.
