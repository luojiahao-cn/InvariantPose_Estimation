# SSL 算法实验参数与设置总结

本文档总结了磁定位项目中与 SSL (Successive Semidefinite Lifting) 算法相关的实验参数、默认配置及调用方式。

## 1. 算法概述

SSL (Successive Semidefinite Lifting) 算法用于在单位球 $S^2$ 上进行主轴估计 (Principal-axis estimation)。

**优化目标：**
$$ \min_{r \in S^2} f(r) = \sum_k \alpha_k ( r^T A^k r - s_k )^2  +  \beta \| \bar{B}^T r - \bar{b}^T u \|^2 $$

## 2. 参数设置来源

参数主要配置在 `utils/run_single_experiment.m` 中，部分默认参数定义在 `Functions/estimate_principal_axis_SSL.m` 的函数头中。

### 2.1 核心计算参数 (`opts_r` / `prob_r`)

这些参数在 `run_single_experiment.m` 中构造，并传递给 `compute_principal_axis_prob` 和 `estimate_principal_axis_SSL`。

| 参数名 | 值 / 设置 | 描述 |
| :--- | :--- | :--- |
| `opts_r.u` | `[0; 0; 1]` | 载体坐标系下的参考主轴方向向量。 |
| `opts_r.r0` | `MatrixExp3(VecToso3(theta_init)) * opts_r.u` | 主轴方向的初始猜测值，基于初始旋转扰动 `theta_init` 计算。 |
| `opts_r.alpha` | `[1, 0]` | 目标函数中高次项的权重向量。这里仅考虑 $k=1$ 的项，即 $(r^T A r - s_1)^2$。 |
| `opts_r.beta` | `1` | 目标函数中线性项 (Linear Term) 的权重，对应 $\| \bar{B}^T r - \bar{b}^T u \|^2$。 |

### 2.2 求解器配置参数 (`opts`)

如果调用时未指定，`estimate_principal_axis_SSL.m` 会使用以下默认值：

| 参数名 | 默认值 | 描述 |
| :--- | :--- | :--- |
| `maxIter` | `500` | 最大迭代次数。 |
| `tol` | `1e-8` | 收敛容差 (Tolerance)。 |
| `proxBeta` | `1e2` | 近端项 (Proximal term) 的权重系数 `beta_prox`。 |
| `verbose` | `0` | 是否打印详细日志 (0: 关闭)。 |

## 3. 调用流程

在 `utils/run_single_experiment.m` 中，SSL 算法通常作为 `estimate_pose_ours` (所提算法) 的后续步骤被调用，用于精细化主轴估计。

```matlab
% 1. 构造配置结构体
opts_r = struct('u', [0;0;1]);
opts_r.r0 = MatrixExp3(VecToso3(theta_init)) * opts_r.u; % 初始猜测
opts_r.alpha = [1, 0];
opts_r.beta = 1;

% 2. 计算问题参数 (prob_r)
% 输入来自 estimate_pose_ours 的统计量 (b_bar, B_bar, A_p, X_opt)
[prob_r, opts_r] = compute_principal_axis_prob(stats_ours.b_bar, stats_ours.B_bar, stats_ours.A_p, stats_ours.X_opt, opts_r);

% 3. 调用 SSL 估计算法
[r_direct_SSL, info_r_SSL] = estimate_principal_axis_SSL(prob_r, opts_r);
```

## 4. 相关文件

- **核心实现**: `code/Functions/estimate_principal_axis_SSL.m`
- **实验调用**: `code/utils/run_single_experiment.m`
- **问题构造**: `code/Functions/compute_principal_axis_prob.m`
