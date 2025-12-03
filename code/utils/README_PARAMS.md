# 实验参数配置说明

## 概述

`exp1_main.m` 已经重构为模块化结构，所有实验参数统一在 `get_experiment_params.m` 函数中配置。

## 参数分类

参数结构体 `params` 包含以下子结构：

### 1. `params.magnet` - 磁铁参数
- `m_pos`: 磁铁位置矩阵 (3×K)，单位：m
- `m_hat`: 磁化方向单位向量 (3×K)，已归一化
- `m_norm`: 磁矩幅值向量 (1×K)，单位：A·m²

### 2. `params.sensor` - 传感器参数
- `d_list`: 传感器位移矩阵 (3×N)，传感器在参考坐标系中的偏移，单位：m

### 3. `params.ground_truth` - 真实姿态参数
- `theta_true`: 真实旋转向量 (3×1)，单位：rad
- `p_true`: 传感器阵列参考点真实位置 (3×1)，单位：m

### 4. `params.uncertainty` - 不确定性参数
- `p_uncertainty`: 位置不确定性
- `r_uncertainty`: 旋转不确定性

### 5. `params.optimization` - 优化算法参数
- `options`: lsqnonlin 优化选项
- `mu`: 所提算法的 mu 参数
- `beta`: 所提算法的 beta 参数

### 6. `params.workspace` - 工作空间约束参数
- `center`: 工作空间中心点 (3×1)
- `radius`: 工作空间半径

### 7. `params.experiment` - 实验设置参数
- `num_experiments`: 实验次数
- `random_seed`: 随机种子

## 如何修改参数

### 方法1：直接编辑配置文件
编辑 `code/utils/get_experiment_params.m` 文件，修改对应参数值。

### 方法2：在主文件中覆盖参数
在 `exp1_main.m` 中调用 `get_experiment_params()` 后，可以覆盖特定参数：

```matlab
params = get_experiment_params();
% 覆盖特定参数
params.magnet.m_norm = [400, 400];
params.optimization.mu = 2e3;
params.experiment.num_experiments = 20;
```

## 代码结构

```
exp1_main.m                    # 主实验文件（简化后）
├── get_experiment_params()    # 参数配置函数
├── generate_magnetic_data()   # 数据生成函数
└── run_single_experiment()    # 单次实验执行函数
```

## 优势

1. **参数集中管理**：所有参数在一个地方配置，易于查找和修改
2. **分类清晰**：参数按功能分类，结构清晰
3. **代码复用**：数据生成和实验执行封装为函数，可在其他实验中使用
4. **易于维护**：模块化设计，便于后续扩展和修改

