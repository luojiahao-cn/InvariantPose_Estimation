# 批量测试使用说明

## 概述

`exp1_batch_main.m` 用于测试算法在多个不同 `p_true` 和 `theta_true` 组合下的恢复能力。这对于评估算法的鲁棒性和收敛域非常有用。

## 使用方法

### 基本使用

直接运行 `exp1_batch_main.m`：

```matlab
exp1_batch_main
```

### 自定义参数

在 `exp1_batch_main.m` 中可以修改以下参数：

```matlab
num_test_points = 100;        % 测试点数量
num_trials_per_point = 5;     % 每个测试点的实验次数
```

### 测试点生成方法

`generate_test_points()` 函数支持三种生成方法：

1. **'random'** (默认): 在工作空间内随机生成位置和旋转
   ```matlab
   test_points = generate_test_points(params, 100, 'random');
   ```

2. **'grid'**: 在工作空间内网格采样（适用于立方数）
   ```matlab
   test_points = generate_test_points(params, 64, 'grid');  % 4×4×4网格
   ```

3. **'uniform_sphere'**: 在球面上均匀采样位置，随机旋转
   ```matlab
   test_points = generate_test_points(params, 100, 'uniform_sphere');
   ```

## 输出结果

### 控制台输出

程序会输出：
- 每个测试点的进度信息
- 每个测试点的统计摘要
- 总体统计结果（均值、标准差、最大值）
- 最佳/最差测试点信息

### 保存的文件

结果会保存到 `results/exp1_batch_results.mat`，包含：
- `batch_results`: 批量实验结果结构体数组
- `test_points`: 测试点信息
- `params`: 实验参数

### 结果结构

`batch_results` 是一个结构体数组，每个元素对应一个测试点：

```matlab
batch_results(i).test_point.p_true      % 测试点的真实位置
batch_results(i).test_point.theta_true   % 测试点的真实旋转向量
batch_results(i).results                 % 该测试点的所有实验结果数组
batch_results(i).summary                 % 该测试点的统计摘要
```

`summary` 包含各算法的统计信息：
- `lm.pos_mean`, `lm.pos_std`, `lm.pos_max`
- `lm.rot_mean`, `lm.rot_std`, `lm.rot_max`
- `elm.*`, `ours.*`, `Rlm.*` (类似结构)

## 示例：分析保存的结果

```matlab
% 加载结果
load('results/exp1_batch_results.mat');

% 分析结果
analyze_batch_results(batch_results);

% 访问特定测试点的结果
point_idx = 1;
fprintf('测试点 %d 的位置误差均值:\n', point_idx);
fprintf('  LM: %.6f m\n', batch_results(point_idx).summary.lm.pos_mean);
fprintf('  Ours: %.6f m\n', batch_results(point_idx).summary.ours.pos_mean);
```

## 性能考虑

- **测试点数量**: 100个测试点 × 5次实验 = 500次优化，可能需要较长时间
- **并行化**: 可以考虑使用 `parfor` 并行执行（需要修改代码）
- **内存**: 大量结果会占用内存，建议定期保存中间结果

## 自定义分析

可以基于保存的结果进行自定义分析：

```matlab
load('results/exp1_batch_results.mat');

% 提取所有位置误差
all_ours_pos_errors = [];
for i = 1:length(batch_results)
    for j = 1:length(batch_results(i).results)
        all_ours_pos_errors(end+1) = batch_results(i).results(j).ours_pos_error;
    end
end

% 绘制误差分布
histogram(all_ours_pos_errors);
xlabel('位置误差 (m)');
ylabel('频次');
title('算法位置误差分布');
```

