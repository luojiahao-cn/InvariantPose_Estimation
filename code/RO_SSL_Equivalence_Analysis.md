# RO 与 SSL 算法等价性推论及分析 (Mathematical Equivalence of RO and SSL)

本文档旨在推导并证明在单位球面 $S^2$ 流形上，**黎曼梯度下降 (Riemannian Optimization, RO)** 与 **连续球面线性化 (Successive Spherical Linearization, SSL)** 两种算法在数学上的等价性，并探讨其实际表现差异的原因。

---

## 1. 黎曼梯度下降 (RO) 的迭代推导

在单位球面 $S^2$ 上优化目标函数 $f(r)$，RO 的典型迭代步骤如下：

1.  **计算欧几里得梯度**：$\nabla f(r)$ (Ambient Gradient)。
2.  **投影到切空间**：得到黎曼梯度 $g_R$。
    $$g_R = \text{Proj}_{T_r S^2}(\nabla f(r)) = (I - rr^T) \nabla f(r) = \nabla f(r) - (r^T \nabla f(r))r$$
3.  **沿负梯度移动**（步长为 $\eta$）：
    $$r_{temp} = r - \eta g_R = r - \eta \left( \nabla f(r) - (r^T \nabla f(r))r \right)$$
    整理得：
    $$r_{temp} = -\eta \nabla f(r) + \left( 1 + \eta r^T \nabla f(r) \right) r$$
4.  **归一化 (Retraction)**：
    $$r_{next} = \frac{r_{temp}}{\|r_{temp}\|}$$

---

## 2. 连续球面线性化 (SSL) 的迭代推导

SSL 的设计初衷是通过近端线性化简化计算，其迭代步骤如下：

1.  **构造更新方向**：引入近端参数 $\beta_{prox}$。
    $$m = -\nabla f(r) + \beta_{prox} r$$
2.  **直接球面投影**：
    $$r_{next} = \frac{m}{\|m\|}$$

---

## 3. 等价性证明 (Equivalence Proof)

比较 RO 的更新方向 $r_{temp}$ 与 SSL 的更新方向 $m$：
*   **RO 方向**: $r_{temp} = \eta \left[ -\nabla f(r) + \left( \frac{1}{\eta} + r^T \nabla f(r) \right) r \right]$
*   **SSL 方向**: $m = -\nabla f(r) + \beta_{prox} r$

由于归一化操作对向量的整体正缩放因子（即 $\eta > 0$）不敏感，因此若要使两者产生的 $r_{next}$ 完全一致，只需满足：
$$\beta_{prox} = \frac{1}{\eta} + r^T \nabla f(r)$$

### 结论
**SSL 算法可以被视为一种特殊的黎曼梯度下降法**。在这种视角下，SSL 隐含地采用了一种随梯度动态变化的步长 $\eta = 1 / (\beta_{prox} - r^T \nabla f(r))$。

---

## 4. 为什么实验中 SSL 表现更好？

尽管数学上可以等价，但在实际磁主轴估计实验中，SSL ($0.037$) 优于 RO ($0.083$)，原因在于两者的**步长策略偏好**不同：

1.  **保守 vs. 贪婪**：
    *   **RO** 通常采用 **BB 步长** 或 **线搜索**，追求每步下降最大化。在磁定位这种复杂的非凸地形中，大步长容易跨过全局最优解，收敛到局部极小值。
    *   **SSL** 通常设定一个固定的、较大的 $\beta_{prox}$（如 $100$），这相当于强制采用了一个**非常小且稳健的步长**。这种“慢步滑行”的策略表现出了更强的全局寻优鲁棒性。
2.  **物理结构的契合度**：
    *   SSL 的形式 $m = -\nabla f(r) + \beta r$ 极像求解特征向量的**幂迭代 (Power Method)**。磁主轴估计的核心本质上就是寻找某个算子的主特征方向。SSL 的迭代动态在物理本质上比通用的流形投影更契合这一结构。
3.  **计算极简性**：
    *   SSL 避开了显式的切空间投影矩阵运算 $(I - rr^T)$，在保留几何约束的同时，通过更简单的代数形式实现了稳健收敛。

---

## 5. 命名意义建议

虽然 SSL 是 RO 的一个特例，但建议保留 **"Successive Spherical Linearization" (SSL)** 这一名称。它准确描述了算法通过“连续线性化 + 球面重投影”这一物理逻辑，且其独特的参数化方式（Proximal-like）在解决特定非凸物理问题时具有显著的实用价值。
