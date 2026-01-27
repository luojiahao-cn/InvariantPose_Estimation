好，这一段本身已经是**论文级别、逻辑非常干净的 LaTeX**了。你现在的目标是：

> 在本地（比如 Markdown + MathJax / Obsidian / Typora / VS Code）先把优化问题“跑通脑子”，方便你后面推 MM 推导和代码实现。

下面我做三件事：

1. **完整转写为 Markdown（保留所有数学公式）**
2. **在关键位置加极少量结构性标题，方便你定位与调试**
3. **不引入任何新的解释，不改变数学含义**

---

## Estimation for Principal Axis (Stage #2)

Under the compact-array approximation, the estimated local gradient tensor $\hat{\mathbf{X}}$ is related to the world-frame gradient tensor $\mathbf{A}(\hat{\mathbf{p}})$ via the orthogonal similarity transform
$$
\hat{\mathbf{X}} \approx \mathbf{R}^\top \mathbf{A}(\hat{\mathbf{p}})\mathbf{R},
\tag{1}
$$
where $\hat{\mathbf{p}}$ is obtained from **Stage #1**.

Instead of aligning full tensors to recover $\mathbf{R}$, we construct **scalar invariants** by probing $\hat{\mathbf{X}}$ along the known local axis $\mathbf{u}$.
For any integer $k \ge 1$, define
$$
s_k \triangleq \mathbf{u}^\top \hat{\mathbf{X}}^k \mathbf{u}.
$$

Using the definition of $\mathbf{r}$ and the similarity relation above, we have
$$
s_k \approx
\mathbf{u}^\top \mathbf{R}^\top \mathbf{A}(\hat{\mathbf{p}})^k \mathbf{R}\mathbf{u}
===================================================================================

\mathbf{r}^\top \mathbf{A}(\hat{\mathbf{p}})^k \mathbf{r},
\quad k \ge 1.
\tag{2}
$$

Therefore, the set ${s_k}$ yields scalar constraints on $\mathbf{r}$ in the world frame, independent of the full rotation $\mathbf{R}$.

---

## Invariant Consistency Cost

In practice, the relations in (2) hold only approximately due to measurement noise and modeling errors.
We therefore enforce these invariant relations in a **least-squares sense** by penalizing the corresponding residuals.

For a given index set $\mathbb{K} \subset \mathbb{N}$, we consider the invariant consistency cost
$$
\sum_{k \in \mathbb{K}} \alpha_k
\left(
\mathbf{r}^\top \mathbf{A}(\hat{\mathbf{p}})^k \mathbf{r} - s_k
\right)^2,
\tag{3}
$$
where $\alpha_k \ge 0$ are weighting coefficients.

---

## First-Order Magnetic Field Consistency

In addition to the gradient-based invariants, we incorporate **first-order magnetic field consistency** to further constrain the axis direction $\mathbf{r}$.

Under the same compact-array approximation, the projected field measurements satisfy
$$
\mathbf{R}^\top \mathbf{B}(\hat{\mathbf{p}}) \approx \bar{\mathbf{B}},
$$
where $\mathbf{B}(\hat{\mathbf{p}})$ denotes the predicted field matrix evaluated at $\hat{\mathbf{p}}$ and $\bar{\mathbf{B}}$ is the corresponding measurement matrix.

Left- and right-multiplying by $\mathbf{u}$ yields
$$
\mathbf{B}(\hat{\mathbf{p}})^\top \mathbf{r}
\approx
\bar{\mathbf{B}}^\top \mathbf{u}.
$$

This relation is likewise enforced in a least-squares sense, leading to the residual
$$
\left|
\mathbf{B}(\hat{\mathbf{p}})^\top \mathbf{r}
--------------------------------------------

\bar{\mathbf{B}}^\top \mathbf{u}
\right|_2^2.
\tag{4}
$$

---

## Joint Optimization Problem

Combining (3) and (4), the estimation of the principal axis $\mathbf{r}$ is formulated as the constrained optimization problem
$$
\begin{aligned}
\min_{\mathbf{r} \in \mathbb{R}^3} \quad
&
\sum_{k \in \mathbb{K}} \alpha_k
\left(
\mathbf{r}^\top \mathbf{A}(\hat{\mathbf{p}})^k \mathbf{r} - s_k
\right)^2
+
\mu
\left|
\mathbf{B}(\hat{\mathbf{p}})^\top \mathbf{r}
--------------------------------------------

\bar{\mathbf{B}}^\top \mathbf{u}
\right|_2^2,
\
\text{s.t.} \quad
&
|\mathbf{r}|_2 = 1,
\end{aligned}
\tag{5}
$$
where $\mu \ge 0$ balances the contributions of the two terms.

---

## Quartic Formulation on the Unit Sphere

Expanding the objective in (5) and discarding constant terms independent of $\mathbf{r}$, we obtain the equivalent formulation
$$
\begin{aligned}
\min_{|\mathbf{r}|*2 = 1} \quad
&
\sum*{k \in \mathbb{K}} \alpha_k
\left(
\mathbf{r}^\top \mathbf{A}^k \mathbf{r}
\right)^2
\
&\quad
+
\mathbf{r}^\top
\left(
\mu \mathbf{G}
--------------

2 \sum_{k \in \mathbb{K}} \alpha_k s_k \mathbf{A}^k
\right)
\mathbf{r}
----------

2\mu, \mathbf{g}^\top \mathbf{r},
\end{aligned}
\tag{6}
$$
where
$$
\mathbf{g} \triangleq \mathbf{B}(\hat{\mathbf{p}})\bar{\mathbf{B}}^\top \mathbf{u},
$$
and we denote $\mathbf{A}(\hat{\mathbf{p}})$ simply as $\mathbf{A}$ for brevity.

Problem (6) is a **quartic minimization over the unit sphere**, with all nonconvexity arising from the squared quadratic forms induced by the gradient invariants.

Rather than treating it as a general fourth-order polynomial optimization, we exploit this composite structure and adopt a **structured majorization–minimization (MM) scheme**.

---

### 下一步我建议你做什么

如果你愿意，下一步我可以直接帮你：

1. **把 (6) 写成标准 MM 可主化形式**

   * 明确 current iterate $\mathbf{r}^{(t)}$
   * 给出每一项的 tight quadratic majorizer
2. **给你一个“可直接抄进论文”的 MM 子问题**
3. **顺手对齐你之前提到的 trust-region / Cayley / 特征值方法**

你只需要告诉我：
👉 **你是更想先推数学，还是先做一个 MATLAB / Python 的可跑版本？**
