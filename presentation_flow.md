# Topic 8 Presentation Flow

## 核心主线

整个 presentation 的主线是：

> 我们要解大规模线性系统 \(Ax=b\)。直接法在大规模稀疏问题中可能太贵，所以介绍三种基础的 stationary iterative methods：Jacobi、Gauss-Seidel 和 SOR。先讲统一框架，再讲三个方法，最后用一个 1D Poisson 例子比较它们的收敛速度。

这三个方法的共同点是都可以写成：

\[
x^{(k+1)} = Gx^{(k)} + c.
\]

它们的区别在于如何分裂矩阵 \(A\)，以及每次迭代使用旧值还是新值。

---

## Slide-by-slide 思路

### 1. Title

介绍题目：

**Stationary Iterative Methods for Solving Linear Systems**

这页不用讲太多。可以直接说：

> Today we will talk about three classical stationary iterative methods: Jacobi, Gauss-Seidel and SOR.

---

### 2. Motivation

这一页回答“为什么需要这些方法”。

很多数值方法最终都会变成解线性系统：

\[
Ax=b.
\]

如果矩阵很大、很稀疏，直接法可能计算量和存储量都很高。迭代法的想法是：从一个初始猜测 \(x^{(0)}\) 开始，不断更新，直到 residual 足够小。

可以讲：

> Direct methods are reliable, but for large sparse systems they may be expensive. Iterative methods improve an initial guess until the residual becomes small.

---

### 3. Where Do These Systems Come From?

这一页给出线性系统的来源。

例子是 1D Poisson equation：

\[
-u''(x)=f(x), \qquad 0<x<1, \qquad u(0)=u(1)=0.
\]

用 centred finite difference 离散后得到：

\[
\frac{-u_{i-1}+2u_i-u_{i+1}}{h^2}=f_i.
\]

这会产生一个三对角稀疏矩阵 \(A\)。

这一页的作用是把 topic 和 Computational Differential Equations 联系起来。

---

### 4. General Stationary Iteration Framework

这一页是理论核心。

把矩阵分裂成：

\[
A=M-N.
\]

那么：

\[
Ax=b
\]

等价于：

\[
Mx=Nx+b.
\]

因此得到迭代格式：

\[
x^{(k+1)}=M^{-1}Nx^{(k)}+M^{-1}b.
\]

写成统一形式：

\[
x^{(k+1)}=Gx^{(k)}+c.
\]

其中：

\[
G=M^{-1}N.
\]

重点讲：

> Jacobi, Gauss-Seidel and SOR are all obtained by choosing different splittings of \(A\).

---

### 5. Convergence Criterion

这一页讲 stationary iteration 什么时候收敛。

设精确解是 \(x^\ast\)，误差是：

\[
e^{(k)}=x^{(k)}-x^\ast.
\]

因为：

\[
x^\ast=Gx^\ast+c,
\]

所以：

\[
e^{(k+1)}=Ge^{(k)}.
\]

于是：

\[
e^{(k)}=G^ke^{(0)}.
\]

核心结论：

\[
\rho(G)<1.
\]

这里 \(\rho(G)\) 是 spectral radius，也就是 \(G\) 的特征值绝对值的最大值。

可以讲：

> If the spectral radius is less than one, powers of \(G\) decay and the error goes to zero.

---

### 6. Matrix Splitting Notation

这一页引入后面三个方法都要用的记号：

\[
A=D+L+U.
\]

其中：

- \(D\)：diagonal part
- \(L\)：strictly lower triangular part
- \(U\)：strictly upper triangular part

三个方法对应：

- Jacobi: \(M=D\)
- Gauss-Seidel: \(M=D+L\)
- SOR: relaxed Gauss-Seidel

这一页是过渡页，目的是让后面的公式更自然。

---

### 7. Jacobi Method

Jacobi 的核心特点：

> 计算 \(x^{(k+1)}\) 的每一个分量时，只使用上一轮的 \(x^{(k)}\)。

公式：

\[
x_i^{(k+1)}
=
\frac{1}{a_{ii}}
\left(
b_i-\sum_{j\neq i}a_{ij}x_j^{(k)}
\right).
\]

矩阵形式：

\[
Dx^{(k+1)}=b-(L+U)x^{(k)}.
\]

所以：

\[
G_J=-D^{-1}(L+U).
\]

讲法：

> Jacobi is very simple: every component is updated using only old values from the previous iteration.

---

### 8. Jacobi: Strengths and Weaknesses

优点：

- 简单
- 容易实现
- 每个分量可以独立更新，所以容易并行

缺点：

- 通常收敛慢
- 需要存 old vector 和 new vector
- 如果 \(\rho(G_J)\geq 1\)，就不会收敛

可以提一个 sufficient condition：

\[
|a_{ii}|>\sum_{j\neq i}|a_{ij}|
\]

也就是 strict diagonal dominance 可以保证 Jacobi 收敛。

---

### 9. Gauss-Seidel Method

Gauss-Seidel 的核心特点：

> 一旦某个新分量算出来，就立刻在后面的计算中使用。

公式：

\[
x_i^{(k+1)}
=
\frac{1}{a_{ii}}
\left(
b_i
-
\sum_{j<i}a_{ij}x_j^{(k+1)}
-
\sum_{j>i}a_{ij}x_j^{(k)}
\right).
\]

矩阵形式：

\[
(D+L)x^{(k+1)}=b-Ux^{(k)}.
\]

所以：

\[
G_{GS}=-(D+L)^{-1}U.
\]

讲法：

> Compared with Jacobi, Gauss-Seidel uses the newest available information immediately.

---

### 10. Jacobi vs Gauss-Seidel

这一页做方法对比。

Jacobi：

- 只用旧值
- 容易并行
- 通常较慢

Gauss-Seidel：

- 使用最新值
- 通常更快
- 顺序依赖更强，不如 Jacobi 容易并行

重点句：

> Gauss-Seidel usually converges faster because it uses fresh information immediately.

---

### 11. SOR: Successive Over-Relaxation

SOR 是 Gauss-Seidel 的加速版本。

公式：

\[
x_i^{(k+1)}
=
(1-\omega)x_i^{(k)}
+
\frac{\omega}{a_{ii}}
\left(
b_i
-
\sum_{j<i}a_{ij}x_j^{(k+1)}
-
\sum_{j>i}a_{ij}x_j^{(k)}
\right).
\]

也可以理解成：

\[
x_i^{(k+1)}
=
(1-\omega)x_i^{(k)}
+
\omega x_{i,GS}^{(k+1)}.
\]

解释 \(\omega\)：

- \(\omega=1\)：Gauss-Seidel
- \(0<\omega<1\)：under-relaxation
- \(1<\omega<2\)：over-relaxation

对于 symmetric positive definite matrix，SOR 在 \(0<\omega<2\) 时收敛。

---

### 12. SOR in Matrix Form

SOR 也可以写成 stationary iteration：

\[
x^{(k+1)}=G_\omega x^{(k)}+c_\omega.
\]

其中：

\[
G_\omega
=
(D+\omega L)^{-1}\left((1-\omega)D-\omega U\right).
\]

这一页不用讲太久。重点是说明 SOR 也属于 stationary iterative methods，所以同样可以用 spectral radius 来分析收敛。

重点句：

> A good \(\omega\) can reduce the spectral radius and speed up convergence, but a bad \(\omega\) may slow down the method.

---

### 13. Choosing the Relaxation Parameter

这一页讲 SOR 最大的问题：怎么选 \(\omega\)。

对于 1D Poisson problem，有公式：

\[
\rho(G_J)=\cos\left(\frac{\pi}{n+1}\right).
\]

最优参数是：

\[
\omega_{\mathrm{opt}}
=
\frac{2}{1+\sqrt{1-\rho(G_J)^2}}
=
\frac{2}{1+\sin\left(\frac{\pi}{n+1}\right)}.
\]

讲法：

> For structured problems, we may know a good choice of \(\omega\). In general, \(\omega\) is problem-dependent and often chosen empirically.

---

### 14. Numerical Experiment

这一页介绍实验设置。

求解：

\[
-u''(x)=1,\qquad 0<x<1,\qquad u(0)=u(1)=0.
\]

使用 \(n=50\) 个 interior grid points。

停止条件：

\[
\frac{\|b-Ax^{(k)}\|_2}{\|b\|_2}<10^{-8}.
\]

这页主要告诉观众实验是公平的：三个方法解的是同一个问题，使用同一个 tolerance。

---

### 15. Numerical Results

这一页是最重要的实验结果。

结果：

| Method | Iterations |
|---|---:|
| Jacobi | 9653 |
| Gauss-Seidel | 4828 |
| SOR, \(\omega=\omega_{\mathrm{opt}}\) | 189 |

可以讲：

> Jacobi is the slowest. Gauss-Seidel is roughly twice as fast. SOR with a good relaxation parameter is much faster.

这页要强调 SOR 的优势来自合适的 \(\omega\)，不是自动保证。

---

### 16. Effect of \(\omega\)

这一页解释 SOR 对 \(\omega\) 的敏感性。

图展示不同 \(\omega\) 对 iteration count 的影响。

可以讲：

> SOR is not automatically faster. If \(\omega\) is poorly chosen, the method may be slow. Near the optimal value, the iteration count drops dramatically.

重点：

- \(\omega\) 太小：类似 under-relaxation，速度慢
- \(\omega\) 接近 optimal：速度最快
- \(\omega\) 太接近 2：可能变慢甚至不稳定

---

### 17. Implementation Sketch

这一页展示代码逻辑。

三个方法代码结构很像：

- Jacobi：所有分量都用 old values
- Gauss-Seidel：前面已经更新过的分量用 new values
- SOR：在 Gauss-Seidel update 上加 relaxation

如果时间紧，这页可以快速带过。

可以说：

> The full reproducible MATLAB code is in `code/stationary_methods_demo.m`.

---

### 18. Practical Use and Limitations

这一页讲现实意义。

这些方法本身不一定是现在最快的 solver，但仍然重要：

- 可以作为简单 baseline
- 可以作为 multigrid smoother
- 可以作为 preconditioner 思想的基础

限制：

- 对困难的大规模问题可能较慢
- 收敛依赖矩阵性质
- SOR 还依赖参数 \(\omega\)

---

### 19. Takeaways

总结五点：

1. Stationary iterations have the form:

   \[
   x^{(k+1)}=Gx^{(k)}+c.
   \]

2. Convergence depends on:

   \[
   \rho(G)<1.
   \]

3. Jacobi is simple and parallel, but often slow.

4. Gauss-Seidel uses new values immediately and is usually faster.

5. SOR can be much faster with a good \(\omega\), but requires tuning.

---

### 20. References

这一页放老师给的参考书。

不用逐本讲。可以说：

> These are the main references suggested in the topic list.

---

### 21. Thank You

结束页。

可以说：

> Thank you. We are happy to take questions.

---

## 15 分钟时间安排

建议分配：

| Part | Time |
|---|---:|
| Motivation + framework | 4 min |
| Jacobi + Gauss-Seidel | 4 min |
| SOR + choosing \(\omega\) | 3 min |
| Numerical experiment/results | 3 min |
| Conclusion | 1 min |

---

## 四人分工建议

如果你们是四个人，可以这样分：

| Speaker | Slides | 内容 |
|---|---|---|
| Speaker 1 | 1--6 | Introduction, motivation, stationary framework, convergence |
| Speaker 2 | 7--10 | Jacobi, Gauss-Seidel, comparison |
| Speaker 3 | 11--13 | SOR, matrix form, choosing \(\omega\) |
| Speaker 4 | 14--21 | Numerical experiment, results, limitations, conclusion |

如果某个人讲得慢，可以把 Slide 17 Implementation Sketch 简略讲，或者直接跳过。

---

## 汇报时最重要的一句话

整场 presentation 可以一直围绕这句话：

> Jacobi, Gauss-Seidel and SOR are all stationary iterations. They differ in how they split \(A\) and how they use old or new information. Their convergence is controlled by the iteration matrix, especially its spectral radius.

只要这个逻辑讲清楚，整个 presentation 就是连贯的。

