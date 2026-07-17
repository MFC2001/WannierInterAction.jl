
## interaction between pairs

$$
p_1 = (i0, lR_1) ,\;
p_2 = (k0, jR_2)
$$

$$
\left\langle i0;jR+R_2 \left| V \right| kR;lR_1 \right\rangle 
= \left\langle p_1 \left| V \right| p_2 \right\rangle 
= \iint \mathrm{d}r \mathrm{d}r' \psi^{}_{p_10}(r) \frac{1}{|r-r'|} \psi^{*}_{p_2R}(r') \\
= \iint \mathrm{d}r \mathrm{d}r' \psi^{}_{p_10}(r) \frac{1}{|r-r'|} \psi^{*}_{p_20}(r'-R) \\
= \iint \mathrm{d}r \mathrm{d}r' \psi^{}_{p_10}(r) \frac{1}{|r-r'-R|} \psi^{*}_{p_20}(r') \\
\lim_{R\rightarrow\infty} \left\langle i0;jR+R_2 \left| V \right| kR;lR_1 \right\rangle 
= \frac{1}{|R|} \int \mathrm{d}r \psi^{}_{p_10}(r) \int \mathrm{d}r'\psi^{*}_{p_20}(r') 
= \frac{1}{|R|} \delta_{il} \delta_{jk} \delta_{R_10} \delta_{R_20}
$$

这里使用了wannier基的局域性、正交性。
只有U-type项涉及到了长程性，故只有U-type项需要头项修正、镜像修正，以及后续计算激子所需的长程修正。

程序中非U-type项也使用了头项修正，但并不会影响收敛的结果。（或许也有加快收敛的作用？）

## Gauss Coulomb potential used for mirror correction of direct term

Under plane-wave basis, the direct terms between wannier basis contain the interaction with all the mirror system. This will cause deviation, expecially for the term near the periodic boundary of system. So we need to correct it.

The direct term contains the interaction with all the mirror system can be written as

$$
U_{i0jR} = \sum_{\mathcal{R}} U'_{i0jR+\mathcal{R}} 
= U'_{i0jR} + \sum_{\mathcal{R} \neq 0} U'_{i0jR+\mathcal{R}} 
\approx U'_{i0jR} + \sum_{\mathcal{R} \neq 0} \phi(\mathcal{R} + R + R_j - R_i) 
= U'_{i0jR} + U^{c}_{i0jR} \\
U^{c}_{i0jR} = \sum_{\mathcal{R} \neq 0} \phi(\mathcal{R} + R + R_j - R_i) 
= \sum_{\mathcal{R}} \phi(\mathcal{R} + R + R_j - R_i) - \phi(R + R_j - R_i)
$$

Here we let $\mathcal{R}$ represents all the mirror system. We require the size of system is big enough, so that each $U'_{i0jR+\mathcal{R}}$ where $\mathcal{R} \neq 0$ approximates to $1/r$ and Gauss Coulomb potential $\phi(r)$.

With Fourier transform, we have

$$
\sum_{\mathcal{R}} \phi(\mathcal{R} + R + R_j - R_i) 
= \frac{1}{N} \sum_{k} \sum_{G} \phi(k+G) e^{i(k+G)(R + R_j - R_i)} + \frac{1}{N} \phi(head) \\
$$

These derivation is appicable to all demensions, only if we use the $\phi(k)$ with appropriate dimension.

Generally speaking, if we don't replace divergent term with head term, we won't use head term here too. But we can't process a infinite term. There is a alternative scheme. Under plane-wave basis, the divergent term for 3D is $\frac{1}{k^2}$, or for 2D is $\frac{1}{k}$, which are the same as gauss potential. So we can use the same head term here.

$$
\sum_{\mathcal{R}} \phi(\mathcal{R} + R + R_j - R_i) 
= \frac{1}{N} \sum_{k} \sum_{G} \phi(k+G) e^{i(k+G)(R + R_j - R_i)} + v(head) \\
$$

进行镜像修正的一个优势为：倒空间中库伦势的发散行为都是相同的，故可以认为其互相抵消，从而完全不需要考虑所谓积分值近似的误差。
只需要保证系统足够大，从而使得周期边界边界附近真实的wannier基之间的相互作用趋于选定的库伦势（这里即为高斯势）即可。

## 3D

omit $\frac{e^2}{4\pi\epsilon_0}$

$$
\phi(r) = \frac{erf(\alpha |r|)}{|r|} \\
\phi(k) = \frac{4\pi}{\Omega} \frac{1}{k^2} e^{-\frac{k^2}{4\alpha^2}} \\
\phi(head) = \frac{4\pi}{\Omega} N\Omega \frac{\alpha}{2\pi^{\frac{3}{2}}} erf(\frac{q_{sz}}{2\alpha}) 
\rightarrow \frac{4\pi}{\Omega} N\Omega \frac{q_{sz}}{2\pi^2} 
,\; q_{sz} = \left( \frac{6\pi^2}{N\Omega} \right)^{\frac{1}{3}} \\
$$

$$
U_{i0jR} = U'_{i0jR} + U^{c}_{i0jR} \\
U^{c}_{i0jR} = \frac{1}{N} \sum_{k} \sum_{G} \phi(k+G) e^{i(k+G)(R + R_j - R_i)} + \frac{1}{N} \phi(head) \\
$$


## 2D

实际测试表明，Keldysh势趋于 $1/r$ 的行为不太好，但可以确定的一点是，通常的基于平面波基组的计算中，超胞的大小最多能保证 Keldysh势以1e-3的精度趋于 $1/r$，并且这是很极限的情况。
故两个方案：Gauss2D和KeldyshGauss，更推荐使用后者。

$$
U_{i0jR} = U'_{i0jR} + U^{c}_{i0jR} \\
U^{c}_{i0jR} = \sum_{\mathcal{R}} \phi(\mathcal{R} + R + R_j - R_i) 
$$

$$
\phi(r) = \frac{\pi}{2\tilde{\epsilon}r_0} \left[ H_0(\frac{r}{r_0}) - Y_0(\frac{r}{r_0}) \right]\\
\phi(k) = \frac{2\pi}{\tilde{\epsilon}A} \frac{1}{k(1+r_0k)} \\
\phi(head) = \frac{4\pi}{\Omega} N\Omega \frac{\alpha}{2\pi^{\frac{3}{2}}} erf(\frac{q_{sz}}{2\alpha}) 
\rightarrow \frac{4\pi}{\Omega} N\Omega \frac{q_{sz}}{2\pi^2} 
,\; q_{sz} = \left( \frac{6\pi^2}{N\Omega} \right)^{\frac{1}{3}} \\
$$


## Keldysh

$$
\nabla^2 \phi(r) = - 4\pi \delta(r) \Rightarrow
\nabla^2 \phi(r) = - 4\pi \delta(r) - 4\pi \alpha_{2D} \nabla_{\rho}^{2} \phi(\rho, z = 0) \delta(z)
$$

$$
\phi(r) = \frac{1}{(2\pi)^3} \int \mathrm{d}k \phi(k) e^{ikr}
$$

$$
-\frac{1}{(2\pi)^3} \int \mathrm{d}k (q^2 + k_z^2) \phi(k) e^{ikr} 
= - 4\pi \frac{1}{(2\pi)^3} \int \mathrm{d}k e^{ikr} + 4\pi \alpha_{2D} q^{2} \frac{1}{(2\pi)^3} 
\int \mathrm{d}q \mathrm{d}k_z' \phi(q, k_z') e^{iq\rho} \frac{1}{2\pi} \int \mathrm{d}k_z e^{ik_zz}
$$

$$
(q^2 + k_z^2) \phi(k)
= 4\pi - 4\pi \alpha_{2D} q^{2} \int \frac{\mathrm{d}k_z'}{2\pi} \phi(q, k_z') 
= 4\pi - 4\pi \alpha_{2D} q^{2} \phi_{2D}(q) \\
\phi_{2D}(q)
= \left[ 4\pi - 4\pi \alpha_{2D} q^{2} \phi_{2D}(q) \right] \int_{-\infty}^{+\infty} \frac{\mathrm{d}k_z}{2\pi} \frac{1}{q^2 + k_z^2} \\
= \left[ 4\pi - 4\pi \alpha_{2D} q^{2} \phi_{2D}(q) \right] \frac{1}{2\pi q} \left. \arctan(\frac{k_z}{q}) \right|_{-\infty}^{+\infty} \\
= \left[ 4\pi - 4\pi \alpha_{2D} q^{2} \phi_{2D}(q) \right] \frac{1}{2 q} \\
$$

$$
\phi_{2D}(q) = \frac{2\pi}{|q|(1+r_0|q|)} ,\; r_0 = 2\pi \alpha_{2D}
$$

in real space 

$$
\phi_{2D}(\rho) = \frac{1}{(2\pi)^2} \int \mathrm{d}q \phi_{2D}(q) e^{iq\rho} \\
= \int_{0}^{+\infty} \mathrm{d}q \frac{1}{1+r_0q} \frac{1}{2\pi} \int_{-\pi}^{\pi} \mathrm{d}\theta e^{iq\rho\sin\theta} \\
= \int_{0}^{+\infty} \mathrm{d}q \frac{J_{0}(q\rho)}{1+r_0q} \\
= \frac{\pi}{2r_0} \left[ H_{0}(\frac{\rho}{r_0}) - Y_{0}(\frac{\rho}{r_0})\right]
$$


if Gauss2D

$$
\delta(r) \rightarrow \left( \frac{\alpha}{\sqrt{\pi}} \right)^{2} e^{-\alpha^2r^2} \delta(z)
= \frac{1}{(2\pi)^3} \int \mathrm{d}q e^{-\frac{q^2}{4\alpha^2}} e^{iq\rho} \int \mathrm{d}k_z e^{ik_zz}
$$

$$
(q^2 + k_z^2) \phi(k)
= 4\pi e^{-\frac{q^2}{4\alpha^2}} - 4\pi \alpha_{2D} q^{2} \phi_{2D}(q) \\
\phi_{2D}(q)
= \frac{2\pi e^{-\frac{q^2}{4\alpha^2}}}{|q|(1+r_0|q|)} ,\; r_0 = 2\pi \alpha_{2D}
$$

in real space 

$$
\phi_{2D}(\rho) 
= \int_{0}^{+\infty} \mathrm{d}q \frac{J_{0}(q\rho) e^{-\frac{q^2}{4\alpha^2}}}{1+r_0q} \\
= \frac{1}{r_0} \int_{0}^{+\infty} \mathrm{d}x \frac{J_{0}(x) e^{-\frac{x^2}{4\alpha^2\rho^2}}}{\frac{\rho}{r_0}+x} \\
$$


if Gauss

$$
\delta(r) \rightarrow \left( \frac{\alpha}{\sqrt{\pi}} \right)^{2} e^{-\alpha^2\rho^2}
\frac{\alpha_z}{\sqrt{\pi}} e^{-\alpha_z^2z^2}
= \frac{1}{(2\pi)^3} \int \mathrm{d}q \mathrm{d}k_z 
e^{-\frac{q^2}{4\alpha^2}} e^{-\frac{k_z^2}{4\alpha_z^2}} e^{iq\rho} e^{ik_zz}
$$

$$
(q^2 + k_z^2) \phi(k)
= 4\pi e^{-\frac{q^2}{4\alpha^2}} e^{-\frac{k_z^2}{4\alpha_z^2}} - 4\pi \alpha_{2D} q^{2} \phi_{2D}(q) \\
$$

$$
\phi_{2D}(q)
=  4\pi e^{-\frac{q^2}{4\alpha^2}} \int_{-\infty}^{+\infty} \frac{\mathrm{d}k_z}{2\pi} 
\frac{e^{-\frac{k_z^2}{4\alpha_z^2}}}{q^2 + k_z^2} 
-4\pi \alpha_{2D} q^{2} \phi_{2D}(q) \int_{-\infty}^{+\infty} \frac{\mathrm{d}k_z}{2\pi} \frac{1}{q^2 + k_z^2} \\
= \frac{2\pi}{q} e^{-\frac{q^2}{4\alpha^2}} erfcx(\frac{q}{2\alpha_z}) - r_0 q \phi_{2D}(q) \\
\phi_{2D}(q) = \frac{2\pi e^{-\frac{q^2}{4\alpha^2}} erfcx(\frac{q}{2\alpha_z})}{|q|(1+r_0|q|)} ,\; r_0 = 2\pi \alpha_{2D}
$$

in real space 

$$
\phi_{2D}(\rho) 
= \int_{0}^{+\infty} \mathrm{d}q \frac{J_{0}(q\rho) e^{-\frac{q^2}{4\alpha^2}} erfcx(\frac{q}{2\alpha_z})}{1+r_0q} \\
$$
