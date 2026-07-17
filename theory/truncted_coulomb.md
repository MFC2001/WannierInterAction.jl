
例如对于二维系统，使用二维截断库伦势，可以在真空层不太厚的情况下获得很好的收敛结果。

理论上来说，计算低维系统时，任何涉及到库伦势的计算（包括DFT、节点函数、GW等）都需要使用截断库伦势，并且要注意系统应该在元胞非周期方向的正中间（这里需要注意软件中默认非周期方向的定义）。

# Truncted Coulomb potential

In order to weaken the effect of the mirror charge in the aperiodic direction, we use a truncated Coulomb potential $\widetilde{v}(r)$ instead of the bare 3D Coulomb potential $v(r)$: (omit $\frac{e^2}{\epsilon_0}$)

$$
\widetilde{v}(r) = \left\{
\begin{aligned}
&\frac{1}{4\pi |r|} &if\;r \in \mathcal{D} \\
&0           &if\;r \notin \mathcal{D}
\end{aligned}
\right.
$$

where $\mathcal{D}$ is the screening region of the system. So we can get modified Fourier integral

$$
\widetilde{v}(k) = \int \mathrm{d}r \widetilde{v}(r) e^{ikr} = \int_{\mathcal{D}} \mathrm{d}r v(r) e^{ikr}
$$

## 3D

$$
v(k) = \frac{1}{4\pi} \int \mathrm{d}r \frac{1}{|r|} e^{ikr} 
= \frac{1}{k^2}
$$

$$
v(k, G) = \frac{1}{|k+G|^2}
$$

head

$$
v(0, 0) \rightarrow \frac{N\Omega}{(2\pi)^3} \int_{k \in miniFBZ} \mathrm{d}k \frac{1}{k^2}
$$

## 2D

若不使用截断库伦势：

From the figure.2 of [paper](https://wvpn.ustc.edu.cn/https/77726476706e69737468656265737421e0e2438f69227b5330079bab/en/content/articlelanding/2021/nr/d1nr04896a), 
关于介电函数，实际上计算的是包含真空层的元胞的3D体系，真空层的贡献会稀释实际二维体系的贡献，导致最终结果为真空层越厚，宏观量越趋于真空。

平面波基组下的相互作用参数计算应使用截断库伦势，满足的经典势为静电力学的Keldysh势，适用的宏观量为屏蔽长度 $r_0$ 或者二维极化率 $\alpha_{2D}$。
宏观量的计算依赖于RPA+不截断的库伦势计算出的三维面内方向宏观介电常数，然后使用公式[16](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.84.085406)。

BerkeleyGW 中定义z方向受限。

$$
\widetilde{v}^{2D}(k) 
= \frac{1}{4\pi} \int_{-R}^{R}\mathrm{d}z \int_{-\infty}^{+\infty} \mathrm{d}y \int_{-\infty}^{+\infty} \mathrm{d}x 
\frac{e^{i(k_xx+k_yy+k_zz)}}{\sqrt{x^2+y^2+z^2}} \\
= \frac{1}{k^2} \left[ 1 + e^{-k_{\parallel}R} \frac{k_z}{k_{\parallel}} \sin(k_zR) - e^{-k_{\parallel}R} \cos(k_zR) \right]
,\;k_{\parallel} =\sqrt{k_x^2+k_y^2} > 0 \\
$$

$$
\widetilde{v}^{2D}(k, G) = \frac{1}{|k+G|^2} \left[ 1 - e^{-|k+G|_{\parallel}R} \cos(G_zR) \right]
,\;G_{z} = \frac{n\pi}{R} \\
$$

BerkeleyGW 中的 average_w 参数指出，会对W使用一个模型介电函数来获得头项。

head

bare

$$
\pi q_{sz}^2 = \frac{(2\pi)^2}{NS} \Rightarrow q_{sz} = \sqrt{\frac{4\pi}{NS}} \\
v(q) = \frac{1}{q^2} \left( 1 - e^{-|q|R} \right) \\
\frac{NS}{(2\pi)^2} \int_{|q| < q_{sz}} \mathrm{d}^2 q v(q) 
= \frac{NS}{2\pi} \int_{0}^{q_{sz}} \frac{1 - e^{-|q|R}}{q} \mathrm{d}q 
$$

screened

$$
\epsilon(q) = 1 + \gamma q^2 v(q) e^{-\alpha|q|} \\
\frac{NS}{(2\pi)^2} \int_{|q| < q_{sz}} \mathrm{d}^2 q \frac{v(q)}{\epsilon(q)} 
= \frac{NS}{2\pi} \int_{0}^{q_{sz}} \frac{1 - e^{-|q|R}}{q\left[ 1 + 4\pi \gamma \left( 1 - e^{-|q|R} \right) e^{-\alpha|q|} \right]} \mathrm{d}q 
$$

In BerkeleyGW, set $\alpha$ as 0, then $\gamma = \frac{\epsilon(q_0) - 1}{q_0^2 v(q_0)}$

这个介电函数的公式可能只适用于静态情况，详见[Fig.1](http://dx.doi.org/10.1103/PhysRevB.73.233103)。

## 1D

BerkeleyGW 中定义xy方向受限。

(in a cylinder, from equation [12-13](http://dx.doi.org/10.1103/PhysRevB.73.205119))
$$
\widetilde{v}^{1D}(k) 
= \frac{1}{4\pi} \int_{\sqrt{x^2+y^2} < R} \mathrm{d}x \mathrm{d}y \int_{-\infty}^{+\infty} \mathrm{d}z 
\frac{e^{i(k_xx+k_yy+k_zz)}}{\sqrt{x^2+y^2+z^2}} \\
= \frac{1}{4\pi} \int_{0}^{R} \mathrm{d}r \int_{0}^{2\pi} \mathrm{d}\theta \int_{-\infty}^{+\infty} \mathrm{d}z 
\frac{e^{ik_zz}e^{ik_{\perp} r \sin\theta}}{\sqrt{r^2+z^2}} \\
= \int_{0}^{R} \mathrm{d}r 
\frac{1}{2} \int_{-\infty}^{+\infty} \mathrm{d}z \frac{e^{ik_zz}}{\sqrt{r^2+z^2}}
\frac{1}{2\pi} \int_{0}^{2\pi} \mathrm{d}\theta e^{ik_{\perp} r \sin\theta} \\
= \int_{0}^{R} \mathrm{d}r r K_0(|k_z|r) J_0(k_{\perp} r) \\
= \frac{1}{k^2} \left[ 1 + k_{\perp} R J_1(k_{\perp} R) K_{0}(|k_z|R) - |k_z|R J_{0}(k_{\perp} R) K_{1}(|k_z|R) \right] 
$$

$$
\widetilde{v}^{1D}(k, G) 
= \frac{1}{|k+G|^2} \left[ 1 + G_{\perp} R J_1(G_{\perp} R) K_{0}(|k_z+G_z|R) - |k_z+G_z|R J_{0}(G_{\perp} R) K_{1}(|k_z+G_z|R) \right] 
$$

head

$$
\widetilde{v}^{1D}(k, 0) = \frac{1}{k^2} \left[ 1 - |k|R K_{1}(|k|R) \right] \\
q_{sz} = \frac{\pi}{NR_z} \\
\frac{NR_z}{2\pi} \int_{-q_{sz}}^{q_{sz}} \widetilde{v}^{1D}(k, 0) \mathrm{d}k 
= \frac{NR_z}{\pi} \int_{0}^{q_{sz}} \frac{1}{k^2} \left[ 1 - kR K_{1}(kR) \right] \mathrm{d}k 
$$

Here $K$ is the Modified Bessel function of the second kind.

(in the unitcell, the method used by BerkeleyGW)
$$
\widetilde{v}^{1D}(k) = \frac{1}{4\pi} 
\int_{-\frac{R_x}{2}}^{\frac{R_x}{2}} \int_{-\frac{R_y}{2}}^{\frac{R_y}{2}} \mathrm{d}x \mathrm{d}y 
\int_{-\infty}^{+\infty} \mathrm{d}z \frac{e^{i(k_xx+k_yy+k_zz)}}{\sqrt{x^2+y^2+z^2}} \\
= \frac{1}{4\pi} 
\int_{-\frac{R_x}{2}}^{\frac{R_x}{2}} \int_{-\frac{R_y}{2}}^{\frac{R_y}{2}} \mathrm{d}x \mathrm{d}y e^{i(k_xx+k_yy)} 
\int_{-\infty}^{+\infty} \mathrm{d}z \frac{e^{ik_zz}}{\sqrt{x^2+y^2+z^2}} \\
= \frac{1}{2\pi} 
\int_{-\frac{R_x}{2}}^{\frac{R_x}{2}} \int_{-\frac{R_y}{2}}^{\frac{R_y}{2}} 
e^{i(k_xx+k_yy)} K_0(|k_z|\sqrt{x^2+y^2}) \mathrm{d}x \mathrm{d}y \\
$$

$$
\widetilde{v}^{1D}(k, G) = \frac{1}{2\pi} 
\int_{-\frac{R_x}{2}}^{\frac{R_x}{2}} \int_{-\frac{R_y}{2}}^{\frac{R_y}{2}} 
e^{i(G_xx+G_yy)} K_0(|k_z + G_z|\sqrt{x^2+y^2}) \mathrm{d}x \mathrm{d}y \\
$$

$$
f^{1D}(k, G) = \frac{|k+G|^2}{2\pi} 
\int_{-\frac{R_x}{2}}^{\frac{R_x}{2}} \int_{-\frac{R_y}{2}}^{\frac{R_y}{2}} 
e^{i(G_xx+G_yy)} K_0(|k_z + G_z|\sqrt{x^2+y^2}) \mathrm{d}x \mathrm{d}y \\
= \frac{2|k+G|^2}{\pi} 
\int_{0}^{\frac{R_x}{2}} \int_{0}^{\frac{R_y}{2}} \cos(|G_x|x) \cos(|G_y|y)
K_0(|k_z + G_z|\sqrt{x^2+y^2}) \mathrm{d}x \mathrm{d}y \\
= \frac{2|k+G|^2}{\pi} \frac{R_x}{2n_i} \frac{R_y}{2n_j} \sum_{i=1}^{n_i} \sum_{j=1}^{n_j}
\cos(\frac{|G_x|R_x}{2\pi}\frac{\pi(i-0.5)}{n_i}) \cos(\frac{|G_y|R_y}{2\pi}\frac{\pi(j-0.5)}{n_j}) 
K_0(|k_z + G_z|\sqrt{x_i^2+y_j^2}) ,\;
x_i = \frac{R_x(i-0.5)}{2n_i} ,\; y_j = \frac{R_y(j-0.5)}{2n_j} \\
$$

type-Ⅱ DCT use a grif with half shift.

$$
\widetilde{G}_x = \frac{|G_x|R_x}{2\pi} ,\; \widetilde{G}_y = \frac{|G_y|R_y}{2\pi} \\
f_{ij}^{kGz} = K_0(|k_z + G_z|\sqrt{x_i^2+y_j^2}) ,\; x = \frac{R_x(i-0.5)}{2n_i} ,\; y = \frac{R_y(j-0.5)}{2n_j} \\
f^{1D}(k, G) = \frac{|k+G|^2}{8\pi} \frac{R_xR_y}{n_in_j} \text{REDFT10}[f^{kGz}]_{\widetilde{G}_x\widetilde{G}_y} \\
$$

head

$$
\widetilde{v}^{1D}(k, 0) = \frac{2}{\pi} 
\int_{0}^{\frac{R_x}{2}} \int_{0}^{\frac{R_y}{2}} 
K_0(|k|\sqrt{x^2+y^2}) \mathrm{d}x \mathrm{d}y \\
q_{sz} = \frac{\pi}{NR_z} \\
\frac{NR_z}{2\pi} \int_{-q_{sz}}^{q_{sz}} \widetilde{v}^{1D}(k, 0) \mathrm{d}k 
= \frac{2NR_z}{\pi^2} \int_{0}^{q_{sz}} \mathrm{d}k 
\int_{0}^{\frac{R_x}{2}} \mathrm{d}x 
\int_{0}^{\frac{R_y}{2}} \mathrm{d}y 
K_0(k\sqrt{x^2+y^2}) 
$$

$$
R = min(\frac{R_x}{2}, \frac{R_y}{2}) \\
\frac{1}{4} \int_{r < R} K_0(k\sqrt{x^2+y^2}) \mathrm{d}x \mathrm{d}y 
= \frac{\pi}{2} \int_{0}^{R} K_0(kr) r \mathrm{d}r 
= \frac{\pi}{2} \frac{1}{k^2} \left[ 1 - kR K_{1}(kR) \right]
$$

screened

$$
\epsilon(q) = 1 + \gamma q^2 v(q) \\
$$

## 0D

(spherical truncation)
$$
\widetilde{v}^{0D}(k) = \frac{1}{4\pi} \int_{|r| < R} \mathrm{d}^3r \frac{e^{i k \cdot r}}{|r|} \\
= \frac{1}{4\pi} \int_{0}^{R} \mathrm{d}r \int_{0}^{\pi} \mathrm{d}\theta \int_{0}^{2\pi} \mathrm{d}\varphi 
r^2 \sin\theta \frac{e^{i|k|r \cos\theta }}{r} \\
= \frac{1}{2} \int_{0}^{R} \mathrm{d}r \int_{-1}^{1} \mathrm{d}x
r e^{i|k|r x} \\
= \frac{1}{|k|} \int_{0}^{R} \sin(|k|r) \mathrm{d}r \\
= \frac{1}{|k|} \int_{0}^{R} \sin(|k|r) \mathrm{d}r \\
= \frac{1}{k^2} \left[ 1 - \cos\left( |k| R\right) \right]
$$

$$
\widetilde{v}^{0D}(k, G) 
= \frac{1}{G^2} \left[ 1 - \cos\left( |G| R\right) \right] ,\; 
\widetilde{v}^{0D}(k, 0) = \frac{R^2}{2} 
$$


(box truncation)
$$
\widetilde{v}^{0D}(k) = \frac{1}{4\pi} 
\int_{-\frac{R_x}{2}}^{\frac{R_x}{2}} \int_{-\frac{R_y}{2}}^{\frac{R_y}{2}} \int_{-\frac{R_z}{2}}^{\frac{R_z}{2}} 
\frac{e^{i k \cdot r}}{|r|} \mathrm{d}^3r \\
$$

$$
\widetilde{v}^{0D}(k, G) = \frac{1}{4\pi} 
\int_{-\frac{R_x}{2}}^{\frac{R_x}{2}} \int_{-\frac{R_y}{2}}^{\frac{R_y}{2}} \int_{-\frac{R_z}{2}}^{\frac{R_z}{2}} 
\frac{e^{i G \cdot r}}{|r|} \mathrm{d}^3r \\
$$

$$
f^{0D}(k, G) = \frac{|G|^2}{4\pi} 
\int_{-\frac{R_x}{2}}^{\frac{R_x}{2}} \int_{-\frac{R_y}{2}}^{\frac{R_y}{2}} \int_{-\frac{R_z}{2}}^{\frac{R_z}{2}} 
\frac{e^{i G \cdot r}}{|r|} \mathrm{d}^3r \\
= \frac{2|G|^2}{\pi} 
\int_{0}^{\frac{R_x}{2}} \int_{0}^{\frac{R_y}{2}} \int_{0}^{\frac{R_z}{2}} 
\cos(|G_x|x) \cos(|G_y|y) \cos(|G_z|z) \frac{1}{\sqrt{x^2 + y^2 + z^2}} \mathrm{d}x \mathrm{d}y \mathrm{d}z \\
= \frac{2|G|^2}{\pi} \frac{R_x}{2n_i} \frac{R_y}{2n_j} \frac{R_z}{2n_k}
\sum_{i = 1}^{n_i} \sum_{j = 1}^{n_j} \sum_{k = 1}^{n_k} 
\cos(\frac{|G_x|R_x}{2\pi}\frac{\pi(i-0.5)}{n_i}) 
\cos(\frac{|G_y|R_y}{2\pi}\frac{\pi(j-0.5)}{n_j}) 
\cos(\frac{|G_z|R_z}{2\pi}\frac{\pi(k-0.5)}{n_k}) 
\frac{1}{\sqrt{x_i^2 + y_j^2 + z_k^2}} ,\;
x_i = \frac{R_x(i-0.5)}{2n_i} ,\; y_j = \frac{R_y(j-0.5)}{2n_j} ,\; z_k = \frac{R_z(k-0.5)}{2n_k} \\
$$

$$
f_{ijk} = \frac{1}{\sqrt{x_i^2 + y_j^2 + z_k^2}} ,\;
x_i = \frac{R_x(i-0.5)}{2n_i} ,\; y_j = \frac{R_y(j-0.5)}{2n_j} ,\; z_k = \frac{R_z(k-0.5)}{2n_k} \\
f^{0D}(k, G) = \frac{|G|^2}{32\pi} \frac{R_xR_yR_z}{n_in_jn_k} \text{REDFT10}[f]_{\widetilde{G}_x\widetilde{G}_y\widetilde{G}_z} \\
$$

# DCT

for `FFTW.dct`

$$
DCT[f]_k = \sqrt{\frac{2}{n}} \sum_{i}^{n} \cos(k\frac{\pi(i-0.5)}{n}) f_i \\
DCT[f]_0 = \sqrt{\frac{1}{n}} \sum_{i}^{n} f_i
$$

it's more convenient to use `FFTW.r2r(A, FFTW.REDFT10)`

$$
\text{REDFT10}[f]_k = 2 \sum_{i}^{n} \cos(k\frac{\pi(i-0.5)}{n}) f_i \\
$$