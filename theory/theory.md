# Coulomb interaction between wannier orbitals




## Hermitian

相互作用算符应为一个厄米算符

$$
\hat{I} = 
\sum_{ijkl} V_{ijkl} \hat{C}_{i}^{\dagger} \hat{C}_{j}^{\dagger} \hat{C}_{k}^{} \hat{C}_{l}^{} \\
\hat{I} = \hat{I}^{\dagger} =
\sum_{ijkl} V_{ijkl}^{*} \hat{C}_{l}^{\dagger} \hat{C}_{k}^{\dagger} \hat{C}_{j}^{} \hat{C}_{i}^{} 
\Rightarrow 
V_{ijkl}^{} = V_{lkji}^{*}
$$

这等价于要求

$$
V(r, r') = V^{*}(r, r') \Rightarrow \\
\iint \psi_{i}^{*}(r) \psi_{j}^{*}(r') V^{}(r, r') \psi_{k}^{}(r') \psi_{l}^{}(r) \mathrm{d}r \mathrm{d}r'
= \iint \psi_{i}^{*}(r) \psi_{j}^{*}(r') V^{*}(r, r') \psi_{k}^{}(r') \psi_{l}^{}(r) \mathrm{d}r \mathrm{d}r' \\
= \left[ \iint \psi_{l}^{*}(r) \psi_{k}^{*}(r') V^{}(r, r') \psi_{j}^{}(r') \psi_{i}^{}(r) \mathrm{d}r \mathrm{d}r' \right]^{*}
$$

屏蔽库伦势类似。

For U-type interaction, we have 

$$
V_{ij}^{*} = V_{ijji}^{} = V_{ijji}^{*} = V_{ij}^{*} 
\Rightarrow \text{Im} V_{ij} = 0
$$

For pair-type interaction, we have 

$$
p_1 = (i,l,R_1); \; p_2 = (k, j, R_2) \\
V_{p_10;p_2R} = \left\langle i0;jR+R_2 \left| V \right| kR;lR_1 \right\rangle 
= \left\langle lR_1;kR \left| V \right| jR+R_2;i0 \right\rangle 
= \left\langle l0;kR-R_1 \left| V \right| jR+R_2-R_1;i-R_1 \right\rangle 
= V_{p_1^*0;p_2^*R+R_2-R_1} \\
p_1^* = (l,i,-R_1); \; p_2^* = (j, k, -R_2) \\
$$

## wospin

时间反演对称性下，wannier函数一般是实函数，无屏蔽以及静态($\omega = 0$)的结果一般为实数。

## wispin

计算 $\rho$ 的时候，分别计算不同自旋的部分，然后求和即可。
故若 $\psi_{ik}$ 具有确定的自旋，可以直接忽略涉及到不同自旋态对应的 $\rho$。

## divergence of $\rho_k$

$\lim_{q \rightarrow 0} \rho_k(G=0)$ 的对角元发散，非对角元不连续，但依赖于$q$趋于0的方向。
与这些项相关的部分为头项和翼项，它们有各自的处理方式。

实际计算中，显式的使用有限值替换 $\rho_k(G=0)$，具体为舍去 $1/|q|$, 剩余部分正常计算。
在最终求和的步骤中，显式地处理掉头项和翼项。

头项：使用积分值替换 $1/|q^2|$，截断情况下略有不同。
翼项：事实上，翼项在任何时候都应设为0，因为其正比于 $q$, 在 $\Gamma$ 点附近取平均自然得到0。(参考自BerkeleyGW-Manual)

正常计算：直接舍去翼项；用户可额外提供negq0的介电矩阵：显式地计入翼项，从结果中可以明确观察到发散性被抵消。

关于negq0的设计逻辑：实际上完全可以用户手动对介电函数取平均然后输入；倘若用户需要在 $\Gamma$ 点附近取多个q点平均值，可以手动对介电矩阵取平均后输入（确保G基组的顺序相同）。

不同类型的相互作用参数的处理细节：
- 对于q0都忽略掉相位因子$e^{iqR}$;
- 头项，U额外计入 head - conj(ρR_left_G0) * ρR_right_G0, 其余额外减去 conj(ρR_left_G0) * ρR_right_G0。
- 翼项，without negq0时，直接删去；with negq0时，显式计入。

# symmetry

关于利用对称性减少对q点求和的方案不可行。
原因在于本质上是要对 $\left\langle iR_1;jR_2 \left| W(V) \right| kR_3;lR_4 \right\rangle$ 整体应用对称性操作，然后寻找等价q点，
但对称操作只会得到实空间相互作用参数之间的映射，而不是参数计算表达式内部q点之间的映射。
或许可以尝试利用这一点减少计算的相互作用参数的数目，但实现起来可能会很复杂，故先搁置。
并且实际上理想的wannier模型也只是继承了晶体的部分对称性，而不是完整的对称群。

随记：在Souza-Marzari-Vanderbilt的原始论文（Phys. Rev. B 65, 035109 (2001)）中明确指出，解纠缠过程可以并且应该在每个k点选择子空间时，强制要求：
对于每个晶体对称操作 {R|T}，如果它将k点映射到k'点，那么在该对称操作下，k点的子空间应该映射到k'点的子空间。




# practical computation

## restriction

介电函数和wannier函数可能来自两个不同的计算任务，二者匹配的最小限制为：
1. k网格，最大局域化要求波函数定义在均匀的MP网格上，不可以存在偏移，介电矩阵一般计算一个无偏移的MP网格上的(gamma点被小q替换)。
   对k网格的限制要求k+q仍在k网格上，故两个MP网格需大小相等，这也是物理上系统大小相同的要求。
2. 如何让离散求和结果接近热力学极限情况下的积分，关键在于q=0的项：根据从BerkeleyGW出发的做法，翼项是伪发散，头项是真发散；
   小量q对应的介电矩阵定义良好，且翼项的发散性会相互抵消，但是介电矩阵的头项不发散，故这里仍需要使用有限积分计算头项。
   By the way, RESCPACK的直接计算q=0的方法一定有问题，其翼项手动除去了发散，导致其局域场效应考虑不完全，介电常数计算有问题，从而导致头项也有问题。
3. 自旋设置相同（目前仅支持标量版本的介电函数，也就是介电函数与自旋无关，BerkeleyGW和yambo都是如此）；
4. 来自于同一个晶体结构（晶体结构相关量仅有 $\Omega$，程序可从介电函数处或波函数处读取，但或许可将确保这一点的责任留给用户）；
5. 一般来说，介电矩阵的截断能比wannier函数低，但实际计算中涉及到的为电荷密度截断能与介电函数截断能，并且电荷密度Glist的生成过程保证最低等于介电函数；

实际计算中，winegq0：正常计算；wonegq0: 程序中会自动将翼项删去。

## tip

1. 插值得到的wavek_smallq影响很小，预计可以直接使用其邻近k点的波函数，但是插值很快，并非计算瓶颈，可以留开关；测试结果表明是否启用误差很小。
2. 头项发散性：连续积分替换，不同维度选用不同的积分形式以及库伦势；仅U需要手动处理头项，其余类型参数中头项趋于0(pair的计算验证了这一点)。
3. 翼项发散性：数值计算验证得到，3D且wospin情况下，发散性体现在异常虚部，且 $\pm q$ 符号相反，故可直接舍去；更严格而言，可以选择额外引入 $-q$ 的计算。
4. 



## head

### 3D

三维球积分

$$
\frac{4}{3} \pi q_{sz}^{3} = \frac{(2\pi)^3}{N\Omega} \Rightarrow q_{sz} = \left( \frac{6\pi^2}{N\Omega} \right)^{\frac{1}{3}} \\
v(k=0) = \frac{N\Omega}{(2\pi)^3} \int_{k<q_{sz}} k^2 \sin \theta \mathrm{d}k \mathrm{d}\theta \mathrm{d}\varphi \frac{1}{k^2} 
= \frac{N\Omega q_{sz} }{2 \pi^2} \\
\frac{1}{|q_0|^2} \rightarrow \frac{N\Omega q_{sz} }{2 \pi^2}
$$

$$
head = \frac{e^2}{N\Omega\epsilon_0} e^{iq_0(R_1-R_3)} 
\rho^{*}(0, R_4-R_1, l, i, q_0) \rho^{}(0, R_2-R_3, j, k, q_0) \\
head\_correction = \frac{e^2}{N\Omega\epsilon_0} e^{iq_0(R_1-R_3)} \left[ 
   \frac{N\Omega q_{sz} }{2 \pi^2} - \frac{1}{q_0^2}
\right] q_0^2 \rho^{*}(0, R_4-R_1, l, i, q_0) \rho^{}(0, R_2-R_3, j, k, q_0) \\
= \frac{e^2}{N\Omega\epsilon_0} e^{iq_0(R_1-R_3)} \left[ 
   \frac{N\Omega q_{sz} }{2 \pi^2} q_0^2 - 1
\right] \rho^{*}(0, R_4-R_1, l, i, q_0) \rho^{}(0, R_2-R_3, j, k, q_0)\\
$$

### 2D

层间距为 $2R$

$$
f(q, G) = 1 - e^{-(q+G)_{\parallel}R} \cos(G_{z}R)
$$

$$
\pi q_{sz}^{2} = \frac{(2\pi)^2}{N A} \Rightarrow q_{sz} = \left( \frac{4\pi}{NA} \right)^{\frac{1}{2}} \\
$$

二维圆积分

$$
v(k=0) = \frac{NA}{(2\pi)^2} \int_{k<q_{sz}} k \mathrm{d}k \mathrm{d}\theta \frac{1}{k^2} \left[ 1 - e^{-kR} \right] \\
= \frac{NA}{2\pi} \int_{0}^{q_{sz}} \mathrm{d}k \frac{1}{k} \left[ 1 - e^{-kR} \right] \\
\approx \frac{NA}{2\pi} \int_{0}^{q_{sz}} \mathrm{d}k R \\
= \frac{N\Omega q_{sz}}{4\pi} \\
$$

三维柱积分

$$
v(k=0) = \frac{NA}{(2\pi)^2} \frac{2R}{2\pi} \int_{k_{\parallel}<q_{sz}} k_{\parallel} \mathrm{d}k_{\parallel} \mathrm{d}\theta \int_{-\frac{\pi}{2R}}^{\frac{\pi}{2R}} \mathrm{d}k_z \frac{1}{k^2} \left[ 1 + e^{-k_{\parallel}R} \frac{k_z}{k_{\parallel}} \sin(k_zR) - e^{-k_{\parallel}R} \cos(k_zR) \right] \\
= \frac{N\Omega}{(2\pi)^2} \int_{0}^{q_{sz}} k_{\parallel} \mathrm{d}k_{\parallel} \int_{-\frac{\pi}{2R}}^{\frac{\pi}{2R}} \mathrm{d}k_z \frac{1}{k^2} \left[ 1 + e^{-k_{\parallel}R} \frac{k_z}{k_{\parallel}} \sin(k_zR) - e^{-k_{\parallel}R} \cos(k_zR) \right] \\
$$

$$
\int_{0}^{q_{sz}} k_{\parallel} \mathrm{d}k_{\parallel} \int_{-\frac{\pi}{2R}}^{\frac{\pi}{2R}} \mathrm{d}k_z \frac{1}{k^2} 
= 2q_{sz} \arctan\frac{\pi}{2R q_{sz}} + \frac{\pi}{2R} \ln(\frac{4R^2 q_{sz}^2}{\pi^2}+1) \\
\int_{0}^{q_{sz}} k_{\parallel} \mathrm{d}k_{\parallel} \int_{-\frac{\pi}{2R}}^{\frac{\pi}{2R}} \mathrm{d}k_z \frac{1}{k^2} e^{-k_{\parallel}R} \frac{k_z}{k_{\parallel}} \sin(k_zR) 
= \frac{2}{R} \int_{0}^{q_{sz}R} e^{-y} \mathrm{d}y \int_{0}^{\frac{\pi}{2}} \mathrm{d}x \frac{x \sin x}{x^2+y^2} \\
-\int_{0}^{q_{sz}} k_{\parallel} \mathrm{d}k_{\parallel} \int_{-\frac{\pi}{2R}}^{\frac{\pi}{2R}} \mathrm{d}k_z \frac{1}{k^2} e^{-k_{\parallel}R} \cos(k_zR) 
= - \frac{2}{R} \int_{0}^{q_{sz}R} y e^{-y} \mathrm{d}y \int_{0}^{\frac{\pi}{2}} \mathrm{d}x \frac{\cos x}{x^2 + y^2} \\
v(k=0) = \frac{N\Omega}{4\pi R} \left\{2 \alpha \arctan \frac{1}{2\alpha} + \frac{1}{2} \ln(4 \alpha^2+1) + \frac{2}{\pi} \int_{0}^{\pi \alpha} \mathrm{d}y \int_{0}^{\frac{\pi}{2}} \mathrm{d}x \frac{e^{-y}}{x^2+y^2} \left[ x \sin x - y \cos x \right] \right\}_{\alpha = \frac{q_{sz}R}{\pi}} \\
\approx \frac{N\Omega}{4\pi R} \left\{2 \alpha \arctan \frac{1}{2\alpha} + \frac{1}{2} \ln(4 \alpha^2+1) + 0.164619 e^{-3.16133a} - 0.164619 \right\}_{\alpha = \frac{q_{sz}R}{\pi}} \\
$$

$$
\frac{1}{|q_0|^2} \rightarrow \frac{N\Omega q_{sz} }{2 \pi^2}
$$

$$
head = \frac{e^2}{N\Omega\epsilon_0} e^{iq_0(R_1-R_3)} 
\rho^{*}(0, R_4-R_1, l, i, q_0) \rho^{}(0, R_2-R_3, j, k, q_0) \\
head\_correction = \frac{e^2}{N\Omega\epsilon_0} e^{iq_0(R_1-R_3)} \left[ 
   \frac{N\Omega q_{sz} }{2 \pi^2} - \frac{1}{q_0^2}
\right] q_0^2 \rho^{*}(0, R_4-R_1, l, i, q_0) \rho^{}(0, R_2-R_3, j, k, q_0) \\
= \frac{e^2}{N\Omega\epsilon_0} e^{iq_0(R_1-R_3)} \left[ 
   \frac{N\Omega q_{sz} }{2 \pi^2} q_0^2 - 1
\right] \rho^{*}(0, R_4-R_1, l, i, q_0) \rho^{}(0, R_2-R_3, j, k, q_0)\\
$$

