## wave on a kgrid with a smallq shift 

The wannier function in real space, which is the result of finite system, is the same as infinite case.

In infinite case:

$$
\left\langle r \middle| kG \right\rangle = \frac{1}{\sqrt{\Omega}} e^{i(k+G)r} \\
$$

$$
\left| ik \right\rangle = \sum_{R} e^{ikR} \left| iR \right\rangle = \sum_{G} C_{G}^{ik} \left| kG \right\rangle \\
\sum_{R} e^{ikR} \left\langle k'G' \middle| iR \right\rangle 
= \sum_{G} C_{G}^{ik} \left\langle k'G' \middle| kG \right\rangle \\
\sum_{R} e^{ikR} \left\langle k'G' \middle| iR \right\rangle 
= \sum_{R} e^{i(k-k')R} \left\langle k'G' \middle| i0 \right\rangle 
= \frac{(2\pi)^3}{\Omega} \delta(k-k') \left\langle k'G' \middle| i0 \right\rangle \\
\sum_{G} C_{G}^{ik} \left\langle k'G' \middle| kG \right\rangle 
= \frac{(2\pi)^3}{\Omega} \delta(k-k') C_{G'}^{ik'} \\
C_{G}^{ik} = \left\langle kG \middle| i0 \right\rangle 
= \frac{1}{\sqrt{\Omega}} \int \mathrm{d}r e^{-i(k+G)r} \psi_{i0}(r) \\
$$

## calculation

### Interpolation

定义 $C^{i}(k+G) = C^{ik}_{G}$ 后可使用插值方法。

$$
C^{i}(q) = \frac{1}{\sqrt{\Omega}} \int \mathrm{d}r e^{-iqr} \psi_{i0}(r) \\
$$

由于 $\psi_{i0}(r)$ 是局域的，故 $C^{i}(q)$ 一般是光滑的。
插值法极快，与权重求和方法不是一个量级。

注意，这里相当于是对函数 $\psi_{i0}(r)$ 进行的傅里叶变换，故即使不同k点对应了不同的单电子态，$\left| ik \right\rangle$ 也没有任意相位。
也就是说，需保证 $\left| ik \right\rangle$ 的展开系数与最大局域化结果中的完全一致，无需人为设置相位规范。


### 权重求和

由于wannier函数指数衰减，故前述无穷积分限可以限定为 $N\Omega$ 范围：

$$
\tilde{C}_{G}^{ik} = \frac{1}{\sqrt{\Omega}} \int_{N\Omega} \mathrm{d}r e^{-i(k+G)r} \psi_{i0}(r) \\
$$

这里 $\tilde{C}_{G}^{ik}$ 强调 $k$ 可取任意值。

在有限大情况，DFT的结果为

$$
\left\langle r \middle| kG \right\rangle = \frac{1}{\sqrt{N\Omega}} e^{i(k+G)r}
$$

$$
\left| iR \right\rangle = \frac{1}{\sqrt{N}} \sum_{k} e^{-ikR} \left| ik \right\rangle 
= \frac{1}{\sqrt{N}} \sum_{k} e^{-ikR} \sum_{G} C_{G}^{ik} \left| kG \right\rangle \\
\psi_{i0}(r) = \frac{1}{N\sqrt{\Omega}} \sum_{kG} C_{G}^{ik} e^{i(k+G)r}
$$

注意 $C_{G}^{ik}$ 仅定义在MP网格上的 $k$。

$$
\tilde{C}_{G}^{ik} = \frac{1}{N\Omega} \int_{N\Omega} \mathrm{d}r e^{-i(k+G)r} \sum_{k'G'} C_{G'}^{ik'} e^{i(k'+G')r} 
= \sum_{k'G'} C_{G'}^{ik'} \frac{1}{N\Omega} \int_{N\Omega} \mathrm{d}r e^{i(k'+G'-k-G)r} 
= \sum_{k'} \sum_{G'} P(k'+G'-k-G) C_{G'}^{ik'} \\
$$

$$
P(k) = \frac{1}{N\Omega} \int_{N\Omega} \mathrm{d}r e^{ikr} 
= \frac{1}{N} \int_{N} \mathrm{d}x\mathrm{d}y\mathrm{d}z e^{i2\pi(k_xx+k_yy+k_zz)} 
= \frac{1}{N_x} \int_{-\frac{N_x}{2}}^{\frac{N_x}{2}} \mathrm{d}x e^{i2\pi k_xx} ...
= \frac{e^{i\pi k_xN_x} - e^{-i\pi k_xN_x}}{i2\pi k_xN_x} ...
= \frac{sin(\pi k_xN_x)}{\pi k_xN_x} ... \\
= sinc(k_xN_x) sinc(k_yN_y) sinc(k_zN_z)
$$

积分限的更改会导致一个依赖于 $G$ 的相位因子，是非全局的，故不可随意设置积分限，其应代表真实的积分区域。
这里通常假定wannier函数局域在0元胞附近，故积分限如此设置。

For a small q

$$
C_{G}^{ik+q} 
= \sum_{k'G'} C_{G'}^{ik'} P(k'-k-q+G'-G) 
= \sum_{k'G'} C_{G'}^{ik'} P(G, G', k', k) \\
$$

real calculation

对 $\Delta G$ 进行截断，

$$
|sinc((\Delta G - q)N)| < \epsilon
$$

### FFT

$$
\mathrm{IFFT}: u_{ik}(r) = \frac{1}{\sqrt{\Omega}} \sum_{G} C_{G}^{ik} e^{iGr} \\
\psi_{i0}(r) = \frac{1}{N} \sum_{k} e^{ikr} u_{ik}(r)
$$

事实上， $\psi_{i0}(r)$ 的周期为 $N\Omega$，或至少需要包含邻近几个元胞，这应该会导致极大的内存占用，所以FFT方法不适用。

$$
\psi'_{ik}(r) = e^{-ikr} \psi_{i0}(r) \\
\mathrm{FFT}: \tilde{C}_{G}^{ik} = \frac{1}{\sqrt{\Omega}} \int_{N\Omega} \mathrm{d}r e^{-iGr} \psi'_{ik}(r) \\
$$

## result

代码中实现了插值方案，但是实际计算表明，是否使用插值后的wave_q0，仅会导致极小的误差，故默认方案为不使用，特殊情况用户可能会需要，例如q0不够小的情况。

