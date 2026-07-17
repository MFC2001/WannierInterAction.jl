## $\epsilon$ on irreducible qpoints

由于实空间介电函数是满足晶体对称性的，故可把不可约q点上的介电矩阵转化为所有q点上的介电矩阵。

Given that the system satisfy $\{R|T\}$, then

$$
\epsilon(r_1,r_2, \omega) = \frac{1}{N\Omega} \sum_{qG} e^{i(q+G_1)r_1} \epsilon_{G_1G_2}(q,\omega) e^{-i(q+G_2)r_2} \\
\epsilon(Rr_1+T,Rr_2+T, \omega) = \epsilon(r_1,r_2, \omega)
$$

$$
\epsilon_{G_1G_2}(q,\omega) 
= \frac{1}{N\Omega} \iint \mathrm{d}r_1 \mathrm{d}r_2 e^{-i(q+G_1)r_1} \epsilon(r_1,r_2, \omega) e^{i(q+G_2)r_2} \\
= \frac{1}{N\Omega} \iint \mathrm{d}r_1 \mathrm{d}r_2 
e^{-i(q+G_1)(Rr_1+T)} \epsilon(Rr_1+T,Rr_2+T, \omega) e^{i(q+G_2)(Rr_2+T)} \\
= e^{i(G_2-G_1)T} \frac{1}{N\Omega} \iint \mathrm{d}r_1 \mathrm{d}r_2 
e^{-iR^{T}(q+G_1)r_1} \epsilon(r_1,r_2, \omega) e^{iR^{T}(q+G_2)r_2} \\
= e^{i(G_2-G_1)T} \epsilon_{R^{T}G_1R^{T}G_2}(R^{T}q,\omega) \\
= e^{i(G_2-G_1)T} \epsilon_{R^{T}G_1+G_{wind}, R^{T}G_2+G_{wind}}(R^{T}q-G_{wind},\omega) \\
$$

$$
\epsilon_{G_1G_2}(q,\omega) e^{i(G_1-G_2)T}
= \epsilon_{R^{T}G_1+G_{wind}, R^{T}G_2+G_{wind}}(R^{T}q-G_{wind},\omega) \\
$$

介电矩阵的逆理应满足相同的关系。

旋转操作不更改模长，故有

$$
\epsilon^{-1}_{G_1G_2}(q,\omega) e^{i(G_1-G_2)T} \frac{|q+G_1|}{|q+G_2|}
= \epsilon^{-1}_{R^{T}G_1+G_{wind}, R^{T}G_2+G_{wind}}(R^{T}q-G_{wind},\omega) \frac{|R^{T}(q+G_1)|}{|R^{T}(q+G_2)|} \\
$$

这说明这些变换对对称化的介电函数也是一样的结果。

### calculation in practice

实际计算中，已知q点介电矩阵的Glist，可求出 $R^{T}q$ 的G-list，二者的Glist顺序一一对应，
于是无需调整介电矩阵元位置便可直接得到 $R^{T}q$ 对应的矩阵元（仅包含一相位因子）。

## Deduce

### polarization

polarization function satisfy dyson equation

$$
\chi_{} = \chi_{0} + \chi_{0} v \chi_{} \\
\chi_{}^{-1} = \chi_{0}^{-1} - v \\
$$

under RPA, we have

$$
\chi_{G_1G_2}^{0}(q, \omega) = \frac{2}{\Omega} \sum_{vck} (f_{vk} - f_{ck+q}) 
\frac{\left\langle vk \right| e^{-i(q+G_1)r} \left| ck+q \right\rangle \left\langle ck+q \right| e^{i(q+G_2)r} \left| vk \right\rangle}{\omega - (E_{ck+q} - E_{vk}) + i\eta}
$$

if don't have spin degeneration, then sum all bands with spin:

$$
\chi_{G_1G_2}^{0}(q, \omega) = \frac{1}{\Omega} \sum_{vck} (f_{vk} - f_{ck+q}) 
\frac{\left\langle vk \right| e^{-i(q+G_1)r} \left| ck+q \right\rangle \left\langle ck+q \right| e^{i(q+G_2)r} \left| vk \right\rangle}{\omega - (E_{ck+q} - E_{vk}) + i\eta}
$$

polarization is a hermitian matrix if $\eta = 0$

$$
\chi_{G_1G_2}^{0*}(q, \omega) = \chi_{G_2G_1}^{0}(q, \omega)
$$

但是对于介电函数情况可能更复杂一点，但是一般静态的对称化的介电矩阵是一个厄密矩阵。

For $q \rightarrow 0$

$$
\chi_{00}^{0}(q, \omega) \propto \sum_{\alpha\beta} A_{\alpha\beta} q_{\alpha}q_{\beta} \\
\chi_{G0}^{0}(q, \omega) \propto \sum_{\alpha} B_{\alpha} q_{\alpha} \\
\chi_{0G}^{0}(q, \omega) \propto \sum_{\alpha} C_{\alpha} q_{\alpha} \\
$$

### dielectric function

$$
\epsilon = 1 - v \chi_{0} \\
\epsilon^{-1} = 1 + v \chi_{}
$$

$$
\epsilon_{G_1G_2}(q) = 1 - \frac{\chi^{0}_{G_1G_2}}{|q+G_1|^2} \\
\tilde{\epsilon}_{G_1G_2}(q) = 1 - \frac{\chi^{0}_{G_1G_2}}{|q+G_1||q+G_2|}
$$

For $q \rightarrow 0$, 头项不连续且是偶函数，翼项不连续且是奇函数。
对称化的介电矩阵完全是有限良定义的。

一般来说，逆介电矩阵也满足同样的性质。

如果涉及到截断库伦势，那么在对称化的时候要把截断因子也开方对称化，以及可以参考[Table 2](https://www.sciencedirect.com/science/article/pii/S0010465511003912)得知截断库伦势下介电函数和屏蔽库伦势在 $q \rightarrow 0$ 的极限行为。

### symmetrize dielectric function

$$
\epsilon W = V \\
V_{GG'} = V(q+G) \delta_{GG'} = \frac{\delta_{GG'}}{|q+G|^2} \\
\epsilon_{GG'} = \delta_{GG'} - \frac{\chi_{GG'}}{|q+G|^2}
$$

$$
\tilde{\epsilon}_{GG'} = \delta_{GG'} - \frac{\chi_{GG'}}{|q+G||q+G'|} = \epsilon_{GG'} \frac{|q+G|}{|q+G'|} \\
\delta_{GG'} = \sum_{G''} \epsilon^{-1}_{GG''} \epsilon_{G''G'} 
= \sum_{G''} \epsilon^{-1}_{GG''} \frac{|q+G|}{|q+G''|} \epsilon_{G''G'} \frac{|q+G''|}{|q+G'|} 
= \sum_{G''} \tilde{\epsilon}^{-1}_{GG''} \tilde{\epsilon}_{G''G'} \\
\tilde{\epsilon}^{-1}_{GG'} = \epsilon^{-1}_{GG'} \frac{|q+G|}{|q+G'|}
$$

$$
\sum_{G''} \epsilon_{GG''} W_{G''G'} = V_{GG'} \Rightarrow
\sum_{G''} \tilde{\epsilon}_{GG''} W_{G''G'}|q+G''| = \frac{\delta_{GG'}}{|q+G|} \Rightarrow
\sum_{G''} \delta_{G'''G''} W_{G''G'}|q+G''| = \sum_{G} \tilde{\epsilon}^{-1}_{G'''G} \frac{\delta_{GG'}}{|q+G|} \Rightarrow
W_{GG'} = \frac{\tilde{\epsilon}^{-1}_{GG'}}{|q+G||q+G'|}
$$

$$
\tilde{\epsilon} = \left[\begin{matrix}
    |q+G_0| & 0 & \dots \\
    0 & |q+G_1| & \dots \\
    \vdots & \vdots & \ddots
\end{matrix}\right] 
\epsilon
\left[\begin{matrix}
    \frac{1}{|q+G_0|} & 0 & \dots \\
    0 & \frac{1}{|q+G_1|} & \dots \\
    \vdots & \vdots & \ddots
\end{matrix}\right] 
= V^{-\frac{1}{2}} \epsilon V^{\frac{1}{2}} \\
\epsilon W = V \Rightarrow
\tilde{\epsilon} V^{-\frac{1}{2}} W = V^{\frac{1}{2}} 
\Rightarrow W = V^{\frac{1}{2}} \tilde{\epsilon}^{-1} V^{\frac{1}{2}}
$$

关于这一点，BerkeleyGW的介电矩阵是没有经过任何对称化操作的；
与此相反，yambo输出的介电函数是完全对称化的，包括截断因子也一起对称化。

## RESPACK

For q=0

redefine

$$
\tilde{\epsilon}_{00} = \epsilon_{00} = |q_0|^2 \left[ 1 - \lim_{q \rightarrow 0} \frac{\chi_{00}(q)}{|q|^2} \right] = 0 \\
\tilde{\epsilon}_{0G'} = - \frac{\chi_{0G'}}{|G'|} = \epsilon_{0G'} \frac{|q_0|^2}{|G'|} \\
\tilde{\epsilon}_{G0} = - \frac{\chi_{G0}}{|G|} = \epsilon_{G0} |G| \\
\tilde{\epsilon} = \left[\begin{matrix}
    |q_0|^2 & 0 & \dots \\
    0 & |G_1| & \dots \\
    \vdots & \vdots & \ddots
\end{matrix}\right] 
\epsilon
\left[\begin{matrix}
    1 & 0 & \dots \\
    0 & \frac{1}{|G_1|} & \dots \\
    \vdots & \vdots & \ddots
\end{matrix}\right] 
= S \epsilon T \\
$$

$$
\tilde{\epsilon}^{-1} = T^{-1} \epsilon^{-1} S^{-1} \\
\tilde{\epsilon}^{-1}_{00} = \frac{\epsilon^{-1}_{00}}{|q_0|^2} \rightarrow \infty
$$

也就是说，RESPACK的方法一定有问题。
即使解析地获取了 $\epsilon_{00}$, 由于理论上 $\tilde{\epsilon}^{-1}_{00}$ 是发散的，仍然无法解析地获取 $\epsilon^{-1}_{00}$, 更不用说进一步地使用 $\epsilon^{-1}_{00}$ 进行头项修正了。

$$
S^{-1}T = V \\
\epsilon W = V \Rightarrow
\tilde{\epsilon} T^{-1} W = SV = T \Rightarrow
W = T \tilde{\epsilon}^{-1} T
$$

不过这么变换确实可以消除掉翼项的发散性。
也就是说，倘若使用一个小量q进行计算，该方法可以将所有的发散性集中到头项上，然后再使用小量q的模计算出 $\epsilon^{-1}_{00}$。

但是对于实际计算，既然使用了小量q，那么就无需特殊处理 q=0 的发散性，可以统一地使用第一种变换形式。
也就是说，这种针对 q=0 的变换形式基本没有实用价值。


## inverse

$$
\epsilon = \left[\begin{matrix}
    a & \frac{b}{q_0^2} \\
    c & d 
\end{matrix}\right]
$$

$$
S = d - \frac{c a^{-1} b}{q_0^2} \\
$$

$$
\epsilon^{-1} = \left[\begin{matrix}
    \epsilon_{00}^{-1} + \epsilon_{00}^{-1} \epsilon_{0*} S^{-1} \epsilon_{*0} \epsilon_{00}^{-1} 
    & -\epsilon_{00}^{-1} \epsilon_{0*} S^{-1} \\
    -S^{-1} \epsilon_{*0} \epsilon_{00}^{-1} & S^{-1} 
\end{matrix}\right]
$$
