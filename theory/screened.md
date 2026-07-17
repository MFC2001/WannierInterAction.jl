In a system with size of $N\Omega$ and periodic boundary condition, under plane basis, we have screened interaction

$$
W(r_1,r_2,\omega) = \frac{e^2}{N\Omega\epsilon_0} \sum_{qG_1G_2} e^{i(q+G_1)r_1} \frac{\epsilon_{G_1G_2}^{-1}(q,\omega)}{|q+G_2|^2} e^{-i(q+G_2)r_2}
$$

symmetrize $G_1$ and $G_2$

$$
\tilde{\epsilon}_{G_1G_2}^{-1}(q,\omega) = \epsilon_{G_1G_2}^{-1}(q,\omega) \frac{|q+G_1|}{|q+G_2|}
$$

$$
W(r_1,r_2,\omega) = \frac{e^2}{N\Omega\epsilon_0} \sum_{qG_1G_2} e^{i(q+G_1)r_1} \frac{\tilde{\epsilon}_{G_1G_2}^{-1}(q,\omega)}{|q+G_1||q+G_2|} e^{-i(q+G_2)r_2}
$$

Generally, we have

$$
W(r_1,r_2,\omega) = \frac{1}{N} \sum_{qG_1G_2} e^{i(q+G_1)r_1} 
\sqrt{V(q, G_1)} \tilde{\epsilon}_{G_1G_2}^{-1}(q,\omega) \sqrt{V(q, G_2)} e^{-i(q+G_2)r_2}
$$

## screened interaction

$$
\left\langle iR_1;jR_2 \left| W(\omega) \right| kR_3;lR_4 \right\rangle 
= \frac{1}{N} \sum_{q} e^{iq(R_1-R_3)} \sum_{G_1G_2} 
\tilde{\epsilon}_{G_1G_2}^{-1}(q,\omega) 
\rho_{V}^{*}(G_1, R_4-R_1, l, i, q) \rho_{V}^{}(G_2, R_2-R_3, j, k, q) \\
$$

For $q=0$ and $G = 0$

$$
\lim_{q \rightarrow 0} \tilde{\epsilon}_{00}^{-1}(q,\omega) V(q, 0)
\frac{1}{N} \sum_{k_1} e^{ik_1(R_1-R_4)} \left\langle u_{ik_1 + q} \middle| u_{lk_1} \right\rangle 
\frac{1}{N} \sum_{k_2} e^{ik_2(R_2-R_3)} \left\langle u_{jk_2} \middle| u_{kk_2 + q} \right\rangle \\
$$

### head term

$$
\lim_{q \rightarrow 0} \tilde{\epsilon}_{00}^{-1}(q,\omega) V(q, 0) \delta_{il} \delta_{jk} \delta_{R_1R_4} \delta_{R_2R_3}
$$

$$
\lim_{q \rightarrow 0} \frac{\epsilon_{00}^{-1}(q)}{|q|^2}
= \lim_{q \rightarrow 0} \frac{\epsilon_{00}^{-1}(0)}{|q|^2} 
+ \frac{1}{2} \sum_{\alpha\beta} \partial_{\alpha} \partial_{\beta} \epsilon_{00}^{-1}(0) \frac{q_{\alpha}q_{\beta}}{|q|^2}
$$

For low dimension case, 目前的做法参考BerkeleyGW，使用模型介电函数在miniFBZ中进行数值积分，基本算是考虑了0阶项。
但是三维的0阶项目前不会算，先忽略，或许可以参考wannier90的有限差分方法，再取方向平均。

### multipole correction

For 3D, $\lim_{q \rightarrow 0} \tilde{\epsilon}_{00}^{-1}(q,\omega)$ 趋于一个常数；
for low dimension, $\lim_{q \rightarrow 0} \tilde{\epsilon}_{00}^{-1}(q,\omega)$ 趋于1。

Multipole correction for 3D is

$$
\tilde{\epsilon}_{00}^{-1}(0,\omega) \frac{e^2}{\Omega\epsilon_0} \left[ 
\frac{1}{3} \sum_{\alpha} A_{\alpha}^{*}(R_4 - R_1, l, i) A_{\alpha}^{}(R_2 - R_3, j, k)
+ \delta_{il} \delta_{R_1R_4} \frac{1}{6} C^{}(R_2-R_3, j, k)
+ \delta_{jk} \delta_{R_2R_3} \frac{1}{6} C^{*}(R_4-R_1, l, i) \right] \\
$$

取某个方向上的结果计算起来较为复杂，如果需要多方向平均，可以直接对介电常数做平均。

### wing term

$$
wing\_left = \lim_{q \rightarrow 0} \sum_{G \neq 0} 
\tilde{\epsilon}_{G0}^{-1}(q,\omega) \rho_{V}^{*}(G, R_4-R_1, l, i, q)
\sqrt{V(q, 0)} \frac{1}{N} \sum_{k_2} e^{ik_2(R_2-R_3)} \left\langle u_{jk_2} \middle| u_{kk_2 + q} \right\rangle \\
= \lim_{q \rightarrow 0} wl(\omega, R_4 - R_1, l, i, q)
\sqrt{V(q, 0)} \frac{1}{N} \sum_{k_2} e^{ik_2(R_2-R_3)} 
\left[ \delta_{jk} + \sum_{\alpha} \left\langle u_{jk_2} \middle| \partial_{\alpha} u_{kk_2} \right\rangle q_{\alpha} \right] \\
= \delta_{jk} \delta_{R_2R_3} \lim_{q \rightarrow 0} wl(\omega, R_4 - R_1, l, i, q) \sqrt{V(q, 0)} 
+ \lim_{q \rightarrow 0} wl(\omega, R_4 - R_1, l, i, q) 
\sum_{\alpha} A_{\alpha}(R_2 - R_3, j, k) q_{\alpha} \sqrt{V(q, 0)} \\
= \delta_{jk} \delta_{R_2R_3} \lim_{q \rightarrow 0} wl(\omega, R_4 - R_1, l, i, \hat{q}) \sqrt{V(q, 0)} 
+ \delta_{jk} \delta_{R_2R_3} \lim_{q \rightarrow 0} \sum_{\alpha} \partial_{\alpha} wl(\omega, R_4 - R_1, l, i, q) q_{\alpha} \sqrt{V(q, 0)} 
+ \lim_{q \rightarrow 0} wl(\omega, R_4 - R_1, l, i, \hat{q}) 
\sum_{\alpha} A_{\alpha}(R_2 - R_3, j, k) q_{\alpha} \sqrt{V(q, 0)} \\
wl(\omega, R_4 - R_1, l, i, q) 
= \sum_{G \neq 0} \tilde{\epsilon}_{G0}^{-1}(q,\omega) \rho_{V}^{*}(G, R_4-R_1, l, i, q)
$$

类似有

$$
wing\_right = \lim_{q \rightarrow 0} \sum_{G \neq 0} 
\tilde{\epsilon}_{0G}^{-1}(q,\omega) 
\sqrt{V(q, 0)} \frac{1}{N} \sum_{k_1} e^{ik_1(R_1-R_4)} \left\langle u_{ik_1 + q} \middle| u_{lk_1} \right\rangle
\rho_{V}^{}(G, R_2-R_3, j, k, q) \\
= \lim_{q \rightarrow 0} \sqrt{V(q, 0)} \frac{1}{N} \sum_{k_1} e^{ik_1(R_1-R_4)} 
\left[\delta_{il} + \sum_{\alpha} \left\langle \partial_{\alpha} u_{ik_1} \middle| u_{lk_1} \right\rangle q_{\alpha} \right]
wr(\omega, R_2 - R_3, j, k, q) \\
= \delta_{il} \delta_{R_1R_4} \lim_{q \rightarrow 0} wr(\omega, R_2 - R_3, j, k, q) \sqrt{V(q, 0)} 
+ \lim_{q \rightarrow 0} wr(\omega, R_2 - R_3, j, k, q) 
\sum_{\alpha} A_{\alpha}^{*}(R_4 - R_1, l, i) q_{\alpha} \sqrt{V(q, 0)} \\
= \delta_{il} \delta_{R_1R_4} \lim_{q \rightarrow 0} wr(\omega, R_2 - R_3, j, k, \hat{q}) \sqrt{V(q, 0)} 
+ \delta_{il} \delta_{R_1R_4} \lim_{q \rightarrow 0} \sum_{\alpha} \partial_{\alpha} wr(\omega, R_2 - R_3, j, k, q) q_{\alpha} \sqrt{V(q, 0)} 
+ \lim_{q \rightarrow 0} wr(\omega, R_2 - R_3, j, k, \hat{q}) 
\sum_{\alpha} A_{\alpha}^{*}(R_4 - R_1, l, i) q_{\alpha} \sqrt{V(q, 0)} \\
wr(\omega, R_2 - R_3, j, k, q) 
= \sum_{G \neq 0} \tilde{\epsilon}_{0G}^{-1}(q,\omega) \rho_{V}^{}(G, R_2-R_3, j, k, q)
$$

在 $q\rightarrow0$ 时, 从相反的方向趋于0介电函数的翼项反号，于是方向平均消去第一项发散项。
并且与多极修正类似，只有3D的翼项存在第二三项的非零贡献。
代入三维库伦势得到

$$
wing\_left 
= \sqrt{\frac{e^2}{\Omega\epsilon_0}} \left[ 
\delta_{jk} \delta_{R_2R_3} \sum_{\alpha} \partial_{\alpha} wl(\omega, R_4 - R_1, l, i, q) \hat{q}_{\alpha}
+ wl(\omega, R_4 - R_1, l, i, \hat{q}) \sum_{\alpha} A_{\alpha}(R_2 - R_3, j, k) \hat{q}_{\alpha} \right] \\
$$

$$
wing\_right 
= \sqrt{\frac{e^2}{\Omega\epsilon_0}} \left[ 
\delta_{il} \delta_{R_1R_4} \sum_{\alpha} \partial_{\alpha} wr(\omega, R_2 - R_3, j, k, q) \hat{q}_{\alpha} 
+ wr(\omega, R_2 - R_3, j, k, \hat{q}) \sum_{\alpha} A_{\alpha}^{*}(R_4 - R_1, l, i) \hat{q}_{\alpha} \right]\\
$$


第一项不会算，目前只计算第二项。

$$
wing\_left = \sqrt{\frac{e^2}{\Omega\epsilon_0}} wl(\omega, R_4 - R_1, l, i, \hat{q}) 
\sum_{\alpha} A_{\alpha}(R_2 - R_3, j, k) \hat{q}_{\alpha} \\
wl(\omega, R_4 - R_1, l, i, \hat{q}) =
\sum_{G \neq 0} \tilde{\epsilon}_{G0}^{-1}(\hat{q},\omega) \rho_{V}^{*}(G, R_4-R_1, l, i, 0)
$$

$$
wing\_right = \sqrt{\frac{e^2}{\Omega\epsilon_0}} wr(\omega, R_2 - R_3, j, k, \hat{q}) 
\sum_{\alpha} A_{\alpha}^{*}(R_4 - R_1, l, i) \hat{q}_{\alpha} \\
wr(\omega, R_2 - R_3, j, k, \hat{q}) 
= \sum_{G \neq 0} \tilde{\epsilon}_{0G}^{-1}(\hat{q},\omega) \rho_{V}^{}(G, R_2-R_3, j, k, 0)
$$

如果单独对wing项取方向平均

$$
wing\_left = \sqrt{\frac{e^2}{\Omega\epsilon_0}} wl(\omega, R_4 - R_1, l, i, \hat{q}) 
\frac{1}{3} \sum_{\alpha} A_{\alpha}(R_2 - R_3, j, k) \\
$$

$$
wing\_right = \sqrt{\frac{e^2}{\Omega\epsilon_0}} wr(\omega, R_2 - R_3, j, k, \hat{q}) 
\frac{1}{3} \sum_{\alpha} A_{\alpha}^{*}(R_4 - R_1, l, i) \\
$$

其实，3D的介电函数的wing可以写作

$$
wing = \sum_{\alpha} w_{\alpha} q_{\alpha}
$$

所以方向平均的结果为

$$
\frac{1}{3} \sum_{\alpha} w_{\alpha} A_{\alpha}
$$

但是这需要提取出计算介电函数用到的偶极矩，太麻烦。
但是一个折中方案是：输入三个方向的介电函数，头项直接平均计算；翼项不取方向平均，内部对三个方向做求和。


## Hermitian

使用yambo+磁性体系计算的结果中，直接相互作用总是带有较大的虚部。
目前并未发现对激子结果有何影响。

对虚部误差的影响的因素：
QE的noinv参数：默认应用了时间反演对称性，如果体系没有时间反演对称性，则应手动设置该参数为false；
yambo的库伦势的RandGvec：本julia包仅支持替换头项为有限积分值，yambo中对应的应设置该参数为1；
（但H-MLG的测试表明仅可以轻微减小虚部）
数值计算原因：可能就是因为磁性的存在放大了数值误差从而导致厄密性的偏移。

yambo计算介电函数的一个特点，不同的q点使用相同的Glist，而不是按照截断能逐q点设置，这一定会引入误差。
可能是这个设计导致了磁性体系结果中过大的虚部。
