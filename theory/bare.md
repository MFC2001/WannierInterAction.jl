In a system with size of $N\Omega$ and periodic boundary condition, under plane basis, we have interaction

$$
V(r_1,r_2) = \frac{1}{N} \sum_{qG} e^{i(q+G)r_1} V(q, G) e^{-i(q+G)r_2}
$$
For 3D,

$$
V(q, G) = \frac{e^2}{\Omega\epsilon_0} \frac{1}{|q+G|^2}
$$

## bare interaction

$$
\left\langle iR_1;jR_2 \left| V \right| kR_3;lR_4 \right\rangle 
= \frac{1}{N^2} \sum_{\{k\}} e^{ik_1R_1} e^{ik_2R_2} e^{-ik_3R_3} e^{-ik_4R_4} 
\left\langle ik_1;jk_2 \left| V \right| kk_3;lk_4 \right\rangle \\
= \frac{1}{N^2} \sum_{\{k\}} e^{ik_1R_1} e^{ik_2R_2} e^{-ik_3R_3} e^{-ik_4R_4} 
\iint \mathrm{d}r_1 \mathrm{d}r_2 \psi_{ik_1}^{*}(r_1) \psi_{jk_2}^{*}(r_2) V(r_1, r_2) \psi_{kk_3}^{}(r_2) \psi_{lk_4}^{}(r_1) \\
= \frac{1}{N^3} \sum_{qG} V(q, G) \sum_{\{k\}} e^{ik_1R_1} e^{ik_2R_2} e^{-ik_3R_3} e^{-ik_4R_4} 
\left\langle ik_1 \left| e^{i(q+G)r} \right| lk_4 \right\rangle 
\left\langle jk_2 \left| e^{-i(q+G)r} \right| kk_3 \right\rangle \\
$$

$$
\left\langle ik_1 \left| e^{i(q+G)r} \right| jk_2 \right\rangle \propto \delta_{q,k_1-k_2}^{2\pi}
$$

$$
\left\langle iR_1;jR_2 \left| V \right| kR_3;lR_4 \right\rangle 
= \frac{1}{N^3} \sum_{qG} V(q, G) \sum_{k_2k_4} e^{i(k_4+q)R_1} e^{ik_2R_2} e^{-i(k_2+q)R_3} e^{-ik_4R_4} 
\left\langle ik_4 + q \left| e^{i(q+G)r} \right| lk_4 \right\rangle 
\left\langle jk_2 \left| e^{-i(q+G)r} \right| kk_2 + q \right\rangle \\
= \frac{1}{N^3} \sum_{qG} V(q, G) e^{iq(R_1-R_3)} 
\sum_{k_1k_2} e^{ik_1(R_1-R_4)} e^{ik_2(R_2-R_3)} 
\left\langle ik_1 + q \left| e^{i(q+G)r} \right| lk_1 \right\rangle 
\left\langle jk_2 \left| e^{-i(q+G)r} \right| kk_2 + q \right\rangle \\
= \frac{1}{N} \sum_{qG} e^{iq(R_1-R_3)} V(q, G)
\frac{1}{N} \sum_{k_1} e^{ik_1(R_1-R_4)} \left\langle ik_1 + q \left| e^{i(q+G)r} \right| lk_1 \right\rangle 
\frac{1}{N} \sum_{k_2} e^{ik_2(R_2-R_3)} \left\langle jk_2 \left| e^{-i(q+G)r} \right| kk_2 + q \right\rangle \\
$$

$$
\rho(k, G, i, j, q) 
= \left\langle ik \left| e^{-i(q+G)r} \right| jk + q \right\rangle \\
\rho(G, R, i, j, q) = \frac{1}{N} \sum_{k} e^{ikR} \rho(k, G, i, j, q) 
$$

$$
\rho_{V}(G, R, i, j, q) = \sqrt{V(q, G)} \rho^{}(G, R, i, j, q) 
$$

$$
\left\langle iR_1;jR_2 \left| V \right| kR_3;lR_4 \right\rangle 
= \frac{1}{N} \sum_{q} e^{iq(R_1-R_3)} 
\sum_{G} \rho_{V}^{*}(G, R_4-R_1, l, i, q) \rho_{V}^{}(G, R_2-R_3, j, k, q) 
$$

不同维度的截断库伦势直接替换 $V(q, G)$ 即可。

### divergent term

For $q=0$ and $G = 0$

$$
\lim_{q \rightarrow 0}
V(q, 0) \rho^{*}(0, R_4-R_1, l, i, q) \rho^{}(0, R_2-R_3, j, k, q) \\
= \lim_{q \rightarrow 0} V(q, 0)
\frac{1}{N} \sum_{k_1} e^{ik_1(R_1-R_4)} \left\langle u_{ik_1 + q} \middle| u_{lk_1} \right\rangle 
\frac{1}{N} \sum_{k_2} e^{ik_2(R_2-R_3)} \left\langle u_{jk_2} \middle| u_{kk_2 + q} \right\rangle \\
= \lim_{q \rightarrow 0} 
\frac{1}{N^2} \sum_{k_1k_2} e^{ik_1(R_1-R_4)} e^{ik_2(R_2-R_3)} V(q, 0) 
\left[\delta_{il} + \sum_{\alpha} \left\langle \partial_{\alpha} u_{ik_1} \middle| u_{lk_1} \right\rangle q_{\alpha} 
+ \frac{1}{2} \sum_{\alpha\beta} 
\left\langle \partial_{\alpha} \partial_{\beta} u_{ik_1} \middle| u_{lk_1} \right\rangle q_{\alpha} q_{\beta} \right] 
\left[\delta_{jk} + \sum_{\alpha} \left\langle u_{jk_2} \middle| \partial_{\alpha} u_{kk_2} \right\rangle q_{\alpha} 
+ \frac{1}{2} \sum_{\alpha\beta} 
\left\langle u_{jk_2} \middle| \partial_{\alpha} \partial_{\beta} u_{kk_2} \right\rangle q_{\alpha} q_{\beta} \right] \\
= \lim_{q \rightarrow 0} \left[
V(q, 0) \delta_{il} \delta_{jk} \delta_{R_1R_4} \delta_{R_2R_3} \right. \\
+ \delta_{il} \delta_{R_1R_4} \frac{1}{N} \sum_{k_2} e^{ik_2(R_2-R_3)} \sum_{\alpha} \left\langle u_{jk_2} \middle| \partial_{\alpha} u_{kk_2} \right\rangle q_{\alpha} V(q, 0) \\
+ \delta_{jk} \delta_{R_2R_3} \frac{1}{N} \sum_{k_1} e^{ik_1(R_1-R_4)} \sum_{\alpha} \left\langle \partial_{\alpha} u_{ik_1} \middle| u_{lk_1} \right\rangle q_{\alpha} V(q, 0) \\
+ \frac{1}{N^2} \sum_{k_1k_2} e^{ik_1(R_1-R_4)} e^{ik_2(R_2-R_3)} \sum_{\alpha\beta} \left\langle \partial_{\alpha} u_{ik_1} \middle| u_{lk_1} \right\rangle \left\langle u_{jk_2} \middle| \partial_{\beta} u_{kk_2} \right\rangle q_{\alpha} q_{\beta} V(q, 0) \\
+ \delta_{il} \delta_{R_1R_4} \frac{1}{N} \sum_{k_2} e^{ik_2(R_2-R_3)} \frac{1}{2} \sum_{\alpha\beta} \left\langle u_{jk_2} \middle| \partial_{\alpha} \partial_{\beta} u_{kk_2} \right\rangle q_{\alpha} q_{\beta} V(q, 0) \\
+ \delta_{jk} \delta_{R_2R_3} \frac{1}{N} \sum_{k_1} e^{ik_1(R_1-R_4)} \frac{1}{2} \sum_{\alpha\beta} \left\langle \partial_{\alpha} \partial_{\beta} u_{ik_1} \middle| u_{lk_1} \right\rangle \left. q_{\alpha} q_{\beta} V(q, 0) \right] \\
$$

不同维度下库伦势的发散行为不同，分别为-2(3D)、-1(2D)、ln(1D)、0(0D)，所以这里只保留至二阶。
$V(q, 0)$ 在 $q=0$ 附近为偶函数，故取方向平均可消去第二三项：

$$
\lim_{q \rightarrow 0} \left[
V(q, 0) \delta_{il} \delta_{jk} \delta_{R_1R_4} \delta_{R_2R_3} \right. \\
+ \frac{1}{N^2} \sum_{k_1k_2} e^{ik_1(R_1-R_4)} e^{ik_2(R_2-R_3)} \sum_{\alpha\beta} \left\langle \partial_{\alpha} u_{ik_1} \middle| u_{lk_1} \right\rangle \left\langle u_{jk_2} \middle| \partial_{\beta} u_{kk_2} \right\rangle q_{\alpha} q_{\beta} V(q, 0) \\
+ \delta_{il} \delta_{R_1R_4} \frac{1}{N} \sum_{k_2} e^{ik_2(R_2-R_3)} \frac{1}{2} \sum_{\alpha\beta} \left\langle u_{jk_2} \middle| \partial_{\alpha} \partial_{\beta} u_{kk_2} \right\rangle q_{\alpha} q_{\beta} V(q, 0) \\
+ \delta_{jk} \delta_{R_2R_3} \frac{1}{N} \sum_{k_1} e^{ik_1(R_1-R_4)} \frac{1}{2} \sum_{\alpha\beta} \left\langle \partial_{\alpha} \partial_{\beta} u_{ik_1} \middle| u_{lk_1} \right\rangle \left. q_{\alpha} q_{\beta} V(q, 0) \right] \\
$$

其中第一项为头项，在任意维度下都存在，通过使用有限积分替换的方式消除发散；后三项为多级修正，只在3D情况下存在。
但是有一点需要注意，这里略掉了系数 $\frac{1}{N}$, 最终导致第一项数值积分 $I_{head}$ 的贡献正比于 $\frac{1}{N^{\frac{1}{3}}}$，后三项正比于 $\frac{1}{N}$ 。

### multipole correction for 3D

代入 $V(q, 0) = \frac{e^2}{\Omega\epsilon_0} \frac{1}{q^2}$ 得到3D情况下的多级修正为

$$
\frac{e^2}{\Omega\epsilon_0} \left[ 
\frac{1}{N^2} \sum_{k_1k_2} e^{ik_1(R_1-R_4)} e^{ik_2(R_2-R_3)} \sum_{\alpha\beta} 
\left\langle \partial_{\alpha} u_{ik_1} \middle| u_{lk_1} \right\rangle 
\left\langle u_{jk_2} \middle| \partial_{\beta} u_{kk_2} \right\rangle 
\hat{q}_{\alpha} \hat{q}_{\beta} \right. \\
+ \delta_{il} \delta_{R_1R_4} \frac{1}{N} \sum_{k_2} e^{ik_2(R_2-R_3)} \frac{1}{2} \sum_{\alpha\beta} 
\left\langle u_{jk_2} \middle| \partial_{\alpha} \partial_{\beta} u_{kk_2} \right\rangle \hat{q}_{\alpha} \hat{q}_{\beta} \\
+ \left. \delta_{jk} \delta_{R_2R_3} \frac{1}{N} \sum_{k_1} e^{ik_1(R_1-R_4)} \frac{1}{2} \sum_{\alpha\beta} 
\left\langle \partial_{\alpha} \partial_{\beta} u_{ik_1} \middle| u_{lk_1} \right\rangle \hat{q}_{\alpha} \hat{q}_{\beta} \right] \\
$$

取方向平均得

$$
\frac{e^2}{\Omega\epsilon_0} \left[ 
\frac{1}{N^2} \sum_{k_1k_2} e^{ik_1(R_1-R_4)} e^{ik_2(R_2-R_3)} \frac{1}{3} \sum_{\alpha} 
\left\langle \partial_{\alpha} u_{ik_1} \middle| u_{lk_1} \right\rangle 
\left\langle u_{jk_2} \middle| \partial_{\alpha} u_{kk_2} \right\rangle \right. \\
+ \delta_{il} \delta_{R_1R_4} \frac{1}{N} \sum_{k_2} e^{ik_2(R_2-R_3)} \frac{1}{6} \sum_{\alpha} 
\left\langle u_{jk_2} \middle| \partial_{\alpha}^{2} u_{kk_2} \right\rangle \\
+ \left. \delta_{jk} \delta_{R_2R_3} \frac{1}{N} \sum_{k_1} e^{ik_1(R_1-R_4)} \frac{1}{6} \sum_{\alpha} 
\left\langle \partial_{\alpha}^{2} u_{ik_1} \middle| u_{lk_1} \right\rangle \right] \\
$$

simplify：

$$
\frac{1}{N} \sum_{k} e^{ikR} \left\langle u_{ik} \middle| \partial_{\alpha} u_{jk} \right\rangle 
= \frac{1}{N} \sum_{k} e^{ikR} A_{\alpha}^{}(k, i, j) 
= A_{\alpha}^{}(R, i, j)
$$

$$
\frac{1}{N} \sum_{k} e^{ikR} \left\langle u_{ik} \middle| \nabla^{2} u_{jk} \right\rangle 
= \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k e^{ikR} \left\langle u_{ik} \middle| \nabla^{2} u_{jk} \right\rangle \\
= \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k e^{ikR} \left[ \nabla \cdot \left\langle u_{ik} \middle| \nabla u_{jk} \right\rangle - \sum_{\alpha} \left\langle \partial_{\alpha} u_{ik} \middle| \partial_{\alpha} u_{jk} \right\rangle \right] \\
= \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k \left[ \nabla \cdot \left( e^{ikR} \left\langle u_{ik} \middle| \nabla u_{jk} \right\rangle \right) - iR \cdot \left( e^{ikR} \left\langle u_{ik} \middle| \nabla u_{jk} \right\rangle \right) \right]
- \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k e^{ikR}\sum_{\alpha} \left\langle \partial_{\alpha} u_{ik} \middle| \partial_{\alpha} u_{jk} \right\rangle \\
= - i\sum_{\alpha} R_{\alpha} \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k e^{ikR} \left\langle u_{ik} \middle| \partial_{\alpha} u_{jk} \right\rangle 
- \sum_{\alpha} \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k e^{ikR} \left\langle \partial_{\alpha} u_{ik} \middle| \partial_{\alpha} u_{jk} \right\rangle \\
= - i\sum_{\alpha} R_{\alpha} \frac{1}{N} \sum_{k} e^{ikR} \left\langle u_{ik} \middle| \partial_{\alpha} u_{jk} \right\rangle 
- \sum_{\alpha} \frac{1}{N} \sum_{k} e^{ikR} \left\langle \partial_{\alpha} u_{ik} \middle| \partial_{\alpha} u_{jk} \right\rangle \\
= - i R \cdot A(R, i, j) - \frac{1}{N} \sum_{k} e^{ikR} B(k, i, j) \\
= - i R \cdot A(R, i, j) - B(R, i, j) \\
= C(R, i, j)
$$

多极修正可以写作

$$
\frac{e^2}{\Omega\epsilon_0} \left[ 
\frac{1}{3} \sum_{\alpha} A_{\alpha}^{*}(R_4 - R_1, l, i) A_{\alpha}^{}(R_2 - R_3, j, k)
+ \delta_{il} \delta_{R_1R_4} \frac{1}{6} C^{}(R_2-R_3, j, k)
+ \delta_{jk} \delta_{R_2R_3} \frac{1}{6} C^{*}(R_4-R_1, l, i) \right] \\
$$

计算方法为参考wannier90计算重叠矩阵的方法，通过搜索shell进行有限差分，涉及到的波函数内积则通过对齐基组直接点积计算，也可以使用FFT处理。
但重叠矩阵并不需要实时计算，.chk文件中有wannier基下的重叠矩阵。

$$
\sum_{b} w_{b} b_{\alpha} b_{\beta} = \delta_{\alpha\beta}
$$

$$
\nabla f(k) = \sum_{b} w_{b} b \left[ f(k+b) - f(k) \right] \\
|\nabla f(k)|^2 = \sum_{b} w_{b} \left[ f(k+b) - f(k) \right]^2
$$

可以把 $f(k+b) = f(k) + \nabla f(k) \cdot b$ 代入进行验证。

$$
A_{}^{}(k, i, j) = \left\langle u_{ik} \middle| \partial_{\alpha} u_{jk} \right\rangle
= \sum_{b} w_{b} b \left[ \left\langle u_{ik} \middle| u_{jk+b} \right\rangle - \delta_{ij} \right]
= \sum_{b} w_{b} b \left[ M_{ij}^{kb} - \delta_{ij} \right]
$$

$$
B_{}^{}(k, i, j) = \sum_{\alpha} \left\langle \partial_{\alpha} u_{ik} \middle| \partial_{\alpha} u_{jk} \right\rangle
= \sum_{b} w_{b} \left[ \left\langle u_{ik+b} \right| - \left\langle u_{ik} \right| \right] \left[ \left| u_{jk+b} \right\rangle - \left| u_{jk} \right\rangle \right]
= \sum_{b} w_{b} \left[ 2\delta_{ij} - \left\langle u_{ik} \middle| u_{jk+b} \right\rangle - \left\langle u_{ik+b} \middle| u_{jk} \right\rangle \right]
= \sum_{b} w_{b} \left[ 2\delta_{ij} - M_{ij}^{kb} - M_{ji}^{kb*} \right]
$$



如果要计算某个方向上的结果，有

$$
\frac{1}{N} \sum_{k} e^{ikR} \left\langle u_{ik} \middle| \partial_{\alpha} \partial_{\beta} u_{jk} \right\rangle 
= \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k e^{ikR} \left\langle u_{ik} \middle| \partial_{\alpha} \partial_{\beta} u_{jk} \right\rangle \\
= \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k e^{ikR} \left[ \partial_{\alpha} \left\langle u_{ik} \middle| \partial_{\beta} u_{jk} \right\rangle - \left\langle \partial_{\alpha} u_{ik} \middle| \partial_{\beta} u_{jk} \right\rangle \right] \\
= \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k \left[ \partial_{\alpha} \left( e^{ikR} \left\langle u_{ik} \middle| \partial_{\beta} u_{jk} \right\rangle \right) - iR_{\alpha} \left( e^{ikR} \left\langle u_{ik} \middle| \partial_{\beta} u_{jk} \right\rangle \right) \right]
- \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k e^{ikR} \left\langle \partial_{\alpha} u_{ik} \middle| \partial_{\beta} u_{jk} \right\rangle \\
= - i R_{\alpha} \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k e^{ikR} \left\langle u_{ik} \middle| \partial_{\beta} u_{jk} \right\rangle 
- \frac{\Omega}{(2\pi)^3} \int \mathrm{d}^{3}k e^{ikR} \left\langle \partial_{\alpha} u_{ik} \middle| \partial_{\beta} u_{jk} \right\rangle \\
= - i R_{\alpha} \frac{1}{N} \sum_{k} e^{ikR} \left\langle u_{ik} \middle| \partial_{\beta} u_{jk} \right\rangle 
- \frac{1}{N} \sum_{k} e^{ikR} \left\langle \partial_{\alpha} u_{ik} \middle| \partial_{\beta} u_{jk} \right\rangle \\
= - i R_{\alpha} A_{\beta}(R, i, j) - \frac{1}{N} \sum_{k} e^{ikR} B_{\alpha\beta}(k, i, j) \\
= - i R_{\alpha} A_{\beta}(R, i, j) - B_{\alpha\beta}(R, i, j) \\
= C_{\alpha\beta}(R, i, j)
$$

$$
\left\langle \partial_{\alpha} u_{ik} \middle| \partial_{\beta} u_{jk} \right\rangle
= \sum_{bb'} w_{b}w_{b'} b_{\alpha} b'_{\beta}
\left[ \left\langle u_{ik+b} \right| - \left\langle u_{ik} \right| \right] \left[ \left| u_{jk+b'} \right\rangle - \left| u_{jk} \right\rangle \right]
$$

这一项的计算需要引入重叠矩阵以外的项。

$$
\left| u_{ik} \right\rangle = \frac{1}{\sqrt{N}} \sum_{R} e^{ik(R-r)} \left| iR \right\rangle \\
\left| \partial_{\alpha} u_{ik} \right\rangle = i\frac{1}{\sqrt{N}} \sum_{R} (R-r)_{\alpha} e^{ik(R-r)} \left| iR \right\rangle \\
$$

$$
\left\langle u_{ik} \middle| \partial_{\alpha} u_{jk} \right\rangle 
= i\frac{1}{N} \sum_{RR'} \left\langle iR \right| e^{-ik(R-r)} (R'-r)_{\alpha} e^{ik(R'-r)} \left| jR' \right\rangle \\
= i\frac{1}{N} \sum_{RR'} e^{ik(R' - R)} \left\langle iR \right| (R'-r)_{\alpha}  \left| jR' \right\rangle \\
= i \sum_{R} e^{ikR} \left\langle i0 \right| (R - r)_{\alpha}  \left| jR \right\rangle 
\approx -i\tau_{i\alpha} \delta_{ij} \\
$$

### test for multipole correction

使用LiF体系进行测试：
1. 镜像修正不收敛，说明U的长程行为极有可能偏离了1/r；
2. 使用相同kgrid和原始数据进行计算，三态激子依旧存在0.2eV左右的能级偏移，说明这个误差与mpc无关。

多极修正有可能是对的，镜像修正的头项处理有些问题，最好把相位因子也展开，可以得到类似于多极修正的结果，二者或许可以相互抵消！
但是，LiF的测试也表明即使完整计入多极修正，结果仍然存在偏移，说明激子能级偏差中mpc至少不是主要贡献。

## notation

### U

$$
U_{i0,jR} = \left\langle i0;jR \left| V \right| jR;i0 \right\rangle 
= \frac{1}{N} \sum_{q} e^{-iqR} 
\sum_{G} \rho_{V}^{*}(G, 0, i, i, q) \rho_{V}^{}(G, 0, j, j, q) 
$$

### ijij

$$
J^{1}_{i0,jR} = \left\langle i0;jR \left| V \right| i0;jR \right\rangle 
= \frac{1}{N} \sum_{q} \sum_{G}
\rho_{V}^{*}(G, R, j, i, q) \rho_{V}^{}(G, R, j, i, q) 
$$

### iijj

$$
J^{2}_{i0,jR} = \left\langle i0;i0 \left| V \right| jR;jR \right\rangle 
= \frac{1}{N} \sum_{q} e^{-iqR} 
\sum_{G} \rho_{V}^{*}(G, R, j, i, q) \rho_{V}^{}(G, -R, i, j, q) 
$$

### pair

$$
(p_1, 0) = (i0, lR_1) \\
(p_2, R) = (kR, jR+R_2)
$$

$$
\left\langle i0;jR+R_2 \left| V \right| kR;lR_1 \right\rangle 
= \left\langle p_1, 0\left| V \right| p_2, R \right\rangle 
= \frac{1}{N} \sum_{q} e^{-iqR} 
\sum_{G} \rho_{V}^{*}(G, R_1, l, i, q) \rho_{V}^{}(G, R_2, j, k, q) 
$$
