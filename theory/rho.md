## regulation $nk \rightarrow ik$

$$
\left| nk \right\rangle = \sum_{G} C_{Gn}^{k} \left| kG \right\rangle
$$

$$
\left| ik \right\rangle = \sum_{n} U_{ni}^{k} \left| nk \right\rangle \\
\left| ik \right\rangle = \sum_{mn} U_{ni}^{k} U_{mn}^{k(dis)} \left| mk \right\rangle \\
$$

$$
C_{Gi}^{k} = \sum_{mn} \left[ U_{ni}^{k} U_{mn}^{k(dis)} \right] C_{Gm}^{k} \\
\left| ik \right\rangle = \sum_{G} C_{Gi}^{k} \left| kG \right\rangle 
$$

## $\rho_{k} (k, G, i, j, q)$

$$
\rho_{k} (k, G, i, j, q) 
= \frac{1}{|q+G|N} \left\langle ik \left| e^{-i(q+G)r} \right| jk + q \right\rangle 
$$

$$
\rho_{k} (k, G, i, j, -q) 
= \frac{1}{|-q+G|N} \left\langle ik \left| e^{-i(-q+G)r} \right| jk - q \right\rangle 
= \frac{1}{|q-G|N} \left\langle i(k-q)+q \left| e^{i(q-G)r} \right| j(k - q) \right\rangle 
= \left[ \frac{1}{|q-G|N} \left\langle j(k - q) \left| e^{-i(q-G)r} \right| i(k-q)+q \right\rangle \right]^{*}
= \rho_{k} (k-q, -G, j, i, q)^{*}
$$

## $\rho_{R} (G, R, i, j, q)$

$$
\rho_{R} (G, R, i, j, q) 
= \sum_{k} e^{ikR} \rho_{k} (k, G, i, j, q) 
$$

$$
\rho_{R} (G, R, i, j, -q) 
= \sum_{k} e^{ikR} \rho_{k} (k, G, i, j, -q) 
= \sum_{k} e^{ikR} \rho_{k} (k-q, -G, j, i, q)^{*}
= \left[ \sum_{k} e^{-ikR} \rho_{k} (k-q, -G, j, i, q) \right]^{*}
?= \left[ e^{-iqR} \sum_{k} e^{-ikR} \rho_{k} (k, -G, j, i, q) \right]^{*}
= e^{iqR} \rho_{R} (-G, -R, j, i, q)^{*} 
$$

<!-- $$
\rho_{R} (0, 0, i, i, q) 
= \sum_{k} \rho_{k} (k, 0, i, i, q) 
= \sum_{k} \rho_{k} (k, 0, i, i, q) 
= \frac{1}{|q|N} \sum_{k} \left\langle ik \left| e^{-iqr} \right| ik + q \right\rangle 
$$ -->

if the shift of MP grid is $k_{s} = \frac{1}{2} [0(1), 0(1), 0(1)]$ then

$$
\rho_{R} (G, R+R_{N}, i, j, q) 
= \sum_{k} e^{ik(R+R_{N})} \rho_{k} (k, G, i, j, q) 
= e^{i k_{s} R_{N}} \rho_{R} (G, R, i, j, q) 
$$

## wispin

$$
\rho_{k} (k, G, i, j, q) 
= \sum_{\sigma} \frac{1}{|q+G|N} \int_{N\Omega} \mathrm{d}r \psi_{ik}^{*}(r,\sigma) \psi_{jk+q}^{}(r,\sigma) e^{-i(q+G)r} 
$$

## calculation

### FFT

$$
\rho_{k} (k, G, i, j, q) 
= \frac{1}{|q+G|N} \int_{N\Omega} \mathrm{d}r \psi_{ik}^{*}(r) e^{-i(q+G)r} \psi_{jk+q}^{}(r)
= \frac{1}{|q+G|N} \int_{N\Omega} \mathrm{d}r u_{ik}^{*}(r) e^{-iGr} u_{jk+q}^{}(r)
= \frac{1}{|q+G|} \int_{\Omega} \mathrm{d}r u_{ik}^{*}(r) e^{-iGr} u_{jk+q}^{}(r)
$$

$$
\psi_{ik+q-G_{wind}} = e^{i(k+q-G_{wind})r} u_{ik+q-G_{wind}}(r)
= \psi_{ik+q} = e^{i(k+q)r} u_{ik+q}(r) \Rightarrow
u_{ik+q}(r) = e^{-iG_{wind}r} u_{ik+q-G_{wind}}(r)
$$

$$
\rho_{k} (k, G, i, j, q) 
= \frac{1}{|q+G|} \int_{\Omega} \mathrm{d}r u_{ik}^{*}(r) u_{jk+q-G_{wind}}(r) e^{-iG''r} 
,\; G'' = G+G_{wind} \\
$$

On a FFT grid

$$
u_{ik}^{*}(r_g) = \frac{N_g}{\sqrt{N\Omega}} \mathrm{IFFT}\{C_{i}^{k}(G)\} 
= \frac{1}{\sqrt{N\Omega}} \mathrm{BFFT}\{C_{i}^{k}(G)\} 
= \frac{1}{\sqrt{N\Omega}} f_{ik}(r_g) \\
\rho_{k} (k, G, i, j, q) 
= \frac{1}{|q+G|} \frac{\Omega}{N_g} \sum_{grid} u_{ik}^{*}(r_g) u_{jk+q-G_{wind}}^{}(r_g) e^{-iG''r_g} 
= \frac{1}{|q+G|} \frac{1}{NN_g} \sum_{grid} f_{ik}^{*}(r_g) f_{jk+q-G_{wind}}^{}(r_g) e^{-iG''r_g} 
$$

$$
= \frac{1}{|q+G|} \frac{1}{NN_g} \mathrm{FFT}[f_{ik}^{*} * f_{jk+q-G_{wind}}](G'') \\
$$

FFT方法更快，与倒空间系数卷积方法耗时不是一个量级。

### 倒空间系数卷积

$$
\rho_{k} (k, G, i, j, q) 
= \frac{1}{|q+G|N} \left\langle ik \left| e^{-i(q+G)r} \right| jk + q \right\rangle 
= \frac{1}{|q+G|N} \sum_{G_1G_2} C_{G_1i}^{k*} C_{G_2j}^{k+q} \delta_{G_1+G, G_2}
= \frac{1}{|q+G|N} \sum_{G'} C_{G'i}^{k*} C_{G'+Gj}^{k+q} 
$$

$$
\psi_{ik+q-G_{wind}} = \sum_{G} C_{Gi}^{k+q-G_{wind}} \left| k+q-G_{wind}+G \right\rangle 
= \psi_{ik+q} = \sum_{G} C_{Gi}^{k+q} \left| k+q+G \right\rangle \Rightarrow
C_{G+G_{wind}i}^{k+q-G_{wind}} = C_{Gi}^{k+q}
$$

$$
\rho_{k} (k, G, i, j, q) 
= \frac{1}{|q+G|N} \sum_{G'} C_{G'i}^{k*} C_{G'+G+G_{wind}j}^{k+q-G_{wind}} 
= \frac{1}{|q+G|N} \sum_{G'} C_{G'i}^{k*} C_{G'+G''j}^{k+q-G_{wind}} 
,\; G'' = G+G_{wind}
$$
