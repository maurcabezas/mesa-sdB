# Welcome to the MESA-sdB repository! 
This is a testing repository for binary/ using the MESA version R10390. The main purpose of this repository is to determine the angular momentum loss during the RLOF in systems where the donor is a red giant.

## Why is this important?

Determining the angular momentum loss during RLOF is important for understanding the evolution of stars and the processes that govern their dynamic. The current and typicall treatment of angular momentum loss during RLOF is made using in MESA fixed mass-loss fractions (Tauris & van den Heuvel,2006), and is not sufficient, where the actual angular momentum loss of the three most common cases is not objectively known (Vos et al. 2020).

![binary_interaction_MESA](https://user-images.githubusercontent.com/8492276/229100181-cf5664e3-430b-48b3-bd63-52c3ce0b2a73.png)

## How we can improve it?

By editing and creating equations in the src/run_binary_extras.f in our working directory inlist_std_sdB/ wehere we expect to make a constrained angular momentum loss model applicable for mass loss on the RGB.

## Edition:

1) lund_std_sdB/src/run_binary_extras:

  - https://github.com/maurcabezas/mesa-sdB/blob/90a03d7e851218a7c5895df9b4cfda7a88517742/inlist_std_sdB/src/run_binary_extras.f#L52

With this line we are able to have the new Jdot which depends of our new parameter e, the new routine is called "jdot_ml_new" and in the history will appear as *other_jdot_ml*.

  - https://github.com/maurcabezas/mesa-sdB/blob/4085193e38366c57d75ffcbbbce968752b8e9f04/inlist_std_sdB/src/run_binary_extras.f#L67-L105

This is how looks our new subroutine, there we declared the new variables in order to use the new equation to calculate the new Jdot (Bobrick,A).

First we define the $\dot{J}_{\beta,min}$;

```math
\dot{J}_{\beta,\text{min}}= \dot{M}_{\text{RG},L1}\frac{M_{\text{RG}}^2}{(M_{\text{RG}}+M_{\text{comp}})^2}\frac{2\pi}{p} a^{2}\left(1-e^2\right)^{\frac{1}{2}}
```

And, the $\dot{J}_{\beta, max}$;

```math
\dot{J}_{\beta, max}= \left[\left(\frac{M_{{RG}}}{M_{\mathrm{RG}}+M_{\text {comp}}}\right)^{2}+\left(0.500-0.227 \cdot \log _{10} \frac{M_{\mathrm{RG}}}{M_{\mathrm{MS}}}\right)^{2}\right]\cdot\frac{2\pi}{p} a^{2}\left(1-e^2\right)^{\frac{1}{2}}\cdot{\dot{M}_{RG}}
```
by using the variables **f1** and **f2** in order to create the new $\dot{J}_{\beta}$;

```math
{J}_{\beta, tot}=(1-\xi) \cdot J_{\beta,min} + \xi\cdot J_{\beta,max}
```
where $0 \leq \xi \leq 1$ be the degree of mixture, represented by the new parameter **e_deg_mix** which will be read from the header. And, the new $\dot{J}_{\beta}$ will be the formula represented by the variable:

  - https://github.com/maurcabezas/mesa-sdB/blob/4085193e38366c57d75ffcbbbce968752b8e9f04/inlist_std_sdB/src/run_binary_extras.f#L102-L103


It can be interesting 'save' the values of $\dot{J}_{\beta, min}$ and $\dot{J}_{\beta, max}$, to do this we add both 'values' in the history by using the subroutines **how_many_extra_binary_history_columns**

  - https://github.com/maurcabezas/mesa-sdB/blob/6649d6c3f94af76fe729e7840eb006f75041a9d0/inlist_std_sdB/src/run_binary_extras.f#L108-L112

And, **data_for_extra_binary_history_columns**

  - https://github.com/maurcabezas/mesa-sdB/blob/6649d6c3f94af76fe729e7840eb006f75041a9d0/inlist_std_sdB/src/run_binary_extras.f#L114-L148

There, *names(n)* are the names of the new variables that we can find in *binary_history.data*. Is important to note we are adding the new equation inside the subroutine **'other_jdot_ml**, with this we are sure the new definition of the angular momentum will affect,for example, the calculation of the period. If we add the equation in the rountine **data_for_extra_binary_history_columns** we will only calculate the values, but it will not taken into account for internal calculations.

2) inlist_std_sdB/inslit_project

  - https://github.com/maurcabezas/mesa-sdB/blob/90a03d7e851218a7c5895df9b4cfda7a88517742/inlist_std_sdB/inlist_project#L16

We must set the 'other' *jdot_* as **.true.**, and choose our parameter $\xi$.

- https://github.com/maurcabezas/mesa-sdB/blob/e0bb62393e96e02bd217981c8ec1046c5e145654/inlist_std_sdB/inlist_project#L15-L16


