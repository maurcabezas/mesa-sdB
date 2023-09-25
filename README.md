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

  - https://github.com/maurcabezas/mesa-sdB/blob/90a03d7e851218a7c5895df9b4cfda7a88517742/inlist_std_sdB/src/run_binary_extras.f#L69-L105

This is how looks our new subroutine, there we declared the new variables in order to use the new equation to calculate the new Jdot (Bobrick,A).

First we define the $\dot{J}_{\beta,min}$;

```math
\dot{J}_{\beta,\text{min}}= \dot{M}_{\text{RG},L1}\frac{M_{\text{RG}}^2}{(M_{\text{RG}}+M_{\text{comp}})^2}\frac{2\pi}{p} a^{2}\left(1-e^2\right)^{\frac{1}{2}}
```

*And, the $\dot{J}_{\beta, max}$;*

$$\dot{J}_{\beta, max}= \left[\left(\frac{M_{{RG}}}{M_{\mathrm{RG}}+M_{\text {comp}}}\right)^{2}+\left(0.500-0.227 \cdot \log _{10} \frac{M_{\mathrm{RG}}}{M_{\mathrm{MS}}}\right)^{2}\right]\cdot\frac{2\pi}{p} a^{2}\left(1-e^2\right)^{\frac{1}{2}}\cdot{\dot{M}_{RG}}$$

**The Cauchy-Schwarz Inequality**
$$\dot{J}_{\beta, max}= \left[\left(\frac{M_{{RG}}}{M_{\mathrm{RG}}+M_{\text {comp}}}\right)^{2}+\left(0.500-0.227 \cdot \log _{10} \frac{M_{\mathrm{RG}}}{M_{\mathrm{MS}}}\right)^{2}\right]\cdot\frac{2\pi}{p} a^{2}\left(1-e^2\right)^{\frac{1}{2}}\cdot{\dot{M}_{RG}}$$

Here, we say how many variables we wanna add to *binary_history*. In our case, for now, we are adding 5 variables (or more later).

  - https://github.com/maurcabezas/mesa-sdB/blob/90a03d7e851218a7c5895df9b4cfda7a88517742/inlist_std_sdB/src/run_binary_extras.f#L143-L186

This is a important part in the *run_bianry_extras.f*. Here we 'calculate' or define the variable which we are passing to the *history_binary*. 
There, *names(n)* are the names of the new variables that we can find in *binary_history.data*. For now, the variable *vals(3)* is the factor between both *jdot_beta*, and the varibales *vals(4)* and *vals(5)*, are the new *jdot_ml* 

2) inlist_std_sdB/inslit_project

  - https://github.com/maurcabezas/mesa-sdB/blob/90a03d7e851218a7c5895df9b4cfda7a88517742/inlist_std_sdB/inlist_project#L15-L16

We must set the 'other' *jdot_* as **.true.**


## NOTES

  - I think is neccesary include the routine of *jdot_min* even if is the default equation, because if we wanna do some change in te *jdot_beta_min* we need to use the pointer **b %** in order to tell mesa that must consider the recalculated *jdot_ml* for any extra internal calculation in the code. 

  - For now I am trying to do it by editing binary/other/binary_jdot.f90, if succes I will updated.

  - Probably, is not neccesary define three times *jdot_alpha* and *jdot_delta*. I am trying ti define them al global variables. 
  
  - Feel free to make any change.
