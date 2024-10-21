  ---
# urOMT (unbalanced regularized optimal mass transport)
  ---

# Introduction

## Mathematical Formulation
This project works on an unbalanced version of regularized optimal mass transport (urOMT) problem by adding a new relative source variable and its indicator function into the formulation. Specaifically, we deal with a problem as follows:

Given the initial mass distribution function $\rho_0(x)\geqslant0$ and the final one $\rho_1(x)\geqslant0$ defined on a bounded region $\Omega\subseteq\mathbb{R}^3$, one solves

$$\underset{\rho,v,r}{\text{min}}\quad \int_0^T\int_{\Omega}\big(\left\lVert v(t,x)\right\rVert^2\rho(t,x) + \alpha \chi(t,x)r(t,x)^2\rho(t,x)\big)dx dt $$
subject to

$$\frac{\partial\rho}{\partial t} + \nabla\cdot(\rho v) = \sigma\Delta\rho + \chi\rho r, $$

$$\rho(0,x) = \rho_0(x), \quad\rho(T,x) = \rho_1(x)$$

where a temporal dimension $t\in[0,T]$ is added to the transport process. In the above expression, $\rho(t,x)$ is the dynamic density function; $v(t,x)$ is the velocity field defining the flow from $\rho_0$ to $\rho_1$; $r(t,x)$ is the relative source variable; $\chi(t,x)$ is the given indicator function of $r(t,x)$ which takes values either 0 or 1 to constrain $r$ to a certain spatial and temporal location; constant $\sigma>0$ is the diffusion coefficient; $\alpha>0$ is the weighting parameter of the source term in the cost functional.

## Previous Work
This work can be considered as the extension of the regularized optimal mass transport (rOMT) problem with code available at https://github.com/xinan-nancy-chen/rOMT and https://github.com/xinan-nancy-chen/rOMT_spdup and with papers available at

> -- <cite>[Visualizing fluid flows via regularized optimal mass transport with applications to neuroscience][1]</cite>,

> -- <cite>[Cerebral amyloid angiopathy is associated with glymphatic transport reduction and time-delayed solute drainage along the neck arteries][2]</cite>.

## Publications Using This Code
1. The numerical method of the current urOMT problem is described and demonstrated with two samples in 

> -- <cite>[Unbalanced regularized optimal mass transport with applications to fluid flows in the brain][3]</cite>

and please cite this paper if you would like to use this code.

2. The first pratical application of the urOMT model is to the glymphatic system measured by rat brain DCE-MRI with cerebral amyloid angiopathy or Alzheimer's disease, which is now in press at iScience with title "Divergent brain solute clearance in rat models of cerebral amyloid angiopathy and Alzheimer’s disease". In this work, we utilized the urOMT model to reveal different fluid clearance patterns in the aforementioned two neuro-diseases.

3. The second practical application of the urOMT mode is to human breast cancer DCE-MRI to visualize and quantify the fluid flows within tumors and evaluate the theraputic responses. Preliminary results were reported at SPIE Medcial Imaging 2024 Conference:

> -- <cite>[Characterizing fluid flows in breast tumor DCE-MRI studies using unbalanced regularized optimal mass transport methods][4]</cite>,

AAPM 2023 Conference and ISMRM 2024 Conference. The latest manuscript has been sent to Magnetic Resonance in Medicine (Magn Reson Med).

[1]: https://link.springer.com/article/10.1007/s10915-023-02337-9
[2]: https://www.nature.com/articles/s43587-022-00181-4
[3]: https://www.nature.com/articles/s41598-023-50874-y
[4]: https://www.spiedigitallibrary.org/conference-proceedings-of-spie/12926/129261T/Characterizing-fluid-flows-in-breast-tumor-DCE-MRI-studies-using/10.1117/12.3005382.short

## Contact
Contact Xinan Chen at chenx7@mskcc.org for questions.




# Sample cases for demonstration

In this section, we use three datasets to briefly show how our method can quantify and visualize the dynamic fluid flows:

(A). **Synthetic geometric data based on Gaussian sphere.**

*Related publication: <cite>[Unbalanced regularized optimal mass transport with applications to fluid flows in the brain][3]</cite>*

(B). **Rat brain DCE-MRI in the glymphatic system.**

*Related publications: <cite>[Unbalanced regularized optimal mass transport with applications to fluid flows in the brain][3]</cite> and "Divergent brain solute clearance in rat models of cerebral amyloid angiopathy and Alzheimer’s disease"*

(C). **Breast cancer DCE-MRI.**

*Related publications: <cite>[Characterizing fluid flows in breast tumor DCE-MRI studies using unbalanced regularized optimal mass transport methods][4]</cite> and "An optimal mass transport model for the analysis of DCE-MRI and its application to breast cancer treatment response"*

For each demonatration, a driver script was created specifically with name ``driver_***.m``. One can run this driver script only and the results will pop up automatically with default parameters.

## (A) Gaussian Spheres

We created five successive Gaussian-based spheres of size $50\times50\times50$ to test our urOMT algorihtm on. One can go to script ``create_gaussian.m`` to check how we created the synthetic data if interested. The five input images are plotted in 3D rendering as follows.

<p align="center">
<img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/images_1_tind_1.png" width="170" /><img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/images_2_tind_2.png" width="170" /><img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/images_3_tind_3.png" width="170" /><img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/images_4_tind_4.png" width="170" /><img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/images_5_tind_5.png" width="170" /><img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/scalebar.png" width="40" />
</p>

To run urOMT algorithm and its Eulerian and Lagrangian post-processing, one can run script ``driver_gauss.m`` with default parameters, which takes about 3 hours and 20 minutes on a cluster with 3 CPUs. The visualized results will pop up automatically.

Here we plot the returned Eulerian results, the <em>time-averaged Eulerian speed map</em> (left) and the <em>time-averaged Eulerian relative-source map</em> (right) as follows.

<p align="center">
<img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/EULA_set001_031623/Speed/gauss1_EulAveSpeed_E01_05.png" width="300" /> <img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/EULA_set001_031623/Speed/gauss1_EulAveR_E01_05.png" width="300" />
</p>

We also plot the returned Lagrangian results, the <em>pathlines</em> (top-left), the <em>velocity flux vectors</em> (top-right), the <em>speed-lines</em> (bottom-left) and the <em>Péclet-lines</em> (bottom-right) as follows.

<p align="center">
<img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/LPPA_set001_031623/Pathlines/gauss1_LagPathlines_E01_05.png" width="300" /> <img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/LPPA_set001_031623/Pathlines/gauss1_LagFluxVector_E01_05.png" width="300" />
</p>
<p align="center">
<img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/LPPA_set001_031623/Pathlines/gauss1_LagSpdlines_E01_05.png" width="300" />
<img src="Gauss/gauss1/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_beta_5000_alpha_9000_smooth0_dtri1_rreinit1_pcg60/LPPA_set001_031623/Pathlines/gauss1_LagPelines_E01_05.png" width="300" />
</p>

## (B) Rat Brain DCE-MRI

We also test out method on real Dynamic Contrast Enhanced MRI (DCE-MRI) data which is from a 3-month-old healthy rat brain. The follows is the sequence of the 15 input images shown in 3D rendering from right-lateral view (colorbar = 'jet', limits = [0,300]).

<p align="center">
<img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_7_tind_1.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_9_tind_2.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_11_tind_3.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_13_tind_4.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_15_tind_5.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_17_tind_6.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_19_tind_7.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_21_tind_8.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_23_tind_9.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_25_tind_10.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_27_tind_11.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_29_tind_12.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_31_tind_13.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_33_tind_14.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/images_35_tind_15.png" width="170" />
</p>

One can use script ``driver_RatBrain.m`` with default parameters to run urOMT algorithm and its post-processings, which takes about 9 hours on a cluster with 3 CPUs. The visualized results will pop up automatically.

The Eulerian results are as follows, the <em>time-averaged Eulerian speed maps</em> (top row) and <em>time-averaged Eulerian speed maps</em> (bottom row).

<p align="center">
<img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveSpeed_E07_11.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveSpeed_E11_15.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveSpeed_E15_19.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveSpeed_E19_23.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveSpeed_E23_27.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveSpeed_E27_31.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveSpeed_E31_35.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/scalebar1.png" width="35" />
</p>
<p align="center">
<img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveR_E07_11.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveR_E11_15.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveR_E15_19.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveR_E19_23.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveR_E23_27.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveR_E27_31.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/C1217_EulAveR_E31_35.png" width="170" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_031623/Speed/scalebar2.png" width="45" />
</p>

The Lagrangian results are as follows, the <em>pathlines</em> (top-left), the <em>velocity flux vectors</em> (top-right), the <em>speed-lines</em> (bottom-left) and the <em>Péclet-lines</em> (bottom-right).

<p align="center">
<img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/LPPA_set001_031623/Pathlines/C1217_LagPathlines_E07_35.png" width="350" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/LPPA_set001_031623/Pathlines/C1217_LagFluxVector_E07_35.png" width="350" />
</p>
<p align="center">
<img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/LPPA_set001_031623/Pathlines/C1217_LagSpdlines_E07_35.png" width="350" /><img src="RatBrainsCAA3M/C1217/diff_2e3_tj_2_dt_0.4_nt_10_ti_7_tf_33_beta_50_alpha_10000_smooth1_dtri1_rreinit0_pcg60/LPPA_set001_031623/Pathlines/C1217_LagPelines_E07_35.png" width="350" />
</p>

Note: One can alternatively use softeware VisIt 3.3.1 (open source) to visualize the Lagrangian results.

| ![pipeline.png](pipeline.png) | 
|:--:| 
| *Pipeline of urOMT applied to the Glymphatic System in Rat Brains* |


## (C) Breast Cancer DCE-MRI

To apply urOMT to cancer DCE-MRI data, we demonstrate our method on a breast cancer case which is kept anonymous. One can use script ``driver_BC.m`` with default parameters to run urOMT algorithm and its post-processings, which takes about 1.5 hours on a cluster with 3 CPUs. The visualized results will pop up automatically.

Below is the 3D rendering of the breast tumor; the orange region is the region of interest to be analyzed by urOMT.

<p align="center">
<img src="BreastTumor/BC46_01/tumor3Drendering.png" width="1000" />
</p>

The first row of the following figure is the sequence of the 12 input images (concentration of contrast agent) shown in 2D slice (colorbar = 'grayscale'). The four below rows are returned Eulerian dynamic metrics.

<p align="center">
<img src="BreastTumor/BC46_01/diff_2e3_tj_2_dt_0.4_nt_10_ti_5_tf_25_beta_1000_alpha_30000_smooth1_dtri1_rreinit0_pcg60/EULA_set001_102124/Speed/BC46_01_Eul_dynamic.png" width="1000" />
</p>

The Lagrangian results are as follows, the <em>pathlines</em> (left), the <em>velocity flux vectors</em> (middle) and <em>speed-lines</em> (right).

<p align="center">
<img src="BreastTumor/BC46_01/diff_2e3_tj_2_dt_0.4_nt_10_ti_5_tf_25_beta_1000_alpha_30000_smooth1_dtri1_rreinit0_pcg60/LPPA_set001_102124/Pathlines/BC46_01_LagPathlines_E05_27.png" width="350" /><img src="BreastTumor/BC46_01/diff_2e3_tj_2_dt_0.4_nt_10_ti_5_tf_25_beta_1000_alpha_30000_smooth1_dtri1_rreinit0_pcg60/LPPA_set001_102124/Pathlines/BC46_01_LagFluxVector_E05_27.png" width="350" /><img src="BreastTumor/BC46_01/diff_2e3_tj_2_dt_0.4_nt_10_ti_5_tf_25_beta_1000_alpha_30000_smooth1_dtri1_rreinit0_pcg60/LPPA_set001_102124/Pathlines/BC46_01_LagSpdlines_E05_27.png" width="350" />
</p>

Note: One can alternatively use softeware VisIt 3.3.1 (open source) to visualize the Lagrangian results.


