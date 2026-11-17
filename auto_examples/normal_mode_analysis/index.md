<a id="sphx-glr-auto-examples-normal-mode-analysis"></a>

# Normal Mode Analysis

pylimer-tools provides a normal mode analysis (NMA) implementation to predict the loss and storage modulus of polymer networks
as published in Gusev and Bernhard [[GB24](../../acknowledgements.md#id7)].
This is useful for understanding the viscoelastic properties of the network.
The corresponding class is [`NormalModeAnalyzer`](../../api/pylimer_tools_cpp.NormalModeAnalyzer.md#pylimer_tools_cpp.NormalModeAnalyzer).

See [KPM vs. Full Eigenvalue Normal Mode Analysis for PDMS Networks](plot_kpm_vs_exact_pdms.md) for a large-scale comparison of the Kernel Polynomial Method (KPM)
against full eigenvalue decomposition on realistic PDMS networks generated via Monte Carlo with
on-demand entanglement sampling.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the NormalModeAnalyzer in pylimer-tools to predict the loss and storage modulus of polymer networks using exact eigenvalue solvers or the Kernel Polynomial Method (KPM).">![](auto_examples/normal_mode_analysis/images/thumb/sphx_glr_plot_normal_mode_analysis_thumb.png)

[Normal Mode Analysis (NMA)](plot_normal_mode_analysis.md)

  <div class="sphx-glr-thumbnail-title">Normal Mode Analysis (NMA)</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Normal Mode Analysis (NMA) in NormalModeAnalyzer predicts the linear viscoelastic storage modulus G&#x27;(\\omega) and loss modulus G&#x27;&#x27;(\\omega) of entangled polymer networks from their bead-spring connectivity graph (gusev_molecular_2024). Computing the full eigenvalue spectrum scales as \\mathcal{O}(N_{\\text{atoms}}^3), which quickly becomes prohibitive for realistic representative volume elements (RVEs). The Kernel Polynomial Method (KPM) approximates the spectral density using Chebyshev moments and stochastic trace estimation in \\mathcal{O}(N_{\\text{atoms}}) time.">![](auto_examples/normal_mode_analysis/images/thumb/sphx_glr_plot_kpm_vs_exact_pdms_thumb.png)

[KPM vs. Full Eigenvalue Normal Mode Analysis for PDMS Networks](plot_kpm_vs_exact_pdms.md)

  <div class="sphx-glr-thumbnail-title">KPM vs. Full Eigenvalue Normal Mode Analysis for PDMS Networks</div>
</div>
<!-- thumbnail-parent-div-close --></div>
