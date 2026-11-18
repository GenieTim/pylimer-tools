Normal Mode Analysis
--------------------

pylimer-tools provides a normal mode analysis (NMA) implementation to predict the loss and storage modulus of polymer networks 
as published in :cite:t:`gusev_molecular_2024`. 
This is useful for understanding the viscoelastic properties of the network.
The corresponding class is :class:`~pylimer_tools_cpp.NormalModeAnalyzer`.

See :doc:`/auto_examples/normal_mode_analysis/plot_kpm_vs_exact_pdms` for a large-scale comparison of the Kernel Polynomial Method (KPM)
against full eigenvalue decomposition on realistic PDMS networks generated via Monte Carlo with
on-demand entanglement sampling.
