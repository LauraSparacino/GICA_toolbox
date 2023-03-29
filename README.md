# GICA_toolbox
Granger Causality (GC), Granger Isolation (GI) and Granger Autonomy (GA) for the assessment of bivariate causal and non-causal interactions in linear autoregressive processes.

Concepts of Granger causality (GC) and Granger autonomy (GA) are central to assess the dynamics of coupled physiologic processes.
While causality measures have been already proposed and largely applied in time and frequency domains,
measures quantifying self-dependencies are still limited to the time-domain formulation and lack of a clear spectral representation.
In this toolbox, we embed into the classical linear parametric framework for computing GC from a driver random process X to a target process Y
a measure of Granger Isolation (GI) quantifying the part of the dynamics of Y not originating from X,
and a new spectral measure of GA assessing frequency-specific patterns of self-dependencies in Y.
The framework is formulated in a way such that the full-frequency integration of the spectral GC, GI and GA measures returns the corresponding time-domain measures.
The measures are applied to representative time series of mean arterial pressure and cerebral blood flow velocity obtained in a healthy subject.
