# Baseline-subtraction


<p align="center">
    <a href="https://ruhr-uni-bochum.sciebo.de/s/1sMpYb00J5v4j3B" target="_blank">
        <img align="center" width = "200" alt="download" src="/assets/download_logo1.png"/>
    </a>
</p>

## Overview
Subtract baseline for XRD libraries based on asymmetric least squares smoothing. In most case, the baseline problems in instrumental methods are characterized by a smooth
baseline and a superimposed signal that carries the analytical information: a series of peaks that are either all positive or all negative. We combine a smoother
with asymmetric weighting of deviations from the (smooth) trend get an effective
baseline estimator. It is easy to use, fast and keeps the analytical peak signal intact. No prior information about peak shapes or baseline (polynomial) is needed
by the method. The performance is illustrated by simulation and applications to
real data.

### Screenshots

<p align = "center">
  <img align = "center" width = "600" src = "/assets/image1.gif"/>
</p>


## Python dependencies
The list of required python packages is contained in the [requirements.txt](requirements.txt) file in this repository. After install the required dependencies, run `main_baseline_subtraction.py` to execute the program.
