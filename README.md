# SIC-Based_Slotted-ALOHA

## Overview
In this repository, we provide MATLAB implementations corresponding to our paper on Spatio-Temporal SIC-Based Slotted ALOHA with Multi-Antenna Reception.

In this work, we analyze the sum rate of a spatio-temporal successive interference cancellation (SIC)-based two-device slotted ALOHA system with real-time feedback and multi-antenna reception over fading channels.
In the decoding process, spatio-temporal SIC is applied jointly across time slots and antenna elements.

The probability density functions of fading channels are approximated using a mixture gamma distribution.
Based on this approximation, the SIC-based slotted ALOHA with real-time feedback is modeled as a Markov process, from which we derive its exact analytical expression for throughput.

In the theoretical analysis, the maximum sum rate is achieved by simultaneously optimizing the transmission probability and coding rate of devices.
The numerical results illustrate the performance gains achieved by employing multiple antennas in the slotted ALOHA scheme over general fading channels.

## A quick Demo
Please run the ```simulation/Sim_2user_Multi_anntenna_SA_MG.m```, ```simulation/Sim_MultiAnt_OptR.m```,  ```analysis/optimization.m``` and ```analysis/Ana_Multi-antenna_SA_MG.m```.
To run the programs in this folder, the following packages need to be installed:

* **matlab2tikz** (for exporting figures to Overleaf)
* **Optimization Toolbox** (for optimization tasks)
* **Statistics and Machine Learning Toolbox**
* **Symbolic Math Toolbox**

From the MATLAB top menu, go to the "Home" tab, click on "Add-Ons", and select "Get Add-Ons".
Then, search for the desired package in the search bar and install it.

## Contributing to the project
Any pull requests or issues are welcome.

## Citations
Please consider citing our paper in your publications if the project helps your research. BibTeX reference is as follows.
```
@ARTICLE{Takahashi2025Sum,
  author={Takahashi, Yuhei and 
          Fukui, Daiki and Song, Guanghui and Kimura, Tomotaka and Liu, Zilong and Jun Cheng}, 
  title={Spatio-Temporal {SIC}-Based Slotted {ALOHA}
with Multi-Antenna Reception}, 
   url={https://github.com/yuhei-takahashi/Multi-Antenna_SA/blob/main/SpatioTemporal_SIC_SA.pdf}
}
```
