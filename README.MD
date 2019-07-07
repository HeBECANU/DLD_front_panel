# Delay Line Detector Front Panel
***Bryce Henson, [Jacob A. Ross](https://github.com/GroundhogState),Roman Khakimov, Rob Dall***   
Matlab based GUI for He* BEC experiment basic data visualization and analysis




## Install
``` 
git clone --recurse-submodules -j8 https://github.com/brycehenson/DLD_front_panel.git
```
then to update 
```
git submodule update --init --recursive --remote --merge
```



| ![The Gui](/pics/gui.png "Fig1") | 
|:--:| 
 **Figure1**- The Gui  |
 
 | ![Two Dimensional Histogram](/pics/2d_hist.png "Fig2") | 
|:--:| 
 **Figure2**- Two dimensional histogram of some test data.   |
 
 | ![Count Rate](/pics/1d_hist.png "Fig3") | 
|:--:| 
 **Figure3**- A one dimensional hisogram of the counts rate as a functiuon of time.  |
 
 

## To Do
- [ ] Unify import code from Core_BEC_Analysis analysis
- [x] smoothed histograms (as kernel density estimators)
- [x] use figure names
- [x] Fix different blur on YT 2d plot
- [x] Padd FFT
- [ ] fix error when try to import single empty file
- [x] Nice pictures
- [ ] Document the features here
- [x] move functions from file exchange into folders


## Contributions
This project would not have been possible without the many open source tools that it is based on.
* ***Ander Biguri*** [Perceptually uniform colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
* ***Jan*** [FileTime](https://au.mathworks.com/matlabcentral/fileexchange/24671-filetime)
* ***Daniel Eaton***    [sfigure](https://au.mathworks.com/matlabcentral/fileexchange/8919-smart-silent-figure)


