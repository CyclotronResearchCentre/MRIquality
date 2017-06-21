## How to install the toolbox

It is recommended to clone the MRIquality repository and keep your version regularly up-to-date. 
This can be fiddly if the toolbox has to be copied into the SPM/toolbox directory.
We provide a simple way to make the toolbox available while keeping it in a folder outside the SPM folder.
Just follow the instructions below!

### Clone the repository

```git clone https://github.com/CyclotronResearchCentre/MRIquality.git```

### Redirection script

- Copy the ```MRIquality/install/MRIquality``` directory (containing the ```tbx_cfg_mriq_redirect.m``` file) to ```<path-to-your-spm>/toolbox```
- Add the path to the full implementation of the MRIquality toolbox (the one you just cloned) to your Maltab Path.

### Very first test

- Start Matlab and the SPM (fMRI) user interface,      
- start the Batch GUI (Batch button in the SPM Menu),     
- check whether you can access the MRI Quality tools: ```SPM>Tools>MRI Quality```
