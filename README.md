# MRIquality
MRIquality is a collection of tools for routine quality assurance (QA) and post-acquisition MRI data quality control (QC). It currently includes mainly tools for EPI data (including fMRI and DWI acquisitions).

## EPI tools

### Quality assurance

Tools for EPI quality assurance, including SNR and stability estimates.

### Quality control

Quality control tools for EPI images (including fMRI and diffusion imaging).

#### Sequential check

Sequential display of EPI series to detect spikes, artefacts and head movements. Suitable for fMRI time series, resting-state,
diffusion-weighted images.

#### Spike check (and correct)

Automated procedure to detect spikes and correct for them when possible, applying a simple correction strategy to reduce the spurious variance introduced by spikes in an fMRI time series.
