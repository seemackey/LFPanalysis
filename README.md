# LFPanalysis

Shell script “ExtractCSDmuaLoopv2” loops through file paths and runs the “EphysExtractFxn.” EphysExtractFxn loads continuous recordings of ephys data, filters it, and epochs it via module_cnt05 and module_trig01. Then the data are present in ExtractCSDmuaLoopv2 workspace, where outliers are discarded (“reject artifacts” and “deleteoutliers”).

Other general purpose functions are here for analyzing oscillations using the Morlet wavelet method, such as snr_wavelet08, rayleigh, circmean, etc. 
written mostly by Peter Lakatos

Recently added functions that extract eye movement data (2/2023) written mostly by Peter 
