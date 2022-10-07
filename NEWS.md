## Changes in NanoTube 1.3.7
- In addition to ruv::RUVIII, the RUVSeq::RUVg method can now be used for data
  normalization. Options have also been added to allow tuning of the parameters
  for these methods. More details are provided in the vignette.
- A csv or txt file containing a design matrix can be input as a 'sampleTab'
  in processNanostringData(). This facilitates easier differential expression 
  analysis with more complex models, and an example has been added to the 
  vignette.
- Various other improvements to vignette.

## Changes in NanoTube 1.3.6
- NanoTube can now process zipped and tarred (.zip or .tar) directories, as
  well as gzipped (.gz) RCC files, such as those downloaded from GEO in many
  cases.

## Changes in NanoTube 1.3.5
- Corrected a bug that caused NanoTube not to recognize reporters labeled as 
  "Endogenous1", "Endogenous2", etc. as Endogenous.

## Changes in NanoTube 1.1.2
- When processing a folder containing RCC files, processNanostringData() 
  now ignores filenames not ending in "RCC" (case-insensitive).
- Corrected a bug that in read_cpdb_sourceDBs() and read_cpdb_tab() that
  sometimes led to the data being misread.

## NanoTube 0.99.0
- Pre-release version
