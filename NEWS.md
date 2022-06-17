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
