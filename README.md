# PenicillinFinder
## Tool to report diagnostic ions of drug-haptenated peptides
When benzylpenicillin-haptenated peptides are subjected to mass spectrometry, the haptenation causes distinctive fragment ions at 160.0432, 217.0647 and 335.1066 m/z in the tandem mass spectra. PenicillinFinder aims to help researchers confirm the presence of haptenated peptides among mass spectrometry search results by reporting the presence or absence of these diagnostic ions for each peptide-spectrum match. While this was the primary goal, the program can also be used to search for other ion masses so it may be used to help identify other drug haptenations or biological post-translational modifications that produce suitable diagnostic ions.

PenicillinFinder was written in the C programming language and has been compiled and run on Ubuntu 18.04, and also compiled on Linux for Windows using MinGW. It was designed to take input mzIdentML and mgf files exported from PEAKS Studio XPro with Bruker TimsTOF or SCIEX data. It may also work for these input file formats derived from other search engines but this has not been tested.

#### Compilation examples:
Linux:

```
cc penicillinFinder.c -lm -o penicillinFinder.o 
```
Linux to run on Windows:
```
x86_64-w64-mingw32-gcc penicillinFinder.c -lm -o penicillinFinder.exe
```

#### Usage:
First place the mzIdentML file and its associated mgfs together in a folder. You only need to specify the location of the mzIdentML file, but PenicillinFinder will look for the mgfs in the same location.

Usage example:
```
/path/to/penicillinFinder.o -p /path/to/peptides_1_1_0.mzid [-b] [-i 123.000,456.000] [-o myOutput] &> myLog.log
```

|Required Input:||
|---|---|
|-p|Specify peptide search results (peptides\_1\_1\_0.mzid from PEAKS, in the same folder as the associated mgf files)|

|Optional:||
---|---
-b|Specify mode is Bruker to include ion mobility information
-i|Specify diagnostic ions in comma-separated list (default is '160.0432,217.0647,335.1066')
-o|Specify name for output csv file (default is 'output')
-h|Print help

**Output**: csv file showing details of PSMs with Y/N columns for each diagnostic ion

#### Citation:
If you use this tool, please cite our paper [manuscript in preparation]
