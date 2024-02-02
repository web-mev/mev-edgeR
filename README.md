# mev-edgeR

This repository contains a WebMeV-compatible tool for performing differential expression analysis using Bioconductor's edgeR tool. 

Note that the default usage assumes you are performing a simple contrast, without consideration for more complex designs. 

The outputs are:
- A tab-delimited file of the differential expression results merged with the normalized counts
- A tab-delimited file of just the normalized counts.

The concatenation of the differential expression results and counts is for convenience with the WebMeV frontend interface since it avoids pulling data from two different files.

---

### To run external of WebMeV:

Either:
- build the Docker image using the contents of the `docker/` folder (e.g. `docker build -t myuser/edgeR:v1 .`) 
- pull the docker image from the GitHub container repository (see https://github.com/web-mev/mev-edgeR/pkgs/container/mev-edger)

To run, change into the directory containing your count matrix you wish to run differential expression on. Then:
```
docker run -it -v $PWD:/work <IMAGE> Rscript /opt/software/edgeR.R \
    /work/<path to raw/integer counts> \
    <base/control condition samples as CSV-string> \
    <experimental condition samples as CSV-string> \
    <base/control condition name> \
    <experimental condition name>
```
Note that we mounted your current directory to `/work` inside the Docker container. Hence, your file is relative to `/work`.

The call to the script assumes the following:
- The input file of expression counts is tab-delimited format and contains only integer entries
- The samples in either the control or experimental groups are given as comma-delimited strings and correspond to the column names contained in the raw count matrix. As an example, if the base/control samples are A, B, and C, specify: `"A,B,C"`. The wrapping quotations are not necessary unless the sample names contain whitespace.
- The condition names are regular strings used to help with naming the output file.
- The output files will be written to the same directory where the input/raw counts file is located.

The choice of samples in each group should be a proper subset of the samples represented in the samples in the count matrix; you are *not* required to subset the count matrix to include only those samples involved in the contrast of interest.