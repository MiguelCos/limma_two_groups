# `limma`-based DE analysis for protein expression data between two groups

The script would take one expression data file and an annotation file and execute differential expression analysis, to evaluate which proteins are affected between two groups.

## How to use:

###  Input requirements:

- 1 File named `input_limma.txt`. It should be tab separated and should contain normalized protein expression data. The first column should be named `ID` and should contain protein IDs (i.e. Uniprot or Symbol). Every other column should correspond to an experimental sample / individual, and should contain the corresponding protein expression values. These column should have a distinctive name.

- 1 File named `annotation.txt` that would map the names of your experimental samples / individuals to each experimental group. The first column should be called `Sample_ID` and should contain exactly the same names/sample codes given in the columns of the expression data file. The second column should be named `Group` and should contain the group that corresponds to each individual / sample.

### Step-by-step:

1. Download this repository in your local computer and initialize it as an R Project in RStudio.
2. Place the two input files in the `Data` folder.
3. Open the `multigroup_limma.R` script. 
4. Answer the questions on line `6` and line `10` and `14` of the script.
5. Define the paired contrasts that you are interested to evaluate, by modifiying `contrast.matrix` input at line `60`. (Check below how to set up your contrast matrix)
6. Click `Source` in the top right corner of the script.

The script will extract the experimental design information from the annotation file and match it to the columns in the expression data file. Then it will perform the DE analysis using `limma` to detect which proteins are differentially expressed between the twop groups. Using the first group in your annotation file as a baseline.


### Modifying your contrast matrix (line `74`)

The contrast matrix would define which groups should be compared against each other (the direction of the contrast).

In line `74` we define the contrast matrix for our comparisons such as this example:

```
contrast.matrix <- makeContrasts(Tumor-normal, levels=design)
```

In this example, we are testing for 2 groups (Tumor vs normal).

By defining `Tumor-normal` in the contrast matrix, the fold-changes would be positive for proteins increased in Tumor.

### Output  

The script will generate an `Output` folder containing: 

- An HTML summary report showing:
  - The summary of the grouping variable (number of samples per group)
  - The number of proteins significantly affected between the groups.
  - A set of boxplots showing the mean values of protein abundance for the top protein hits (those with lowest p-values) for each group.
  - Volcano plots for the comparison.
  - A tabular file with the output statistics after `limma` analysis, including the p-values after testing for each protein.


