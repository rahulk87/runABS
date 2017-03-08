# runABS
A collection of R scripts and functions to run ABSOLUTE on pipeline mutation summary outputs


## Generating a standard package summary

### Metrics file
1. Check that none of the tumor samples have sequencing statistics below 90

### geneCN
1. Go through each and every sample and check whether *amplifications* or *homozygous deletions* were called correctly and ensure their prsense in the sample *geneCN.txt* and visible in the copy number corresponding plot
  - our methods allow for manual curation, thus, inspect and confirm all regions marked as *homozygous deletions*.  Some of these may need to be changed to a *loss* if the signal lacks strength
2. FACETS can miss calls depending on how the genome is segmented.  These parameters change the sensitivity and size of CNV detected.  If there are regions that should be considered having a CNV please discuss it with the group as to ascertain what parameters should be modified to improve our CNV sensitivity.

### Mutation Summary
1. From the summary excel file, save the **MUTATION_SUMMARY** sheet to a .csv file, we suggest *muts.csv*
2. Move the FACETS *\*.cncf.txt* files from each sample to the same directory as *muts.csv*
3. Run the script *Make_ABS_muts_Input.R* to generate per sample input files for ABSOLUTE
4. Run *Make_SegFromFacets.R* to create the copy number input files for ABSOLUTE
5. Run *RunAbsolute.R*
  1. If you run this script in the same R session as the previous two, *sampleList* object is preserved and will run the ABSOLUTE on a loop
  2. Run **Step 1** and **Step 2** of *RunAbsolute.R*
  3. An `Output/` directory containing *absolute_results/* will be created
  4. Select ABSOLUTE solutions for each case.
    - In `Output/absolute_results/` open the *\*.PP-calls_tab.txt* file.  
    - Create a new first column labeled "solution"
    - Open the sample *.pdf* from *Output/absolute_results/* .  
      + Best solution should have an **SSNV multiploidy peak between 1.0 or 2.0** and **CCF plot should have a green peak around 1.0**
    - If you are undecided between two solutions, pick the lower number solution and consult with a colleague
  5.Once solutions are set run **Step 3** in *RunAbsolute.R*
6. Run *CombineAbsMaf.R*
  1. Change the MAF directory to where the results from ABSOLUTE are present. This is typically *Output/absolute_results/Output/reviewed/SEG_MAF*
  2. Asign to *results.directory* the location of the *muts.csv* file
  3. This is all annotated in the code.  If the hotspots are already annotated you do not need to redo this step.  Currently, **REDO the LOH**
7. Run *Mutation_Figure_Patient.R*
  1. Comments and suggestions throughout the code
  2. Will need to make cosmetic chages to the figures
    - In the CCF heatmap, surround mutations with yellow circles with a yellow box indicating these are clonal.  Those with LOH need to have a white slash through them
    - Make all **samples bold**, *mutations italicized* and everything should be in Arial Font
8. Run *make mutSigReport*
