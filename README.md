# Estimating the Micro-Indel Mutation Rate in *Plasmodium falciparum*

This repository contains information from my S3 Master's thesis project at the University of Montpellier, titled:\
**"Estimating the micro-indel mutation rate in *Plasmodium falciparum* using genomes from mutation accumulation experiments."**

*P. falciparum* has a unique genome that is eighty percent AT rich, providing adequate opportunities for mutations. The underlying mutations aid in the selection of drug-resistant strains in an environment where drugs are prevalent. I hypothesized that Asian strains have a higher rate of micro-indel mutations as compared to the African strains of *P.falciparum*. To test the hypothesis, I used data from Claessens et al. (2014) who generated six clone trees of P. falciparum strains that belong to four geographically distinct regions. Hamilton et al (2016) used the same data to calculate SNP mutation rate in the six strains. I extended the experiment to calculate the micro-indel mutation rate, which could now be appropriately calculated after the publication of twelve newly assembled PacBio reference genomes for *P. falciparum* strains, which also include the strains that were used for clone tree generation. I created a GATK-based pipeline to detect, micro-indels in non-3D7 strains.

The project was carried out at LPHI (Laboratory of Host-Pathogen Interactions) under the supervision of Dr. Antoine Classens (PI) and Marc Antoine Guery, from September 2020 to February 2021.

------------------------------------------------------------------------

## Table of Contents

-   [Repository Structure](#repository-structure)
    -   [Pipeline](#pipeline)
    -   [Test Scripts](#test-scripts)
    -   [Internship Results](#internship-results)
-   [Additional Resources](#additional-resources)
-   [Important Note](#important-note)

------------------------------------------------------------------------

## Repository Structure {#repository-structure}

This repository contains three main folders:

### Pipeline {#pipeline}

-   Contains the main scripts used in the analysis pipeline.
-   The schematic below illustrates the major steps in the pipeline:

![](pipeline_schematic.png)

For more details about the pipeline and each step, see: `Internship_results/S3_report_Aakanksha_Singh.pdf`

------------------------------------------------------------------------

### Test Scripts {#test-scripts}

-   Includes test scripts and their documentation.\
-   Provides insight into how I coded during the project.\
-   Explore this folder to see the scripts and the corresponding documentation.

------------------------------------------------------------------------

### Internship Results {#internship-results}

-   Contains the full Master's thesis report.\
-   Includes the presentation used to summarize the project for quick reference.\
-   Useful for reviewing both my writing and presentation skills.

------------------------------------------------------------------------

## Additional Resources {#additional-resources}

For an overview of the scripts and tools used in the analysis, see: `summary_of_the_update_meetings.pptx`

This presentation was used to update my supervisors during project meetings. It also served as a way for me to keep track of what I was testing and developing. Looking back, it's a valuable resource for understanding how the project was organized.

------------------------------------------------------------------------

## Important Note {#important-note}

The project continued after the internship ended, and a manuscript draft is currently in preparation.\
Because of this, not all results and data generated during the project can be shared in this repository.

The raw files used in the analysis can be found in the following study:

-   Claessens, A., Hamilton, W. L., Kekre, M., Otto, T. D., Faizullabhoy, A., Rayner, J. C., & Kwiatkowski, D. (2014). Generation of Antigenic Diversity in Plasmodium falciparum by Structured Rearrangement of Var Genes During Mitosis. PLoS Genetics, 10(12), e1004812. <https://doi.org/10.1371/journal.pgen.1004812>

-   Hamilton, W. L., Claessens, A., Otto, T. D., Kekre, M., Fairhurst, R. M., Rayner, J. C., & Kwiatkowski, D. (2016). Extreme mutation bias and high AT content inPlasmodium falciparum. Nucleic Acids Research, gkw1259. <https://doi.org/10.1093/nar/gkw1259>
