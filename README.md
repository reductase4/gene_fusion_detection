Detection pipeline for TNL-type fusion genes in green plants.

Version: v1.0

Author: Jiang Qian


**<left>Description</left>**
---------------
This project aims to detect gene fusions of TNL-type resistance genes. We modified the <a href="https://bitbucket.org/yaanlpc/rgaugury/src/master/">RGAugury</a> project for identification of resistance gene analogs and their domain-related genes. We designed a novel detection pipeline for TNL-type gene fusion referring to <a href="https://github.com/zhangcj2022/GriffinDetector">GriffinDetector</a>. 

Advantages of the pipeline:
1. The pipeline is feasible to all green plants. We provided two basic databases of candidate fused and parental genes covering major lineages of green plants, and they could be flexibly updated by additional data from customized species beneficial for detection. 
2. Considering the conserved domains of proteins, we limited detection of gene fusion events to the two candidate databases, making parental filter not only based on protein sequence similarity. 
3. The concept of out-group extended to all lineages diverged before the lineage of focused species, which makes the filter of gene fission more reliable.

Li, P., Quan, X., Jia, G., Xiao, J., Cloutier, S., & You, F. M. (2016). RGAugury: a pipeline for genome-wide prediction of resistance gene analogs (RGAs) in plants. BMC genomics, 17(1), 852. https://doi.org/10.1186/s12864-016-3197-x

Zhou, Y., Zhang, C., Zhang, L., Ye, Q., Liu, N., Wang, M., Long, G., Fan, W., Long, M., & Wing, R. A. (2022). Gene fusion as an important mechanism to generate new genes in the genus Oryza. Genome biology, 23(1), 130. https://doi.org/10.1186/s13059-022-02696-w


**<left>Contents</left>**
---------------
<a href="#Installation">1. Installation</a>
<a href="#Quick_start">2. Quick start</a>
<a href="#Workflow">3. Workflow</a>

<a name="Installation">1. Installation</a>

Install all softwares and set up for environment variables according to <a href="https://bitbucket.org/yaanlpc/rgaugury/wiki/Home">RGAugury</a>. Make sure RGAugury correctly working.

Set the  environment variables to the path of RGAugury_pipeline in scripts folder of this project.


<a name="Quick_start">2. Quick start</a>

The main script is gene_fusion_pipeline.py. 

usage: gene_fusion_pipeline.py [-h] -d DIR -s SPECIES -o ORDER [-e EXPRESSION]


BASIC OPTIONS:

**-h, --help** \<verbose help of all options>

**-d, --dir**  \<dir of genome files>

**-s, --species** \<target species> Focused species list.

**-o, --order** \<species order> Order of focused species and pre-stored 45 species.

**-e, --expression** \<expression dir> Any type of gene expression(Optional).

<a name="Workflow">3. Workflow</a>

**Step1. RGAugury**

1) run RGAugury
   
   _run_RGAugury.sh_

   _RGAugury_pipeline:_ This is a modified version of RGAugury

2) Statistics of results from RGAugury
   
   _extract_RGAugury_results.sh_
   

**Step2. extract protein sequences**

extract TNL/NL/TX/TN peps

_extract_target_info:_
Contains all scripts used to extract proteins and domains identified from RGAugury.

Note: 

_extract_target_info.py_ for all; 

_extract_target_info_v2.0.py_ for LRRfilter version

**Step3. detection of fusion genes**

_blastp_filter_two_steps.py_

**Step4. extract expression data of fused and parental genes (Optional)**

_extract_gene_expression.py_


**Test data**
---------
with expression data

    python scripts/gene_fusion_pipeline.py -d test_data/genomes -s test_data/species.list -o test_data/species_order.txt -e test_data/trans_data

without expression data

    python scripts/gene_fusion_pipeline.py -d test_data/genomes -s test_data/species.list -o test_data/species_order.txt


##END
