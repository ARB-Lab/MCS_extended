MCS extended
============

An extended and generalized MCS framework

2020/06/15

Philipp Schneider, Axel von Kamp, Steffen Klamt

Added Features:
---------------

1.  Definition of multiple target regions (**T**·**r** ≤ **t**) and multiple desired
    regions (**D**·**r** ≤ **d**).

2.  Specification of **gene** and **reaction deletion** and **addition**
    candidates and **individual cost factors** for each intervention.

3.  Fast computation of **gene-based MCS** using **GPR associations** and novel
    **compression** techniques.

Software Requirements:
----------------------

1.  MATLAB 2016b® or later

2.  IBM ILOG® CPLEX® 12.7, 12.8, 12.9 or 12.10 (Make sure that your CPLEX® and MATLAB® versions are compatible)

3.  CellNetAnalyzer2020.1

4.  Set up *CellNetAnalyzer* to access the CPLEX-Matlab-API (replace the CPLEX-default paths in startcna.m with the paths to your CPLEX installation as described by *CellNetAnalyzer* manual)
    
Getting Started:
----------------------
1. Download this project to your computer (see release page https://github.com/ARB-Lab/MCS_extended/releases) and extract all files.
2. Start MATLAB.
3. Navigate to the main directory of your *CellNetAnalyzer* Toolbox. By executing **pwd** you can verify that the CNA main directory is now your current woking directory.
4. Start *CellNetAnalyzer* by executing **startcna** or **startcna(1)**  (silent start).
5. *Either* navigate to the path of the script in this project that you want to execute *or* add all folders from this project to your MATLAB® path by a right click on the project folder and "Add to path" -> "Selected folders and subfolders"
6. Run script

Script Files:
-------------

1. cofeeding_example/**cofeeding_example.m**

   Demonstrates how CNAMCSEnumerator2 can be used to compute strain designs
   for single subtrate feeding and substrate co-feeding from a single MCS setup.

2. GPR_example/**GPR-example.m** 

   Demonstrates how CNAgeneMCSEnumerator2 computes gene-based MCS using
   GPR associations and advanced compression routines.

3. e_coli/**benchmark.m**

   Computes and characterizes gene-MCS for the production of 2,3-BDO in *E. coli* in
   a core (ECC2) and a genome-scale (iML1515) setup. The scipt benchmarks the runtime reduction 
   achieved through applying the novel GPR rule compression. The computed MCS for the genome-scale
   setup (scenario 1) also serve as a reference for scripts 4-6.
   Results are saved to a .mat file in the working directory and the MCS
   characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because 
   of a return statement after MCS computation. To activate it, remove this return statement.
   By changin variables in the script to e.g. model='full', options.compression_GPR or options.preproc_compression, 
   different setups and compression routines can be used.

   model='ECC2': Computation from an E. coli core model with ca. 500 reactions
   model='full': Computation from a genome-scale E. coli model (iML1515) with ca. 2700 reactions

4. e_coli/**desired2.m**

   Computes and characterizes genome-scale gene-MCS for the production of
   2,3-BDO in *E.coli* from a similar setup as in scenario 1 (benchmark.m). A second
   desired region is added to the setup to ensure that strain designs support higher ATP
   maintanance rates. The results of scenario 2 are saved to a .mat file in the working directory and
   the MCS characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because 
   of a return statement after MCS computation. To activate it, remove this return statement.

5. e_coli/**des1tar2.m**

   Computes and characterizes genome-scale gene-MCS for the production of
   2,3-BDO in *E.coli.* A second target region is added to scenario 1 to compute, at the same time, 
   single substrate and co-feeding strategies using glucose, acetate and glycerol. Therefore, 
   the supply reactions for glucose, acetate and glycerol are specified as addition candidates. 
   Results of scenario 3 are saved to a .mat file in the working directory and the MCS
   characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because 
   of a return statement after MCS computation. To activate it, remove this return statement.

6. e_coli/**des2tar2.m**

   Computes and characterizes genome-scale gene-MCS for the production of
   2,3-BDO in *E. coli.* In addition to the changes in scenario 3, a second desired
   region is added to demand the support of higher ATP maintanance rates. This
   setup shows that a combination of multiple target and desired regions is
   possible and generates again qualitatively new solutions. Results of 
   scenario 4 are saved to a .mat file in the working directory and the MCS
   characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because 
   of a return statement after MCS computation. To activate it, remove this return statement.
   
7. pseudomonas/**pseudomonas_benchmark.m**

    Computes and characterizes gene-MCS for the production of 2,3-BDO in *Pseudomonas putida* from Glucose. 
    By changing options.compression_GPR or options.preproc_compression, GPR- and/or network compression can
    be activated or deactivated. Results are saved to a .mat file in the working directory and the MCS
    characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because of a
    return statement after MCS computation. To activate it, remove this return statement.

8. pseudomonas/**pseudomonas_desired2.m**

   Computes and characterizes genome-scale gene-MCS for the production of
   2,3-BDO in *P. putida* from a similar setup as in pseudomonas_benchmark. A second
   desired region is added to the setup to ensure that strain designs support higher ATP
   maintanance rates. The results of this computation are saved to a .mat file in the working directory and
   the MCS characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because 
   of a return statement after MCS computation. To activate it, remove this return statement.

9. pseudomonas/**pseudomonas_des1tar2.m**

   Computes and characterizes genome-scale gene-MCS for the production of
   2,3-BDO in *P. putida* A second target region is added to the benchmark setup to compute, at the same time, 
   single substrate and co-feeding strategies using glucose, acetate and glycerol. Therefore, 
   the supply reactions for glucose, acetate and glycerol are specified as addition candidates. 
   Results of this computation are saved to a .mat file in the working directory and the MCS
   characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because 
   of a return statement after MCS computation. To activate it, remove this return statement.
   
10. pseudomonas/**pseudomonas_des2tar2.m**

   Computes and characterizes genome-scale gene-MCS for the production of
   2,3-BDO in *P. putida* In addition to the changes in pseudomonas_des1tar2, a second desired
   region is added to demand the support of higher ATP maintanance rates. This
   setup shows that a combination of multiple target and desired regions is
   possible and generates again qualitatively new solutions. Results of 
   this computation are saved to a .mat file in the working directory and the MCS
   characterization/ranking is saved as a .tsv table.
   
11. yeast/**yeast_benchmark.m**

    Computes and characterizes gene-MCS for the production of 2,3-BDO in *Saccharomyces cerevisiae* from Glucose. 
    By changing options.compression_GPR or options.preproc_compression, GPR- and/or network compression can
    be activated or deactivated. Results are saved to a .mat file in the working directory and the MCS
    characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because of a
    return statement after MCS computation. To activate it, remove this return statement.

12. yeast/**yeast_desired2.m**

   Computes and characterizes genome-scale gene-MCS for the production of
   2,3-BDO in *S. cerevisiae* from a similar setup as in yeast_benchmark. A second
   desired region is added to the setup to ensure that strain designs support higher ATP
   maintanance rates. The results of this computation are saved to a .mat file in the working directory and
   the MCS characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because 
   of a return statement after MCS computation. To activate it, remove this return statement.

13. yeast/**yeast_des1tar2.m**

   Computes and characterizes genome-scale gene-MCS for the production of
   2,3-BDO in *S. cerevisiae* A second target region is added to the benchmark setup to compute, at the same time, 
   single substrate and co-feeding strategies using glucose and acetate. Therefore, 
   the supply reactions for glucose and acetate are specified as addition candidates. 
   Results of this computation are saved to a .mat file in the working directory and the MCS
   characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because 
   of a return statement after MCS computation. To activate it, remove this return statement.
   
14. yeast/**yeast_des2tar2.m**

   Computes and characterizes genome-scale gene-MCS for the production of
   2,3-BDO in *S. cerevisiae* In addition to the changes in yeast_des1tar2, a second desired
   region is added to demand the support of higher ATP maintanance rates. This
   setup shows that a combination of multiple target and desired regions is
   possible and generates again qualitatively new solutions. Results of 
   this computation are saved to a .mat file in the working directory and the MCS
   characterization/ranking is saved as a .tsv table. By default MCS characterization is skipped because 
   of a return statement after MCS computation. To activate it, remove this return statement.
   
15. synthetic_lethals/**synthetic_lethals_iML.m**
   
   Computes synthetic lethals in the genome scale *E. coli* model iML1515 up to the size of 4 gene knockouts. 
   The results of both functions are checked for consistency and the runtimes of both methods are returned in 
   a command window output.

Minor functions required for scripts:
-------------

16. functions/**verify_mcs.m** 

17. functions/**cell2csv.m**

18. functions/**relev_indc_and_mdf_Param.m**

19. functions/**start_parallel_pool_on_SLURM_node.m**

20. functions/**block_non_standard_products**

21. functions/**compare_mcs_sets.m**

22. functions/**text2num_mcs.m**

23. functions/**verify_lethals.m**

Model files:
-------------

24. e_coli/**benchmark_iJOcore.mat** - E. coli core model required for script (3)

25. e_coli/**iML1515.mat** - genome scale E. coli model required for scripts (3-6)

26. pseudomonas/**iJN746.mat** - genome scale P. putida model required for scripts (7-10)

27. yeast/**yeastGEM.xml** - genome scale S. cerevisiae model required for scripts (11-14)

28. yeast/**yeast_BiGGmetDictionary.csv** - Dictionary to replace the species identifiers from yeastGEM with Bigg identifiers

29. yeast/**yeast_BiGGrxnDictionary.csv** - Dictionary to replace the reaction identifiers from yeastGEM with Bigg identifiers


New (API) functions - included in the most recent *CellNetAnalyzer* toolbox and not available in this repository:
-------------

* **CNAgeneMCSEnumerator2**

   Function wrapper for CNAMCSEnumerator2 that allows the computation of
   gene-MCS using GPR association and GPR-rule compression, multiple target and desired regions and gene- and reaction
   deletions and additions with individual intervention cost factors.

* **CNAMCSEnumerator2**

   MCS computation with multiple target and desired regions, reaction additions
   and deletions and individual cost factors.

* **CNAgenerateGPRrules.m** 

   Translates GPR-rules provided in text form into a gene-protein-reaction
   mapping.

* **CNAintegrateGPRrules.m**

    Extends a metabolic network model with genes and GPR rules represented by
    pseudoreaction and pseudometabolites. Uses mapping generated by
    CNAgenerateGPRrules.m

* **CNAcharacterizeGeneMCS.m**

    Characterizes and ranks geneMCS by different criteria, such as product
    yield, ability to grow, implementation effort
    
* **testRegionFeas.m** 

    Tests if a model/mutant has feasible steady state flux vectors in a flux space spanned by a set of constraints (**V**·**r** ≤ **v**).

Remarks:
--------

-   If a fast but incomplete iterative MCS computation/search is preferred over
    a full MCS enumeration, set the parameter "enum_method" (e.g. in scripts
    3-6) from 2 to 1 and set a time or solution limit.
