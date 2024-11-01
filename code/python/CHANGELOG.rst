Changes
=======

Version 0.2.9, 2024/11/01
-------------------------

  - update: STAR/Salmon pipeline integrates variant calling

Version 0.2.8, 2023/06/29
-------------------------

  - update: kallisto pipeline now integrates file finder

Version 0.2.7, 2023/04/28
-------------------------

  - bugfix: run SAMtools in shell environment for deduplication

Version 0.2.6, 2023/04/27
-------------------------

  - bugfix: run SAMtools as joined string

Version 0.2.5, 2023/04/27
-------------------------

  - bugfix: typo in deduplication

Version 0.2.4, 2023/04/18
-------------------------

  - collect_trimmed_files() now finds single trimmed reads without specific number

Version 0.2.3, 2023/04/18
-------------------------

  - find_files, can now find single reads without a read num pattern

Version 0.2.2, 2023/04/13
-------------------------

  - Postrun: made organisms flexible, configurable

Version 0.2.1, 2023/04/11
-------------------------

  - STAR: support for multi-valued quantMode

Version 0.2.0, 2023/04/11
-------------------------

  - Salmon using transcriptome aligned file when possible

Version 0.1.9, 2023/04/07
-------------------------

  - added more STAR options
  - added Salmon options

Version 0.1.8, 2023/04/07
-------------------------

  - single end files support

Version 0.1.7, 2023/02/28
-------------------------

  - fixed STAR pipeline deindent

Version 0.1.6, 2023/02/09
-------------------------

  - configurable STAR indexing
  - more options for STAR
  - GFF files can be used in STAR

Version 0.1.5, 2023/02/09
-------------------------

  - outsamoptions in gs_prepare fixed

Version 0.1.4, 2023/01/24
-------------------------

  - SLURM Array throttling added
  - user can specify post run output directory
  - fixed outSAMattributes option

Version 0.1.3, 2023/01/24
-------------------------

  - Refactoring: renaming of Python modules
  - Refactoring: indexing code moved to index module
  - Post-run step adds TPM/Count extraction and MultiQC

Version 0.1.2, 2023/01/23
-------------------------

  - STAR/Salmon indexing is a separate step of gs_submit

Version 0.1.1, 2023/01/23
-------------------------

  - log file names are now specific to array and index names
  - subprocess.run() for salmon_quant step

Version 0.1.0, 2023/01/17
-------------------------

fixes

Version 0.0.9, 2023/01/17
-------------------------

fixes

Version 0.0.8, 2023/01/09
-------------------------

Run Kallisto jobs from gs_submit

Version 0.0.7, 2023/01/06
-------------------------

RNA seq job submitted as array now

Version 0.0.6, 2023/01/05
-------------------------

gs_prepare as separate command tool

Version 0.0.5, 2023/01/05
-------------------------

relative import in run_star_salmon.py

Version 0.0.4, 2023/01/05
-------------------------

configurable outSAMattributes

Version 0.0.3, 2022/12/15
-------------------------

Fix: deduplication

Version 0.0.2, 2022/12/15
-------------------------

Fix: gs_submit, doesn't need GS_HOME path anymore

Version 0.0.1, 2022/12/09
-------------------------

Initial PyPI version
