# SARS-CoV-2 full spike BA.1 library deep mutational scanning validations

In this project we validate SARS-CoV-2 full spike deep muttaional scanning for BA.1 library. 

Original BA.1 full spike deep mutational scanning data can be found [here](https://dms-vep.github.io/SARS-CoV-2_Delta_spike_DMS/).

These scripts preform the following:

- [virus_titers_for_mAb_escape.ipynb](virus_titers_for_mAb_escape.ipynb) and [virus_titers_functional_mutants.ipynb](virus_titers_functional_mutants.ipynb) calculate titers for pseudotyped viruses used either for antibody escape or spike function analysis. 
-[functional_virus_titer_fold_change.ipynb](functional_virus_titer_fold_change.ipynb) plots fold chance in virus titers compared to unmutated spike
- [spike_neutralization.ipynb](spike_neutralization.ipynb) runs neutralization analysis for spike mutants against different monoclonal antibodies
- [VSVG_neutralization.ipynb](VSVG_neutralization.ipynb) runs neutralization analysis for VSV-G or unmutated BA.1 spike pseudotyped lentivirus against different monoclonal antibodies
- [Lycov1404_yeast_lenti_dms_comparison.ipynb](Lycov1404_yeast_lenti_dms_comparison.ipynb) compares deep mutational scanning data for full spike lentivirus-based DMS and [yeast-based RBD DMS](https://www.biorxiv.org/content/10.1101/2022.09.20.508745v1). 

[plasmid_maps](./plasmid_maps) subfolder contains plasmid maps for splike and VSV-G expression plasmids used for making pseudoviruses used in performing neutralization ploted in this analysis.  

To run the analysis using snakemake pipeline type:

```
snakemake --cores 1
```
After the run is finished analysis results are found in [results](./results) folder.
