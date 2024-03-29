## Complete map of SARS-CoV-2 RBD mutations that escape the monoclonal antibody LY-CoV555 and its cocktail with LY-CoV016

For experimental background, see our paper **[here](https://doi.org/10.1016/j.xcrm.2021.100255)**.

### What data are shown here?
We are showing mutations to the SARS-CoV-2 RBD that escape antibody binding as measured using mutational antigenic profiling. Raw data is available [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_MAP_LY-CoV555/blob/main/results/supp_data/LY_cocktail_raw_data.csv).
The drop-down menus can be used to select the escape-mutation maps for each antibody.

When you click on sites, they will be highlighted on the protein structure of the ACE2-bound RBD ([PDB 6M0J](https://www.rcsb.org/structure/6M0J), or antibody-bound RBD structures ([PDB 7KMG](https://www.rcsb.org/structure/7KMG) and [PDB 7C01](https://www.rcsb.org/structure/7C01)).

At the site level you can visualize one of two quantities:

 - *total escape* is the sum of the escape from all mutations at a site.

 - *max escape* is the magnitude of the largest-effect escape mutation at each site.

At the mutation level, the height of each letter is proportional to the extent to which that amino-acid mutation escapes antibody binding.
You can color the logo plot letters in four ways:

 - *escape color ACE2 bind* means color letters according to how that mutation affects ACE2 binding as measured in our prior deep mutational scanning ([Starr et al. (2020)](https://doi.org/10.1016/j.cell.2020.08.012), with yellow meaning highly deleterious, and brown meaning neutral or beneficial for ACE2 binding. See Figure 1 in the preprint linked above for color scale.
 
 - *escape color RBD expr* means color letters according to how that mutation affects RBD expression as measured in [Starr et al. (2020)](https://doi.org/10.1016/j.cell.2020.08.012).

 - *escape color gray* means color all letters gray.

 - *escape color func group* means color letters by functional group.
