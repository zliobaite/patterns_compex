# Patterns of Competitive Exclusion in the Mammalian Fossil Record, Support material

### Original article:
Galbrun, E., Hermansen, J.S., Žliobaitė, I. (2023). *Patterns of Competitive Exclusion in the Mammalian Fossil Record.* In: Casanovas-Vilar, I., van den Hoek Ostende, L.W., Janis, C.M., Saarinen, J. (eds) *Evolution of Cenozoic Land Mammal Faunas and Ecosystems.* Vertebrate Paleobiology and Paleoanthropology. Springer, Cham. <https://doi.org/10.1007/978-3-031-17491-9_9>

Presented at the [Society of Vertebrate Paleontology (SVP)](https://vertpaleo.org/) [81st annual meeting](https://vertpaleo.org/wp-content/uploads/2021/10/SVP_2021_VirtualBook_final.pdf) (November 1-5, 2021)

#### Abstract

Due to recent common ancestry, species belonging to the same genus are expected to be more similar with respect to their phenotype, and hence exhibit less niche divergence than species belonging to different genera. Consequently, congeneric species are expected to compete intensely for resources, and therefore to be segregated in space. Yet, despite the longstanding history of this hypothesis of congeneric competitive exclusion, empirical evidence in support of it is at best limited. Here, we analyze co-occurrence patterns of species that belong to the same genera in the mammalian fossil record kept in the NOW database, considering separately Europe during the Neogene, and North America during the Oligocene--Neogene. We assess co-occurrence patterns in comparison to baselines where competitive exclusion is  obfuscated through randomization. We find that congeneric species occur together notably less than would be expected at random, with large herbivores being more segregated than large carnivorans and small mammals.

## List of contents

- The `scripts` folder contains the Python scripts for carrying out the computational experiments
- The `times` folder contains files defining the time bins used in the experiments
- The `db_dumps` folder can be used to store NOW database dump as csv files
- The `localities` folder containing lists of fossil localities
- The `docs` folder containing code for interactive plotly figures of fossil localities in space and time, as well as supporting svg figures, prepared using the scripts
- The `rnd_xps.zip` archive contains raw results from the randomization experiments with one thousand repetitions for each of the three null-models on the NOW database dump downloaded on November 25, 2020
- The `manuscript.pdf` file contains the main text
- The `appendix.pdf` file contains supplementary information about the data and an exhaustive report of computational experiments
- The `SVP_presentation.mp4` file contains the recording of the presentation given at the SVP annual meeting
- The `SVP_slides.pdf` file contains the slides of the presentation given at the SVP annual meeting
- The `tikz_fig.tex` LaTeX file to compile tikz figures
- The `README.md` file, this file

## Running the computational analysis

1. Download a dump of the [NOW database](https://nowdatabase.org/) and save it as `./db_dumps/NOW_latest_public.csv`:
    1. Go to <https://nowdatabase.org/now/database/>, then *Enter Database without login -->*
    2. Click on *Locality* in the left-hand side panel, then *Export*
    3. Select *include species lists* and set *Field separator* to *comma*, then click on *All NOW localities*
    4. When the dump is ready, select *Save as CSV* and save the file as `NOW_latest_public.csv` in the `db_dumps` folder 

2. Go to the `scripts` folder and run the `filter_fossils.py` script to prepare the data:

    ```
     python filter_fossils.py
    ```
    
    This will extract the subsets for the different ecological groups and continents, and save them as separate files in a `prepared_data` folder
    
3. Run the randomization experiments with the `run_rnd.py` script. For instance, to carry out the experiments for North America, large mammals, 1000 data copies randomized with the Curveball algorithm, computing the number of genera with co-occurring species:

    ```
     python run_rnd.py  --continent NA --group L --rnd_nb 1000 --null_model CB --measure gen
    ```
    
     This will generate randomized datasets, compute the co-occurrence statistics from these datasets as well as from the original dataset, and save the raw results to a `xps_rnd` folder.

4. Produce figures to visualize the results with the `run_plots.py` script. For instance, to produce figures for the previously run experiments for North America, large mammals:

    ```
     python run_plots.py  --continent NA --group L --null_model CB --measure gen
    ```
    
    This will generate PDF figures from the raw results stored in the `xps_rnd` folder and save them to a `figs` folder.
    Tikz figures can be compiled using the `tikz_fig.tex` LaTeX file, changing the input command to the desired tikz figure source. 
    
## SVG figures

Scalable vector graphics versions with species labels of Appendix Figures 1--6, i.e. distribution of species occurrences between different orders and genera over time, are available online:

- [Europe, large herbivores](https://zliobaite.github.io/patterns_compex/bbl-genO_EU-L.svg)
- [Europe, small mammals](https://zliobaite.github.io/patterns_compex/bbl-genO_EU-S.svg)
- [Europe, large carnivorans](https://zliobaite.github.io/patterns_compex/bbl-genO_EU-C.svg)
- [North America, large herbivores](https://zliobaite.github.io/patterns_compex/bbl-genO_NA-L.svg)
- [North America, small mammals](https://zliobaite.github.io/patterns_compex/bbl-genO_NA-S.svg)
- [North America, large carnivorans](https://zliobaite.github.io/patterns_compex/bbl-genO_NA-C.svg)
