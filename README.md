# 2020_Grant_Thiery_C3S511_Task7.4

Scripts used in the lake temperature - precipitation coupling analysis for the C3S_511 thematic assessment task 7.4.

## To install
Python scripts should be run on python3.
Python environment available here as "environment.yml"

## For users
This repository includes the processing and plotting scripts.

"main.py" controls output by options commented with "<< SELECT >>"
* Processing options for:
  * ERA5 vs ERA5-Land dataset choice (flag_ds)
  * Pearson vs Spearman correlation (flag_corrtst)
  * Process and plot vs plot only (flag_proc)
  * Variables to assess (flag_coup)
* Plotting options:
  * Global vs regional (flag_subset)
  * Wind vector overlay (flag_wind)
  * Projection choice (flag_proj)
  * General plotting options for figure size, panel spacing, colorbar, boundaries
  
"funcs.py" called by "main.py" for processing and plotting.

## Versions
Version 0.1.0 - August 2020  


## License
This project is licensed under the MIT License. See also the [LICENSE](https://github.com/VUB-HYDR/2020_Grant_etal_NGEO/blob/master/LICENSE.md) file.
