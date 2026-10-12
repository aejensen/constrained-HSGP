# Leaf light-response data

`leaf1818_raw.csv` contains the 16 records for sample 1818 from Rogers,
Ely, and Serbin (2019), *Leaf photosynthetic parameters: Quantum yield,
convexity, respiration, gross CO2 assimilation rate and raw gas exchange
data, Utqiagvik (Barrow), Alaska, 2016*, NGEE Arctic.
The dataset DOI is [10.5440/1482338](https://doi.org/10.5440/1482338).

The records come from the authors' [released raw data](https://github.com/TESTgroup-BNL/Rogers_etal_NGEEArctic_LightResponse/blob/master/input_data/NGEE-Arctic_Light_Response_Barrow_2016_Raw_Data.csv).
The column names and values are retained. `Sample_ID` identifies the leaf,
`PARi` is incident irradiance in micromoles of quanta per square meter per
second, and `Photo` is net CO2 assimilation in micromoles per square meter
per second. `Time` records measurement time and `QC` records manual flags.

`read_leaf_observations()` applies the source study's filters for flags,
conductance, intercellular CO2, and leaf temperature, then excludes the
dark measurement. This leaves the 13 observations used in the manuscript,
sorted by irradiance. The time order is retained for predictive checks.

The accompanying study is Rogers, Serbin, Ely, and Wullschleger (2019),
*New Phytologist* 223, 167–179, [doi:10.1111/nph.15750](https://doi.org/10.1111/nph.15750).
