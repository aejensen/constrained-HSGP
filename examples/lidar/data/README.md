# Lidar observations

`lidar.csv` contains all 221 released observations, in their original order.
`range_m` is recorded distance in meters and `logratio` is the recorded log
ratio of received laser signals. The example models the negative log ratio.
No observations have been excluded, averaged further, or simulated.

The source is `data/lidar.txt.gz` in **SemiPar 1.0-4.2**, by Matt Wand and
maintained in that release by Billy Aung Myint. The package is archived on
CRAN. Its [data documentation](https://search.r-project.org/CRAN/refmans/SemiPar/html/lidar.html)
credits Sigrist, M. (ed., 1994), *Air Monitoring by Spectroscopic Techniques*,
Wiley, and cites Ruppert, Wand, and Carroll (2003), *Semiparametric Regression*.
The original package is available from the
[CRAN source archive](https://cran.r-project.org/src/contrib/Archive/SemiPar/SemiPar_1.0-4.2.tar.gz).

The CSV renames `range` to `range_m` and retains all recorded values.

The source package declares **GPL (>= 2)**. That declaration and the upstream
author information are preserved in [SemiPar-DESCRIPTION](SemiPar-DESCRIPTION).
The bundled observations are redistributed with that source-package notice.
No separate data-specific license is stated in the inspected package
description or lidar documentation. A copy of the GNU General Public License,
version 2, is provided in [GPL-2.txt](GPL-2.txt). This notice preserves the
upstream package terms and does not claim new ownership of the measurements.
