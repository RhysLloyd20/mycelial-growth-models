# Mycelial growth models

Code accompanying:

> Lloyd, R., Boddy, L., Windsor, F., Grieneisen, V., Christofides, S. Mycelial growth models and their relevance to cord-forming fungi. *In preparation.*

Reimplementations of published mycelial growth models, accompanying the above paper. Each script is seeded to reproduce the corresponding figure in the paper (remove or change the seed to explore different outputs).

## Repository contents

### population-level/

| File | Model |
|---|---|
| `figure_2.R` | Edelstein-Keshet & Ermentrout (1989) |
| `figure_3.R` | Davidson (1998) |
| `figure_4.R` | Regalado et al. (1996) |

### individual-based/

| File | Model |
|---|---|
| `figure_5.nlogo` | Ermentrout & Edelstein-Keshet (1993) |
| `figure_6.nlogo` | Boswell et al. (2007) |
| `figure_9.nlogo` | Hopkins & Boswell (2012) |

## Requirements

- [R](https://www.r-project.org/) (v4.4.3) with packages: `deSolve`, `fields`, `ggplot2`
- [NetLogo](https://ccl.northwestern.edu/netlogo/) (v6.4.0)

## References

Boswell, G.P. et al. (2007). Bull. Math. Biol. 69, 605–634. doi:[10.1007/s11538-005-9056-6](https://doi.org/10.1007/s11538-005-9056-6)

Davidson, F.A. (1998). J. Theor. Biol. 195, 281–292. doi:[10.1006/jtbi.1998.0739](https://doi.org/10.1006/jtbi.1998.0739)

Edelstein-Keshet, L. & Ermentrout, B. (1989). SIAM J. Appl. Math. 49, 1136–1157. doi:[10.1137/0149068](https://doi.org/10.1137/0149068)

Ermentrout, G.B. & Edelstein-Keshet, L. (1993). J. Theor. Biol. 160, 97–133. doi:[10.1006/jtbi.1993.1007](https://doi.org/10.1006/jtbi.1993.1007)

Hopkins, S.M. (2011). A hybrid mathematical model of fungal mycelia: Tropisms, polarised growth and application to colony competition. PhD Thesis, University of Glamorgan.

Hopkins, S.M. & Boswell, G.P. (2012). Fungal Ecol. 5, 124–136. doi:[10.1016/j.funeco.2011.06.006](https://doi.org/10.1016/j.funeco.2011.06.006)

Nychka, D., Furrer, R., Paige, J., Sain, S. (2021). fields: Tools for spatial data. R package. doi:[10.5065/D6W957CT](https://doi.org/10.5065/D6W957CT)

R Core Team (2025). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. [https://www.R-project.org/](https://www.R-project.org/)

Regalado, C.M. et al. (1996). Mycol. Res. 100, 1473–1480. doi:[10.1016/s0953-7562(96)80080-3](https://doi.org/10.1016/s0953-7562(96)80080-3)

Soetaert, K., Petzoldt, T., Setzer, R.W. (2010). Solving Differential Equations in R: Package deSolve. J. Stat. Softw. 33(9), 1–25. doi:[10.18637/jss.v033.i09](https://doi.org/10.18637/jss.v033.i09)

Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4. [https://ggplot2.tidyverse.org](https://ggplot2.tidyverse.org)

Wilensky, U. (1999). NetLogo. Center for Connected Learning and Computer-Based Modeling, Northwestern University. [https://ccl.northwestern.edu/netlogo/](https://ccl.northwestern.edu/netlogo/)
