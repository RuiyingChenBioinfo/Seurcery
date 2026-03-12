# Seurcery

<img src="Public/Logo-Seurcery.png" align="right" width="115" height="135">

An R package designed for further analysis and visualization based on utilities of R package Seurat.

Please refer to [Tutorial html](Public/Seurcery_v0.8.6_tutorial.html) for detailed tutorials.

<details> <summary> Note about package development </summary>

Seurcery is under active development, so you may run into errors and small typos. We welcome users to write GitHub issues to report bugs, ask for help, and to request potential enhancements. 

</details>

## Installation

In R:

If dependencies (ggplot2, patchwork, rlang) are not installed:
```
install.packages("remotes")
remotes::install_github("RuiyingChenBioinfo/Seurcery@master", dependencies = TRUE)
```
Else:
```
install.packages("remotes")
remotes::install_github("RuiyingChenBioinfo/Seurcery@master", dependencies = FALSE)
```

## Citation
If you use Seurcery, please cite:

Chen, R., Wang, Y., & Li, T. (2026). Seurcery: An R package designed for further analysis and visualization based on utilities of R package Seurat (v0.8.6). Zenodo. https://doi.org/10.5281/zenodo.18983549

## Contact
* Ruiying Chen (chenruiying21@mails.ucas.ac.cn)
* Yixiao Wang (wangyixiao21@mails.ucas.ac.cn)
