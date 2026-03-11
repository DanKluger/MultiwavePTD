# MultiwavePTD
Reproducibility repository for the paper titled "M-estimation under Two-Phase Multiwave Sampling with Applications to Prediction-Powered Inference" by Dan M. Kluger and Stephen Bates. It gives the code and datasets that were used to run the experiments and generate the figures in the paper.


Please direct any questions or bug reports to dkluger@mit.edu. 

## Preferred citation

If using software from this repository, please cite the corresponding paper.

> Dan M. Kluger and Stephen Bates (2026). M-estimation under Two-Phase Multiwave Sampling with Applications to Prediction-Powered Inference. 	arXiv:2602.16933 [stat.ME]  [https://arxiv.org/abs/2602.16933](https://arxiv.org/abs/2602.16933)


## Dataset citations

If using the processed datasets, please cite their original source rather than this github repository.

The housing price, forest cover and tree cover datasets came from: 
> Rolf, E., Proctor, J., Carleton, T., Bolliger, I., Shankar, V., Ishihara, M., Recht, B., and
Hsiang, S. (2021a). A generalizable and accessible approach to machine learning with
global satellite imagery. Nature Communications, 12(1):4392.
>
> 
> Rolf, E., Proctor, J., Carleton, T., Bolliger, I., Shankar, V., Ishihara, M., Recht, B., and
Hsiang, S. (2021b). A generalizable and accessible approach to machine learning with
global satellite imagery. https://www.codeocean.com/capsule/6456296/tree/v2.


The AlphaFold dataset had origins and was later processed in:
> Bludau, I., Willems, S., Zeng, W.-F., Strauss, M. T., Hansen, F. M., Tanzer, M. C., Karayel,
O., Schulman, B. A., and Mann, M. (2022). The structural context of posttranslational
modifications at a proteome-wide scale. PLoS biology, 20(5):e3001636.
>
> 
> Angelopoulos, A. N., Bates, S., Fannjiang, C., Jordan, M. I., and Zrnic, T. (2023b).
Prediction-powered inference: Data sets. 10.5281/zenodo.8397451.

Note that for each of these datasets, [additional processing steps](https://github.com/DanKluger/Predict-Then-Debias_Bootstrap/blob/main/Datasets/ReadAndProcessData.Rmd) were taken before uploading the data to this repository as part of the following paper:

> Dan M. Kluger, Kerri Lu, Tijana Zrnic, Sherrie Wang, and Stephen Bates (2025). Prediction-Powered Inference with Imputed Covariates and Nonuniform Sampling. 2501.18577 [stat.ME]  [https://arxiv.org/abs/2501.18577](https://arxiv.org/abs/2501.18577)
