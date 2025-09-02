Overview

Hackathon hosted by the OpenGeoHub Foundation as part of the 2025 OpenGeoHub Summer School, designed to explore the potential of cutting-edge remote sensing embeddings for forest monitoring. This challenge invites students, early-career researchers, and GIS/EO practitioners to develop machine learning models that predict forest canopy height from Google DeepMind’s AlphaEarth Foundations (AEF) embeddings, derived from Sentinel-1, Sentinel-2, and Landsat imagery.

Participants will compete to benchmark the accuracy and generalizability of AEF embeddings against GEDI-derived canopy height measurements, gaining hands-on experience with real-world data, open-source tools, and environmental AI workflows.

The winning team will present their results at the 2025 OpenGeoHub Summer School—a non-profit event bringing together developers, scientists, and institutions committed to open environmental data and software. The 2025 edition is co-organized with EO Council, Netherlands Space Office, and SURF, and includes live demos, workshops, panel discussions, and social gatherings.

Hackathons are supported by the EU-funded Open-Earth-Monitor project, under the Horizon Europe research and innovation programme (grant agreement No. 101059548).

Start
4 days ago
Close
3 days to go
Description

The hackathon dataset contains ~50.000 quality screened GEDI L2A shots distributed across southern Netherlands and northern Belgium. For each 25m footprint we provide a single response, rh98, in centimeters (cm). The data extraction and stratification are documented and presented in the lecture Introduction of cloud-native vector format: hands-on in Python environment by Yu-Feng Ho.The acronym rh stands for relative height and measures the height difference between the detected ground return and the point at which 98% (hence rh98) of the cumulative lidar-waveform energy is reached. This metric is widely used as a proxy for dominant tree height because it captures the upper canopy while trimming most noise.

Fig 1: Distribution of GEDI Sample points in Benelux

Fig 2: Vanilla model CHM over Wageningen

The predictor variables you will be using to measure the canopy height are the 64 Google DeepMind AlphaEarth Foundations (AEF) remote-sensing embeddings; values are in the [-1, +1] interval. They capture multi-season Sentinel-1, Sentinel-2 + Landsat spectral-texture context at a spatial resolution of 10m. Attached is a user-friendly description of the predictors, their production and usage.
Evaluation

The result will be evaluated by Root Mean Square Error (RMSE). This formula calculates the square root of the average of the squared differences between the predicted values ŷ and the actual values y, providing a measure of the differences in the units of the response variable.

Citation

Carmelo Bonannella, Yu-Feng Ho. Canopy height mapping using Google `AEF embeddings. https://kaggle.com/competitions/canopy-height-mapping-using-google-aef-embeddings, 2025. Kaggle.

Dataset Description

The hackathon dataset contains ~50.000 quality screened GEDI L2A shots distributed across southern Netherlands and northern Belgium. For each 25m footprint we provide a single response, rh98, in centimeters (cm).
Files
train.csv

The training set contains ~39,500 points collected over southern Netherlands and northern Belgium, where each sample includes a single response variable, spatial coordinates (see below) and predictors.
test.csv

the test set contains ~17,000 of the whole 56k points overlay and includes points available for internal validation of your models. It will help you fine-tune and assess the accuracy of your predictions during the development phase. In additional, the hidden dataset will be used that represents the final test of your model's effectiveness in real-world scenarios. Your goal is to predict these hidden values accurately. Height values for the hidden set come from trees surveyed in the field with traditional methods. Field-measured dominant heights are typically 1–2 m higher than GEDI RH98. If you apply a bias / linear offset to compensate, please:

    document the method (e.g. constant shift, regression, eco-type specific, …);
    state the exact offset value(s);
    cite the evidence or rationale you used.

Include this note with your submission.
Columns

    rh98 - height difference between the detected ground return and the point at which 98%of the cumulative lidar-waveform energy is reached
    lon - longtitude of the sample point in EPSG:4326
    lat - latitude of the sample point in EPSG:4326
    cv_{x}_google.embedding_nl - value of google embedding vector number x, from 1 - 64

Files

2 files
Size

73.89 MB
Type

csv
License
