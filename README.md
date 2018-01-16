# LINCS target prediction

A python pipeline for predicting drug targets using the NIH LINCS project's L1000 dataset.

- **LINCS_data_prep.ipynb** : A notebook for cleaning raw LINCS L1000 mRNA data and preparing it for feature construction and classification.
- **STRING_interaction_data_prep.ipynb** : A notebook for cleaning protein-protein interaction data from the String 4.0 database.
- **data_profiling.ipynb** : A notebook for examining properties of the complete L1000 dataset to determine which data we can use.
- **feature_construction.ipynb** : A notebook for using the cleaned L1000 data to construct the actual features used for classification.
- **classification.ipynb** : A notebook for constructing the training set and training the classifier.
- **prediction.ipynb** : A notebook for using our trained classifier to predict targets for new drugs.
