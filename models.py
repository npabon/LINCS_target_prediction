
import numpy as np
import scipy
import itertools
from sklearn.ensemble import RandomForestClassifier

class LincsRandomForestClassifier(object):
    
    "WE ASSUME THE DATA IS GROUPED BY CELL LINE AND HAS 4 FEATURES PER CELL LINE"
   
    def __init__(self, n_cells_per_forest, n_estimators_per_forest=10, max_depth=None, max_features="auto", random_state=1):
        self.n_cells_per_forest = n_cells_per_forest
        self.n_estimators_per_forest = n_estimators_per_forest
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        
    def fit(self, X, y):
        '''
        Train several random forests, each one on a different
        subset of cells. Store forests in a dictionary called
        self.forests.
        '''
        # make sure we have enough data to work with
        min_num_cells = self.get_min_num_cells(X)
        assert min_num_cells >= self.n_cells_per_forest, "Too much missing data for n_cells_per_forest = %s. (Some samples only tested in %d cells)" % \
                                                         (self.n_cells_per_forest, min_num_cells)
        
        # generate cell subsets for training
        # ASSUMES 4 FEATURES PER CELL
        total_num_cells = int(X.shape[1] / 4) # THIS IS HARDCODED IN
        cell_subsets = itertools.combinations(np.arange(total_num_cells), self.n_cells_per_forest)
        
        # initialize dictionary to hold the forests
        self.forests = {}
        
        # train forest on each subset
        for cell_subset in cell_subsets:
            
            # find samples that have complete data from the cell subset
            cell_subset_idx = np.array([ 4*i + np.array([0, 1, 2, 3])for i in cell_subset ]).reshape(1,-1)[0].astype(int)
            cell_subset_data = X[:,cell_subset_idx]
            bad_sample_idx = np.isnan(cell_subset_data).any(axis=1)
            good_samples = cell_subset_data[~bad_sample_idx]
            good_labels = y[~bad_sample_idx]
            
            # train and store a RF classifier on this training subset
            # print('Growing forest for cell subset: %s' % str(cell_subset))
            forest = RandomForestClassifier(criterion='gini',
                                            n_estimators=self.n_estimators_per_forest,
                                            max_depth=self.max_depth,
                                            max_features=self.max_features,
                                            random_state=self.random_state,
                                            n_jobs=-1)
            forest.fit(good_samples, good_labels)
            self.forests[cell_subset] = forest            

        
    def get_min_num_cells(self, X):
        '''
        Calculate the minimum number of cells any sample has data for
        ASSUMES 4 FEATURES PER CELL LINE
        '''
        X_not_missing = ~np.isnan(X)
        num_cells_not_missing = np.count_nonzero(X_not_missing, axis=1) / 4
        min_num_cells = np.min(num_cells_not_missing)
        return min_num_cells
    
    def predict_proba(self, X):
        '''
        Return the class probabilities label OF ONE SINGLE SAMPLE FOR FUCKS SAKE
        '''
        # figure out which cell lines we have data for
        non_nan_idx = np.where(np.isnan(X) == False)[0]
        good_cells = (non_nan_idx[np.where(non_nan_idx/4%1 == 0)[0]] / 4).astype(int)
        
        # select appropriate forests and predict
        cell_subsets = itertools.combinations(good_cells, self.n_cells_per_forest)
        tree_predictions_ = []
        for cell_subset in cell_subsets:
            # extract appropriate data
            cell_subset_idx = np.array([ 4*i + np.array([0, 1, 2, 3])for i in cell_subset ]).reshape(1,-1)[0].astype(int)
            cell_subset_data = X[cell_subset_idx].reshape(1,-1) 
            # extract appropriate forest and make prediction
            forest = self.forests[cell_subset]
            tree_predictions = [ tree.predict(cell_subset_data) for tree in forest.estimators_ ]
            tree_predictions_.append(tree_predictions)
        
        # majority vote of all the trees in all the forests
        results = np.array(tree_predictions_).flatten()
        proba = results.sum() / len(results)
        return np.array([1.-proba, proba])
    
    def predict(self, X):
        '''
        Return the predicted class label OF ONE SINGLE SAMPLE FOR FUCKS SAKE
        '''
        class_probabilities = self.predict_proba(X)
        return np.argmax(class_probabilities)
    
    def predict_proba_(self, X):
        '''
        for a multidimentional X
        '''
        proba_ = np.array([ self.predict_proba(x) for x in X ])
        return proba_
    
    def predict_(self, X):
        '''
        for a multidimentional X
        '''
        predicted_classes = np.array([ self.predict(x) for x in X ])
        return predicted_classes
