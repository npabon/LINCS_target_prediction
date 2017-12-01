import numpy as np

class Perceptron(object):
    """Perceptron classifier.
    
    Parameters **(the same for all members of the class)
    ------------
    eta : float
      Learning rate (between 0.0 and 1.0)
    n_iter : int
      Passes over training dataset (epochs)
    random_state : int
      Random number generator seed for random weight initialization
    
    
    Attributes **(specific to each class member)
    ------------
    w_ : 1d-array
      Weights after fitting
    errors_ : list
      Number of misclassifications (updates) in each epoch
    """
    
    def __init__(self, eta=0.01, n_iter=50, random_state=1):
        self.eta = eta
        self.n_iter = n_iter
        self.random_state = random_state
        
    def fit(self, X, y):
        """Fit training data.
        
        Parameters
        ------------
        X : {array-like} shape = [n_samples, n_features]
          Training samples, where n_samples is the number
          of samples and n_features is the number of features
        y : array-like, shape = [n_samples]
          Target values.
          
        Returns
        ---------
        self : object
        """
        
        # draw random samples from a normal distribution
        rgen = np.random.RandomState(self.random_state)
        # loc = mean, scale = st.dev., size = number of samples
        self.w_ = rgen.normal(loc=0.0, scale = 0.01, size=1 + X.shape[1])
        self.errors_ = []
        
        for _ in range(self.n_iter):
            errors = 0
            for xi, target in zip(X,y): # pairs each sample with its label
                update = self.eta * (target - self.predict(xi))
                self.w_[0] += update # bias update
                self.w_[1:] += update * xi
                errors += int(update != 0) # number of misclassifications
            self.errors_.append(errors)
        return self
    
    def net_input(self, X):
        """Calculate net input to perceptron"""
        return np.dot(X, self.w_[1:]) + self.w_[0]
    
    def predict(self, X):
        """Return predicted class label (in a 1-element array) after unit step"""
        return np.where(self.net_input(X) >= 0.0, 1, -1)
        
        
        
        
        
        
        
        
        
        
        