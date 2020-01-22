# MATLAB script for LDA Training and Simulation
A MATLAB script used for training and testing a machine learning model, using a LDA algorithm (provided by Dr. Xiaorong Zhang of the ICE Lab at SFSU). 

- MATLAB script that reads data from an Excel file and simulates a machine learning model
- The main MATLAB script is named "Training_Simulation.m"
- Predicts the effectiveness of the ML model
   - Extracts features from the raw data
   - Splits the data into 5 sets of training and testing datasets
   - Cross validates results by testing each set of data and then averaging them for an ACA
- Generates weights for the trained LDA model and saves them in 2 text files (named Wg.txt and Cg.txt)
   - These weights can be loaded into our Tiva Code to make LDA predictions
- Currently set to train a ML model to recognize 7 gestures