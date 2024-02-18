"""
parameters file consists of all the parameters need to build the RNN.

"""

# parameters needed for one hot encoding
"""
each amino acid represented in aas is given a number between 0 and 20.
do not modify the sequence in aas as it will change the identity of the amino acids.
"""
aas = ['C', 'A', 'G', 'R', 'H', 'I', 'L', 'K', 'D', 'F', 'P', 'S', 'M', 'Y', 'N', 'E', 'Q', 'W', 'V', 'T']


# RNN parameters ######
random_seed = 2012
lstm_size = 64  # can be 16 32 64 128...
lstm_layers = 1
batch_size = 1
learning_rate = 0.001  # can be 0.01, 0.001, 0.0001...
