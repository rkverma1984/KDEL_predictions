# coding=utf-8
"""
Background script to create RNN.

"""

import os
import pickle

import numpy as np
import tensorflow as tf

np.warnings.filterwarnings('ignore')


def get_batches(x, y, b_size = 5):
    n_batches = len(x) // b_size
    x, y = x[:n_batches * b_size], y[:n_batches * b_size]
    for ii in range(0, len(x), b_size):
        yield x[ii:ii + b_size], y[ii:ii + b_size]


def mean_squared_error(y, pred):
    # mean squared error calculation
    return np.mean((y - pred) ** 2)


def dump_file(filename, data):
    f = open(filename, 'wb')
    pickle.dump(data, f)
    f.close()


def graph_initialization(lstm_size = 64, lstm_layers = 1, batch_size = 1, learning_rate = 0.001, epochs = 500,
                         whole_x = None,
                         whole_y = None, train_index = None, test_index = None, ckpt_dir = None, ckpt = None,
                         caltype = None,
                         cell = 'LSTM', sequences_for_prediction = None, num2aa = '', out_file_name = None,
                         fixed_seed = '2012', gpu_id = '', cv = True):
    """
    :param cv: define if to do cross validation is to be done.
    :param gpu_id: which gpu to use. used to set CUDA_VISIBLE_DEVICES.
    :param cell: type of RNN to be used. Currently accepted RNN includes "GRU" or "LSTM"
                    LSTM refers to Standard LSTM
                    GRU refers to Gated Recurrent Unit which is a modification over LSTM
    :param caltype: if the model is needed to be built from scratch use "model" else use "data" to use previously
    saved checkpoints
    :param ckpt: prefix of the checkpoint file
    :param ckpt_dir: directory where checkpoints are to be saved
    :param test_index: testing data
    :param train_index: training data
    :param whole_y: Tg-response values from the input xlsx file
    :param whole_x: sequences in one-hot encoding
    :param epochs: one epoch refers to one forward pass and one backward pass of all the training examples
    :param learning_rate: relation with new weights: new_weight = existing_weight — learning_rate * gradient
    :param batch_size: the number of training examples in one forward/backward pass. The higher the batch size,
    the more memory space you'll need.
    :param lstm_layers: no of layers to be used
    :type  lstm_size:  no of hidden state of an LSTM unit
    :param fixed_seed: random number seed
    :param out_file_name: name of the pickle file name
    :param num2aa: reverse dictionary containing one-hot encoding as keys and amino-acid identities as values
    :param sequences_for_prediction:
    """
    
    # config = tf.ConfigProto(device_count = {'GPU': int(gpu_id)})
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '4'
    
    if gpu_id:
        print('\nrunning on gpu_device =', "/gpu:%s\n" % gpu_id)
    else:
        print('\nrunning on cpu\n')
    
    if cell == 'GRU' or cell == 'LSTM':
        
        # print('random seed:', fixed_seed)
        
        np.random.seed(fixed_seed)
        tf.set_random_seed(fixed_seed)
        
        graph = tf.Graph()
        # Add nodes to the grp
        with graph.as_default():
            inputs_ = tf.placeholder(tf.int32, [None, None], name='inputs')
            labels_ = tf.placeholder(tf.float32, [None, None], name='labels')
            keep_prob = tf.placeholder(tf.float32, name='keep_prob')
        
        with graph.as_default():
            # Your basic LSTM/GRU cell (choose one below)
            
            lstm = ''
            if cell == 'GRU':
                lstm = tf.contrib.rnn.GRUCell(lstm_size)
            elif cell == 'LSTM':
                lstm = tf.contrib.rnn.BasicLSTMCell(lstm_size)
            
            # Add dropout to the cell
            drop = tf.contrib.rnn.DropoutWrapper(lstm, output_keep_prob=keep_prob, seed=fixed_seed)
            
            # Stack up multiple LSTM layers, for deep learning
            cell = tf.contrib.rnn.MultiRNNCell([drop] * lstm_layers)
            
            # Getting an initial state of all zeros
            initial_state = cell.zero_state(batch_size, tf.float32)
        
        with graph.as_default():
            # one hot encoding
            one_hot_code = tf.one_hot(inputs_, 20)
            outputs, final_state = tf.nn.dynamic_rnn(cell, one_hot_code, initial_state=initial_state)
        
        with graph.as_default():
            # get the prediction from the final output
            # using RELU as the nonlinear function
            # initialize the weight matrix using truncated normal distribution
            
            predictions = tf.contrib.layers.fully_connected(outputs[:, -1], 1, activation_fn=tf.nn.relu,
                                                            weights_initializer=tf.truncated_normal_initializer(
                                                                stddev=0.1, seed=fixed_seed))
            
            # minimize the mean squared error
            cost = tf.losses.mean_squared_error(labels_, predictions)
            
            optimizer = tf.train.AdamOptimizer(learning_rate).minimize(cost)
        
        with graph.as_default():
            accuracy = tf.losses.mean_squared_error(labels_, predictions)
        
        with graph.as_default():
            saver = tf.train.Saver()
        
        if caltype == 'data':
            result_all = {}
            with tf.Session(graph=graph) as sess:
                saver = tf.train.import_meta_graph(os.path.join(ckpt_dir, ckpt + '.meta'))
                saver.restore(sess, tf.train.latest_checkpoint(ckpt_dir))
                
                val_state = sess.run(cell.zero_state(batch_size, tf.float32))
                for x in sequences_for_prediction:
                    # print(x)
                    feed = {
                        inputs_: [x],
                        keep_prob: 1,
                        initial_state: val_state
                    }
                    prediction, _ = sess.run([predictions, final_state], feed_dict=feed)
                    
                    seq = ''.join([num2aa[s] for s in x])
                    result_all[seq] = prediction
                
                sorted_result = sorted(result_all.items(), key=lambda pair: (pair[1], pair[0]), reverse=True)
                
                o_file_h = open(out_file_name.replace('.pkl', '.txt'), 'w')
                for a in sorted_result:
                    o_file_h.write("%s\t%.5f\n" % (a[0], a[1][0][0]))
                dump_file(out_file_name, sorted_result)
        
        elif caltype == 'model':
            iteration = 1
            
            # print(accuracy)
            with tf.Session(graph=graph) as sess:
                
                # training...
                sess.run(tf.global_variables_initializer())
                
                for e in range(epochs):
                    # print('epoch', e)
                    state = sess.run(initial_state)
                    
                    for x, y in get_batches(whole_x[train_index], whole_y[train_index], batch_size):
                        # print(x, y)
                        feed = {
                            inputs_: x,
                            labels_: y[:, None],
                            keep_prob: 0.8,
                            initial_state: state
                        }
                        loss, _, _ = sess.run([cost, final_state, optimizer], feed_dict=feed)
                        # print('loss', loss)
                        # output some training error every 1000 steps
                        if iteration % 1000 == 0:
                            print("Epoch: {}/{}".format(e, epochs),
                                  "Iteration: {}".format(iteration),
                                  "Train loss: {:.3f}".format(loss))
                        
                        iteration += 1
                
                # predicting...
                data_y = []
                predict_y = []
                
                val_state = sess.run(cell.zero_state(batch_size, tf.float32))
                for x, y in get_batches(whole_x[test_index], whole_y[test_index], batch_size):
                    feed = {
                        inputs_: x,
                        keep_prob: 1,
                        initial_state: val_state
                    }
                    prediction, _ = sess.run([predictions, final_state], feed_dict=feed)
                    
                    data_y.append(y)
                    predict_y.append(prediction)
                
                data_y = np.ravel(data_y)
                predict_y = np.ravel(predict_y)
                
                # output the score by R squared or mean_squared_error
                mse = mean_squared_error(data_y, predict_y)
                print('MEAN Squared Error = %.3f' % mse)
                
                if cv == True:
                    saver.save(sess, os.path.join(ckpt_dir, ckpt))
                else:
                    pass
                
                return mse

