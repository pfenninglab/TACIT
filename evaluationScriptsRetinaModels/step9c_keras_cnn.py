from step10_roc_evaluations import get_filtered_input
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Conv2D, MaxPooling2D, Flatten, BatchNormalization
from tensorflow.keras.metrics import AUC, TruePositives, FalsePositives, TrueNegatives, FalseNegatives, Precision, Recall
from tensorflow.keras.regularizers import l2
from tensorflow.keras.optimizers import Adam, SGD, Nadam
from tensorflow.keras.models import load_model
from sklearn.metrics import f1_score, recall_score, precision_score
from sklearn.model_selection import KFold
import tensorflow as tf
import tensorflow.keras.backend as K
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
import argparse, sys, os
from tqdm import tqdm
os.system("nvidia-smi")
print(tf.__version__)
print(tf.keras.__version__)

# set seed
from numpy.random import seed
seed(1)
from tensorflow.random import set_seed
set_seed(2)

def f1_weighted(true, pred): #shapes (batch, 4)

  #for metrics include these two lines, for loss, don't include them
  #these are meant to round 'pred' to exactly zeros and ones
  #predLabels = K.argmax(pred, axis=-1)
  #pred = K.one_hot(predLabels, 4)


  ground_positives = K.sum(true, axis=0) + K.epsilon()       # = TP + FN
  pred_positives = K.sum(pred, axis=0) + K.epsilon()         # = TP + FP
  true_positives = K.sum(true * pred, axis=0) + K.epsilon()  # = TP
      #all with shape (4,)

  precision = true_positives / pred_positives
  recall = true_positives / ground_positives
      #both = 1 if ground_positives == 0 or pred_positives == 0
      #shape (4,)

  f1 = 2 * (precision * recall) / (precision + recall + K.epsilon())
      #still with shape (4,)

  weighted_f1 = f1 * ground_positives / K.sum(ground_positives)
  weighted_f1 = K.sum(weighted_f1)


  return 1 - weighted_f1 #for metrics, return only 'weighted_f1'



class Metrics(tf.keras.callbacks.Callback):
  def __init__(self, validation, human_validation, evals):
    super(Metrics, self).__init__()
    self.validation = validation
    self.human_validation = human_validation
    (self.eval1x, self.eval1y, self.eval2x, self.eval2y, self.eval3x, self.eval3y, self.eval4x, self.eval4y, self.eval5x, self.eval5y, self.eval6x, self.eval6y, self.eval7x, self.eval7y) = evals


  def on_train_begin(self, logs={}):
    self.val_f1s = []
    self.val_recalls = []
    self.val_precisions = []
    self.human_val_f1s = []

  def on_epoch_end(self, epoch, logs={}):

    eval1score = np.sum(np.asarray(self.model.predict(self.eval1x)).round().flatten() == self.eval1y)/len(self.eval1y)
    eval2score = np.sum(np.asarray(self.model.predict(self.eval2x)).round().flatten() == self.eval2y)/len(self.eval2y)
    eval3score = np.sum(np.asarray(self.model.predict(self.eval3x)).round().flatten() == self.eval3y)/len(self.eval3y)
    eval4score = np.sum(np.asarray(self.model.predict(self.eval4x)).round().flatten() == self.eval4y)/len(self.eval4y)
    eval5score = np.sum(np.asarray(self.model.predict(self.eval5x)).round().flatten() == self.eval5y)/len(self.eval5y)
    eval6score = np.sum(np.asarray(self.model.predict(self.eval6x)).round().flatten() == self.eval6y)/len(self.eval6y)
    eval7score = np.sum(np.asarray(self.model.predict(self.eval7x)).round().flatten() == self.eval7y)/len(self.eval7y)

    val_targ = self.validation[1]
    val_predict = (np.asarray(self.model.predict(self.validation[0]))).round()

    val_f1 = f1_score(val_targ, val_predict)
    #val_recall = recall_score(val_targ, val_predict)
    #val_precision = precision_score(val_targ, val_predict)

    human_val_predict = (np.asarray(self.model.predict(self.human_validation[0]))).round()
    human_val_f1 = f1_score(self.human_validation[1], human_val_predict)

    self.val_f1s.append(round(val_f1, 6))
    #self.val_recalls.append(round(val_recall, 6))
    #self.val_precisions.append(round(val_precision, 6))
    self.human_val_f1s.append(round(human_val_f1))

    print(f' — val_f1: {val_f1}, - human_val_f1: {human_val_f1}, - eval1: {eval1score}, - eval2: {eval2score}, -eval3: {eval3score}, - eval4: {eval4score}, - eval5: {eval5score}, - eval6: {eval6score}, - eval7: {eval7score} ')

class CyclicLR(tf.keras.callbacks.Callback):
    def __init__(self,base_lr, max_lr, step_size, base_m, max_m, cyclical_momentum):
      self.base_lr = base_lr
      self.max_lr = max_lr
      self.base_m = base_m
      self.max_m = max_m
      self.cyclical_momentum = cyclical_momentum
      self.step_size = step_size
      self.clr_iterations = 0.
      self.cm_iterations = 0.
      self.trn_iterations = 0.
      self.history = {}

    def clr(self):
      cycle = np.floor(1+self.clr_iterations/(2*self.step_size))
      if cycle == 2:
        x = np.abs(self.clr_iterations/self.step_size - 2*cycle + 1)
        return self.base_lr-(self.base_lr-self.base_lr/100)*np.maximum(0,(1-x))
      else:
        x = np.abs(self.clr_iterations/self.step_size - 2*cycle + 1)
        return self.base_lr + (self.max_lr-self.base_lr)*np.maximum(0,(1-x))

    def cm(self):
      cycle = np.floor(1+self.clr_iterations/(2*self.step_size))
      if cycle == 2:
        x = np.abs(self.clr_iterations/self.step_size - 2*cycle + 1)
        return self.max_m
      else:
        x = np.abs(self.clr_iterations/self.step_size - 2*cycle + 1)
        return self.max_m - (self.max_m-self.base_m)*np.maximum(0,(1-x))

    def on_train_begin(self, logs={}):
      logs = logs or {}
      if self.clr_iterations == 0:
        K.set_value(self.model.optimizer.lr, self.base_lr)
      else:
        K.set_value(self.model.optimizer.lr, self.clr())
      if self.cyclical_momentum == True:
        if self.clr_iterations == 0:
          K.set_value(self.model.optimizer.momentum, self.cm())
        else:
          K.set_value(self.model.optimizer.momentum, self.cm())

    def on_batch_begin(self, batch, logs=None):
      logs = logs or {}
      self.trn_iterations += 1
      self.clr_iterations += 1
      self.history.setdefault('lr', []).append(K.get_value(self.model.optimizer.lr))
      self.history.setdefault('iterations', []).append(self.trn_iterations)
      if self.cyclical_momentum == True:
        self.history.setdefault('momentum', []).append(K.get_value(self.model.optimizer.momentum))
      for k, v in logs.items():
        self.history.setdefault(k, []).append(v)
      K.set_value(self.model.optimizer.lr, self.clr())
      if self.cyclical_momentum == True:
        K.set_value(self.model.optimizer.momentum, self.cm())

def onehot_seq(seq):
    letter_to_index =  {'A':0, 'a':0,
                        'C':1, 'c':1,
                        'G':2, 'g':2,
                        'T':3, 't':3}
    to_return = np.zeros((len(seq),4), dtype='int8')
    for idx,letter in enumerate(seq):
      if letter in letter_to_index:
        to_return[idx,letter_to_index[letter]] = 1
    return to_return

def encode_sequence(fasta_pos, fasta_neg, shuffleOff = True):
    x_pos = np.array([onehot_seq(seq) for seq in tqdm(SeqIO.parse(fasta_pos, "fasta"))] +
    [onehot_seq(seq.reverse_complement()) for seq in tqdm(SeqIO.parse(fasta_pos, "fasta"))])
    x_neg = np.array([onehot_seq(seq) for seq in tqdm(SeqIO.parse(fasta_neg, "fasta"))] +
    [onehot_seq(seq.reverse_complement()) for seq in tqdm(SeqIO.parse(fasta_neg, "fasta"))])
    # concatenate positives and negatives
    print(f'There are {x_pos.shape[0]} positives and {x_neg.shape[0]} negatives.')
    x = np.expand_dims(np.concatenate((x_pos, x_neg)), axis=3)
    y = np.concatenate((np.ones(len(x_pos)),np.zeros(len(x_neg))))
    # need to shuffle order of training set for validation splitting last
    if not shuffleOff:
      indices = np.arange(y.shape[0])
      np.random.shuffle(indices)
      x = x[indices,:]
      y = y[indices]
    return x, y

def macro_f1(y, y_hat, thresh=0.5):
    """Compute the macro F1-score on a batch of observations (average F1 across labels)
    Args:
        y (int32 Tensor): labels array of shape (BATCH_SIZE, N_LABELS)
        y_hat (float32 Tensor): probability matrix from forward propagation of shape (BATCH_SIZE, N_LABELS)
        thresh: probability value above which we predict positive

    Returns:
        macro_f1 (scalar Tensor): value of macro F1 for the batch
    """
    y_pred = tf.cast(tf.greater(y_hat, thresh), tf.float32)
    tp = tf.cast(tf.math.count_nonzero(y_pred * y, axis=0), tf.float32)
    fp = tf.cast(tf.math.count_nonzero(y_pred * (1 - y), axis=0), tf.float32)
    fn = tf.cast(tf.math.count_nonzero((1 - y_pred) * y, axis=0), tf.float32)
    f1 = 2*tp / (2*tp + fn + fp + 1e-16)
    macro_f1 = tf.reduce_mean(f1, axis=-1)
    return macro_f1

def get_model(input_shape, hp):
    model = Sequential()
    # Convolutional Layers
    for i in range(len(hp["sizes"])):
      model.add(Conv2D(filters = hp["sizes"][i][0], kernel_size = (11,4 if i == 0 else 1),
                      activation='relu', kernel_initializer='he_normal',
                      bias_initializer='he_normal', kernel_regularizer = l2(l=hp["l2_reg"]),
                      input_shape=input_shape))
      model.add(Dropout(rate = hp["sizes"][i][1]))
    # Pool
    model.add(MaxPooling2D(pool_size=(26,1), strides=26))
    model.add(Flatten())
    # Linear Layers
    for _ in range(1):
      model.add(Dense(units = hp["dense_filters"], activation = 'relu',
                      kernel_initializer='he_normal', bias_initializer='he_normal',
                      kernel_regularizer = l2(l=hp["l2_reg"])
                      #kernel_constraint = tf.keras.constraints.MaxNorm(max_value=4, axis=0)
                      )
      )
      model.add(Dropout(rate = hp["dropout"]))

    # output layer
    model.add(Dense(units = 1, activation = 'sigmoid',
                  kernel_initializer='he_normal', bias_initializer='he_normal',
                  kernel_regularizer = l2(l=hp["l2_reg"])))
    myoptimizer = SGD(lr=hp["base_lr"], momentum=hp["max_m"])
    #myoptimizer = Nadam(lr=hp["base_lr"])
    model.compile(optimizer=myoptimizer,
                  loss="binary_crossentropy",
                  #loss=f1_weighted,
                  metrics=[
                          TruePositives(name='TP'),FalsePositives(name='FP'),
                          TrueNegatives(name='TN'), FalseNegatives(name='FN'),
                          AUC(name='auroc', curve='ROC'), AUC(name='auprc', curve='PR')
                          ]
                  )
    model.summary()
    return model

def train_model_clr(x_train, y_train, x_valid, y_valid, x_human_valid, y_human_valid, evals, hp):
    # dealing w/ class imbalance
    total = y_train.shape[0]
    weight_for_0 = (1 / np.sum(y_train==0))*(total)/2.0
    weight_for_1 = (1 / np.sum(y_train==1))*(total)/2.0
    class_weight = {0: weight_for_0, 1: weight_for_1}
    # An epoch is calculated by dividing the number of training images by the batchsize
    iterPerEpoch = y_train.shape[0] / hp["batch"] # unit is iter
    # number of training iterations per half cycle.
    # Authors suggest setting step_size = (2-8) x (training iterations in epoch)
    iterations = list(range(0,round(y_train.shape[0]/hp["batch"]*hp["epoch"])+1))
    step_size = len(iterations)/(hp["n_cycles"])
    #
    # set cyclic learning rate
    scheduler =  CyclicLR(base_lr=hp["base_lr"],
                max_lr=hp["max_lr"],
                step_size=step_size,
                max_m=hp["max_m"],
                base_m=hp["base_m"],
                cyclical_momentum=True)
    '''
    scheduler = tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss',
                                                factor=0.1,
                                                patience=3,
                                                verbose=0,
                                                mode='auto',
                                                min_delta=0.0001,
                                                cooldown=0,
                                                min_lr=hp["base_lr"],
                                              )
    '''
    if hp["model"]:
      model = load_model(hp["model"])
    else:
      model = get_model(x_train.shape[1:], hp)
    model_checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(filepath=hp["data"]+"weights.h5",
                                                                   save_weights_only=False,
                                                                   monitor='val_loss',
                                                                   mode='min',
                                                                   save_best_only=True)
    if hp["k"]:
      kfold = KFold(n_splits=hp["k"], shuffle=True)
      fold_no = 1
      acc_per_fold, loss_per_fold = list(), list()
      for train, test in kfold.split(x_train, y_train):
        print('------------------------------------------------------------------------')
        print(f'Training for fold {fold_no} ...')
        hist = model.fit(x_train[train],
                        y_train[train],
                        batch_size = hp["batch"],
                        epochs = hp["epoch"],
                        verbose = 1,
                        class_weight = class_weight,
                        validation_data=(x_valid, y_valid),
                        callbacks = [scheduler,
                                    model_checkpoint_callback,
                                    Metrics(validation=(x_valid, y_valid),
                                            human_validation=(x_human_valid, y_human_valid),
                                            evals=evals)])
        scores = model.evaluate(x_train[test], y_train[test], verbose=0)
        print(f'Score for fold {fold_no}: {model.metrics_names[0]} of {scores[0]}; {model.metrics_names[5]} of {scores[5]}; {model.metrics_names[6]} of {scores[6]}')
        acc_per_fold.append(scores[1] * 100)
        loss_per_fold.append(scores[0])
        fold_no += 1
      return model, scheduler, hist, acc_per_fold, loss_per_fold
    else:
      hist = model.fit(x_train,
                        y_train,
                        batch_size = hp["batch"],
                        epochs = hp["epoch"],
                        verbose = 1,
                        class_weight = class_weight,
                        validation_data=(x_valid, y_valid),
                        callbacks = [scheduler,
                                    model_checkpoint_callback,
                                    Metrics(validation=(x_valid, y_valid),
                                            human_validation=(x_human_valid, y_human_valid),
                                            evals=evals)])
      return model, scheduler, hist, [], []

def load_data(hp, load=True):
  TRAIN_FASTA_POS=hp["data"]+"/pos_TRAINING.fa"
  TRAIN_FASTA_NEG=hp["data"]+"/neg_TRAINING.fa"
  VAL_FASTA_POS=hp["data"]+"/pos_VALIDATION.fa"
  VAL_FASTA_NEG=hp["data"]+"/neg_VALIDATION.fa"
  if (os.path.isfile(hp["data"]+"train_fasta.npy")):
    # load
    print("loading training data...")
    x_train = np.load(hp["data"]+"train_fasta.npy", allow_pickle=True)
    y_train = np.load(hp["data"]+"train_labels.npy", allow_pickle=True)
    print("loading validation data...")
    x_valid = np.load(hp["data"]+"val_fasta.npy", allow_pickle=True)
    y_valid = np.load(hp["data"]+"val_labels.npy", allow_pickle=True)
  else:
    # encode and save
    print("loading training data...")
    (x_train, y_train) = encode_sequence(TRAIN_FASTA_POS, TRAIN_FASTA_NEG, shuffleOff = False)
    print("loading validation data...")
    (x_valid, y_valid) = encode_sequence(VAL_FASTA_POS, VAL_FASTA_NEG, shuffleOff = True)
    np.save(hp["data"]+"train_fasta.npy", x_train, allow_pickle=True)
    np.save(hp["data"]+"train_labels.npy", y_train, allow_pickle=True)
    np.save(hp["data"]+"val_fasta.npy", x_valid, allow_pickle=True)
    np.save(hp["data"]+"val_labels.npy", y_valid, allow_pickle=True)
  return x_train, y_train, x_valid, y_valid

def main(hp):
  model_name = hp["name"] if hp["name"] else "tmp.h5"
  # load data
  x_train, y_train, x_valid, y_valid = load_data(hp)
  human_valid_path ="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/human4/"
  mouse_valid_path = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/"
  x_human_valid, y_human_valid = encode_sequence(human_valid_path+"pos_VALIDATION.fa", human_valid_path+"neg_VALIDATION.fa", shuffleOff = False)
  x_mouse_valid, y_mouse_valid = encode_sequence(mouse_valid_path+"pos_VALIDATION.fa", mouse_valid_path+"neg_VALIDATION.fa", shuffleOff = False)
  mouse_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/"
  human_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/"
  eval1 = mouse_data+"mouse_specific_mm10.fa"
  eval2 = mouse_data+"mouse_specific_hg38.fa"
  eval3 = human_data+"human_specific_hg38.fa"
  eval4 = human_data+"human_specific_mm10.fa"
  eval5 = mouse_data+"mm10_ret_nooverlap_liver.fa"
  eval6 = mouse_data+"mm10_ret_overlap_liver.fa"
  eval7 = mouse_data+"mm10_ret_noTSS_filtered_500bp_neg_overlap_liver.fa"

  eval1x, eval1y = get_filtered_input(eval1)
  eval2x, _ = get_filtered_input(eval2)
  eval2y = np.zeros(len(eval2x))
  eval3x, eval3y = get_filtered_input(eval3)
  eval4x, _ = get_filtered_input(eval4)
  eval4y = np.zeros(len(eval4x))
  eval5x, eval5y = get_filtered_input(eval5)
  eval6x, eval6y = get_filtered_input(eval6)
  eval7x, _ = get_filtered_input(eval7)
  eval7y = np.zeros(len(eval7x))
  evals = (eval1x, eval1y, eval2x, eval2y, eval3x, eval3y, eval4x, eval4y, eval5x, eval5y, eval6x, eval6y, eval7x, eval7y)
  K.clear_session()
  # train model
  model, clr, hist, acc_per_fold, loss_per_fold = train_model_clr(x_train, y_train,
                                                                  x_valid, y_valid,
                                                                  x_human_valid, y_human_valid,
                                                                  evals,
                                                                  hp)
  model.save(hp["data"]+ model_name)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-d", "--data", type=str, help="Directory of split input sequences")
  parser.add_argument("-b", "--batch", type=int, help="Batch size")
  parser.add_argument("-e", "--epoch", type=int, help="Number of epochs")
  parser.add_argument("-k", "--kfold", type=int, help="Number of folds in k-fold cross validation")
  parser.add_argument("-m", "--model", type=str, help="Pre-trained TensorFlow model")
  parser.add_argument("-n", "--name", type=str, help="Output model name")


  options, args = parser.parse_known_args()

  if (len(sys.argv)==1):
    parser.print_help(sys.stderr)
    sys.exit(1)
  elif (options.data is None or
        options.batch is None or
        options.epoch is None or
        (options.kfold and options.kfold < 1)
        ):
    parser.print_help(sys.stderr)
    sys.exit(1)
  else:
    # lr 1e-3 - 0.01
    hp = {
          "data" : options.data, "batch" : options.batch, "epoch" : options.epoch,
          "base_lr" : 1e-5, "max_lr" : 0.1, "base_m": 0.875, "max_m": 0.99,
          "l2_reg" : 1e-6, "n_cycles" : 2.5,
          "sizes" : [(256, 0.4), (128, 0.3), (64, 0.2)], "dense_filters" : 300,
          "dropout" : 0.3,
          "k" : options.kfold,
          "model" : options.model, "name" : options.name
          }
    main(hp)
