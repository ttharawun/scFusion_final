# -*- coding: utf-8 -*-
"""
Safe Prediction Script (Huang Wenjian, adjusted)
"""

from keras.models import Sequential
from keras.layers import Embedding, Dropout, Bidirectional, Flatten, Dense, LSTM, TimeDistributed, Activation
from keras.callbacks import ModelCheckpoint, CSVLogger
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D
from tensorflow.keras.optimizers import Adam
import numpy as np
from tensorflow.keras.utils import to_categorical
from Model1 import Cla_LSTM
import os
import sys
import shutil
import random

np.random.seed(1122)

# === Args
final_output_path = sys.argv[1]         # Output file to write
weightfile = sys.argv[2]                # Initial weights
prefix = sys.argv[3] if len(sys.argv) == 4 else ""

# === Derived paths
findgang = final_output_path.rfind('/')
filedir = final_output_path[:findgang+1]
codedir = os.path.dirname(sys.argv[0]) + "/"

# === Data loading
Good_for_Tra = np.load(os.path.join(filedir, prefix + 'Reads.npy'))
Good_for_Tra_rev = np.load(os.path.join(filedir, prefix + 'Reads_rev.npy'))
Tst_x = np.squeeze(Good_for_Tra)
Tst_x_rev = np.squeeze(Good_for_Tra_rev)

# === Build model and load weights from a temp copy
model = Cla_LSTM()

if os.path.isfile(weightfile):
    print(f"🔄 Loading weights from: {weightfile}")
    tmp_weight = f"/tmp/tmp_model_weight_{random.randint(1000,9999)}.h5"
    shutil.copy(weightfile, tmp_weight)
    model.load_weights(tmp_weight)

# === Predict
print("🚀 Predicting probabilities...")
batch_size = 500
Prob = model.predict(Tst_x, batch_size)
Prob_rev = model.predict(Tst_x_rev, batch_size)
AveProb = (Prob[:, 0] + Prob_rev[:, 0]) / 2

# === Write output to temp then move
tmp_output = f"/tmp/tmp_output_{random.randint(1000,9999)}.txt"
with open(tmp_output, 'w') as outfile:
    for i in range(len(AveProb)):
        outfile.write(f"{Prob[i,0]}\t{Prob_rev[i,0]}\t{AveProb[i]}\n")

shutil.copy(tmp_output, final_output_path)
print(f"✅ Output written to: {final_output_path}")
