import os
import sys
import time
import random
import shutil
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Embedding, Bidirectional, LSTM, Dense
from tensorflow.keras.callbacks import ModelCheckpoint, CSVLogger
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.utils import to_categorical

def Cla_LSTM_debug():
    input_layer = Input(shape=(61,))
    embedded = Embedding(input_dim=5, output_dim=5)(input_layer)
    x = Bidirectional(LSTM(32))(embedded)
    x = Dense(32, activation='relu')(x)
    output = Dense(2, activation='softmax')(x)
    return Model(inputs=input_layer, outputs=output)

if __name__ == '__main__':
    np.random.seed(1122)
    random.seed(1122)
    tf.random.set_seed(1122)

    # === Parse args
    npydir = sys.argv[1]
    weightfile = sys.argv[2]
    ignore1 = sys.argv[3]
    ignore2 = sys.argv[4]

    # === Load data
    Good_for_Tra = np.load(os.path.join(npydir, 'Good_for_Tra.npy'))
    Simu_for_Tra = np.load(os.path.join(npydir, 'Simu_for_Tra.npy'))
    Good_for_Tst = np.load(os.path.join(npydir, 'Good_for_Tst.npy'))
    Simu_for_Tst = np.load(os.path.join(npydir, 'Simu_for_Tst.npy'))

    Tra_x = np.squeeze(np.concatenate((Good_for_Tra, Simu_for_Tra), axis=0))
    Tra_y = np.concatenate((np.zeros((Good_for_Tra.shape[0], 1)),
                            np.ones((Simu_for_Tra.shape[0], 1))), axis=0)
    Tst_x = np.squeeze(np.concatenate((Good_for_Tst, Simu_for_Tst), axis=0))
    Tst_y = np.concatenate((np.zeros((Good_for_Tst.shape[0], 1)),
                            np.ones((Simu_for_Tst.shape[0], 1))), axis=0)

    Tra_y = to_categorical(Tra_y)
    Tst_y = to_categorical(Tst_y)

    # Shuffle
    idx = np.random.permutation(Tra_x.shape[0])
    Tra_x, Tra_y = Tra_x[idx], Tra_y[idx]
    idx = np.random.permutation(Tst_x.shape[0])
    Tst_x, Tst_y = Tst_x[idx], Tst_y[idx]

    # === Build model
    model = Cla_LSTM_debug()
    model.compile(loss='binary_crossentropy', optimizer=Adam(learning_rate=0.001), metrics=['accuracy'])

    # === Load weights
    if os.path.isfile(weightfile):
        print(f"🔄 Loading initial weights from: {weightfile}", flush=True)
        model.load_weights(weightfile)

    # === Use /tmp for safe saving
    tmp_h5 = f"/tmp/DebugWeight-{random.randint(1000,9999)}.hdf5"
    tmp_log = f"/tmp/debug_log_{random.randint(1000,9999)}.csv"

    model_checkpoint = ModelCheckpoint(filepath=tmp_h5, verbose=1,
                                       monitor='val_loss', save_best_only=True)
    csv_logger = CSVLogger(tmp_log, append=True, separator=';')

    # === Train
    print("🚀 Training on 100 samples (debug mode)", flush=True)
    model.fit(
        x=Tra_x[:100], y=Tra_y[:100],
        batch_size=16,
        epochs=2,
        validation_data=(Tst_x[:50], Tst_y[:50]),
        verbose=1,
        shuffle=True,
        callbacks=[model_checkpoint, csv_logger]
    )

    # === Copy results to npydir
    final_h5 = os.path.join(npydir, "DebugWeight-001.hdf5")
    final_log = os.path.join(npydir, "debug_log.csv")

    print(f"📁 Copying results to: {npydir}", flush=True)
    shutil.copy(tmp_h5, final_h5)
    shutil.copy(tmp_log, final_log)

    print("✅ Training complete. All weights and logs saved.")
