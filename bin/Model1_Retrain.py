import os
import sys
import time
import random
import shutil
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Embedding, Bidirectional, LSTM, Dense, Dropout, Activation
from tensorflow.keras.callbacks import ModelCheckpoint, CSVLogger
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.utils import to_categorical
from glob import glob

def Cla_LSTM_full():
    input_layer = Input(shape=(61,))
    x = Embedding(5, 5, input_length=61)(input_layer)
    x = Bidirectional(LSTM(32, return_sequences=True))(x)
    x = Dropout(0.5)(x)
    x = Bidirectional(LSTM(64, return_sequences=True))(x)
    x = Dropout(0.5)(x)
    x = Bidirectional(LSTM(128, return_sequences=True))(x)
    x = Dropout(0.5)(x)
    x = Bidirectional(LSTM(256))(x)
    x = Dense(256)(x)
    x = Dense(2)(x)
    output = Activation('softmax')(x)
    return Model(inputs=input_layer, outputs=output)

if __name__ == '__main__':
    np.random.seed(1122)
    random.seed(1122)
    tf.random.set_seed(1122)

    npydir = sys.argv[1]
    weightfile = sys.argv[2]
    epochoutdir = sys.argv[3]
    itere = int(sys.argv[4])

    # Load data
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

    idx = np.random.permutation(Tra_x.shape[0])
    Tra_x, Tra_y = Tra_x[idx], Tra_y[idx]
    idx = np.random.permutation(Tst_x.shape[0])
    Tst_x, Tst_y = Tst_x[idx], Tst_y[idx]

    # Build model
    model = Cla_LSTM_full()
    model.compile(loss='binary_crossentropy', optimizer=Adam(learning_rate=0.0001), metrics=['accuracy'])

    # Load initial weights if given
    if os.path.isfile(weightfile):
        print(f"🔄 Loading initial weights from: {weightfile}")
        model.load_weights(weightfile)

    # Save files to /tmp first
    tmp_prefix = f"/tmp/retrain_{random.randint(1000,9999)}"
    tmp_h5_epoch = tmp_prefix + "-{epoch:03d}.hdf5"
    tmp_h5_final = tmp_prefix + "-final.hdf5"
    tmp_log = tmp_prefix + "-log.csv"

    model_checkpoint_epoch = ModelCheckpoint(tmp_h5_epoch, monitor='val_loss', save_best_only=False, verbose=1)
    model_checkpoint_best = ModelCheckpoint(tmp_h5_final, monitor='val_loss', save_best_only=True, verbose=1)
    csv_logger = CSVLogger(tmp_log, append=True, separator=';')

    print(f"🚀 Training full model for {itere} epochs...")
    model.fit(
        x=Tra_x, y=Tra_y,
        batch_size=64,
        epochs=itere,
        validation_data=(Tst_x, Tst_y),
        verbose=1,
        shuffle=True,
        callbacks=[model_checkpoint_epoch, model_checkpoint_best, csv_logger]
    )

    print(f"📦 Copying files to: {epochoutdir}")
    os.makedirs(epochoutdir, exist_ok=True)

    # Copy per-epoch weights
    for file in glob(tmp_prefix + "-*.hdf5"):
        shutil.copy(file, epochoutdir)

    shutil.copy(tmp_log, os.path.join(epochoutdir, "log.csv"))

    print("✅ Training complete. All files saved to:", epochoutdir)
