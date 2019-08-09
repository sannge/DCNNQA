import config
import subprocess
import os
import numpy as np
import scipy.stats
import random

import tensorflow as tf

import load_data
import model
import argparse

CONF = config.load_config()
FLAGS = None

def predict(sess, maps_placeholder, is_training, logits, filename):

    # mapping protein
    print('# Scoring '+filename)
    mapFilename = CONF.TEMP_PATH + 'map_'+str(os.getpid())+ '_pred.bin'
    subprocess.call(
      [CONF.MAP_GENERATOR_PATH, "--mode", "map", "-i", filename, "--native", "-m", "24", "-v", "0.8", "-o",
      mapFilename])
    if not os.path.exists(mapFilename):
        print('# Mapping failed, ignoring protein')
        return None
    predDataset = load_data.read_data_set(mapFilename)
    os.remove(mapFilename)


    preds = []
    # compute prediction res by res
    for i in range(predDataset.num_res):
        f_map = np.reshape(predDataset.maps[i], (1, model.GRID_VOXELS * model.NB_TYPE))
        feed_dict = {maps_placeholder: f_map, is_training: False}
        pred = sess.run(logits,feed_dict=feed_dict)
        preds.append(pred)
        outline='RES {:4d} {:c} {:5.4f}'.format(predDataset.meta[i][0], predDataset.meta[i][1], pred)
        print(outline)
        #print(predDataset.meta[i][0]+)
        #print(pred)

def main():
    sess = tf.Session()
    print('Restore existing model: %s' % CONF.MODEL_PATH)
    saver = tf.train.import_meta_graph(CONF.MODEL_PATH + '.meta')
    saver.restore(sess, CONF.MODEL_PATH)

    graph = tf.get_default_graph()

    # getting placeholder for input data and output
    maps_placeholder = graph.get_tensor_by_name('main_input:0')
    is_training = graph.get_tensor_by_name('is_training:0')
    logits = graph.get_tensor_by_name("main_output:0")

    if FLAGS.structure != None :
        predict(sess, maps_placeholder, is_training, logits, FLAGS.structure)
    if FLAGS.directory != None :
        for filename in os.listdir(FLAGS.directory):
            predict(sess, maps_placeholder, is_training, logits, FLAGS.directory+'/'+filename)






if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-d',
    '--directory',
    type=str,
    help='Path to the validation data'
  )
  parser.add_argument(
    '-s',
    '--structure',
    type=str,
    help='Path to the structure to score (in pdb format)'
  ) 
  FLAGS = parser.parse_args()
  main()


