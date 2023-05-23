import h5py
import numpy as np
import os

np.random.seed(31337)

OUTPUT_FILE_NAME = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data.h5')

def gen_array(name, shape, func=np.random.rand):
    arr = func(*shape)
    with h5py.File(OUTPUT_FILE_NAME, 'a') as file:
        file.create_dataset(name, data=arr)
    print(f'Generated {name}')

try:
    os.remove(OUTPUT_FILE_NAME)
except OSError:
    pass

gen_array('small_1d_1', (1, 1))
gen_array('small_1d_2', (1, 2))
gen_array('small_1d_3', (1, 4))
gen_array('small_1d_4', (1, 8))

gen_array('large_1d_1', (1, 1024))
gen_array('large_1d_2', (1, 4096))

gen_array('small_2d_1', (1, 1))
gen_array('small_2d_2', (2, 2))
gen_array('small_2d_3', (4, 4))
gen_array('small_2d_4', (8, 8))

gen_array('large_2d_1', (512, 512))
gen_array('large_2d_2', (1024, 2048))

gen_array('small_5d_1', (1, 1, 1, 1, 1))
gen_array('small_5d_2', (2, 2, 2, 2, 2))
gen_array('small_5d_3', (4, 4, 4, 4, 4))
gen_array('small_5d_4', (8, 8, 8, 8, 8))

gen_array('large_5d_1', (32, 32, 8, 16, 16))
gen_array('large_5d_2', (4, 8, 16, 32, 64))
