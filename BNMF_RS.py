import scipy.io as sio



#x = sio.loadmat('init_data.mat').get('X')
#R = sio.loadmat('init_data.mat').get('R')
#M = sio.loadmat('init_data.mat').get('M')
#test_vals = sio.loadmat('init_data.mat').get('test_vals')
#train_vals = sio.loadmat('init_data.mat').get('train_vals')
#users = sio.loadmat('init_data.mat').get('users')
#folder_id = 45

I=50
(W,K) = np.shape(R)


a_tm = 10*np.ones((W, I))
b_tm = np.ones((W, I))  
a_ve = 1*np.ones((I, K))
b_ve = np.ones((I, K))

g = gnmf_vb_poisson_mult_fast(R, R, a_tm, b_tm, a_ve, b_ve,tie_a_tm='free' ,tie_b_tm='free' ,tie_a_ve='free' ,tie_b_ve='free' , UPDATE = 30, EPOCH=100,print_period = 30)

