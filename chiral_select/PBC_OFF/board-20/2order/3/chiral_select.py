import numpy as np
import math
import copy
import matplotlib.pyplot as plt

def initiate_matrix(dim, l_ratio, r_ratio):
    # 0 for achiral, 1 for left, 2 for right
    num_row = dim[0]
    num_col = dim[1]
    
    total_num = num_row * num_col
    num_left = int(total_num * l_ratio)
    num_right = int(total_num * r_ratio)
    
    matrix = np.zeros(dim)
    positions = range(total_num)
    np.random.shuffle(positions)
    for i, pos in enumerate(positions[: num_left+num_right]):
        row = pos / num_col
        col = pos - row * num_col
        if i < num_left:
            matrix[row, col] = 1
        else:
            matrix[row, col] = 2
    return matrix

def calc_nn_env(matrix, pos, PBC):
    # return (num_achiral, num_left, num_right)
    nn_env = np.zeros(3)
    num_row = matrix.shape[0]
    num_col = matrix.shape[1]
    row = pos[0]
    col = pos[1]
    if PBC == True:
        nn_env[matrix[(row+num_row+1)%num_row, col]] += 1
        nn_env[matrix[(row+num_row-1)%num_row, col]] += 1
        nn_env[matrix[row, (col+num_col+1)%num_col]] += 1
        nn_env[matrix[row, (col+num_col-1)%num_col]] += 1
    elif PBC == False:
        if row+1 < num_row:
            nn_env[matrix[row+1, col]] += 1
        if row-1 >= 0:
            nn_env[matrix[row-1, col]] += 1
        if col+1 < num_col:
            nn_env[matrix[row, col+1]] += 1
        if col-1 >= 0:
            nn_env[matrix[row, col-1]] += 1
    return nn_env

def transition_one(matrix, pos, prob_params, PBC):
    # return list of [pos, next_state, rate]
    transition = []
    k0 = prob_params['k0']
    k1 = prob_params['k1']
    k2 = prob_params['k2']
    gamma = prob_params['gamma']
    
    if matrix[pos[0], pos[1]] == 1 or matrix[pos[0], pos[1]] == 2:
        transition.append([pos, 0, gamma])
        return transition
    if matrix[pos[0], pos[1]] == 0:
        nn_achiral, nn_left, nn_right = calc_nn_env(matrix, pos, PBC)
        if nn_left == 0:
            transition.append([pos, 1, k0])
        elif nn_left == 1:
            transition.append([pos, 1, k0+k1])
        elif nn_left >= 2:
            transition.append([pos, 1, k0+k1+k2])
        if nn_right == 0:
            transition.append([pos, 2, k0])
        elif nn_right == 1:
            transition.append([pos, 2, k0+k1])
        elif nn_right >= 2:
            transition.append([pos, 2, k0+k1+k2])
    return transition


def transition_all(matrix, prob_params, PBC):
    transition = []
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            transition.extend(transition_one(matrix, [row, col], prob_params, PBC))
    return transition

def jump(matrix, prob_params, PBC):
    # return [pos, state, time]
    transition = transition_all(matrix, prob_params, PBC)
    rates = [i[2] for i in transition]
    cum_rates = copy.copy(rates)
    for i in range(1, len(rates)):
        cum_rates[i] += cum_rates[i-1]
    rn = np.random.rand() * cum_rates[-1]
    selected_trans = np.searchsorted(cum_rates, rn) # binary search for rn in cum_rates
    time = prob_params['k0'] / cum_rates[-1]
    return transition[selected_trans][0], transition[selected_trans][1], time

def print_MC(step, matrix, time):
    num_left = sum(sum(matrix == 1))
    num_right = sum(sum(matrix == 2))
    num_achiral = matrix.shape[0]*matrix.shape[1] - num_left - num_right
    with open('MC.log', 'a') as fp:
        fp.write('%4d%16.9f%10d%10d%10d\n'%(step, time, num_achiral, num_left, num_right))
    np.savetxt('%04d.conf'%step, matrix, fmt='%d')
    
    # fig, ax = plt.subplots()
    # ax.imshow(matrix, cmap=plt.cm.gray, interpolation='none')
    # ax.set_title('%04d'%step)
    # plt.show()
    # plt.savefig('%04d.png'%step, bbox_inches = 0, dpi = 400)
    # plt.close()

        

def MC(matrix, prob_params, PBC, steps, out_step):
    with open('MC.log', 'w') as fp:
        fp.write('step\ttime\tachiral\tleft\tright\n')
    cum_time= 0.0
    for i in range(steps):
        (pos, state, time) = jump(matrix, prob_params, PBC)
        matrix[pos[0], pos[1]] = state
        cum_time += time
        if i % out_step == 0:
            print_MC(i/out_step, matrix, cum_time)
    return matrix, cum_time

def main():
    dim = [20, 20]
    l_ratio, r_ratio = 0.0, 0.0
    prob_params = {'k0': 0.1, 'k1': 0.0, 'k2': 20, 'gamma': 0.1}
    PBC = False
    steps = 40000
    out_step = 100
    matrix = initiate_matrix(dim, l_ratio, r_ratio)
    matrix, cum_time = MC(matrix, prob_params, PBC, steps, out_step)
    
if __name__ == '__main__':
    main()
