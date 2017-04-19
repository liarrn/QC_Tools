import numpy as np
import matplotlib.pyplot as plt

def plotBoard(matrix, name):
    fig, ax = plt.subplots()
    ax.imshow(matrix, cmap=plt.cm.gray, interpolation='none')
    ax.set_title(name)
    # plt.show()
    plt.savefig(name+'.png', bbox_inches = 0, dpi = 400)
    plt.close()
    
if __name__ == '__main__':
    prefix = './PBC_OFF/'
    conf_range = range(90)
    for i in conf_range:
        matrix = np.loadtxt(prefix+'%04d.conf'%i)
        plotBoard(matrix, '%04d.conf'%i)
    # filename = './PBC_OFF/0088.conf'
    # matrix = np.loadtxt(filename)
    # plotBoard(matrix)