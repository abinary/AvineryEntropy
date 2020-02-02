
import numpy, glob, pickle, scipy, scipy.io

folder = "d:\\Files\\Workspace\\Sosnick\\201807\\nug2\\"
filepaths = glob.glob(folder + "*.pkl")
data = []

for f in filepaths:
    print(f)

    with open(f, 'rb') as fp:
        data = pickle.load(fp, encoding='latin1')
        #numpy.savetxt(f + ".txt", data)
        scipy.io.savemat(f + ".mat", dict(data=data))
        


