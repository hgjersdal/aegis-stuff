import random
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


class track:
    """Simulated track, starting in (0,0,0)"""
    x = 0.0
    y = 0.0
    z = 0.0
    fs = 0.0

    def __init__(self, mevs):
        """
        Init a helix to random dip-angle, charge and phi.
        """
        self.ebeam = mevs/1000.0
        # Random angle smaller than 10deg
        lmbda_loss = random.random() * 0.17 * 2
        self.lmbda = math.pi * 0.5 - lmbda_loss
        self.phi = random.random() * math.pi * 2
        self.chrge = 1
        if (random.random() > 0.5):
            self.chrge = -1
        self.R = self.ebeam * math.sin(self.lmbda) / 0.3

    def propagate(self, s):
        """
        Propagate track by s meters. 
        Incremental change in s from last posittion.
        """
        self.fs += s
        phi1 = self.phi + self.chrge * self.fs * math.cos(self.lmbda) / self.R
        self.x = self.R * (math.cos(phi1) - math.cos(self.phi))
        self.y = self.R * (math.sin(phi1) - math.sin(self.phi))
        self.z = self.fs * math.sin(self.lmbda)

    def find_z(self, to_z):
        """
        Dumb way of propagating to given z_position
        """
        while(True):
            z_diff = to_z - self.z
            if(math.fabs(z_diff < 1e-5)):
                return()
            if(self.z > 1e-4):
                self.propagate(z_diff)
            elif(z_diff > 0):
                self.propagate(1e-6)
            else:
                self.propagate(-1e-6)


def plot_track(mevs, lmbda=False):
    trck = track(mevs)
    if(lmbda):
        trck.lmbda = lmbda
    xes = []
    ys = []
    zs = []
    counter = 0
    while trck.z < 0.08 and (counter < 10000):
        xes.append(trck.x * 1e6)
        ys.append(trck.y * 1e6)
        zs.append(trck.z * 1e6)
        # trck.find_z(i * 0.01)
        trck.propagate(0.001)
        counter += 1

    print("Counter ", counter)

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(xes, ys, zs)
    plt.show()


def make_track_file(fname, nTracks, mevs, planes=[0, 0.005, 0.04, 0.08]):
    """
    Make a file that can be read by the track fitter.

    one line is the y-value for a track at z-positions defined by 'planes',
    separated by ';'
    """
    with open(fname, "w") as f:
        for i in range(nTracks):
            trck = track(mevs)
            for z_pos in planes:
                trck.find_z(z_pos)
                f.write(str(int(trck.y * 1e6)) + ";")  # print y-val in um
            f.write("\n")


plot_track(50, math.pi/100)
# make_track_file("/tmp/test_out.csv", 10000, 50)
