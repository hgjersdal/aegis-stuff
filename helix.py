import random
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy
import subprocess


class track:
    """Simulated track, starting in (0,0,0)"""
    # width and hright of detector planes
    max_x = 0.09
    min_x = 0.0
    max_y = 0.05
    min_y = 0.0

    def __init__(self, mevs=0.0):
        """
    Init a helix to random dip-angle, charge and phi.
        """
        if(mevs < 10.0):
            self.ebeam = self.rand_mevs()
        else:
            self.ebeam = mevs
        # Dip angle is random, within 90degrees.
        # In reality it should be 180, so estimated acceptance is 2x.
        # self.lmbda = random.random() * math.pi * 0.5
        self.lmbda = math.pi * 0.5 - math.acos(random.random())
        self.phi = random.random() * math.pi * 2
        # pi + or -, 50-50 chance.
        self.chrge = 1
        if (random.random() > 0.5):
            self.chrge = -1
        # Tracks start out sonewhere on plane 0
        self.x = random.random() * self.max_x
        self.y_pos()
        # self.y = random.random() * self.max_y
        self.z = 0.0

    def y_pos(self):
        """ Find random y position, fringe like"""
        pitch = 150e-6
        period = pitch/2.0
        while(True):
            p = random.random()
            self.y = random.random() * self.max_y
            limit = 0.5 + 0.5 * math.sin((math.pi * self.y)/period)
            if (limit < p):
                return()

    def rand_mevs(self):
        while(True):
            selector = random.random()
            position = random.random()
            if(position < selector):
                return(50 + position * 450) 

    def scatter(self):
        """ Scattering in detector planes, 50 um thick. Highland formula. """
        pion_mass = 139.0
        radlength = 50e-6/0.1  # 50um / 10cm
        beta = self.ebeam/math.sqrt(self.ebeam * self.ebeam + pion_mass * pion_mass)
        theta = 13.6/(beta * self.ebeam) * math.sqrt(radlength) * (1.0 + 0.038 * math.log(radlength))
        self.lmbda += numpy.random.normal(0, theta)
        self.phi += numpy.random.normal(0, theta)

    def propagate_destructive(self, s, largeR=False):
        self.R = self.ebeam * math.cos(self.lmbda) / (0.3 * 1000)
        if(largeR):
            self.R = 100000
        phi1 = self.phi + (self.chrge * s * math.cos(self.lmbda) / self.R)
        self.x += self.R * (math.cos(phi1) - math.cos(self.phi))
        self.y += self.R * (math.sin(phi1) - math.sin(self.phi))
        self.z += s * math.sin(self.lmbda)
        self.phi = phi1

    def find_z(self, to_z, largeR=False):
        """
        Dumb way of propagating to given z_position
        Return true if within bounds, false if outside
        """
        count = 0
        while(True):
            if(count > 1000):
                # print("Took too long to propagate!!!")
                return(False)
            count += 1

            z_diff = to_z - self.z
            if(math.fabs(z_diff < 1e-5)):  # Close enough, stop looking
                if((self.max_x > self.x > self.min_x) and
                   (self.max_y > self.y > self.min_y)):
                    return(True)
                else:
                    return(False)
            if(self.z > 1e-4):
                self.propagate_destructive(z_diff, largeR)
            elif(z_diff > 0):
                self.propagate_destructive(1e-6, largeR)
            else:
                self.propagate_destructive(-1e-6, largeR)


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
        trck.propagate_destructive(0.001)
        counter += 1

    print("Counter ", counter)

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(zs, xes, ys)
    plt.show()


def plot_track_planes(mevs,  planes=[0, 0.005, 0.04, 0.08], largeR=False):
    """ Propagation with scattering in planes """
    xes = []
    ys = []
    zs = []

    trck = track(mevs)
    for z_pos in planes:
        if(trck.find_z(z_pos, largeR)):
            xes.append(trck.x * 1e6)
            ys.append(trck.y * 1e6)
            zs.append(trck.z * 1e6)
        else:
            print("Out of bounds!")
            break

    if(len(xes) == 4):
        print("In bounds!!!")
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(xes, ys, zs)
    plt.show()


def make_track_scatter(fname, nTracks, mevs=0.0, planes=[0, 0.005, 0.02, 0.04], largeR=False):
    """
    Make a file that can be read by the track fitter.

    one line is the y-value for a track at z-positions defined by 'planes',
    separated by ';'

    Scattering in detector planes. 25um/sqrt(12) measurement resolution.

    First plane is true decay position, the rest are simulated measurements
    """
    with open(fname, "w") as f:
        accepted_tracks = 0
        lambdas = []
        phis = []
        for i in range(nTracks):
            ys = []
            trck = track(mevs)
            lambdas.append(trck.lmbda)
            phis.append(trck.phi)
            for z_pos in planes:
                if(trck.find_z(z_pos, largeR)):
                    val = trck.y * 1e6
                    # Measurement uncertainties for all planes except the first
                    if(z_pos > 0.0001):
                        val += numpy.random.normal(0, 25.0/math.sqrt(12.0))
                    ys.append(val)
                    trck.scatter()
                else:
                    break
            if(len(ys) == 4):
                accepted_tracks += 1
                for ii in range(len(planes)):
                    f.write(str(int(ys[ii])) + ";")  # print y-val in um
                f.write("\n")
        print("Accepted: ", accepted_tracks)
        # plt.scatter(lambdas, phis)
        # plt.show()


# plot_track(50, lmbda=math.pi/100.0)

print("simulating")
make_track_scatter("/tmp/test_out.csv", 100000, mevs=0, largeR=True)
print("cmopiling")
subprocess.run(['make', 'simple'], cwd='/home/haavagj/git/eigen-track-fitter')
print("running")
subprocess.run(['./simple'], cwd='/home/haavagj/git/eigen-track-fitter')
