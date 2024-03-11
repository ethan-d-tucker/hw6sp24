# region imports
import numpy as np
import math
from scipy.optimize import fsolve
import random as rnd
# endregion

# region class definitions
class Fluid:
    # region constructor
    def __init__(self, mu=0.00089, rho=1000):
        """
        Initializes the Fluid with its properties.
        :param mu: dynamic viscosity in Pa*s
        :param rho: density in kg/m^3
        """
        self.mu = mu  # dynamic viscosity in Pa*s
        self.rho = rho  # density in kg/m^3
        self.nu = mu / rho  # kinematic viscosity in m^2/s
    # endregion

class Node:
    # region constructor
    def __init__(self, Name='a', Pipes=[], ExtFlow=0):
        """
        Initializes the Node in the pipe network.
        :param Name: the name of the node
        :param Pipes: a list of pipes connected to this node
        :param ExtFlow: any external flow into (+) or out (-) of this node in L/s
        """
        self.name = Name
        self.pipes = Pipes  # a list of pipes connected to this node
        self.extFlow = ExtFlow  # external flow into (+) or out (-) of this node in L/s
    # endregion

    # region methods/functions
    def getNetFlowRate(self):
        """
        Calculates the net flow rate into this node in m^3/s.
        :return: Net flow rate into the node as a float.
        """
        Qtot = self.extFlow / 1000  # convert L/s to m^3/s and count the external flow first
        for p in self.pipes:
            Qtot += p.getFlowIntoNode(self.name) / 1000  # convert L/s to m^3/s
        return Qtot
    # endregion

class Loop:
    # region constructor
    def __init__(self, Name='A', Pipes=[]):
        """
        Initializes the Loop in the pipe network.
        :param Name: the name of the loop
        :param Pipes: a list/array of pipes in this loop
        """
        self.name = Name
        self.pipes = Pipes
    # endregion

    # region methods/functions
    def getLoopHeadLoss(self):
        """
        Calculates the net head loss as we traverse around the loop, in m of fluid.
        :return: Net head loss as a float.
        """
        deltaP = 0  # initialize to zero
        startNode = self.pipes[0].startNode  # begin at the start node of the first pipe
        for p in self.pipes:
            phl = p.getFlowHeadLoss(startNode)
            deltaP += phl
            startNode = p.endNode if startNode != p.endNode else p.startNode  # move to the next node
        return deltaP
    # endregion

class Pipe:
    # region constructor
    def __init__(self, Start='A', End='B', L=100, D=200, r=0.00025, fluid=Fluid()):
        """
        Initializes the Pipe in the pipe network.
        :param Start: the start node (alphabetically)
        :param End: the end node (alphabetically)
        :param L: the length of the pipe in meters
        :param D: the diameter of the pipe in millimeters
        :param r: the roughness of the pipe in meters
        :param fluid: the fluid flowing through the pipe
        """
        self.startNode = min(Start, End)  # start node (alphabetically)
        self.endNode = max(Start, End)  # end node (alphabetically)
        self.length = L
        self.diameter = D / 1000.0  # diameter in m
        self.r = r
        self.fluid = fluid
        self.A = math.pi * (self.diameter / 2) ** 2  # cross-sectional area
        self.Q = 10  # initial guess for flow rate in L/s
        self.vel = self.V()  # velocity in m/s
        self.reynolds = self.Re()  # Reynolds number
    # endregion

    # region methods/functions
    def V(self):
        """
        Calculate average velocity in the pipe for volumetric flow Q.
        :return: Average velocity as a float.
        """
        self.vel = (self.Q / 1000) / self.A  # convert L/s to m^3/s and calculate velocity
        return self.vel

    def Re(self):
        """
        Calculate the Reynolds number under current conditions.
        :return: Reynolds number as a float.
        """
        self.reynolds = (self.fluid.rho * self.vel * self.diameter) / self.fluid.mu
        return self.reynolds

    def FrictionFactor(self):
        """
        Calculates the friction factor for the pipe flow based on flow regime.
        :return: The Darcy friction factor as a float.
        """
        Re = self.Re()
        rr = self.r / self.diameter

        def CB():
            cb = lambda f: 1 / (f ** 0.5) + 2.0 * np.log10(rr / 3.7 + 2.51 / (Re * f ** 0.5))
            result = fsolve(cb, [0.01])
            return result[0]

        def lam():
            return 64 / Re

        if Re >= 4000:
            return CB()
        elif Re <= 2000:
            return lam()
        else:
            CBff = CB()
            Lamff = lam()
            mean = Lamff + ((Re - 2000) / (4000 - 2000)) * (CBff - Lamff)
            sig = 0.2 * mean
            return rnd.normalvariate(mean, sig)

    def frictionHeadLoss(self):
        """
        Calculate head loss through a section of pipe in m of fluid using the Darcy-Weisbach equation.
        :return: Head loss as a float.
        """
        g = 9.81  # m/s^2
        ff = self.FrictionFactor()
        hl = (ff * self.length / self.diameter) * (self.vel ** 2 / (2 * g))
        return hl

    def getFlowHeadLoss(self, s):
        """
        Calculate the signed head loss for the pipe.
        :param s: the node we're starting with in a traversal of the pipe
        :return: Signed headloss through the pipe in m of fluid.
        """
        nTraverse = 1 if s == self.startNode else -1
        nFlow = 1 if self.Q >= 0 else -1
        return nTraverse * nFlow * self.frictionHeadLoss()

    def Name(self):
        """
        Gets the pipe name.
        :return: Pipe name as a string.
        """
        return f'{self.startNode}-{self.endNode}'

    def oContainsNode(self, node):
        """
        Checks if the pipe is connected to the node.
        :param node: Node to check against.
        :return: True if node is connected, False otherwise.
        """
        return self.startNode == node or self.endNode == node

    def printPipeFlowRate(self):
        """
        Prints the flow rate in the pipe.
        """
        print(f'The flow in segment {self.Name()} is {self.Q:0.4f} m^3/s')

    def getFlowIntoNode(self, n):
        """
        Determines the flow rate into the node.
        :param n: Node to check flow rate for.
        :return: Flow rate into node as a float.
        """
        if n == self.startNode:
            return -self.Q
        return self.Q
    # endregion

class PipeNetwork:
    # region constructor
    def __init__(self, Pipes=[], Loops=[], Nodes=[], fluid=Fluid()):
        """
        Initializes the PipeNetwork with pipes, nodes, loops, and fluid.
        :param Pipes: List of Pipe objects in the network.
        :param Loops: List of Loop objects in the network.
        :param Nodes: List of Node objects in the network.
        :param fluid: Fluid object representing the fluid in the network.
        """
        self.loops = Loops
        self.nodes = Nodes
        self.fluid = fluid
        self.pipes = Pipes
    # endregion

    # region methods/functions
    def findFlowRates(self):
        """
        Analyzes the pipe network and finds the flow rates in each pipe.
        :return: List of flow rates in the pipes.
        """
        N = len(self.nodes) + len(self.loops)
        Q0 = np.full(N, 10)  # initial guess

        def fn(q):
            for i in range(len(self.pipes)):
                self.pipes[i].Q = q[i]  # update flow rate in pipes
            L = self.getNodeFlowRates()  # net flow rates at nodes
            L += self.getLoopHeadLosses()  # net head losses in loops
            return L

        FR = fsolve(fn, Q0)
        return FR

    def getNodeFlowRates(self):
        """
        Calculates net flow rates at nodes.
        :return: List of net flow rates at nodes.
        """
        qNet = [n.getNetFlowRate() for n in self.nodes]
        return qNet

    def getLoopHeadLosses(self):
        """
        Calculates net head losses in loops.
        :return: List of net head losses in loops.
        """
        lhl = [l.getLoopHeadLoss() for l in self.loops]
        return lhl

    def getPipe(self, name):
        """
        Retrieves a pipe object by its name.
        :param name: Name of the pipe to retrieve.
        :return: Pipe object.
        """
        for p in self.pipes:
            if name == p.Name():
                return p

    def getNodePipes(self, node):
        """
        Retrieves a list of pipe objects that are connected to the node.
        :param node: Node to check for connected pipes.
        :return: List of connected Pipe objects.
        """
        l = []
        for p in self.pipes:
            if p.oContainsNode(node):
                l.append(p)
        return l

    def nodeBuilt(self, node):
        """
        Checks if a node object has already been constructed.
        :param node: Name of the node to check.
        :return: True if the node exists, False otherwise.
        """
        for n in self.nodes:
            if n.name == node:
                return True
        return False

    def getNode(self, name):
        """
        Retrieves one of the node objects by name.
        :param name: Name of the node to retrieve.
        :return: Node object.
        """
        for n in self.nodes:
            if n.name == name:
                return n

    def buildNodes(self):
        """
        Automatically creates the node objects based on the pipe ends.
        """
        for p in self.pipes:
            if not self.nodeBuilt(p.startNode):
                self.nodes.append(Node(p.startNode, self.getNodePipes(p.startNode)))
            if not self.nodeBuilt(p.endNode):
                self.nodes.append(Node(p.endNode, self.getNodePipes(p.endNode)))

    def printPipeFlowRates(self):
        """
        Prints the flow rate in each pipe in the network.
        """
        for p in self.pipes:
            p.printPipeFlowRate()

    def printNetNodeFlows(self):
        """
        Prints the net flow into each node in the network.
        """
        for n in self.nodes:
            print(f'net flow into node {n.name} is {n.getNetFlowRate():0.4f} m^3/s')

    def printLoopHeadLoss(self):
        """
        Prints the head loss for each loop in the network.
        """
        for l in self.loops:
            print(f'head loss for loop {l.name} is {l.getLoopHeadLoss():0.4f} m of fluid')
    # endregion

# endregion

# region function definitions
def main():
    '''
    This program analyzes flows in a given pipe network based on the following:
    1. The pipe segments are named by their endpoint node names:  e.g., a-b, b-e, etc.
    2. Flow from the lower letter to the higher letter of a pipe is considered positive.
    3. Pressure decreases in the direction of flow through a pipe.
    4. At each node in the pipe network, mass is conserved.
    5. For any loop in the pipe network, the pressure loss is zero
    Approach to analyzing the pipe network:
    Step 1: build a pipe network object that contains pipe, node, loop and fluid objects
    Step 2: calculate the flow rates in each pipe using fsolve
    Step 3: output results
    Step 4: check results against expected properties of zero head loss around a loop and mass conservation at nodes.
    :return:
    '''
    #instantiate a Fluid object to define the working fluid as water
    water= Fluid()
    roughness = 0.00025  # in meters

    #instantiate a new PipeNetwork object
    PN= PipeNetwork()
    #add Pipe objects to the pipe network (see constructor for Pipe class)
    PN.pipes.append(Pipe('a','b',250, 300, roughness, water))
    PN.pipes.append(Pipe('a','c',100, 200, roughness, water))
    PN.pipes.append(Pipe('b','e',100, 200, roughness, water))
    PN.pipes.append(Pipe('c','d',125, 200, roughness, water))
    PN.pipes.append(Pipe('c','f',100, 150, roughness, water))
    PN.pipes.append(Pipe('d','e',125, 200, roughness, water))
    PN.pipes.append(Pipe('d','g',100, 150, roughness, water))
    PN.pipes.append(Pipe('e','h',100, 150, roughness, water))
    PN.pipes.append(Pipe('f','g',125, 250, roughness, water))
    PN.pipes.append(Pipe('g','h',125, 250, roughness, water))
    #add Node objects to the pipe network by calling buildNodes method of PN object
    PN.buildNodes()

    #update the external flow of certain nodes
    PN.getNode('a').extFlow=60
    PN.getNode('d').extFlow=-30
    PN.getNode('f').extFlow=-15
    PN.getNode('h').extFlow=-15

    #add Loop objects to the pipe network
    PN.loops.append(Loop('A',[PN.getPipe('a-b'), PN.getPipe('b-e'),PN.getPipe('d-e'), PN.getPipe('c-d'), PN.getPipe('a-c')]))
    PN.loops.append(Loop('B',[PN.getPipe('c-d'), PN.getPipe('d-g'),PN.getPipe('f-g'), PN.getPipe('c-f')]))
    PN.loops.append(Loop('C',[PN.getPipe('d-e'), PN.getPipe('e-h'),PN.getPipe('g-h'), PN.getPipe('d-g')]))

    #call the findFlowRates method of the PN (a PipeNetwork object)
    flow_rates = PN.findFlowRates()

    #get output
    PN.printPipeFlowRates()
    print()
    print('Check node flows:')
    PN.printNetNodeFlows()
    print()
    print('Check loop head loss:')
    PN.printLoopHeadLoss()
    #PN.printPipeHeadLosses()
# endregion


# region function calls
if __name__ == "__main__":
    main()
# endregion
