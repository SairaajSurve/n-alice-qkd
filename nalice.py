from qiskit import *
import math
import numpy as np
import hashlib
import hmac
import time
import json
from Crypto.Cipher import AES
from Crypto.Random import get_random_bytes
import Padding
from qiskit.tools.visualization import circuit_drawer
import qiskit
qiskit.__qiskit_version__

simulator = Aer.get_backend('qasm_simulator')

class QChannel():
    def __init__(self, n, delta, c ,numberOfCommunicators, qber):
        # size of INFO BITS
        self.n = n
        
        # delta parameter
        self.delta = delta

        # N
        self.N = math.ceil(8*self.n*(1+self.delta))

        # probability of not measuring by communicators
        self.c = c

        self.numberOfCommunicators = numberOfCommunicators

        self.qber = qber

    def performTrip(self):
        
        # log
        print("Starting Algorithm")
        print(f'{self.n=}')
        print(f'{self.delta=}')
        print(f'{self.N=}')
        print(f'{self.c=}')
        print(f'{self.numberOfCommunicators=}')
        print(f'{self.qber=}')


        # quantum channel
        circuit = QuantumCircuit(self.N, self.N)

        # encoding of initial state
        # 0 = |+> and 1 = |->
        initialStateDecision = np.random.randint(2, size=self.N)

        print("initialStateDecision")
        print(initialStateDecision)

        for i in range(len(initialStateDecision)):
            if initialStateDecision[i]:
                circuit.x(i)
            circuit.h(i)

        circuit.barrier()

        # decisions by each communicator

        allMeasureDecisions = []

        for i in range(self.numberOfCommunicators):
            allMeasureDecisions.append(list(np.random.choice([0, 1], self.N, p=[self.c, 1-self.c])))
        
        print("allMeasureDecisions")
        for amd in allMeasureDecisions:
            print(amd)

        # communicator result

        values = []

        for i in range(self.numberOfCommunicators):
            for j in range(self.N):
                if allMeasureDecisions[i][j]:
                    circuit.measure(j, j)
            result = execute(circuit, backend=simulator, shots=1).result()
            count = result.get_counts(circuit)
            values.append([int(x) for x in list(list(count.keys())[0][::-1])])
            for j in range(len(values[-1])):
                if allMeasureDecisions[i][j]:
                    if  values[-1][j] == 0:
                        circuit.initialize([1, 0], j)
                    else:
                        circuit.initialize([0, 1], j)
            circuit.barrier()

        print("values")
        for v in values:
            print(v)

        # sift bits

        sift = dict()

        for i in range(self.numberOfCommunicators):
            for j in range(i+1, self.numberOfCommunicators):
                sift[(i, j)] = []

        # (p0) None of Bi measured the travel qubit
        # (p1) One of Bi measured the travel qubit but the other did not
        # (p2) Two of Bi measured the travel qubit.

        partiesInvolved = [[] for i in range(self.N)]

        for j in range(self.N):
            for i in range(self.numberOfCommunicators):
                if allMeasureDecisions[i][j] and len(partiesInvolved[j]) < 2:
                    partiesInvolved[j].append(i)

        print("partiesInvolved")
        print(partiesInvolved)

        # get sift bits

        for q in range(self.N):
            if len(partiesInvolved[q]) == 2:
                sift[(partiesInvolved[q][0], partiesInvolved[q][1])].append(q)

        print("sift")
        for k, v in sift.items():
            print(k, v)

        for k, v in sift.items():
            if len(v) < self.n + 1:
                print("Abroted: SIFT of "+str(k)+" is less than required size "+str(self.n + 1)+"(n+1)")
                return None

        cases = []

        for i in range(self.N):
            cases.append(len(partiesInvolved[i]))

        print("cases")
        print(cases)
        
        ctrlA = set()
        ctrlB = set()

        for i in range(len(cases)):
            if cases[i] == 0:
                ctrlA.add(i)
            elif cases[i] == 1:
                ctrlB.add(i)

        print(f'{ctrlA=}')
        print(f'{ctrlB=}')

        # separate siftCheck and info bits

        siftCheck = dict()
        info = dict()

        for k, v in sift.items():
            info[k] = v[:self.n]
            siftCheck[k] = v[self.n:]


        print("siftCheck")
        for k, v in siftCheck.items():
            print(k, v)

        print("info")
        for k, v in info.items():
            print(k, v)

        # decisions of measure basis by channel
        # 0 = Z and 1 = H

        finalMeasurementBasis = np.random.randint(2, size=self.N)

        for i in range(self.N):
            if finalMeasurementBasis[i]:
                circuit.h(i)
            circuit.measure(i, i)

        print("finalMeasurementBasis")
        print(finalMeasurementBasis)

        result2 = execute(circuit, backend=simulator, shots=1).result()
        count2 = result2.get_counts(circuit)
        value2 = [int(x) for x in list(list(count2.keys())[0][::-1])]

        print("value2")
        print(value2)

        # checking error

        qberCtrlA = 0
        qberCtrlB = 0

        qberSift = dict()
        for k, v in siftCheck.items():
            qberSift[k] = 0

        for q in ctrlA:
            if finalMeasurementBasis[q] == 1:
                qberCtrlA += int(value2[q] != initialStateDecision[q])
        
        for q in ctrlB:
            if finalMeasurementBasis[q] == 0:
                qberCtrlB += int(value2[q] != values[partiesInvolved[q][0]][q])

        for k, v in siftCheck.items():
            n1, n2 = k
            for q in v:
                qberSift[k] += int(values[n1][q] != values[n2][q])/len(v)


        qberCtrl = (qberCtrlA+qberCtrlB)/(len(ctrlA)+len(ctrlB))

        print(f'{qberCtrl=}')
        
        print("qberSift")
        for k, v in qberSift.items():
            print(k, v)

        if qberCtrl > self.qber:
            print("Abroted: qberCtrl is greater than permissible amount "+str(self.qber))
            return None
        
        for k, v in qberSift.items():
            if v > self.qber:
                print("Abroted: qberSift is greater than permissible amount "+str(self.qber)+" for key "+(str(k)))
                return None

        print("Algorithm successful")

        infoStrings = dict()

        for k, v in info.items():
            n1, n2 = k
            infoStrings[k] = [values[n1][q] for q in v]

        print("infoStrings")
        for k, v in infoStrings.items():
            print(k, v)

        circuit_drawer(circuit, output='mpl', filename='nalice.png')
        
        return infoStrings

QChannel(n=2, delta=1/2, c=0.5, numberOfCommunicators=3, qber=0.5).performTrip()