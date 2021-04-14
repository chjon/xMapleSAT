import sys
import math
import os

if len(sys.argv) != 3:
	print "Usage: {} <NUM_HOLES> <GENERATE_EXTENSION_VARS>".format(sys.argv[0])
	exit()

n = int(sys.argv[1])
generateExtensionVars = int(sys.argv[2])

def generatePHP(n, generateExtensionVars, outf):
    # for each pigeon clause                                                           
    for i in range(1, n + 1):
        outstr = ""
        for j in range(1, n):
            P_ij = str((i - 1) * (n - 1) + j)
            outstr += P_ij
            outstr += " "
        outstr += "0\n"
        outf.write(outstr)

    # for each hole clause                                                             
    for k in range(1, n):
        for j in range(1, n + 1):
            for i in range(1, j):
                P_ik = str((i - 1) * (n - 1) + k)
                P_jk = str((j - 1) * (n - 1) + k)
                outf.write("-{} -{} 0\n".format(P_ik, P_jk))

    # for each extension variable
    if generateExtensionVars:
        offset = 0
        for k in range(n, 1, -1):
            for i in range(1, k):
                for j in range(1, k - 1):
                    P_ij = str(offset + (i - 1) * (k - 1) + j)
                    P_ik = str(offset + (i - 1) * (k - 1) + k)
                    P_kj = str(offset + (k - 1) * (k - 1) + j)
                    Q_ij = str(offset + (k    ) * (k - 1) + (i - 1) * (k - 2) + j)
                    outf.write("{} -{} 0\n".format(Q_ij, P_ij))
                    outf.write("{} -{} -{} 0\n".format(Q_ij, P_ik, P_kj))
                    outf.write("-{} {} {} 0\n".format(Q_ij, P_ij, P_ik))
                    outf.write("-{} {} {} 0\n".format(Q_ij, P_ij, P_kj))
            offset += k * (k - 1)
    return

if not os.path.exists("instances"):
    os.makedirs("instances")
outf = open("instances/{0}_{1}.cnf".format(n, generateExtensionVars),"w")

nVar = n * (n - 1)
nClause = int(n + (n - 1) * (n - 1 + 1) * (n - 1) / 2)
if generateExtensionVars:
    nxVar = int((n - 1 - 1) * (n - 1) * (n - 1 + 1) / 3)
    nVar += nxVar
    nClause += 4 * nxVar
outf.write("p cnf {0} {1}\n".format(nVar, nClause))
generatePHP(n, generateExtensionVars, outf)

outf.close()

