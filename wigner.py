import sys
from sympy.physics.wigner import wigner_6j
from sympy import N
l1 = int(sys.argv[1])
l2 = int(sys.argv[2])
l3 = int(sys.argv[3])
l4 = int(sys.argv[4])
l5 = int(sys.argv[5])
l6 = int(sys.argv[6])

res = wigner_6j(l1,l2,l3,l4,l5,l6)
f = open('wigner_results.dat','w')
f.write(str(N(res)))
f.close()
