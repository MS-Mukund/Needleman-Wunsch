import matplotlib.pyplot as plt

icc = open('wparl_icc_pro.txt', 'r')
gcc0 = open('wparl_gcc0_pro.txt', 'r')
gcc1 = open('wparl_gcc1_pro.txt', 'r')
gcc2 = open('wparl_gcc2_pro.txt', 'r')
gcc3 = open('wparl_gcc3_pro.txt', 'r')

icc = icc.readlines()
gcc0 = gcc0.readlines()
gcc1 = gcc1.readlines()
gcc2 = gcc2.readlines()
gcc3 = gcc3.readlines()

for i in range(15):
    icc[i] = float(icc[i].strip())
    gcc0[i] = float(gcc0[i].strip())
    gcc1[i] = float(gcc1[i].strip())
    gcc2[i] = float(gcc2[i].strip())
    gcc3[i] = float(gcc3[i].strip())


xaxis = [1000 * i for i in range(len(icc))]
plt.plot( xaxis, icc, label='icc' )
plt.plot( xaxis, gcc0, label='gcc0' )
plt.plot( xaxis, gcc1, label='gcc1' )
plt.plot( xaxis, gcc2, label='gcc2' )
plt.plot( xaxis, gcc3, label='gcc3' )

plt.xlabel('DNA Sequence length')
plt.ylabel('Gflops')

plt.title('Along the anti-diagonals of the matrix')

plt.legend()
plt.show()