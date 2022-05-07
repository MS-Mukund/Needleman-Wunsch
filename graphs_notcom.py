import matplotlib.pyplot as plt

wparl = open('wparl_gcc0_pro.txt', 'r')
opti_256 = open('opti_128_gcc3_pro.txt', 'r')
# opti_128 = open('opti_128_gcc3_pro.txt', 'r')
# opti_256 = open('opti_256_gcc3_pro.txt', 'r')
# naive = open('naive_processed_gcc3.txt', 'r')

wparl = wparl.readlines()
# opti_128 = opti_128.readlines()
opti_256 = opti_256.readlines()
# naive = naive.readlines()

for i in range(30):
    # opti_128[i] = float(opti_128[i].strip())
    opti_256[i] = float(opti_256[i].strip())
    # naive[i] = float(naive[i].strip())

for i in range(16):
    wparl[i] = float(wparl[i].strip())

xaxis = [1000 * i for i in range(len(opti_256))]
# plt.plot( xaxis, opti_128, label='tiling 128' )
plt.plot( xaxis, opti_256, label='tiling 128' )
# plt.plot( xaxis, naive, label='naive' )

x2 = [1000 * i for i in range(len(wparl))]
plt.plot( x2, wparl, label='anti-diagonal' )


plt.xlabel('DNA Sequence length')
plt.ylabel('execution time (ms)')

plt.title('Performance comparison')

plt.legend()
plt.show()