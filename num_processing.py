#process a txt file with numbers and take the average



with open('wparl_gcc2_raw.txt', 'r') as f:
    lines = f.readlines()

    f2 = open('wparl_gcc2_raw2.txt', 'r')
    l2 = f2.readlines()

    arr = []
    asum = 0
    for i, line in enumerate(lines, start=0):
        line2 = float(l2[i].strip())
        print( (float(line) + line2)/2 )
