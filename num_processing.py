#process a txt file with numbers and take the average

with open('wparl_gcc3_raw.txt', 'r') as f:
    lines = f.readlines()

    arr = []
    asum = 0
    for i, line in enumerate(lines, start=1):
        if i % 10 == 0:
            arr.append(asum/10)
            asum = 0
        
        asum += float(line)
        i += 1

    arr.append(asum/10)
    #write to file
    with open('wparl_processed_gcc3.txt', 'w') as f:
        for i in arr:
            f.write(str(i) + '\n')
