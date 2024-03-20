with open("input.txt") as f:
    t = int(f.readline())
    for i in range(t):
        k = int(f.readline())
        line = [int(x) for x in f.readline().split()]
        #print(line)
        locmin, length = 999999999, 0
        answer = []
        for j in line:
            if locmin > j:
                locmin = j
            if locmin >= length+1:
                length += 1
            else:
                answer.append(length)
                locmin = j
                length = 1
        answer.append(length)
        print(len(answer))
        print(' '.join(map(str, answer)))
