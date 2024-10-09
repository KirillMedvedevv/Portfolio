with open("input.txt") as f:
    dict1 = {}
    n = int(f.readline())
    for i in range(n):
        k = int(f.readline())
        sings = f.readline().split()
        for i in sings:
            try:
                dict1[i] += 1
            except KeyError:
                dict1[i] = 1

count = 0
answer = []
for k, v in dict1.items():
    if dict1[k] == n:
        count += 1
        answer.append(k)

answer.sort()
print(count)
print(" ".join(answer))

