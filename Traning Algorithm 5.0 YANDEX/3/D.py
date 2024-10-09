n, k = [int(x) for x in input().split()]
numbers = [int(x) for x in input().split()]
dict1 = {}
logic = 0
for i in numbers:
    try:
        dict1[i] += 1
    except KeyError:
        dict1[i] = 1

for k1 in dict1.keys():
    if dict1[k1] > 1:
        for i in range(n):
            if numbers[i] == k1:
                try:
                    for j in range(1, k+1):
                        if numbers[i] == numbers[i+j]:
                            logic = 1
                            break
                except IndexError:
                    break
        if logic == 1:
            break
    if logic == 1:
        break

if logic == 1:
    print("YES")
else:
    print("NO")
