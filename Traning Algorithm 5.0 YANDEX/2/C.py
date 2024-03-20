with open("input.txt") as f:
    N = int(f.readline())
    list1 = [int(x) for x in f.readline().split()]

max1 = max(list1)
for i in range(N):
    if list1[i] == max1:
        list1[i] = 0
        break

if sum(list1) >= max1:
    print(sum(list1) + max1)
else:
    print(max1 - sum(list1))
