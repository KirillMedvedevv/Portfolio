with open("input.txt") as f:
    N, K = [int(x) for x in f.readline().split()]
    price = [int(x) for x in f.readline().split()]

profit = [0]

i = 0
while 1:
    j = i
    while 1:
        if j >= i+K or j >= N-1:
            break
        j += 1
        profit.append(price[j] - price[i])
    i += 1
    if i >= N-1:
        break

if max(profit) <= 0:
    print(0)
else:
    print(max(profit))
