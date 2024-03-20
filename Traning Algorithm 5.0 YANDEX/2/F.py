n = int(input())
wheel = [int(x) for x in input().split()]
a, b, k = [int(x) for x in input().split()]
maxpwin = max(wheel)
maxwin = 0

for i in range(a, b+1, k):
    m = i//k - 1 if i % k == 0 else i//k
    win = max(wheel[m % n], wheel[-m % n])
    if win == maxpwin:
        maxwin = win
        break
    elif win > maxwin:
        maxwin = win

print(maxwin)
