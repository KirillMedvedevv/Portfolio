n = int(input())
num = list(map(lambda x: int(x) % 2, input().split()))
if num.count(1) % 2 == 1:
    print((n-1)*"+")
else:
    for i in range(n-1):
        if (num[i] == 1 and num[i+1] == 1) or (num[i] == 1 and num[i+1] == 1) or (num[i] == 0 and num[i+1] == 1):
            print(i*"+" + "x" + (n-i-2)*"+")
            break
