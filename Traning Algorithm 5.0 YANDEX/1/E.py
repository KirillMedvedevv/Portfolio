import sys
sys.set_int_max_str_digits(0)
n, k, d = input().split()
n, k, d = int(n), int(k), int(d)
add = -1

if n % k == 0:
    print(n*10**d)
else:
    for j in range(10):
        buffer = 10 * n + j
        if buffer % k == 0:
            add = buffer
            break
    if add == -1:
        print(-1)
    else:
        print(add * (10 ** (d - 1)))
