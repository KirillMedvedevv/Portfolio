number  = 0 #Число пробелов строке
a = 0 #Результат деления на 4
b = 0 # Остаток деления на 4
p = 0
result = 0

with open("input.txt") as f:
    for i in f:
        if p == 0:
            p = 1
        else:
            number = int(i)
            a = number // 4
            b = number % 4
            result += a
            if b == 3:
                result += 2
            else:
                result += b

print(result)