x = int(input())  # Число моих солдат
y = int(input())  # Здоровье казармы
p = int(input())  # Производительность противника
s = 0  # число солдат противника
Count = 0
Answer1 = 999999999
buff = 0
num = 0
loose = 0


def hod(n):  # Обычный ход
    global x, y, p, s, num
    y -= n
    if y <= 0:
        y = 0
    s = s - x + n
    if s < 0:
        s = 0
    x -= s
    if y != 0:
        s = s + p
    num += 1


def analyser(x1, y1, s1, n):
    for i in range(n):
        s1 = s1 - (x1 - y1)
        y1 = 0
        if s1 < 0:
            s1 = 0
        x1 = x1 - s1
        if x1 <= 0:
            return -1
        if x1 >= s1:
            return i+2

    return -1


hod(x)
while 1:
    buff = analyser(x, y, s, 10)
    if s == 0 and y == 0:
        break
    elif x <= 0 or (x == y and p == x):
        loose = -1
        break
    elif (x - p) <= 0:
        hod(y)
    elif buff == -1:
        hod(x - p)
    else:
        count = num
        count += buff #Получаем количетсво ходов за которые можем победить
        if count <= Answer1:
            Answer1 = count
        hod(x - p)

if loose == 0:
    print(min(Answer1, num))
else:
    print(loose)