def cheker(array):
    number = 0
    for i in range(1, n + 3):
        for j in range(1, m + 3):
            if array[i][j] == '1':
                number += 1

    return number


def get_rect(array, x1, y1): #Возвращает координаты прямоугольника и удаляет найденные
    x2, y2 = -1, -1
    x3, y3 = -1, -1
    logi = 0
    for j in range(x1+1, m+3):
        if array[y1][j] == "1" or array[y1][j] == "3":
            x2 = j
            y2 = y1
            break

    for i in range(y2+1, n+3):
        if array[i][x2] == '1' or array[i][x2] == '3':
            x3 = x2
            y3 = i
            break

    if x3 == -1 and y3 ==  -1:
        x3 = x2
        y3 = y2
        logi = 1

    if x2 == -1 and y2 == -1:
        x2 = x1
        y2 = y1
        x3 = x2
        y3 = y2
        logi = 2
    if logi == 0:
        array[y1][x1] = '#' if array[y1][x1] == "1" else array[y1][x1]
        array[y2][x2] = '#' if array[y2][x2] == "1" else array[y2][x2]
        array[y3][x3] = '#' if array[y3][x3] == "1" else array[y3][x3]
        array[y3][x1] = '#' if array[y3][x1] == "1" else array[y3][x1]
    elif logi == 1:
        array[y1][x1] = '#' if array[y1][x1] == "1" else array[y1][x1]
        array[y2][x2] = '#' if array[y2][x2] == "1" else array[y2][x2]
    elif logi == 2:
        array[y1][x1] = '#' if array[y1][x1] == "1" else array[y1][x1]

    return [[y1, x1], [y2, x2], [y3, x3], [y3, x1]]


def builder(coords, p): # Строим прямоугольник из границ
    print(coords)
    for j in range(coords[0][1], coords[1][1]+1):
        for i in range(coords[0][0], coords[3][0]+1):
            massive[i][j] = p

# Считываем входные данные
with open("input.txt") as f:
    massive = []
    n, m = [int(x) for x in f.readline().split()]
    massive.append(['.' for x in range(m+4)])
    massive.append(['.' for x in range(m + 4)])
    for i in range(n):
        line = list(f.readline())
        try:
            line.remove('\n')
        except ValueError:
            pass
        line.append('.')
        line.append('.')
        line.insert(0, '.')
        line.insert(0, '.')
        massive.append(line)
    massive.append(['.' for x in range(m + 4)])
    massive.append(['.' for x in range(m + 4)])

count1 = 0 # Число крестиков
count2 = 0 # Число кольев
count3 = 0
coords1 = []
coords2 = []

# Вбиваем колышки
for i in range(1, n+3):
    lcoords = []
    for j in range(1, m+3):
        if massive[i][j] == '.':
            # Отмечаем границы
            if massive[i][j + 1] == '.' and massive[i + 1][j] == '.' and massive[i + 1][j + 1] != '.':
                massive[i + 1][j + 1] = '1'
                count2 += 1
            elif massive[i][j - 1] == '.' and massive[i + 1][j] == '.' and massive[i + 1][j - 1] != '.':
                massive[i + 1][j - 1] = '1'
                count2 += 1
            elif massive[i][j + 1] == '.' and massive[i - 1][j] == '.' and massive[i - 1][j + 1] != '.':
                massive[i - 1][j + 1] = '1'
                count2 += 1
            elif massive[i][j - 1] == '.' and massive[i - 1][j] == '.' and massive[i - 1][j - 1] != '.':
                massive[i - 1][j - 1] = '1'
                count2 += 1
            elif massive[i][j+1] == '#' and massive[i+1][j] == '#' and massive[i+1][j+1] == '#':
                massive[i][j + 1] = '3'
                massive[i+1][j + 1] = '4'
                count2 += 1
                count3 += 1
                coords2.append([i, j+1])
            elif massive[i][j - 1] == '#' and massive[i + 1][j] == '#' and massive[i + 1][j - 1] == '#':
                massive[i][j - 1] = '3'
                massive[i + 1][j - 1] = '4'
                count3 += 1
                count2 += 1
                coords2.append([i, j-1])
            elif massive[i][j+1] == '#' and massive[i-1][j] == '#' and massive[i-1][j+1] == '#':
                massive[i][j + 1] = '3'
                massive[i - 1][j + 1] = '4'
                count3 += 1
                count2 += 1
                coords2.append([i, j+1])
            elif massive[i][j-1] == '#' and massive[i-1][j] == '#' and massive[i-1][j-1] == '#':
                massive[i][j - 1] = '3'
                massive[i - 1][j - 1] = '4'
                count3 += 1
                count2 += 1
                coords2.append([i, j-1])
        else:
            count1 += 1

logic = 2
if count2 > 8 or (count1 == 1 and count2 == 4):
    logic = 0
elif count2 == 4:
    # Один прямоугольник сразу же разбиваем
    logic = 1
    for i in range(1, n+3):
        parm = 0
        for j in range(1, m+3):
            if massive[i][j] != ".":
                if parm == 0:
                    massive[i][j] = "a"
                    parm = 1
                else:
                    massive[i][j] = "b"
                    logic = 5

    if logic == 5:
        print("YES")
        for i in range(2, n + 2):
            print((' '.join(map(str, massive[i][2:m + 2]))).replace(' ', ''))
    else:
        for j in range(1, m + 3):
            parm = 0
            for i in range(1, n + 3):
                if massive[i][j] != ".":
                    if parm == 0:
                        massive[i][j] = "a"
                        parm = 1
                    else:
                        massive[i][j] = "b"
        print("YES")
        for i in range(2, n + 2):
            print((' '.join(map(str, massive[i][2:m + 2]))).replace(' ', ''))


elif count2 < 8:
    #Построим отражение точки 3
    for k in coords2:
        for j in range(1, m + 3):
            massive[k[0]][k[1]] = '1'
            x1, y1 = 0, 2147483647
            x2, y2 = 0, -1
            for i in range(1, k[0]):
                if massive[i][j] == '1':
                    y1, x1 = i, j
                    break

            for i in range(k[0] + 1, n + 3):
                if massive[i][j] == '1':
                    y2, x2 = i, j
                    break
            if (k[0] < y2) and (k[0] > y1):
                massive[k[0]][x1] = '1'
                count2 += 1
    #Построим отражение точки 4
    for k in coords2:
        if massive[k[0] + 1][k[1]] == '4':
            massive[k[0] + 1][k[1]] = '#'
            for j in range(1, k[1]):
                if (massive[k[0] + 1][j] == '.' and  massive[k[0] + 1][j+1] == '#'):
                    massive[k[0]+1][j+1] = '1'
                    count2 += 1 if massive[k[0] + 1][j+1] != "1" else 0

        elif massive[k[0] - 1][k[1]] == '4':
            massive[k[0] - 1][k[1]] = '#'
            for j in range(k[1]+1, m+3):
                if massive[k[0] - 1][j-1] == '.' and  massive[k[0] + 1][j+1] == '#':
                    massive[k[0]+1][j+1] = '1'
                    count2 += 1 if massive[k[0] - 1][j+1] != "1" else 0

    number1 = cheker(massive)

    if (number1 < 8 and number1 > 4) or number1 > 8:
        logic = 0


if logic == 0:
    print("NO")
elif logic == 1:
    pass
elif logic == 2:
    # Строим прямоугольник
    for i in massive:
        print(i)
    x1, y1 = -1, -1
    for i in range(1, n + 3):
        for j in range(1, m + 3):
            if massive[i][j] == "1":
                x1, y1 = j, i
                break
        if (x1 != -1) and (y1 != -1):
            break
    rectangular = []
    rectangular.append(get_rect(massive, x1, y1))

    builder(rectangular[0], 'a')

    x1, y1 = -1, -1
    for i in range(1, n + 3):
        for j in range(1, m + 3):
            if massive[i][j] == "1":
                x1, y1 = j, i
                break
        if (x1 != -1) and (y1 != -1):
            break
    rectangular = []
    for i in massive:
        print(i)
    rectangular.append(get_rect(massive, x1, y1))
    builder(rectangular[0], 'b')
    for i in massive:
        print(i)

    print("YES")
    for i in range(2, n + 2):
        print((' '.join(map(str, massive[i][2:m + 2]))).replace(' ', ''))
