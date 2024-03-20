with open("input.txt") as f:
    coords = []
    N = int(f.readline())
    for i1 in range(N):
        coords.append([int(x)-1 for x in f.readline().split()])


def centralize(array):  # 0 - строка 1 - столбец
    # Определяем близлежащий столбец
    m = -1
    length = 2147483647
    for j in range(N):
        local = 0
        for k in array:
            local += abs(k[1] - j)
        if local < length or m == -1:
            length = local
            m = j
    line = [0 for k in range(N)]  # Пустая корабельная линия
    # Считаем количество ходов по стягиванию и генерируем корабельную линию
    hod = 0
    for i in array:
        line[i[0]] += 1
        hod += abs(i[1] - m)


    noone = []
    for i in range(N):
        if line[i] > 1:
            for j in range(line[i]-1):
                noone.append(i)

    j = 0
    for i in range(N):
        if line[i] == 0:
            hod += abs(i - noone[j])
            j += 1

    return hod


print(centralize(coords))
