def maxvalue(array, k, p):
    value, im, jm = 0, 0, 0
    for i in range(n):
        if i != k:
            for j in range(m):
                if j != p:
                    if array[i][j] > value:
                        value = array[i][j]
                        im = i
                        jm = j
    return value, im, jm


with open("input.txt") as f:
    persons = []
    n, m = [int(x) for x in f.readline().split()]
    for i1 in range(n):
        persons.append([int(x) for x in f.readline().split()])

# Находим максимальный элемент
maxpower = maxvalue(persons, -1, -1)
# Вычеркнем строку потом столбец
maxpower1 = maxvalue(persons, maxpower[1], -1)
maxpower11 = maxvalue(persons, maxpower[1], maxpower1[2])
# Вычеркнем столбец потом строку
maxpower2 = maxvalue(persons, -1, maxpower[2])
maxpower22 = maxvalue(persons, maxpower2[1], maxpower[2])

if maxpower11[0] < maxpower22[0]:
    print(maxpower[1]+1, maxpower1[2]+1)
else:
    print(maxpower2[1]+1, maxpower[2]+1)
