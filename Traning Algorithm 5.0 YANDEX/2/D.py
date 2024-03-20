doska = [["#" for i in range(10)] for j in range(10)] # Шахматная доска
with open("input.txt") as f:
    w = int(f.readline())
    for i in range(w):
        x, y = [int(x) for x in f.readline().split()]
        y = y
        x = x
        doska[x][y] = '0'

per = 0
for i in range(1, 9):
    for j in range(1, 9):
        if doska[i][j] == "0":
            if doska[i-1][j] == "#":
                per += 1
            if doska[i+1][j] == "#":
                per += 1
            if doska[i][j-1] == "#":
                per += 1
            if doska[i][j+1] == "#":
                per += 1

print(per)
