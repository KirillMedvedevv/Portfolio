doska = [["", "", "", "", "", "", "", ""],
         ["", "", "", "", "", "", "", ""],
         ["", "", "", "", "", "", "", ""],
         ["", "", "", "", "", "", "", ""],
         ["", "", "", "", "", "", "", ""],
         ["", "", "", "", "", "", "", ""],
         ["", "", "", "", "", "", "", ""],
         ["", "", "", "", "", "", "", ""]]


def risov(doska, x, y, a):
    if a == "R":
        i, j = y, x
        while (i < 7):
            i += 1
            if doska[i][j] == "B" or doska[i][j] == "R":
                break
            elif doska[i][j] == "*":
                doska[i][j] = "P"
        i, j = y, x
        while (i > 0):
            i -= 1
            if doska[i][j] == "B" or doska[i][j] == "R":
                break
            elif doska[i][j] == "*":
                doska[i][j] = "P"
        i, j = y, x
        while (j > 0):
            j -= 1
            if doska[i][j] == "B" or doska[i][j] == "R":
                break
            elif doska[i][j] == "*":
                doska[i][j] = "P"
        i, j = y, x
        while (j < 7):
            j += 1
            if doska[i][j] == "B" or doska[i][j] == "R":
                break
            elif doska[i][j] == "*":
                doska[i][j] = "P"
    elif a == "B":
        i, j = y, x
        while(i < 7) and (j < 7):
            i += 1
            j += 1
            if doska[i][j] == "B" or doska[i][j] == "R":
                break
            elif doska[i][j] == "*":
                doska[i][j] = "P"
        i, j = y, x
        while (i > 0) and (j < 7):
            i -= 1
            j += 1
            if doska[i][j] == "B" or doska[i][j] == "R":
                break
            elif doska[i][j] == "*":
                doska[i][j] = "P"
        i, j = y, x
        while (i > 0) and (j > 0):
            i -= 1
            j -= 1
            if doska[i][j] == "B" or doska[i][j]== "R":
                break
            elif doska[i][j] == "*":
                doska[i][j] = "P"
        i, j = y, x
        while (i < 7) and (j > 0):
            i += 1
            j -= 1
            if doska[i][j] == "B" or doska[i][j] == "R":
                break
            elif doska[i][j] == "*":
                doska[i][j] = "P"
    else:
        pass


with open("input.txt") as f:
    a = 0
    for i in f:
        for j in range(8):
            doska[a][j] = i[j]
        a += 1

for i in range(8):
    for j in range(8):
        risov(doska, j, i, doska[i][j])

count = 0

for i in range(8):
    for j in range(8):
        if doska[i][j] == "*":
            count += 1

print(count)