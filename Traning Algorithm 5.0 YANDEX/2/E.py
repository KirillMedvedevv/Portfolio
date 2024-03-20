highmMAX = 0
coordm = -1
dhighpMAX = 0
coordp = -1
cherrym = []
cherryp = []

value = 0
with open("input.txt") as f:
    w = int(f.readline())
    for i in range(w):
        x, y = [int(x) for x in f.readline().split()]
        if (x - y) > 0:
            cherryp.append(i + 1)
            if y > dhighpMAX:
                dhighpMAX = y
                coordp = i+1
            value += x - y
        else:
            cherrym.append(i + 1)
            if x > highmMAX:
                highmMAX = x
                coordm = i+1

number = [0]
string = ""
answer = value
if highmMAX < dhighpMAX:
    number[0] = str(coordp)
    cherryp.remove(coordp)
    answer += dhighpMAX
else:
    number[0] = str(coordm)
    cherrym.remove(coordm)
    answer += highmMAX

string += ' '.join(map(str, (cherryp + number + cherrym)))

print(answer)
print(string)
