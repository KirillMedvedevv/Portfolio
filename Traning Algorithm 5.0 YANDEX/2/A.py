with open("input.txt") as f:
    w = int(f.readline())
    coords = []
    for i in range(w):
        coords.append([int(x) for x in f.readline().split()])

xcoords = [coords[i][0] for i in range(w)]
ycoords = [coords[i][1] for i in range(w)]
x1, y1 = min(xcoords), min(ycoords)
x2, y2 = max(xcoords), max(ycoords)

print(x1, y1, x2, y2)