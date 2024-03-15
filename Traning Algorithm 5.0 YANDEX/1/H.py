from math import pi
L, x1, v1, x2, v2 = input().split()
L, x1, v1, x2, v2 = float(L), float(x1), float(v1), float(x2), float(v2)
R = L/(2.0*pi)
p1, w1, p2, w2 = x1/R, v1/R, x2/R, v2/R
List = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]


def maxroot():
    def root1(k):
        return (2.0 * pi * k - p1 - p2) / (w1 + w2)

    def root2(k):
        return (2.0 * pi * k - p1 + p2) / (w1 - w2)

    k1 = round((p1 + p2) / (2.0 * pi))
    k2 = round((p1 - p2) / (2.0 * pi))

    if (w1 + w2) != 0.0:
        List[0] = root1(k1)
        List[1] = root1(k1 - 1.0)
        List[2] = root1(k1 + 1.0)
        List[3] = root1(k1 + 2.0)
        List[4] = root1(k1 - 2.0)

    if (w1 - w2) != 0.0:
        List[5] = root2(k2 - 1.0)
        List[6] = root2(k2)
        List[7] = root2(k2 + 1.0)
        List[8] = root2(k2 + 2.0)
        List[9] = root2(k2 - 2.0)
    try:
        return min(filter(lambda val: val >= 0.0, List))
    except ValueError:
        return -1


if p1 == p2:
    print("YES")
    print(0)
elif maxroot() == -1:
    print("NO")
else:
    print("YES")
    print(maxroot())
