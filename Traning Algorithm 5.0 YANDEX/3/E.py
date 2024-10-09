with open("input.txt") as f:
    k1 = f.readline()
    list1 = [int(x) for x in f.readline().split()]
    k2 = f.readline()
    list2 = [int(x) for x in f.readline().split()]
    k3 = f.readline()
    list3 = [int(x) for x in f.readline().split()]

numbers = []
dict1 = {}
dict2 = {}
dict3 = {}
for i in list1:
    try:
        dict1[i] += 1
