dict1 = {}
dict2 = {}
stroke1 = input()
stroke2 = input()
len1 = len(stroke1)
len2 = len(stroke2)

for i in range(len1):
    dict1[stroke1[i]] = 0

for i in range(len2):
    dict2[stroke2[i]] = 0

for i in range(len1):
    dict1[stroke1[i]] += 1

for i in range(len2):
    dict2[stroke2[i]] += 1

logic = 1
for k, v in dict1.items():
    try:
        if dict2[k] != v:
            logic = 0
            break
    except KeyError:
        logic = 0
        break

if logic == 1:
    print('YES')
else:
    print("NO")