dict1 = {1: "Monday", 2: "Tuesday", 3: "Wednesday", 4: "Thursday", 5: "Friday", 6: "Saturday", 0: "Sunday"}
dict2 = {"January": 31, "February": 28, "March": 31, "April": 30, "May": 31, "June": 30,
         "July": 31, "August": 31, "September": 30, "October": 31, "November": 30, "December": 31}
Listp = [0, 0, 0, 0, 0, 0, 0]
Lista = [0, 0, 0, 0, 0, 0, 0]
Listx = ["", ""]
with open("input.txt") as f:
    n = int(f.readline())
    year = int(f.readline())
    days = 365
    if ((year % 4 == 0) and year % 100 != 0) or (year % 400 == 0):
        days = 366
        dict2["February"] = 29
    for i in range(1, days + 1, 1):
        Lista[i % 7] += 1
    for j in range(n):
        Listx = [item for item in f.readline().split()]
        cord = 0
        for k, v in dict2.items():
            if k != Listx[1]:
                cord += dict2[k]
            else:
                cord += int(Listx[0])
                break
        Listp[cord % 7] += 1
    day = f.readline().split()


def get_key(d, value):
    for k1, v1 in d.items():
        if v1 == value:
            return k1


if day[0] != dict1[1]:
    code = 7 if day[0] == "Sunday" else get_key(dict1, day[0])
    values = list(dict1.values())
    dict1 = dict(zip(dict1.keys(), values[-(8 - code):] + values[:-(8 - code)]))

for i in range(7):
    Lista[i] = Lista[i] + sum(Listp) - Listp[i]

print(dict1[Lista.index(max(Lista))], dict1[Lista.index(min(Lista))])
