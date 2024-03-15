a11, b11 = input().split(":")
a21, b21 = input().split(":") #Cчёт матчa 2
a11, a21, b11, b21 = int(a11), int(a21), int(b11), int(b21)
c1 = int(input()) # 1 - дома - гости 2 - гости - дома
deltap = 0


def comment(a1, a2, b1, b2, c, delta):
    if (a1 + a2) > (b1 + b2):
        return 0 + delta  #Общий счёт в нашу пользу. Пофиг
    elif (a1+a2) == (b1 + b2): #Общий счёт одинаков
        if c == 2:
            if a1 > b2:
                return 0 + delta
            else:
                return 1 + delta
        elif c == 1:
            if a2 > b1:
                return 0 + delta
            else:
                return 1 + delta
    else: #Разный общий счёт
        #delta1 = b1 + b2 - a1 - a2 #Сравнять счёт
        a2 = a2 + 1
        delta1 = delta+1
        return comment(a1, a2, b1, b2, c, delta1)


print(comment(a11, a21, b11, b21, c1, deltap))
