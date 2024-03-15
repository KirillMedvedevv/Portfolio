P, V = (input().split())
Q, M = (input().split())
P, V, Q, M = int(P), int(V), int(Q), int(M)
N = 0
P2, P1 = P + V, P - V #Область пацана
Q2, Q1 = Q + M, Q - M #Область бабы
N = abs(P1 - P2) + abs(Q1 - Q2) + 2

#Отношения между областями
if (Q1 >= P1) and (Q2 <= P2):
    N = N - abs(Q2 - Q1) - 1
elif (P1 >= Q1) and (P2 <= Q2):
    N = N - abs(P2 - P1) - 1
elif (P2 >= Q1) and (Q2 >= P2):
    N = N - abs(P2 - Q1) - 1
elif (Q2 >= P1) and (P2 >= Q2):
    N = N - abs(Q2 - P1) - 1

print(N)