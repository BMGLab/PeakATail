f = open("/home/user/D/BAMdata/proje/ProjectEMA/emaout/filterdmatrix.mtx", "r")
f.readline()
f.readline()
m1 = 0
m2 = 0
s = 0
for l in f:
    s += 1
    d = l.split()
    if len(d) != 3:
        print(s)
        print("ewsrdfghjk")
    if int(d[0]) > m1:
        m1 = int(d[0])
        
    if int(d[1]) > m2:
        m2 = int(d[1])

print(m1)
print(m2)
print(s) 