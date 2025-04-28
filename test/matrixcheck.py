matrix = open("/home/user/D/BAMdata/proje/ProjectEMA/emaout/posmatrix.mtx", "r")
c = 0
for line in matrix:
    lis = line.split()
    c+=1
    try:
        for i in range(3):
            if type(lis[i]) is int:
                print("hey")
                
    except:
        print(line)
        print(c)