from random import randrange, uniform
a = [int(x) for x in input("enter for matrix A : ").split()]
b = [int(x) for x in input("enter for matrix B : ").split()]

the_file = open("matA.txt", "a")
the_file.write(' '.join(str(e) for e in a)+'\n')
for x in range(1,a[0]+1):
	li = ""
	for y in range(1,a[1]+1):
		li = li + str(uniform(-100, 100)) + " "
	li = li + "\n"
	the_file.write(li)
the_file.close()

the_file = open("matB.txt", "a")
the_file.write(' '.join(str(e) for e in b)+'\n')
for x in range(1,b[0]+1):
	li = ""
	for y in range(1,b[1]+1):
		li = li + str(uniform(-100, 100)) + " "
	li = li + "\n"
	the_file.write(li)
the_file.close()