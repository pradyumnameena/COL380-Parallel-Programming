import random

def write_to_file(file_name,rows,cols):
	file = open(file_name,"w")
	for i in range(rows):
		for j in range(cols):
			file.write(str(10.0*random.random()) + " ")
		if i!=rows-1:
			file.write("\n")
	file.close()

def main():
	n = input()
	file_name1 = "A_" + str(n) + ".txt"
	file_name2 = "B_" + str(n) + ".txt"
	write_to_file(file_name1,n,32)
	write_to_file(file_name2,32,n)

if __name__ == '__main__':
	main()