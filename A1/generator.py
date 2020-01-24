import random

def main():
	n = input()
	file = open("data.txt","w")
	file.write(str(n) + "\n")
	for i in range(n):
		for j in range(n):
			file.write(str(10.0*random.random()) + " ")
		if i!=n-1:
			file.write("\n")
	file.close()

if __name__ == '__main__':
	main()