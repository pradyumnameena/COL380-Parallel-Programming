import random

def generate_data():
	n = input()
	file = open("data.txt","w")
	file.write(str(n) + "\n")
	for i in range(n):
		for j in range(n):
			file.write(str(10000.0*random.random()) + " ")
		if i!=n-1:
			file.write("\n")
	file.close()

def checker():
	return

def main():
	param = input()
	if param==0:
		generate_data()
	else:
		checker()

if __name__ == '__main__':
	main()