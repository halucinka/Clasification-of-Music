import math
import copy

def suma(a, b):
	#a, b su uz logaritmy, len ich treba scitat spravne, ln(0) = 47 :)
	e = 2.71828182845904523536
	if (a == 47):
		return b
	elif (b == 47):
		return a
	else:
		return 	a+math.log(1+e**(b-a))

def mul(a, b):
	#a, b su uz logaritmy, teda mul = scitanie
	if ((a == 47) or (b == 47)):
		return 47
	return a+b

class HMM:
	def __init__(self):
		self.s = [47, 47, 47]
		self.e = [[0 for x in range(-11,12)] for x in range(3)] 
	
		self.Pr_seq=0
		self.pocet_stavov = 3
		self.t = [[47 for x in range(self.pocet_stavov)]for x in range(self.pocet_stavov)]
		self.f = ''
		self.lines = 0

	def nastav(self, s, e, t, f, pocet_stavov):
		# nastavi pociatocne parametre uz v zlogaritmovanom stave
		self.pocet_stavov = pocet_stavov
		self.t = [[47 for x in range(self.pocet_stavov)]for x in range(self.pocet_stavov)]
		self.t = t
		self.e = [[0 for x in range(-11,12)] for x in range(self.pocet_stavov)]
		self.s = [47 for x in range(self.pocet_stavov)]
		self.s = s
		self.e = e
		self.f = f

		for i in range (0, len(self.s)):
			if self.s[i] == 0:
				self.s[i] = 47
			else:
			 self.s[i]=math.log(self.s[i])

		for j in range (0, self.pocet_stavov):
			for i in range (-11, 12):
				if self.e[j][i] == 0:
					self.e[j][i] = 47
				else:
					self.e[j][i]=math.log(self.e[j][i])
		for i in range (0, self.pocet_stavov):
			for j in range (0, self.pocet_stavov):
				if (self.t[i][j] == 0):
					self.t[i][j] = 47
				else:
					self.t[i][j]=math.log(self.t[i][j])

	def ForwardTabulka(self, sekvencia):
		#vo Forward tabulke si pamatame vzajomne pravdepodobnosti, nie podmienene, teda Pr(na i tej pozicii je stav j A sekvencia je seq[0:i])
		seq = [0 for x in range(len(sekvencia))]
		for i in range(0, len(sekvencia)):
			seq[i] = int(sekvencia[i])
		F= [[0 for x in range(self.pocet_stavov)] for x in range(len(seq))] 
		for j in range(0, self.pocet_stavov):
			F[0][j] = mul(self.s[j],self.e[j][seq[0]])
		for i in range (1,len(seq)):
			for j in range (0, self.pocet_stavov):
				temp=[0 for x in range(self.pocet_stavov)]
				for k in range(0, self.pocet_stavov):
						temp[k] =mul(mul(F[i-1][k],self.t[k][j]),self.e[j][seq[i]])
						#generujeme i-te pismenko v stave j, cez vsetky stavy j
				#scitame temp od najmensieho po najvacsie
				temp.sort()
				F[i][j]=1
				for k in range(self.pocet_stavov-1, -1, -1):
					F[i][j]=suma(F[i][j], temp[k])
		self.Pr_seq = 1
		temp = [0 for x in range(self.pocet_stavov)]
		for i in range(0, self.pocet_stavov):
			temp[i] = F[len(seq)-1][i]
		temp.sort()
		for j in range(self.pocet_stavov-1, -1, -1):		
			self.Pr_seq = suma(self.Pr_seq, temp[j])
		
		return F
		

	def BackwardTabulka(self, sekvencia):
		# v B zaciname ratanim prechodovej pravdepodobnosti
		seq = [0 for x in range(len(sekvencia))]
		for i in range(0, len(sekvencia)):
			seq[i] = int(sekvencia[i])
		B = [[1 for x in range(self.pocet_stavov)] for x in range(len(seq))] 
		for j in range (0, self.pocet_stavov):
			B[len(seq)-1][j] = 0 # toto je 0, lebo ln 1=0
		for i in range (len(seq)-2, -1, -1):
			for j in range (0, self.pocet_stavov):
				B[i][j] = 1
				temp=[1 for x in range(self.pocet_stavov)]
				for k in range(0, self.pocet_stavov):
					temp[k] =mul(mul(self.t[j][k],self.e[k][seq[i + 1]]),B[i+1][k])	
				temp.sort()
				for k in range(self.pocet_stavov-1, -1, -1):
					B[i][j] = suma(B[i][j],temp[k])

		return B


	def new_t_e(self):
		noty_rel_cisla = open (self.f, 'r')
		self.lines = 0
		for line in noty_rel_cisla:
			self.lines += 1
		noty_rel_cisla.close()
		noty_rel_cisla = open(self.f, 'r')
		#scitavanie cez vsetky trenovacie sekvencie:
		a=-1
		A = [[[1 for x in range(self.pocet_stavov)] for x in range(self.pocet_stavov)] for x in range(self.lines)]
		E = [[[1 for x in range(-11, 12)] for x in range(self.pocet_stavov)] for x in range(self.lines)]
		for line in noty_rel_cisla:
			a+=1
			seq=line.split()
			F= [[[47 for x in range(self.pocet_stavov)] for x in range(len(seq))] ]
			B= [[47 for x in range(self.pocet_stavov)] for x in range(len(seq))] 
			F = self.ForwardTabulka(seq)
			B = self.BackwardTabulka(seq)
			pravd_stavu = [47 for x in range(self.pocet_stavov)]
			for k in range (0, self.pocet_stavov):
				pravdep_stavu=[47 for x in range(len(seq))]
				for i in range (0, len(seq)): 
					pravd_stavu[k] = suma(mul(F[i][k],B[i][k]), pravd_stavu[k])
			# pre kazdu sekvenciu si budeme pamatat t[k, l] a e[k, b]
		
			#pocitame t (pomocne pole A[a][k][l] vravi, ze sekvencia a ma tolko prechodov zo stavu k do l) 
			for k in range(0, self.pocet_stavov):
				for l in range(0, self.pocet_stavov):
					temp = [47 for x in range(len(seq))]
					for i in range(0, len(seq)-1):
						temp[i] =mul(mul(F[i][k], self.t[k][l]), mul(self.e[l][int(seq[i+1])],B[i+1][l]))
					temp.sort()
					for i in range(len(seq)-1, -1, -1):
						A[a][k][l] = suma(A[a][k][l], temp[i])
					A[a][k][l] = mul(A[a][k][l], -pravd_stavu[k])
			#pocitame e (teda E[a][k][b] - sekvencia a ma taku emisiu v stave k pismenka b)
			for k in range(0, self.pocet_stavov):
				for b in range(-11, 12):
					for i in range(0, len(seq)):
						if (int(seq[i]) == b):
							E[a][k][b] = suma(E[a][k][b], mul(F[i][k],B[i][k]))
					E[a][k][b] = mul(E[a][k][b], -pravd_stavu[k])
		new_s = self.s
		new_t= [[47 for x in range(self.pocet_stavov)] for x in range(self.pocet_stavov)]
		new_e= [[47 for x in range(-11,12)] for x in range(self.pocet_stavov)]
					


		#spocitame new_t tak, ze zosumujeme A[a][k][l] cez vsetky sekvencie a vydelime sumou A[k][l] cez vsetky l
		all_t = [[47 for x in range(self.pocet_stavov)] for x in range(self.pocet_stavov)]
		for k in range(0, self.pocet_stavov):
			for l in range(0, self.pocet_stavov):
				temp = [47 for x in range(self.lines)]
				for u in range(0, self.lines):
					#print('u', u)
					temp[u] = A[u][k][l]
				temp.sort()
				for u in range(0, self.lines):
					all_t[k][l] = suma(all_t[k][l], temp[u])
		for k in range(0, self.pocet_stavov):
			temp = [47 for x in range(self.pocet_stavov)]
			for l in range(0, self.pocet_stavov):
				temp[l] = all_t[k][l]
			temp.sort()
			sucet_cez_l = 47
			for l in range(0, self.pocet_stavov):
				sucet_cez_l = suma(sucet_cez_l, temp[l])
			for l in range(0, self.pocet_stavov):
				new_t[k][l] = mul(all_t[k][l], -sucet_cez_l)

		#spocitame new_e tak, ze zosumujeme E[a][k][b] cez vsetky sekvencie a vydelime sumou E[k][b] cez vsetky b
		all_e = [[47 for x in range(-11, 12)] for x in range(self.pocet_stavov)]
		for k in range(0, self.pocet_stavov):
			for b in range(-11, 12):
				temp = [47 for x in range(self.lines)]
				for u in range(0, self.lines):
					temp[u] = E[u][k][b]
				temp.sort()
				
				for u in range(0, self.lines):
					all_e[k][b] = suma(suma(all_e[k][b], temp[u]), math.log(0.001)) # pripocitavam pseudocount
		for k in range(0, self.pocet_stavov):
			temp = [47 for x in range(-11,12)]
			for b in range(-11, 12):
				temp[b] = all_e[k][b]
			temp.sort()
			sucet_cez_b = 47
			for b in range(-11, 12):
				sucet_cez_b = suma(sucet_cez_b, temp[b])
			for b in range(-11, 12):
				new_e[k][b] = (mul(all_e[k][b], -sucet_cez_b)) 
	
		self.s = copy.deepcopy(new_s)
		self.t = copy.deepcopy(new_t)
		self.e = copy.deepcopy(new_e)

	def natrenuj(self,cykly):
		for i in range (0, cykly): #iteruje t a e
			self.new_t_e()


	def Viterbi(self,sekvencia):
		seq=sekvencia
		V = [[47 for x in range (self.pocet_stavov)] for x in range(len(seq))]
		C = [[-1 for x in range (self.pocet_stavov)] for x in range(len(seq))] # C[i][j]=k , ked do V[i][j] sa prislo zo stavu k (k vysledku to nepotrebujeme)
		maximum=0
		for j in range (0, self.pocet_stavov):
			V[0][j] = mul(self.s[j],self.e[j][int(seq[0])])
		for i in range(1, len(seq)):
			for j in range (0, self.pocet_stavov):
				C[i][j] = -1
				maximum = -9999999999  # hladame najvacsie zaporne cislo (alebo az 0 = vtedy je pravdepodobnost 1)
				for k in range (0, self.pocet_stavov):
					if (mul(mul(V[i-1][k],self.t[k][j]),self.e[j][int(seq[i])])>maximum):
						maximum = mul(mul(V[i-1][k],self.t[k][j]),self.e[j][int(seq[i])])
						V[i][j] = maximum
						C[i][j] = k
		
		#najdeme maximum v poslednom stlpci
		maximum = -999999999
		for j in range(0, self.pocet_stavov):
				if (V[len(seq)-1][j])>maximum:
					maximum=V[len(seq)-1][j]
		return maximum




def vstup(a, b, c): #a, b, c su subory
	# v a je vstup v notach
	# v b su prepisane noty na ciselka
	# v c su rozdiely vysok not
	noty = open(a, 'r')
	abs_cisla = open(b, 'w')
	for line in noty:
		pesnicka=line.split()
		for i in range(0, len(pesnicka)):
			if pesnicka[i] == 'c':
				pesnicka[i] = '0'
			elif pesnicka[i] == 'cis' or pesnicka[i] == 'db':
				pesnicka[i] = '1'
			elif pesnicka[i] == 'd':
				pesnicka[i] = '2'
			elif pesnicka[i] == 'dis' or pesnicka[i] == 'eb':
				pesnicka[i] = '3'
			elif pesnicka[i] == 'e':
				pesnicka[i] = '4'
			elif pesnicka[i] == 'f':
				pesnicka[i] = '5'
			elif pesnicka[i] == 'fis' or pesnicka[i] == 'gb':
				pesnicka[i] = '6'
			elif pesnicka[i] == 'g':
				pesnicka[i] = '7'
			elif pesnicka[i] == 'gis' or pesnicka[i] == 'ab':
				pesnicka[i] = '8'
			elif pesnicka[i] == 'a' :
				pesnicka[i] = '9'
			elif pesnicka[i] == 'ais' or pesnicka[i] == 'hb':
				pesnicka[i] = '10'
			elif pesnicka[i] == 'h':
				pesnicka[i] = '11'
			elif pesnicka[i] == 'c2':
				pesnicka[i] = '12'
			elif pesnicka[i] == 'cis2' or pesnicka[i] == 'db2':
				pesnicka[i] = '13'
			elif pesnicka[i] == 'd2':
				pesnicka[i] = '14'
			elif pesnicka[i] == 'dis2' or pesnicka[i] == 'eb2':
				pesnicka[i] = '15'
			elif pesnicka[i] == 'e2':
				pesnicka[i] = '16'
			elif pesnicka[i] == 'f2':
				pesnicka[i] = '17'
			elif pesnicka[i] == 'fis2' or pesnicka[i] == 'gb2':
				pesnicka[i] = '18'
			elif pesnicka[i] == 'g2':
				pesnicka[i] = '19'
			elif pesnicka[i] == 'gis2' or pesnicka[i] == 'ab2':
				pesnicka[i] = '20'
			elif pesnicka[i] == 'a2':
				pesnicka[i] = '21'
			elif pesnicka[i] == 'ais2' or pesnicka[i] == 'hb2':
				pesnicka[i] = '22'
			abs_cisla.write(pesnicka[i] + ' ')
		abs_cisla.write('\n')

	noty.close()
	abs_cisla.close()
	
	abs_cisla = open(b, 'r')
	rel_cisla = open(c, 'w')
	for line in abs_cisla:
		pesnicka=line.split()
		pesnicka2={}
		for i in range(1, len(pesnicka)):
			pesnicka2[i] = str(int(pesnicka[i])-int(pesnicka[i-1]))
			rel_cisla.write(pesnicka2[i]+' ')
		rel_cisla.write('\n')
	abs_cisla.close()
	rel_cisla.close()

def HMM1(test, vystup):
	print('Vyrabam slovenske hmm1...')
	slovenske_hmm = HMM()
	slovenske_hmm.nastav([0.5, 0.3, 0.2], [[0.04347 for x in range(23)] for x in range(3)], [[0.2, 0.8,0],[0, 0.2,0.8], [0, 0, 1]], 'slovenske_rel_cisla.txt', 3)
	slovenske_hmm.natrenuj(10)
	print('Vyrabam nemecke hmm1...')
	nemecke_hmm = HMM()
	nemecke_hmm.nastav([0.5, 0.3, 0.2], [[0.04347 for x in range(23)] for x in range(3)], [[0.2, 0.8,0],[0, 0.2,0], [0, 0, 1]], 'nemecke_rel_cisla.txt', 3)
	nemecke_hmm.natrenuj(10)
	print('Vyrabam cinske hmm1...')
	china_hmm = HMM()
	china_hmm.nastav([0.5, 0.3, 0.2], [[0.04347 for x in range(23)] for x in range(3)],[[0.2, 0.8,0],[0, 0.2,0], [0, 0, 1]], 'china_rel_cisla.txt', 3)
	china_hmm.natrenuj(10)
	test = open(test, 'r')
	vystup = open (vystup, 'w')
	for line in test:
		seq = line.split()
		a=slovenske_hmm.Viterbi(seq)
		b=nemecke_hmm.Viterbi(seq)
		c=china_hmm.Viterbi(seq)
		if (a>b) and (a>c):
			vystup.write('s\n')
			#print(a)
		elif (b>a) and (b>c):
			vystup.write('n\n')
			#print(b)
		else:
			vystup.write('c\n')
			#print(c)

def HMM2(test, vystup):
	print('Vyrabam slovenske hmm2...')
	slovenske_hmm = HMM()
	slovenske_hmm.nastav([0.5, 0.3, 0.2], [[0.04347 for x in range(23)] for x in range(3)], [[0.2, 0.4,0.4],[0, 0.5,0.5], [0, 0, 1]], 'slovenske_rel_cisla.txt', 3)
	slovenske_hmm.natrenuj(10)
	print('Vyrabam nemecke hmm2...')
	nemecke_hmm = HMM()
	nemecke_hmm.nastav([0.5, 0.3, 0.2], [[0.04347 for x in range(23)] for x in range(3)], [[0.2, 0.4,0.4],[0, 0.5,0.5], [0, 0, 1]], 'nemecke_rel_cisla.txt', 3)
	nemecke_hmm.natrenuj(10)
	print('Vyrabam cinske hmm2...')
	china_hmm = HMM()
	china_hmm.nastav([0.5, 0.3, 0.2], [[0.04347 for x in range(23)] for x in range(3)], [[0.2, 0.4,0.4],[0, 0.5,0.5], [0, 0, 1]], 'china_rel_cisla.txt', 3)
	china_hmm.natrenuj(10)
	test = open(test, 'r')
	vystup = open (vystup, 'w')
	for line in test:
		seq = line.split()
		a=slovenske_hmm.Viterbi(seq)
		b=nemecke_hmm.Viterbi(seq)
		c=china_hmm.Viterbi(seq)
		if (a>b) and (a>c):
			vystup.write('s\n')
			#print(a)
		elif (b>a) and (b>c):
			vystup.write('n\n')
			#print(b)
		else:
			vystup.write('c\n')
			#print(c)

def HMM3(test, vystup):
	print('Vyrabam slovenske hmm3...')
	slovenske_hmm = HMM()
	slovenske_hmm.nastav([0.5, 0.3, 0.2], [[0.04347 for x in range(23)] for x in range(3)], [[0.3, 0.4,0.3],[0.3, 0.3,0.4], [0.3, 0.3, 0.4]], 'slovenske_rel_cisla.txt', 3)
	slovenske_hmm.natrenuj(10)
	print('Vyrabam nemecke hmm3...')
	nemecke_hmm = HMM()
	nemecke_hmm.nastav([0.5, 0.3, 0.2], [[0.04347 for x in range(23)] for x in range(3)], [[0.3, 0.4,0.3],[0.3, 0.3,0.4], [0.3, 0.3, 0.4]], 'nemecke_rel_cisla.txt', 3)
	nemecke_hmm.natrenuj(10)
	print('Vyrabam cinske hmm3...')
	china_hmm = HMM()
	china_hmm.nastav([0.5, 0.3, 0.2], [[0.04347 for x in range(23)] for x in range(3)], [[0.3, 0.4,0.3],[0.3, 0.3,0.4], [0.3, 0.3, 0.4]], 'china_rel_cisla.txt', 3)
	china_hmm.natrenuj(10)
	test = open(test, 'r')
	vystup = open (vystup, 'w')
	for line in test:
		seq = line.split()
		a=slovenske_hmm.Viterbi(seq)
		b=nemecke_hmm.Viterbi(seq)
		c=china_hmm.Viterbi(seq)
		if (a>b) and (a>c):
			vystup.write('s\n')
			#print(a)
		elif (b>a) and (b>c):
			vystup.write('n\n')
			#print(b)
		else:
			vystup.write('c\n')
			#print(c)



def HMM4(test, vystup):
	print('Vyrabam slovenske hmm4...')
	slovenske_hmm = HMM()
	slovenske_hmm.nastav([0.5, 0.2, 0.1, 0.1, 0.1, 0.1], [[0.04347 for x in range(23)] for x in range(6)], [[0.4, 0.2, 0.1, 0.1, 0.1, 0.1],[0, 0.6, 0.4, 0, 0, 0], [0, 0, 0.6, 0.4, 0, 0], [0, 0, 0, 0.6, 0.4, 0], [0,0, 0, 0, 0.6, 0.4], [0, 0, 0, 0, 0, 1]], 'slovenske_rel_cisla.txt', 6)
	slovenske_hmm.natrenuj(10)
	print('Vyrabam nemecke hmm4...')
	nemecke_hmm = HMM()
	nemecke_hmm.nastav([0.5, 0.1, 0.1, 0.1, 0.1, 0.1], [[0.04347 for x in range(23)] for x in range(6)], [[0.4, 0.2, 0.1, 0.1, 0.1, 0.1],[0, 0.6, 0.4, 0, 0, 0], [0, 0, 0.6, 0.4, 0, 0], [0, 0, 0, 0.6, 0.4, 0], [0,0, 0, 0, 0.6, 0.4], [0, 0, 0, 0, 0, 1]], 'slovenske_rel_cisla.txt', 6)
	nemecke_hmm.natrenuj(10)
	print('Vyrabam cinske hmm4...')
	china_hmm = HMM()
	china_hmm.nastav([0.5, 0.1, 0.1, 0.1, 0.1, 0.1], [[0.04347 for x in range(23)] for x in range(6)], [[0.4, 0.2, 0.1, 0.1, 0.1, 0.1],[0, 0.6, 0.4, 0, 0, 0], [0, 0, 0.6, 0.4, 0, 0], [0, 0, 0, 0.6, 0.4, 0], [0,0, 0, 0, 0.6, 0.4], [0, 0, 0, 0, 0, 1]], 'slovenske_rel_cisla.txt', 6)
	china_hmm.natrenuj(10)
	test = open(test, 'r')
	vystup = open (vystup, 'w')
	for line in test:
		seq = line.split()
		a=slovenske_hmm.Viterbi(seq)
		b=nemecke_hmm.Viterbi(seq)
		c=china_hmm.Viterbi(seq)
		if (a>b) and (a>c):
			vystup.write('s\n')
			#print(a)
		elif (b>a) and (b>c):
			vystup.write('n\n')
			#print(b)
		else:
			vystup.write('c\n')
			#print(c)




def HMM5(test, vystup):
	print('Vyrabam slovenske hmm5...')
	slovenske_hmm = HMM()
	slovenske_hmm.nastav([0.4, 0.2, 0.1, 0.1, 0.1, 0.1], [[0.04347 for x in range(23)] for x in range(6)], [[0.4, 0.2, 0.1, 0.1, 0.1, 0.1],[0, 0.4, 0.2, 0.2, 0.1, 0.1], [0, 0, 0.5, 0.3, 0.1, 0.1], [0, 0, 0, 0.5, 0.3, 0.2], [0,0, 0, 0, 0.6, 0.4], [0, 0, 0, 0, 0, 1]], 'slovenske_rel_cisla.txt', 6)
	print(slovenske_hmm.t)
	slovenske_hmm.natrenuj(10)
	print('Vyrabam nemecke hmm5...')
	nemecke_hmm = HMM()
	nemecke_hmm.nastav([0.4, 0.2, 0.1, 0.1, 0.1, 0.1], [[0.04347 for x in range(23)] for x in range(6)], [[0.4, 0.2, 0.1, 0.1, 0.1, 0.1],[0, 0.4, 0.2, 0.2, 0.1, 0.1], [0, 0, 0.5, 0.3, 0.1, 0.1], [0, 0, 0, 0.5, 0.3, 0.2], [0, 0, 0, 0, 0.6, 0.4], [0, 0, 0, 0, 0, 1]], 'nemecke_rel_cisla.txt', 6)
	print('Vyrabam cinske hmm5...')
	china_hmm = HMM()
	china_hmm.nastav([0.4, 0.2, 0.1, 0.1, 0.1, 0.1], [[0.04347 for x in range(23)] for x in range(6)], [[0.4, 0.2,0.1, 0.1, 0.1, 0.1],[0, 0.4,0.2, 0.2, 0.1, 0.1], [0, 0,0.5, 0.3, 0.1, 0.1], [0, 0,0, 0.5, 0.3, 0.2], [0, 0,0, 0, 0.6, 0.4], [0, 0, 0, 0, 0, 1]], 'china_rel_cisla.txt', 6)
	china_hmm.natrenuj(10)
	test = open(test, 'r')
	vystup = open (vystup, 'w')
	for line in test:
		seq = line.split()
		a=slovenske_hmm.Viterbi(seq)
		b=nemecke_hmm.Viterbi(seq)
		c=china_hmm.Viterbi(seq)
		if (a>b) and (a>c):
			vystup.write('s\n')
			#print(a)
		elif (b>a) and (b>c):
			vystup.write('n\n')
			#print(b)
		else:
			vystup.write('c\n')
			#print(c)

def HMM6(test, vystup):
	print('Vyrabam slovenske hmm6...')
	slovenske_hmm = HMM()
	slovenske_hmm.nastav([0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [[0.04347 for x in range(23)] for x in range(6)], [[0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166]], 'slovenske_rel_cisla.txt', 6)
	slovenske_hmm.natrenuj(10)
	print('Vyrabam nemecke hmm6...')
	nemecke_hmm = HMM()
	nemecke_hmm.nastav([0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [[0.04347 for x in range(23)] for x in range(6)], [[0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166]], 'slovenske_rel_cisla.txt', 6)
	nemecke_hmm.natrenuj(10)
	print('Vyrabam cinske hmm6...')
	china_hmm = HMM()
	china_hmm.nastav([0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [[0.04347 for x in range(23)] for x in range(6)], [[0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166], [0.175, 0.166, 0.166, 0.166, 0.166, 0.166]], 'slovenske_rel_cisla.txt', 6)
	test = open(test, 'r')
	vystup = open (vystup, 'w')
	for line in test:
		seq = line.split()
		a=slovenske_hmm.Viterbi(seq)
		b=nemecke_hmm.Viterbi(seq)
		c=china_hmm.Viterbi(seq)
		if (a>b) and (a>c):
			vystup.write('s\n')
			#print(a)
		elif (b>a) and (b>c):
			vystup.write('n\n')
			#print(b)
		else:
			vystup.write('c\n')
			#print(c)




#transformujeme vstup na relativne cisla, teda bereme do uvahy len rozdiely medzi susednymi tonmi v kazdej pesnicke	
vstup('slovenske_noty.txt', 'slovenske_abs_cisla.txt', 'slovenske_rel_cisla.txt')
vstup('nemecke_noty.txt', 'nemecke_abs_cisla.txt', 'nemecke_rel_cisla.txt')
vstup('china_noty.txt', 'china_abs_cisla.txt', 'china_rel_cisla.txt')
vstup('test.txt', 'test_abs.txt', 'test_rel.txt')

HMM1('test_rel.txt', 'vystup1.txt')
HMM2('test_rel.txt', 'vystup2.txt')
HMM3('test_rel.txt', 'vystup3.txt')
HMM4('test_rel.txt', 'vystup4.txt')
HMM5('test_rel.txt', 'vystup5.txt')
HMM6('test_rel.txt', 'vystup6.txt')

