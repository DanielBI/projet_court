import math

atomes_dinteret = ["CA",  "N", "C", "O" , "H"]

 
def parse_pdb(pdbfile):
	dico_residu={}
	liste_residu=[]
	liste_atome=[] 
	num_residu=1
	with open(pdbfile, 'r') as inputfile:
		compteur=0
		dico_atome={} 
		for line in inputfile:
			if line[0:6].strip() == "ATOM" and line[12:16].strip() in atomes_dinteret:
					dico_coord={}
					dico_coord['num_residu'] = line[22:26].strip()
					dico_coord['x'] = float(line[30:38])
					dico_coord['y'] = float(line[38:46])
					dico_coord['z'] = float(line[46:54])
					if line[12:16].strip() not in dico_atome and int(dico_coord['num_residu']) == num_residu:						
						dico_atome[line[12:16].strip()]=dico_coord
						
					else:						
						num_residu=int(line[22:26].strip())
						liste_atome.append(dico_atome)
						print(len(dico_atome))
						dico_atome={}
						dico_atome[line[12:16].strip()]=dico_coord
												
						
	print(liste_atome[1])
	return liste_atome

def calcul_distance(coord1, coord2):
	distance = math.sqrt((coord2['x'] - coord1['x'])**2 + (coord2['y'] - coord1['y'])**2 + (coord2['z'] - coord1['z'])**2)
	return distance


def calcul_Energie(liste_atome):
	liste_energie={}
	q1 = 0.42 * 1.60217662 * 10**-19
	q2 = 0.20 * 1.60217662 * 10**-19
	for i in range(len(liste_atome)):
		rON = calcul_distance(liste_atome[i]['O'],liste_atome[i+4]['N'])
		rCH = calcul_distance(liste_atome[i]['C'],liste_atome[i+4]['H'])
		rOH = calcul_distance(liste_atome[i]['O'],liste_atome[i+4]['H'])
		rCN = calcul_distance(liste_atome[i]['C'],liste_atome[i+4]['N'])
		E = q1*q2*(1/rON + 1/rCH - 1/rOH - 1/rCN) * 332
		print(E)
		liste_energie[i] = E
	return liste_energie

liste=[]
liste=parse_pdb("1bta.pdb")
Energie=calcul_Energie(liste)






            
