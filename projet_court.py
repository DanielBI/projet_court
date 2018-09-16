#dssp -i 1bta.pdb

import math

def parse_pdb(pdbfile):
	dico_residu={}
	Tout_atome={} 
	#num_residu=1
	with open(pdbfile, 'r') as inputfile:
		atome_coord={} 
		for line in inputfile:
			if line[10:13].strip()=='1' and line[0:6].strip() == "ATOM":
				num_residu=int(line[22:26].strip())
			if line[0:6].strip() == "ATOM" and line[12:16].strip() in ["N", "C", "O" , "H"]:
					dico_coord={}
					dico_coord['nom_residu'] = line[17:20].strip()
					dico_coord['chaine'] =  line[21:23].strip()
					dico_coord['num_residu'] = line[22:26].strip()
					dico_coord['x'] = float(line[30:38])
					dico_coord['y'] = float(line[38:46])
					dico_coord['z'] = float(line[46:54])
					if line[12:16].strip() not in atome_coord and int(dico_coord['num_residu']) == num_residu:												
						atome_coord[line[12:16].strip()]=dico_coord
				
						
					else:
						Tout_atome[num_residu]=atome_coord						
						num_residu=int(line[22:26].strip())
						Tout_atome[num_residu]=atome_coord						
						atome_coord={}
						atome_coord[line[12:16].strip()]=dico_coord
	
	Tout_atome[num_residu]=atome_coord
	return Tout_atome

def calcul_distance(coord1, coord2):
	distance = math.sqrt((coord2['x'] - coord1['x'])**2 + (coord2['y'] - coord1['y'])**2 + (coord2['z'] - coord1['z'])**2)
	return distance

def calcul_energie(residu1, residu2):
	q1 = 0.42 
	q2 = 0.20
	if 'H' not in residu2 or 'H' not in residu1:
		E = 0
	else:
		rON = calcul_distance(residu1['O'], residu2['N'])
		rCH = calcul_distance(residu1['C'], residu2['H'])
		rOH = calcul_distance(residu1['O'], residu2['H'])
		rCN = calcul_distance(residu1['C'], residu2['N'])
		E = q1*q2*(1/rON + 1/rCH - 1/rOH - 1/rCN) * 332
	return E

def liaisons_consecutives(dico_liaisons):
	liaisons_consecutives={}
	for key in dico_liaisons:
		if key + 1 in dico_liaisons or key - 1 in dico_liaisons:
			liaisons_consecutives[key]=dico_liaisons[key]
	return liaisons_consecutives
			

def calcul_helices(liste_atome):
	liaisons_alpha={}
	liaisons_310={}
	liaisons_pi={}
	helices_alpha={}
	helices_310={}
	helices_pi={}
	helices_totales={}
	for k in [3, 4, 5]: 
		for key in liste_atome:
			if key + k in liste_atome:
				if liste_atome[key+k]['C']['nom_residu']=='PRO':
					E = 0
				else :  
					E = calcul_energie(liste_atome[key], liste_atome[key+k])
					if E < -0.5 :
						if k == 3 : liaisons_310[key] = E
						elif k == 4 : liaisons_alpha[key] = E 
						elif k == 5 : liaisons_pi[key] = E
	helices_alpha=liaisons_consecutives(liaisons_alpha)
	helices_310=liaisons_consecutives(liaisons_310)
	helices_pi=liaisons_consecutives(liaisons_pi)
	helices_totales['helices_pi']=helices_pi
	helices_totales['helices_310']=helices_310
	helices_totales['helices_alpha']=helices_alpha
	return helices_totales


def calcul_feuillet(liste_atome):
	Energie_feuillet_parallele={}
	Energie_feuillet_antiparallele={}
	feuillets_paralleles={}
	feuillets_antiparalleles={}
	feuillets_total={}
	indice=list(liste_atome.keys())[0]
	for i in range(indice + 1, indice + len(liste_atome)-2):
		for j in range (i+4, indice + len(liste_atome)-2):
			#feuillets parallèles
			P1 = calcul_energie(liste_atome[i-1], liste_atome[j])
			P2 = calcul_energie(liste_atome[j], liste_atome[i+1])
			P3 = calcul_energie(liste_atome[j-1], liste_atome[i])
			P4 = calcul_energie(liste_atome[i], liste_atome[j+1])
			#feuillets anti-parallèles
			A1 = calcul_energie(liste_atome[i], liste_atome[j])
			A2 = calcul_energie(liste_atome[j], liste_atome[i])
			A3 = calcul_energie(liste_atome[i-1], liste_atome[j+1])
			A4 = calcul_energie(liste_atome[i+1], liste_atome[j-1])

			if (P1 < -0.5 and P2 < -0.5) or (P3 < -0.5 and P4 < -0.5):	
				Energie_feuillet_parallele[i]=j
				continue
			if (A1 < -0.5 and A2 < -0.5) or (A3 < -0.5 and A4 < -0.5):
				Energie_feuillet_antiparallele[i]=j
	feuillets_paralleles=liaisons_consecutives(Energie_feuillet_parallele)
	feuillets_antiparalleles=liaisons_consecutives(Energie_feuillet_antiparallele)
	feuillet_total['feuillets_paralleles']=feuillets_paralleles
	feuillet_total['feuillets_antiparalleles']=feuillets_antiparalleles
	return feuillets_paralleles

def affichage(liste_atome, helices_totales, feuillets):
	Helices=['helices_alpha', 'helices_310', 'helices_pi']
	feuillets=['
	print ('RESIDUE', 'AA', 'STRUCTURE', 'BP1')
	for key in liste_atome:
		for type_helice in Helices:
			if key in helices_totales[type_helice]:
				if type_helice=='helices_alpha':
					print(key, liste_atome[key]['C']['nom_residu'], "      H")
				elif type_helice=='helices_310':
					print(key, liste_atome[key]['C']['nom_residu'], "      G")
				elif type_helice=='helices_pi':
					print(key, liste_atome[key]['C']['nom_residu'], "      I")
		for 
			elif key in feuillets:
				print(key, liste_atome[key]['C']['nom_residu'], "       E")
			elif key not in helices_totales[type_helice] and key not in feuillets:
				print(key, liste_atome[key]['C']['nom_residu'], "      ")
			break
					
	return 0
	


liste=[]
#liste=parse_pdb("1uol_H.pdb")
liste=parse_pdb("1uol_H.pdb")
Helices=calcul_helices(liste)
Feuillets = calcul_feuillet(liste)
affichage(liste, Helices, Feuillets)






            
