#dssp -i 1bta.pdb

import math

def parse_pdb(pdbfile):
	dico_residu={}
	liste_atome=[] 
	#num_residu=1
	with open(pdbfile, 'r') as inputfile:
		dico_atome={} 
		for line in inputfile:
			if line[9:12].strip()=='1' and line[0:6].strip() == "ATOM":
				num_residu=int(line[22:26].strip())
			if line[0:6].strip() == "ATOM" and line[12:16].strip() in ["N", "C", "O" , "H"]:
					dico_coord={}
					dico_coord['nom_residu'] = line[17:20].strip()
					dico_coord['chaine'] =  line[21:23].strip()
					dico_coord['num_residu'] = line[22:26].strip()
					dico_coord['x'] = float(line[30:38])
					dico_coord['y'] = float(line[38:46])
					dico_coord['z'] = float(line[46:54])
					if line[12:16].strip() not in dico_atome and int(dico_coord['num_residu']) == num_residu:						
						dico_atome[line[12:16].strip()]=dico_coord
						
					else:						
						num_residu=int(line[22:26].strip())
						liste_atome.append(dico_atome)
						dico_atome={}
						dico_atome[line[12:16].strip()]=dico_coord
	
	liste_atome.append(dico_atome)
	print(liste_atome[16])
	return liste_atome

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
			liaisons_consecutives[key+1]=dico_liaisons[key]
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
		for i in range(len(liste_atome) - k):
			if liste_atome[i+k]['C']['nom_residu']=='PRO':
				E = 0
			else :  
				E = calcul_energie(liste_atome[i], liste_atome[i+k])
				if E < -0.5 :
					if k == 3 : liaisons_310[i] = E
					elif k == 4 : liaisons_alpha[i] = E 
					elif k == 5 : liaisons_pi[i] = E
	for key in liaisons_alpha:
		print(key,liste_atome[key]['C']['nom_residu'])	
	helices_alpha=liaisons_consecutives(liaisons_alpha)
	helices_310=liaisons_consecutives(liaisons_310)
	helices_pi=liaisons_consecutives(liaisons_pi)
	helices_totales['helices_pi']=helices_pi
	helices_totales['helices_310']=helices_310
	helices_totales['helices_alpha']=helices_alpha
	print(helices_totales['helices_alpha'])
	return helices_alpha

def calcul_feuillet(liste_atome):
	Energie_feuillet_parallele={}
	Energie_feuillet_antiparallele={}
	feuillets_paralleles={}
	feuillets_antiparalleles={}
	for i in range(1, len(liste_atome)-1):
		for j in range (i+2, len(liste_atome)-1):
			#feuillets parallèles
			P1 = calcul_energie(liste_atome[i-1], liste_atome[j])
			P2 = calcul_energie(liste_atome[j], liste_atome[i+1])
			P3 = calcul_energie(liste_atome[j-1], liste_atome[i])
			P4 = calcul_energie(liste_atome[i], liste_atome[j+1])
			#feuillets anti-parallèles
			k = len(liste_atome) - i -1
			A1 = calcul_energie(liste_atome[i], liste_atome[k])
			A2 = calcul_energie(liste_atome[k], liste_atome[i])
			A3 = calcul_energie(liste_atome[i-1], liste_atome[k+1])
			A4 = calcul_energie(liste_atome[k+1], liste_atome[i-1])
			if (P1 < -0.5 and P2 < -0.5) or (P3 < -0.5 and P4 < -0.5):
			#	print(i,j, liste_atome[88])
			#	for k in (i, len(liste_atome)-1):
			#		for w in (j, len(liste_atome)):
			#			P1 = calcul_energie(liste_atome[k-1], liste_atome[w])
			#			P2 = calcul_energie(liste_atome[w], liste_atome[k+1])
			#			if P1 and P2 < -0.5:
				Energie_feuillet_parallele[i+1]=j+1
				continue
			if A1 < -0.5 and A2 < -0.5 or A3 < -0.5 and A4 < -0.5:
				Energie_feuillet_antiparallele[i]=k
				continue
	feuillets_paralleles=liaisons_consecutives(Energie_feuillet_parallele)
	feuillets_antiparalleles=liaisons_consecutives(Energie_feuillet_antiparallele)
	print(feuillets_paralleles)
	print(feuillets_antiparalleles)
	return Energie_feuillet_parallele


liste=[]
liste=parse_pdb("1mh1_withH.pdb")
Energie_helice=calcul_helices(liste)
Energie_feuillet = calcul_feuillet(liste)






            
