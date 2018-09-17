#dssp -i 1bta.pdb

import sys
import math
import string

def analyse_pdb(pdbfile):
	"""Cette fonction lit le fichier pdb place en argument, extrait et associe pour chaque residu du fichier, les coordonnees des atomes "N", "C", "O" , "H" etudies .

    Parameters
    ----------
    pdbfile : fichier .pdb
        Ce fichier contient les coordonnees des atomes d une proteine.

    Returns
    -------
    type : dictionnaire
        L objet retourne contient pour chaque residu du fichier, les coordonnees x, y et z des atomes "N", "C", "O" , "H" etudies.
"""
	dico_residu={}
	Tout_atome={}
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
						atome_coord={}
						atome_coord[line[12:16].strip()]=dico_coord
						Tout_atome[num_residu]=atome_coord
	Tout_atome[num_residu]=atome_coord
	return Tout_atome

def calcul_distance(coord1, coord2):
	"""Cette fonction calcule la distance en angström entre deux atomes places en argument.

    Parameters
    ----------
    coord1 : coordonnées d'un atome du pdb
        Coordonnees x, y et z de cet atome
    coord2 : coordonnées d'un autre atome du pdb
        Coordonnees x, y et z de cet atome

    Returns
    -------
      type int
          L objet retourne contient la distance en angström entre deux atomes places en argument.
"""
	distance = math.sqrt((coord2['x'] - coord1['x'])**2 + (coord2['y'] - coord1['y'])**2 + (coord2['z'] - coord1['z'])**2)
	return distance

def calcul_energie(residu1, residu2):
	"""Cette fonction calcule l energie d interaction electrostatique entre deux residus places en argument.

    Parameters
    ----------
    residu1 : type dictionnaire
		Contient les coordonnées x, y et z des atomes "N", "C", "O" , "H" du residu1
    residu2 : type dictionnaire
    	Contient les coordonnées x, y et z des atomes "N", "C", "O" , "H" du residu2

    Returns
    -------
    type int
        L objet retourne contient l energie d interaction electrostatique obtenue entre les deux residus.

"""
	q1 = 0.42
	q2 = 0.20
	if 'H' not in residu2:
		E = 0
	else:
		rON = calcul_distance(residu1['O'], residu2['N'])
		rCH = calcul_distance(residu1['C'], residu2['H'])
		rOH = calcul_distance(residu1['O'], residu2['H'])
		rCN = calcul_distance(residu1['C'], residu2['N'])
		E = q1*q2*(1/rON + 1/rCH - 1/rOH - 1/rCN) * 332
	return E

def liste_feuillets(dico_liaisons):
	"""Cette fonction renvoie la position des feuillets Betas potentiels.

    Parameters
    ----------
    dico_liaisons : type dictionnaire
        Dictionnaire contenant la liste de tous les résidus ayant etabli des liaisons hydrogenes sous forme de parallele bridge ou anti-parallele bridge.

    Returns
    -------
    type dictionnaire
        L objet retourne contient la position des feuillets Betas potentiels ainsi que le numero des residus intervenant.
"""
	feuillets={}
	for key in dico_liaisons:
		if key + 1 in dico_liaisons or key - 1 in dico_liaisons:
			feuillets[key]=dico_liaisons[key]
	return feuillets

def liste_helices(dico_liaisons):
	"""Cette fonction renvoie la position des helices (alpha, pi et 310) potentielles.

    Parameters
    ----------
    dico_liaisons : type dictionnaire
        Dictionnaire contenant la liste de tous les résidus ayant établi des n-turns.

    Returns
    -------
    type dictionnaire
        L objet retourne la position des helices potentielles (helice alpha, pi ou 310).
"""
	helices={}
	for key in dico_liaisons:
		if key in dico_liaisons and key + 1 in dico_liaisons:
			helices[key]=1
			helices[key+1]=1
			helices[key+2]=1
			helices[key+3]=1
			helices[key+4]=1
	return helices

def calcul_helices(liste_atome):
	"""
      Cette fonction identifie les n-turns et renvoie la position des helices potentielles.

    Parameters
    ----------
    liste_atome : type dicitonnaire
        Liste_atome contient les coordonnees des atomes importants pour chaque residu.

    Returns
    -------
    type dictionnaire
        L objet retourne contient les positions des helices potentielles.

"""
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
				E = calcul_energie(liste_atome[key], liste_atome[key+k])
				if E < -0.5 :
					if k == 3 : liaisons_310[key] = E
					elif k == 4 : liaisons_alpha[key] = E
					elif k == 5 : liaisons_pi[key] = E
	helices_alpha=liste_helices(liaisons_alpha)
	helices_310=liste_helices(liaisons_310)
	helices_pi=liste_helices(liaisons_pi)
	helices_totales['helices_pi']=helices_pi
	helices_totales['helices_310']=helices_310
	helices_totales['helices_alpha']=helices_alpha
	return helices_totales


def calcul_feuillet(liste_atome):
	"""Cette fonction renvoie l emplacement des feuillets Betas potentiels a partir de la liste des residus.

    Parameters
    ----------
    liste_atome : type dictionnaire
        Liste_atome contient les coordonnees des atomes importants pour chaque residu.

    Returns
    -------
    type dictionnaire
        L objet retourne contient la position des feuillets Betas potentiels.

"""
	parallele_bridge={}
	antiparallele_bridge={}
	feuillets_paralleles={}
	feuillets_antiparalleles={}
	feuillets_total={}
	indice=list(liste_atome.keys())[0]
	for i in range(indice + 1, indice + len(liste_atome)-2):
		for j in range (i + 2, indice + len(liste_atome)-2):
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
				if i in parallele_bridge:
					BP=[]
					BP.append(parallele_bridge[i])
					BP.append(j)
					parallele_bridge[i]=BP
				elif i not in parallele_bridge:
					parallele_bridge[i]=j
				continue
			if (A1 < -0.5 and A2 < -0.5) or (A3 < -0.5 and A4 < -0.5):
				if i in parallele_bridge:
					BP=[]
					BP.append(antiparallele_bridge[i])
					BP.append(j)
					antiparallele_bridge[i]=BP
				elif i not in parallele_bridge:
					antiparallele_bridge[i]=j
				continue
	print(antiparallele_bridge)
	feuillets_paralleles=liste_feuillets(parallele_bridge)
	feuillets_antiparalleles=liste_feuillets(antiparallele_bridge)
	feuillets_total['feuillets_paralleles']=feuillets_paralleles
	feuillets_total['feuillets_antiparalleles']=feuillets_antiparalleles
	return feuillets_total

def affichage(liste_atome, helices_totales, feuillets):
	"""Cette fonction permet d'afficher, pour chaque residu du fichier pdb, le type de structure qu il serait susceptible de former.

    Parameters
    ----------
    liste_atome : type dictionnaire
        Liste des coordonnees des atomes d interet pour chaque residu.
    helices_totales : type dictionnaire
        Liste de la position des helices potentielles.
    feuillets : type dictionnaire
        Description of parameter `feuillets`.

    Returns
    -------
    type
        Description of returned object.

"""
	Structures=['helices_alpha', 'helices_310', 'helices_pi', 'feuillets_paralleles']
	print ('RESIDUE', '  AA', '  STRUCTURE', 'BP1', 'BP2')
	liste_lower = list(string.ascii_lowercase)
	liste_upper = list(string.ascii_uppercase)
	compteur=0
	for key in liste_atome:
		if key in helices_totales['helices_alpha']:
			print(key, "      ", liste_atome[key]['C']['nom_residu'], "     H")
		elif key in helices_totales['helices_310']:
			print(key, "      ", liste_atome[key]['C']['nom_residu'], "     G")
		elif key in helices_totales['helices_pi']:
			print(key, "      ", liste_atome[key]['C']['nom_residu'], "     I")
		elif key in feuillets['feuillets_paralleles']:
			print(key, "      ", liste_atome[key]['C']['nom_residu'], "     Ep", "    ", feuillets['feuillets_paralleles'][key])
		elif key in feuillets['feuillets_antiparalleles']:
			print(key, "      ", liste_atome[key]['C']['nom_residu'], "     Ea", "    ", feuillets['feuillets_antiparalleles'][key])
		elif key not in helices_totales and key not in feuillets:
			print(key, "      ", liste_atome[key]['C']['nom_residu'], "      ")
	return

liste=[]
if len(sys.argv) == 1:
	print("Veuillez entrer un argument")
elif len(sys.argv) == 2 :
	liste=analyse_pdb(sys.argv[1])
	Helices=calcul_helices(liste)
	Feuillets = calcul_feuillet(liste)
	affichage(liste, Helices, Feuillets)
elif len(sys.argv) > 2:
	print("Vous avez entre trop d arguments")
