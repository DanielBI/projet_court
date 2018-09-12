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
						dico_atome={}
	print(liste_atome[0])
	return liste_atome

#dico_residu[line[22:26].strip()]=dico_atome[line[14:16].strip()]=dico_coord


test=parse_pdb("1bta.pdb")






            
