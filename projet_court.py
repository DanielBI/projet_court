representative_atoms = ["CA",  "N", "C", "O" , "H"]


def parse_pdb(pdbfile):
	dico_residu={}
	liste_residu=[]
	liste_atome=[] 
	with open(pdbfile, 'r') as inputfile:
		compteur=0
		for line in inputfile:
			current_residu = line[22:26].strip()
			dico_atome={} 
			if line[0:6].strip() == "ATOM" and line[12:16].strip() in representative_atoms:
				new_residu = line[22:26].strip()
				if new_residu == current_residu:
					dico_coord={}
					dico_coord['num_residu'] = line[22:26].strip()
					dico_coord['x'] = float(line[30:38])
					dico_coord['y'] = float(line[38:46])
					dico_coord['z'] = float(line[46:54])
					dico_atome[line[12:16].strip()]=dico_coord
					compteur=compteur + 1
					liste_atome.append(dico_atome)
				#liste_atome.append(dico_atome)
				#compteur=compteur+1
				#if compteur > 5:
				#dico_residu[line[25:28].strip()]=liste_atome
				#compteur=0
				#liste_atome=[]
	num_residu=1
	liste1=[]
	num_residu=1
	print (liste_atome[1]['CA'])
	for i in range(len(liste_atome)):
		for atome in representative_atoms:
			print(atome)
			if liste_atome[i][atome][num_residu] == num_residu:
				liste1.append(liste_atome[i][atome][num_residu])
			num_residu = num_residu + 1 
		liste_residu.append(liste1)
		liste1=[]
		num_residu=num_residu+1
	return dico_atome

#dico_residu[line[22:26].strip()]=dico_atome[line[14:16].strip()]=dico_coord


test=parse_pdb("/home/sdv/m2bi/dde_murat/projet_court/1bta.pdb")


#with open("/home/sdv/m2bi/dde_murat/projet_court/1bta.pdb", "r") as pdb_file:
 #   res_count = 0
  #  for line in pdb_file:
   #     if line.startswith("ATOM"):
    #        atom_name = line[12:16].strip()
     #       if atom_name in representative_atoms:
      #          res_count += 1
       #         res_name = line[17:20].strip()
        #        res_num = int(line[22:26])
         #       x = float(line[30:38])
          #      y = float(line[38:46])
           #     z = float(line[46:54])
            #    el = (line[76:78])
                #print(res_name, res_num, atom_name, x, y, z, el)
             #   test=parse_pdb("/home/sdv/m2bi/dde_murat/projet_court/1bta.pdb")
                #print(test)





            
