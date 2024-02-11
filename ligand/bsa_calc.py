from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl


class Bsa():

    def magnitude(self, vector):
        return np.linalg.norm(vector)

    def creating_vector(self, vector,direction,length):
        """
        Функция, создающая векторы
        переменные: vector - вектор, direction - направление ,length - длина шага
        """

        xc = (vector[0] / direction) * length
        yc = (vector[1] / direction) * length
        zc = (vector[2] / direction) * length
        new_vector = np.array([xc, yc, zc])

        return new_vector

    def id_separator(self, structure):
        """
        Функция, разделяющая атомы лиганда и молекулы
        переменные: structure_name - путь к файлу структуры
        """
        atom_count1 =[]
        atom_count2 =[]
        lig_id =[]
        rec_id =[]
        chain_ids =[]

        for model in structure:
            for chain in model:
                for residue in chain:
                    chain_ids.append(residue.get_id())

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id() == chain_ids[1]:
                        for atom in residue:
                            atom_count1.append(atom.name)
                    else:
                        for atom in residue:
                            atom_count2.append(atom.name)

        if len(atom_count1) < len(atom_count2):
            lig_id = chain_ids[1]
            rec_id = chain_ids[0]
        else:
            lig_id = chain_ids[0]
            rec_id = chain_ids[1]
        return lig_id, rec_id

    def check_vectors(self, lig_vector, structure_name, steps, coeff, threshold, chain_ids):
        """
        Функция, проверяющая сообщённые вектора на вшивость
        переменные: lig_vector - вектор для проверки, structure_name - путь к файлу структуры, steps - расстояние (А)
        на которое будет смещаться структура лиганда при проверке, coeff - -1 или 1 показывает направление по которому будет
        смещаться структура лиганда, threshold - расстояние в А, которое определяет отсутсвие ковалентной связи между атомами, 
        нужно для проверки пересечения атомов возвращает массив с проверкой
        """
        parser = PDB.PDBParser()
        check_structure=[]
        coords=[]
        rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))
        for i in range(steps):
            structure = parser.get_structure(str(structure_name), structure_name)
            b = float(1) * (1 + i)
            a = int(coeff) * np.linalg.norm(lig_vector)
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id() == chain_ids[0]:
                            for atom in residue:
                                atom.transform(rotation_matrix, self.creating_vector(lig_vector, a, b))
            coords = self.coord_lig_atoms(structure, self.id_separator(structure)[0])

            for i in range(len(coords[0])): 
                for j in range(len(coords[1])):
                    xvector = coords[0][i] - coords[1][j]
                    if self.magnitude(xvector) > threshold:
                        check_structure.append(1)
                    else:
                        check_structure.append(0)

            coords[0].clear()
            coords[1].clear()

        return check_structure

    def coord_lig_atoms(self, structure, lig_id):
        """
        Функция, находящая все координаты атомов для лиганда и рецептора
        переменные: structure - структура,  lig_id - id цепочки
        """
        lig_coor =[]
        rec_coor =[]

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id() == lig_id:
                        for atom in residue:
                            lig_coor.append(atom.get_coord())
                    else:
                        for atom in residue:
                            rec_coor.append(atom.get_coord())

        return lig_coor, rec_coor

    def ligand_vectors(self, coeff, chain_com):
        """
        Функция, создающая все векторы для лиганда
        переменные:  coeff - -1 или 1 показывает направление по которому будет
        смещаться структура лиганда, chain_com - координаты атомов лиганда 
        """
        lig_vectors=[]
        for i in range(len(chain_com)):
            for j in range(len(chain_com)):
                if i != j:
                    vector = chain_com[i] - chain_com[j]
                    b = float(1) * (1)
                    a = int(coeff) * np.linalg.norm(vector)
                    lig_vectors.append(self.creating_vector(vector, a, b))
        return lig_vectors

    def creating_model(self, ranges, structure_name, coeff, good_vectors, chain_ids, f=0, name_save='Cyclodextrine_ligand_names'):
        """
        Функция, прорисовывающая модель
        переменные: ranges - кол-во повторов итераций, structure_name - путь к файлу структуры, coeff - -1 или 1 показывает направление по которому будет
        смещаться структура лиганд, good_vectors - хорошие векторы, прошедшие проверку,  chain_ids - id молекул
        """
        parser = PDB.PDBParser()
        rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))
        io = PDB.PDBIO()
        for i in range(ranges):
            structure = parser.get_structure(str(structure_name), structure_name)
            b = float(1) * (i + 1)
            vector = good_vectors[f]
            a = int(coeff) * np.linalg.norm(vector)

            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id() == chain_ids[0]:
                            for atom in residue:
                                atom.transform(rotation_matrix, self.creating_vector(vector, a, b)) 

            io.set_structure(structure)


            io.save(f"{self.dir_calc_model}/{name_save}_{str(b*int(coeff))}.pdb")

    def make_good_vectors(self, vectors, structure_name, chains_ids, ncheck=5):

        """
        Функция, находящая все хорошие векторы
        переменные: vectors - полученные векторы, structure_name - путь к файлу структуры, chain_ids - id молекул
        """

        good_vectors=[]
        count = 0

        for vector in vectors:
            count = count + 1
            print(round(count/len(vectors)*100, 2), "%  ", sep=' ', end='', flush=True)
            check_result_pos = self.check_vectors(vector, structure_name, ncheck, 1, self.distance, chains_ids)
            check_result_neg = self.check_vectors(vector, structure_name, ncheck, -1, self.distance, chains_ids)
            if 0 not in check_result_pos and check_result_neg:
                good_vectors.append(vector)

        return good_vectors

    def check_good_vector(self, good_vectors, structure_name, chains_ids, nshift=5):

        """
        Функция, 
        переменные:
        """
        parser = PDB.PDBParser()
        avg_sasa = []
        all_max = []
        count=0
        for f in range(len(good_vectors)): #вместо 2 нужно len(good_vectors)
            count=count+1
            print(round(count/len(good_vectors)*100, 2), " %")
            a = []  # complex
            sasa_a = []
            self.creating_model(nshift, structure_name, 1, good_vectors, chains_ids, f = f, name_save = 'check_good_vectors_'+str(f))
            self.creating_model(nshift, structure_name, -1, good_vectors, chains_ids, f = f, name_save = 'check_good_vectors_'+str(f))
            for i in range(nshift):
                current_name_pos = f'{self.dir_calc_model}/check_good_vectors_'+str(f)+'_' + str(i+1)+'.0' + '.pdb'
                current_name_neg = f'{self.dir_calc_model}/check_good_vectors_'+str(f)+'_' + str(-(i+1))+'.0' + '.pdb'
                structure_pos = parser.get_structure(current_name_pos, current_name_pos)
                structure_neg = parser.get_structure(current_name_neg, current_name_neg)
                a.append(structure_pos)
                a.append(structure_neg)
                os.remove(os.path.join(os.getcwd(),current_name_pos))
                os.remove(os.path.join(os.getcwd(),current_name_neg))
            for struct in a:
                sr = ShrakeRupley()
                sr.compute(struct, level="S")
                sasa_a.append(round(struct.sasa, 2))
            max_sasa=sasa_a[0]
            for i in sasa_a:
                if i>max_sasa:
                    max_sasa=i
            average_sasa=sum(sasa_a)/len(sasa_a)

            all_max.append(max_sasa)
            avg_sasa.append(average_sasa)
            #print(sasa_a)
            #print(average_sasa)
            #print(max_sasa)
        #print(all_max)
        #print(avg_sasa)
        good_index=self.check_sasa_massive(all_max, avg_sasa)
        print(good_index)
        return good_index


    def reject_outliers(self, data, m = 2.):
        data=np.array(data)
        d = np.abs(data - np.median(data))
        mdev = np.median(d)
        s = d/mdev if mdev else np.zeros(len(d))
        return data[s<m]

    def check_sasa_massive(self, max, avg):
        if self.reject_outliers(max).all() == np.array(max).all():
            min_maxv=min(max)
            good_index=max.index(min_maxv)
        else:
            min_avgv=min(avg)
            good_index=avg.index(min_avgv)
        return good_index


    def getting_final_sasa(self, nshift, structure_name, good_vectors, chains_ids, f):
        parser = PDB.PDBParser()
        sasa_final = []
        a = []  # complex

        self.creating_model(nshift, structure_name, -1, good_vectors, chains_ids, f , name_save = 'check_good_vectors_'+str(f))
        self.creating_model(nshift, structure_name, 1, good_vectors, chains_ids, f , name_save = 'check_good_vectors_'+str(f))
        for i in range(nshift):

            
            current_name_neg = f'{self.dir_calc_model}/check_good_vectors_'+str(f)+'_' + str(-(i+1))+'.0' + '.pdb'
            current_name_pos = f'{self.dir_calc_model}/check_good_vectors_'+str(f)+'_' + str(i+1)+'.0' + '.pdb'
            structure_pos = parser.get_structure(current_name_pos, current_name_pos)
            structure_neg = parser.get_structure(current_name_neg, current_name_neg)
            a.append(structure_pos)
            a.append(structure_neg)
            os.remove(os.path.join(os.getcwd(),current_name_pos))
            os.remove(os.path.join(os.getcwd(),current_name_neg))
        for struct in a:
            sr = ShrakeRupley()
            sr.compute(struct, level="S")
            sasa_final.append(round(struct.sasa, 2))
            # print(sasa_final)
        return sasa_final

    def making_graf(self, sasa):
        x = list(range(int((-len(sasa))/2),int((len(sasa))/2)+1))
        x.remove(int(0))
        sasa_pos=sasa[::2]
        sasa_neg=list(reversed(sasa[1::2]))
        sasa_neg.extend(sasa_pos)
        y=sasa_neg
        plt.plot(x,y)
        plt.savefig(f'{self.dir_calc_model}/graph.png')

    def main(self, modal_num, distance, dir_path):

        self.distance = distance
        self.dir_calc_model = dir_path
        
        # np.int = np.int_ can fix problem with numpy
        parser = PDB.PDBParser()
        # Путь к файлу 
        structure_name = ".\Cyclodextrine_ligand_names_good.pdb"
        structure = parser.get_structure(structure_name, structure_name)
        chains_com = self.coord_lig_atoms(structure, self.id_separator(structure)[0])
        chains_ids = self.id_separator(structure)
        print(len(chains_com[0]))
        print(len(chains_com[1]))
        print(chains_ids)
        vectors = self.ligand_vectors(1, chains_com[0])
        print(len(vectors))
        good_vectors=self.make_good_vectors(vectors,structure_name, chains_ids, ncheck=5)
        print(len(good_vectors))
        best_index=self.check_good_vector(good_vectors, structure_name, chains_ids)
        self.creating_model(modal_num, structure_name, 1, good_vectors, chains_ids, f = best_index)
        self.creating_model(modal_num, structure_name, -1, good_vectors, chains_ids, f = best_index)
        final_sasa=self.getting_final_sasa(modal_num, structure_name, good_vectors, chains_ids, f = best_index)
        print(final_sasa)
        self.making_graf(final_sasa)

    
