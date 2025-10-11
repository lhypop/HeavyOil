import pandas as pd
import numpy as np
import MDAnalysis as mda
from tqdm import tqdm

from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


class Interactanalysis():
    def __init__(self, traject, frame_clusters):

        self.traject = traject
        self.frame_clusters = frame_clusters
        self.ONS = self.traject.select_atoms("element O N S")

        self.frames = [frame for frame, _ in frame_clusters.items()]
        self.Hete_C, self.Hete_pairs = self.contacts_within_cutoff()

    def contacts_within_cutoff(self, radius = 5, min_radius = 2.5):

        timeseries = []
        pairs_per_frame = []
        for frame in tqdm(self.frames, desc="cluster heteroatoms analysis"):
            # calculate distances between group_a and group_b
            dist = contacts.distance_array(
                                            self.ONS.positions, self.ONS.positions, 
                                            box=self.traject.trajectory[frame].dimensions, 
                                            backend="openmp")
            # determine which distances <= radius
            triu = np.triu_indices_from(dist, k=1)
            # mask = dist[triu] <= radius
            mask = (dist[triu] >= min_radius) & (dist[triu] <= radius)
            
            #统计接触数
            n_contacts= mask.sum()
            timeseries.append([frame, n_contacts])

            # 保存原子对
            pairs = [(self.ONS[i].index, self.ONS[j].index) for i, j, m in zip(triu[0], triu[1], mask) if m]
            pairs_per_frame.append([frame, pairs])

        return np.array(timeseries), pairs_per_frame
    
    def Heter_statistic_analysis(self):
        
        frame_Heter = {
            "A_heter": [],
            "I_heter": [],
            "O_heter": [],
            "S_heter":[]
        }
        
        for frame in tqdm(self.frames, desc="Heteroatoms analysis"):

            Heter_pair_frame = [item[1] for item in self.Hete_pairs if item[0] == frame][0]

            # --- 1. 整帧杂原子（簇内+簇外） ---
            IOut_resids = [res for cluster in self.frame_clusters[frame] for res in cluster]
            IO_atoms = [atom.index for res in IOut_resids for atom in res.atoms]
            IO_heter = [ pair for pair in Heter_pair_frame
                if pair[0] in IO_atoms and pair[1] in IO_atoms 
            ]
            # print(IO_heter) 
            # print(len(IO_heter))

            # --- 2. 簇内杂原子 ---
            I_heter = []
            for res_set in self.frame_clusters[frame]: 
                res_set = [res for res in res_set] 
                In_atoms = [atom.index for res in res_set for atom in res.atoms]
                Heter_I = [ pair for pair in Heter_pair_frame
                            if pair[0] in In_atoms and pair[1] in In_atoms 
                        ]
                if Heter_I:
                    I_heter.extend(Heter_I)
            # print(I_heter)
            # print(len(I_heter))

            # --- 3. 簇间杂原子 ---
            O_heter = [pair for pair in IO_heter if pair not in I_heter]
            # print(O_heter)
            # print(len(O_heter))
            
            # --- 4. 簇表面杂原子
            sur_heter = []
            sur_set = set()
            for res_set in self.frame_clusters[frame]:
                res_set = [res for res in res_set]
                sur_atoms = [atom.index for res in res_set for atom in res.atoms]

                for pair in Heter_pair_frame:
                    donor, acceptor = pair

                    if (donor in sur_atoms) != (acceptor in sur_atoms) \
                            and (donor, acceptor) not in sur_set \
                            and (acceptor, donor) not in sur_set:
                        
                        sur_heter.append(pair)
                        sur_set.add((donor,acceptor))
                
            # print(sur_heter)
            # print(len(sur_heter))

            frame_Heter["A_heter"].append(len(IO_heter))
            frame_Heter["I_heter"].append(len(I_heter))
            frame_Heter["O_heter"].append(len(O_heter))
            frame_Heter["S_heter"].append(len(sur_heter))

        return pd.DataFrame(frame_Heter)
    
    def Hydro_statistic_analysis(self, select_info):

        frame_Hydro = {
            "A_hydro": [],
            "I_hydro": [],
            "O_hydro": [],
            "S_hydro":[]
        }


        T_hbonds = HydrogenBondAnalysis(universe=self.traject)
        T_hbonds.hydrogens_sel = T_hbonds.guess_hydrogens(select_info)
        T_hbonds.acceptors_sel = T_hbonds.guess_acceptors(select_info)
        T_hbonds.run(verbose=True, start= self.frames[0], stop = self.frames[-1])

        for frame in tqdm(self.frames, desc="cluster hydrogen bonds analysis"):
            
            T_hbonds_array = T_hbonds.results.hbonds
            T_hbonds_frame = T_hbonds_array[T_hbonds_array[:,0].astype(int) == frame]
            

            # --- 1. 整帧氢键（簇内+簇外） ---
            IOut_resids = [res for cluster in self.frame_clusters[frame] for res in cluster]

            IO_atoms = [atom.index for res in IOut_resids for atom in res.atoms]
            IO_hbonds = [hbond for hbond in T_hbonds_frame
                          if int(hbond[1]) in IO_atoms and int(hbond[3]) in IO_atoms]
            # print(IO_hbonds)

            # --- 2. 簇内氢键 ---
            I_hbonds = [] 
            for res_set in self.frame_clusters[frame]: 
                res_set = [res for res in res_set] 

                In_atoms = [atom.index for res in res_set for atom in res.atoms]
                hbonds_In = [hbond for hbond in T_hbonds_frame
                                if int(hbond[1]) in In_atoms and int(hbond[3]) in In_atoms]
                if hbonds_In:
                    I_hbonds.extend(hbonds_In)

            # print(I_hbonds)

            # --- 3. 簇间氢键 ---
            I_hbonds_set = {(int(h[1]), int(h[3])) for h in I_hbonds}
            O_hbonds = [h for h in IO_hbonds if (int(h[1]), int(h[3])) not in I_hbonds_set]
            # print(O_hbonds)

            # --- 4. 簇表面氢键
            sur_hbonds = []
            sur_set = set()
            for res_set in self.frame_clusters[frame]:
                res_set = [res for res in res_set]

                sur_atoms = [atom.index for res in res_set for atom in res.atoms]
                for hbond in T_hbonds_frame:
                    donor, acceptor = int(hbond[1]), int(hbond[3])

                    if (donor in sur_atoms) != (acceptor in sur_atoms) \
                          and (donor, acceptor) not in sur_set \
                          and (acceptor, donor) not in sur_set:
                        
                        sur_hbonds.append(hbond)
                        sur_set.add((donor,acceptor))
    

            frame_Hydro["A_hydro"].append(len(IO_hbonds))
            frame_Hydro["I_hydro"].append(len(I_hbonds))
            frame_Hydro["O_hydro"].append(len(O_hbonds))
            frame_Hydro["S_hydro"].append(len(sur_hbonds))

        return pd.DataFrame(frame_Hydro)