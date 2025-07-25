from collections import defaultdict

class postprocess:    
    def __init__(self,cluster_labels):
        self.labels = cluster_labels
        self.unique_res_clusters = defaultdict(set)


    def _postprocss(self, select):


        self._group_by_labels()
        self.tolnum,self.filtered_clusters = self._filter_single(select)
        self.merge_clusters = self._merge_clusters()
        self.valiate_cluster()

        return self.tolnum,self.merge_clusters

    def _group_by_labels(self):

        # Group residues by cluster ID
        for cluster_id, res in self.labels:
            self.unique_res_clusters[cluster_id].add(res)

    def _filter_single(self,select):

        if select == "part" :
            total_numbers = len([cluster for cluster in 
                    (set(residues) for residues in self.unique_res_clusters.values())])
        else:
            total_numbers = len([cluster for cluster in 
                (set(residues) for residues in self.unique_res_clusters.values()) 
                if len(cluster) > 1])
        filtered_clusters = [cluster for cluster in 
                (set(residues) for residues in self.unique_res_clusters.values()) 
                if len(cluster) > 1]
        
        return total_numbers,filtered_clusters
    
    def _merge_clusters(self):
        """
        Merge clusters that share common residues.
        
        Parameters:
            unique_res_clusters (dict): A dictionary where keys are cluster IDs and values are sets of residues.
            
        Returns:
            list of sets: A list of merged clusters, each represented as a set of residues.
        """
        # Initialize cluster sets with the input dictionary values
        cluster_sets = self.filtered_clusters
        
        changed = True
        while changed:
            changed = False
            merged_clusters = []
            
            # Process all clusters to merge overlapping ones
            while cluster_sets:
                current = cluster_sets.pop()
                merged = False
                
                # Check if the current cluster overlaps with any merged cluster
                for i, merged_set in enumerate(merged_clusters):
                    if current & merged_set:  # If there is any shared residue
                        merged_clusters[i].update(current)  # Merge clusters
                        merged = True
                        changed = True
                        break
                
                # If no overlap found, add the current cluster as a new merged cluster
                if not merged:
                    merged_clusters.append(current)
            
            # Update cluster sets for the next iteration
            cluster_sets = merged_clusters
        
        return cluster_sets
    
    def valiate_cluster(self):

        if sum(len(s) for s in self.merge_clusters) != len(set().union(*self.merge_clusters)):
            print("There are repetitive residue molecules")
            raise

