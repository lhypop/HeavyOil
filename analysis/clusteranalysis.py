import numpy as np
import pandas as pd
from collections import Counter

class ClusterAnalysis:
    """Analyzes molecular cluster properties and statistics."""
    
    def __init__(self,datafile):
        self.datafile = datafile
        
    def get_cluster_statistics(self, summarys: list, cluster_time) -> pd.DataFrame:
        """
        Calculate descriptive statistics for cluster sizes per frame.
        
        Args:
            cluster_time_array: list of lists or 2D array, cluster sizes per frame
            total_number: total number of particles/nodes
            cluster_time: list/array of frame times
        
        Returns:
            DataFrame with descriptive statistics per frame
        """
        statis_cluster = []
        
        for cluster_sizes, time in zip(summarys, cluster_time):
            sizes = list(cluster_sizes)

            # print(f"Total number of molecular in clusters is : {np.sum(sizes)}")
            size_series = pd.Series(sizes, name='cluster_size')
            stats = size_series.describe().to_dict()
            stats['time'] = time
            statis_cluster.append(stats)
        
        return pd.DataFrame(statis_cluster)
    
    def cal_Csize_probability(self, summarys: list):

        flat_list = [x for row in summarys for x in row]
        counts = Counter(flat_list)
        n_frame = len(summarys)

        values = np.array(list(counts.keys()), dtype=float)
        freqs = np.array(list(counts.values()), dtype=float)

        f_counts = freqs / n_frame
        p_counts = values * f_counts / sum(summarys[0])

        sort_idx = np.argsort(values)
        values_sorted = values[sort_idx]
        p_counts_sorted = p_counts[sort_idx]

        print("values_sorted:", values_sorted)
        print("p_counts_sorted:", p_counts_sorted)

    # def cal_radius_of_gyration():



    # def save_outcome(summarys, time, data_file):

    #     lines_list = []
    #     header = f"{'times':<10} {'count':<10} {'mean':<10} {'std':<10} {'min':<10} {'25%':<10} {'50%':<10} {'75%':<10} {'max':<10}\n"
    #     for j in range(len(summarys)):
    #         line = f"{time[j]:<10} {summarys[j]['count']:<10} {summarys[j]['mean']:<10.2f} {summarys[j]['std']:<10.2f} " \
    #                f"{summarys[j]['min']:<10} {summarys[j]['25%']:<10} {summarys[j]['50%']:<10} {summarys[j]['75%']:<10} " \
    #                f"{summarys[j]['max']:<10}\n"
    #         print(line)
    #         lines_list.append(line)

    #     with open(data_file, "w") as f:
    #         f.writelines(lines_list)

    # def plot_MaxSize_and_AverSize():
