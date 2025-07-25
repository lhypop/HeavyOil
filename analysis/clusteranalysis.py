import numpy as np
import pandas as pd

class ClusterAnalysis:
    """Analyzes molecular cluster properties and statistics."""
    
    def __init__(self):
        pass
        
    def get_cluster_statistics(self,merged_clusters: list, tolal_number: int) -> pd.DataFrame:
        """
        Calculate descriptive statistics for cluster sizes.
        
        Returns:
            DataFrame containing:
            - count: Number of clusters
            - mean: Average cluster size
            - std: Standard deviation
            - min/max: Size range
            - 25/50/75%: Quartiles
        """
        # Get actual cluster sizes
        cluster_sizes = [len(cluster) for cluster in merged_clusters]
        
        for i in range(tolal_number - np.sum(cluster_sizes)):
            cluster_sizes.append(1)
        
        # Generate statistics
        size_series = pd.Series(cluster_sizes, name='cluster_size')
        stats = size_series.describe()
        
        # Convert to DataFrame for better formatting
        return pd.DataFrame(stats).transpose()

    def save_outcome(self, summarys, time, data_file):

        lines_list = []
        header = f"{'times':<10} {'count':<10} {'mean':<10} {'std':<10} {'min':<10} {'25%':<10} {'50%':<10} {'75%':<10} {'max':<10}\n"
        lines_list.append(header)
        for j in range(len(summarys)):
            line = f"{time[j]:<10} {summarys[j]['count']:<10} {summarys[j]['mean']:<10.2f} {summarys[j]['std']:<10.2f} " \
                   f"{summarys[j]['min']:<10} {summarys[j]['25%']:<10} {summarys[j]['50%']:<10} {summarys[j]['75%']:<10} " \
                   f"{summarys[j]['max']:<10}\n"
            lines_list.append(line)

        if not data_file.exists():
            data_file.parent.mkdir(exist_ok=True,parents=True)

        with open(data_file, "w") as f:
            f.writelines(lines_list)