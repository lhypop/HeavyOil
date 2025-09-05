import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager

class MyPlot:
    def __init__(self, nrows=1, ncols=1, figsize=(8, 5), dpi=100, sharex=False, sharey=False):
        self.fig, self.axes = plt.subplots(nrows=nrows, ncols=ncols,
                                           figsize=figsize, dpi=dpi,
                                           sharex=sharex, sharey=sharey)
        if isinstance(self.axes, np.ndarray):
            self.axes = self.axes.flatten()
        else:
            self.axes = np.array([self.axes])
        self.ax_index = 0
        self.colors = ["#71c9ce","#ffb6b9","#ff165d","#1f39ab"]

        font_path = "/usr/share/fonts/Arial/arial.ttf"
        if os.path.exists(font_path):
            self.prop = font_manager.FontProperties(fname=font_path, size=8)
        else:
            self.prop = None  # 使用默认字体

    def set_axis(self, ax=0, xlabel=None, ylabel=None, title=None,
                 xlabel_size=12, ylabel_size=12, title_size=14,
                 xlim=None, ylim=None, xticks=None, yticks=None,
                 grid=False, 
                 legend=True,legend_loc='best',flegend_size=10):
        """
        统一设置坐标轴信息

        参数：
        ax : int
            指定第几个子图（axes）进行设置。
        xlabel : str or None
            x轴标签文字。
        ylabel : str or None
            y轴标签文字。
        title : str or None
            图表标题文字。
        xlabel_size : int
            x轴标签字体大小，默认12。
        ylabel_size : int
            y轴标签字体大小，默认12。
        title_size : int
            标题字体大小，默认14。
        xlim : tuple or list or None
            x轴显示范围，例如 (xmin, xmax)。
        ylim : tuple or list or None
            y轴显示范围，例如 (ymin, ymax)。
        xticks : list or np.ndarray or None
            x轴刻度位置。
        yticks : list or np.ndarray or None
            y轴刻度位置。
        grid : bool
            是否显示网格，True 显示，False 不显示。
        legend : bool
            是否显示图例，True 显示，False 不显示。
        legend_loc : str or int
            图例位置，可选值：
            字符串：
                'best', 'upper right', 'upper left', 'lower left', 'lower right',
                'right', 'center left', 'center right', 'lower center',
                'upper center', 'center'
            整数编码：
                0='best', 1='upper right', 2='upper left', 3='lower left', 4='lower right',
                5='right', 6='center left', 7='center right', 8='lower center',
                9='upper center', 10='center'
        """
        ax = self.axes[ax]

        if xlabel is not None:
            ax.set_xlabel(xlabel, fontsize=xlabel_size, fontproperties=self.prop)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=ylabel_size, fontproperties=self.prop)
        if title is not None:
            ax.set_title(title, fontsize=title_size, fontproperties=self.prop)
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        if xticks is not None:
            ax.set_xticks(xticks)
        if yticks is not None:
            ax.set_yticks(yticks)

        ax.tick_params(
            axis="both",      # 同时作用于 x 和 y
            direction="in",   # 刻度朝内
            width=1,          # 刻度线宽度
            length=3,         # 刻度线长度
            top=False,         # 显示上边刻度
            right=False        # 显示右边刻度
        )

        ax.grid(grid)

        if legend:
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(loc = legend_loc, fontsize= flegend_size)

    def bar(self, x_array, y_array, 
            yerr=None,
            label=None,bar_width=0.8,
            color='skyblue', 
            alpha=0.7, 
            capsize=5, error_kw={'elinewidth':2, 'ecolor':'black'},
            ax=0):

        ax = self.axes[ax]
        x = np.asarray(x_array).ravel()
        y = np.asarray(y_array).ravel()

        # 如果传入 yerr，也展开成一维
        if yerr is not None:
            yerr = np.asarray(yerr).ravel()

        if label:
            ax.bar(x, y,
                yerr=yerr,  # 添加误差棒
                label=label,
                color=color, alpha=alpha,
                width=bar_width,
                capsize=capsize, error_kw=error_kw)  # capsize 控制误差棒横线长度
        else:
            ax.bar(x, y,
                yerr=yerr,
                color=color, alpha=alpha,
                width=bar_width,
                capsize=capsize, error_kw=error_kw)


    def step(self, x_array, y_array,
             label = None,
             color='blue', alpha=1.0, lw=2, where="mid", ax=0):

        ax = self.axes[ax]
        x = np.asarray(x_array).ravel()
        y = np.asarray(y_array).ravel()

        if label:
            ax.step(x, y, 
                    label = label,
                    color=color, alpha=alpha, lw=lw, where=where)
        else:
            ax.step(x, y, 
                    color=color, alpha=alpha, lw=lw, where=where)
            
    def scatter(self, x_array, y_array,
            label=None,
            color='blue', alpha=1.0, s=20, ax=0, marker='o'):

        ax = self.axes[ax]
        x = np.asarray(x_array).ravel()
        y = np.asarray(y_array).ravel()

        if label:
            ax.scatter(x, y,
                    label=label,
                    color=color, alpha=alpha, s=s, marker=marker)
        else:
            ax.scatter(x, y,
                    color=color, alpha=alpha, s=s, marker=marker)
            
    def line(self, x_array, y_array,
         label=None,
         color='blue', alpha=1.0, lw=2, ax=0, linestyle='-'):

        ax = self.axes[ax]
        x = np.asarray(x_array).ravel()
        y = np.asarray(y_array).ravel()

        if label:
            ax.plot(x, y, label=label, color=color, alpha=alpha, lw=lw, linestyle=linestyle)
        else:
            ax.plot(x, y, color=color, alpha=alpha, lw=lw, linestyle=linestyle)

    def probility_densit(self, x, label=None,
                            color='blue', fill=True, bw_adjust=0.5, ax=0):
        
        ax = self.axes[ax]
        sns.kdeplot(x.ravel(), fill=fill, bw_adjust=bw_adjust,
                    ax=ax, color=color, label=label)
        

    def show(self, show_if=False, save_path=None):
        self.fig.tight_layout()
        if save_path:
            self.fig.savefig(save_path, dpi=self.fig.dpi)
        if show_if:
            plt.show()

