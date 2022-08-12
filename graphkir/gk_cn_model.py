"""
Module for clustering depths into CN
* CNgroup (proposed)
* KDE
"""
from __future__ import annotations
import json
import numpy as np
from typing import Any
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema


class Dist:
    def __init__(self):
        self.base = 0

    def fit(self, values: list[float]):
        """ Determine the parameters by data """
        raise NotImplementedError

    def assignCN(self, values: list[float]) -> list[int]:
        """ Assign CN by parameters """
        raise NotImplementedError

    def save(self, filename: str):
        """ Save parameters """
        json.dump(self.getParams(), open(filename, "w"))

    def load(self, filename: str) -> Dist:
        """ Load parameters """
        data = json.load(open(filename))
        return self.setParams(data)

    def getParams(self) -> dict:
        """ Save parameters """
        raise NotImplementedError

    @classmethod
    def setParams(cls, data: dict) -> Dist:
        """ Load parameters """
        raise NotImplementedError


class CNgroup(Dist):
    """
    Our method: CN_group
    """
    def __init__(self):  # , assume_base=None):
        # const
        self.bin_num = 500    # discrete the values
        self.max_cn = 6       # CN limits

        # parameters
        self.base = None      # u1
        self.x_max = 1        # for normalize
        self.base_dev = 0.08  # SD1
        self.y0_dev = 1.5     # SD0

        # KIR assumption
        # self.assume_base = assume_base
        # self.assume_3DL3_diploid = True

        # result (saved for plotting)
        self.values = []
        self.likelihood = []

    def getParams(self) -> dict:
        """ Save parameters """
        return {
            'method': "CNgroup",
            'base': self.base,
            'base_dev': self.base_dev,
            'x_max': self.x_max,
            'y0_dev': self.y0_dev,
            'bin_num': self.bin_num,
            'max_cn': self.max_cn,
        }

    @classmethod
    def setParams(cls, data: dict) -> CNgroup:
        """ Load parameters """
        assert data['method'] == "CNgroup"
        self = cls()
        self.base = data['base']
        self.base_dev = data['base_dev']
        self.x_max = data['x_max']
        self.max_cn = data['max_cn']
        self.y0_dev = data['y0_dev']
        self.bin_num = data['bin_num']
        return self

    def fit(self, values: list[float]):
        """
        Fit the data, it will set

        * x_max, base_dev
        * base
        """
        assert self.base is None
        # normalize
        max_depth = max(values) * 1.2
        self.base_dev *= max_depth
        self.x_max = max_depth
        self.values = values

        # discrete (bin_num)
        # Calculate the probility that CN groups fit the data
        density, _ = np.histogram(values, bins=self.bin_num, range=(0, self.x_max))
        likelihood = []
        for base in np.linspace(0, self.x_max, self.bin_num):
            # all probility of cn group across 0 ~ x_max
            cn_group = self.calcCNGroupProb(base)
            # Highest probility in each x
            max_prob = cn_group.max(axis=0)
            # log-probility = depth \cdot the log(highest probility)
            likelihood.append((base, np.sum(np.log(max_prob + 1e-9) * density)))
        self.likelihood = np.array(likelihood)  # n x 2(base, likelihood of the base)

        # Find best fit x = base
        max_point = self.likelihood[np.argmax(self.likelihood[:, 1]), :]
        self.base = max_point[0]

        # special case
        # if self.assume_base and base / self.assume_base > 1.7:  # TODO: maybe more CN?
        #     base /= 2
        #     self.base = base

    def assignCN(self, values: list[float]) -> list[int]:
        """ Assign CN group for each depths """
        assert self.base is not None
        cn_group = self.calcCNGroupProb(self.base)
        cn_max = cn_group.argmax(axis=0)
        space = self.x_max / self.bin_num
        return [cn_max[int(depth / space)] for depth in values]

    def predictCN(self, read_depth: dict[str, float]) -> dict[str, int]:
        """ fit + assign """
        assert self.base is None
        values = list(read_depth.values())
        self.fit(values)
        cn = self.assignCN(values)
        return dict(zip(
            read_depth.keys(),
            cn
        ))

    def calcCNGroupProb(self, base: float):
        """
        Returns:
            ( CN x norm_read_depth(500bin) ) array

            indicate the probility of read_depth belong to the CN
        """
        x = np.linspace(0, self.x_max, self.bin_num)
        cn = np.arange(1, self.max_cn)
        y0 = norm.pdf(x, loc=0, scale=self.base_dev * self.y0_dev)
        y = np.stack([y0, *[norm.pdf(x, loc=base*n, scale=self.base_dev*n) for n in cn]])
        space = self.x_max / self.bin_num  # * space is to make y-sum = 1
        return y * space
    # def predictKIRCN(self, read_depth: dict[str, float]) -> dict[str, int]:
    #     """ fit + assign + 3DL3 calibrated """
    #     predict_cn = self.predictCN(read_depth)
    #     if self.assume_3DL3_diploid and predict_cn["KIR3DL3*BACKBONE"] != 2:
    #         print("WARNING  KIR3DL3 is not diploid, trying to fix it")
    #         assert predict_cn["KIR3DL3*BACKBONE"] == 1
    #         # we only deal with this case
    #         self.base /= 2
    #         predict_cn = dict(zip(
    #             read_depth.keys(),
    #             self.assignCN(list(read_depth.values()))
    #         ))
    #         assert predict_cn["KIR3DL3*BACKBONE"] == 2
    #     return predict_cn

    def plot(self):
        assert self.base is not None

        # loss function of distribution
        fig_loss = px.line(x=self.likelihood[:, 0], y=self.likelihood[:, 1])
        fig_loss.add_vline(x=self.base, line_dash="dash", line_color="gray",
                           annotation_text=f"max value={self.base}")
        fig_loss.update_layout(
            xaxis_title="read-depth",
            yaxis_title="likelihood",
        )

        # histogram of read depth
        fig_dist = make_subplots(specs=[[{"secondary_y": True}]])
        fig_dist.add_trace(go.Histogram(x=list(self.values),
                                        name="Observed",
                                        nbinsx=self.bin_num),
                           secondary_y=True)
        fig_dist.update_layout(
            xaxis_title="read-depth",
            yaxis_showticklabels=False,
            yaxis_title="probility",
        )
        fig_dist.update_yaxes(title_text="Number of read depth", secondary_y=True)

        # distribution
        x = np.linspace(0, self.x_max, self.bin_num)
        y = self.calcCNGroupProb(self.base)
        for n in range(len(y)):
            fig_dist.add_trace(go.Scatter(x=x, y=y[n], name=f"cn={n}"))

        prev_y = 0
        for i, yi in enumerate(np.argmax(y, axis=0)):
            if yi != prev_y:
                fig_dist.add_vline(x=x[i], annotation_text=f"CN {yi-1}-{yi}",
                                   line_dash="dash", line_color="gray")
            prev_y = yi

        return [fig_loss, fig_dist]


class KDEcut(Dist):
    """ CN estimation via KDE """

    def __init__(self):
        self.bandwidth = 0.05  # KDE parameters
        self.points = 100      # discrete
        self.neighbor = 5      # thresholding 5 / 100 = 0.05

        # data
        self.kde = None        # KDE object
        self.x_max = None      # normalize
        self.local_min = []    # thresholding

        # for plot
        self.data = []         # normalized data

    def getParams(self) -> dict:
        """ Save parameters """
        return {
            'method': "KDEcut",
            'x_max': self.x_max,
            'local_min': self.local_min,
            'kde': self.kde.get_params(),
        }

    @classmethod
    def setParams(cls, data: dict) -> KDEcut:
        """ Load parameters """
        assert data['method'] == "KDEcut"
        self = cls()
        self.x_max = data['x_max']
        self.local_min = data['local_min']
        self.kde = KernelDensity().set_params(**data['kde'])
        return self

    def fit(self, values: list[float]):
        """ Fit the data to KDE and save the threshold for cutting local minimal """
        assert self.kde is None
        # normalize to 0 - 1
        self.x_max = np.max(values)
        data = np.array(values)[:, None] / self.x_max
        self.kde = KernelDensity(kernel='gaussian', bandwidth=self.bandwidth).fit(data)

        # cut
        x = np.linspace(0, 1.1, self.points)
        y = self.kde.score_samples(x[:, None])
        self.local_min = x[argrelextrema(y, np.less, order=self.neighbor)[0]]
        print("Threshold", self.local_min)

        # for plot
        self.data = values

    def assignCN(self, values: list[float]) -> list[int]:
        """ Assign CN group for each depths """
        assert self.kde is not None
        x = np.array(values) / self.x_max
        cn = np.searchsorted(self.local_min, x)
        return cn

    def predictCN(self, read_depth: dict[str, float]) -> dict[str, int]:
        """ fit + assign """
        values = list(read_depth.values())
        self.fit(values)
        cn = self.assignCN(values)
        return dict(zip(
            read_depth.keys(),
            cn
        ))

    def plot(self):
        x = np.linspace(0, 1.1, self.points)
        y = self.kde.score_samples(x[:, None])
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(go.Scatter(x=x, y=y, name="kde"))
        for cn, m in enumerate(self.local_min):
            fig.add_vline(x=m, line_width=2, line_dash="dash", annotation_text=f"cn={cn}")

        fig.add_trace(go.Histogram(
            x=np.array(self.data) / self.x_max,
            name="Relative Depth", nbinsx=100, histnorm="probability"), secondary_y=True)
        return [fig]
