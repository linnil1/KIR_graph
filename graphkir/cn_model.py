"""
Module for clustering depths into CN

* CNgroup (proposed)
* KDE
"""
from __future__ import annotations
from typing import Any
import json
import numpy as np

from scipy.stats import norm
from scipy.signal import argrelextrema
from sklearn.neighbors import KernelDensity
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.graph_objects as go

from .utils import NumpyEncoder


class Dist:
    """ Abstract class of CN prediction model """
    def __init__(self):
        pass

    def fit(self, values: list[float]) -> None:
        """ Determine the parameters by data """
        raise NotImplementedError

    def assignCN(self, values: list[float]) -> list[int]:
        """ Assign CN with depth input by parameters """
        raise NotImplementedError

    def save(self, filename: str) -> None:
        """ Save parameters """
        with open(filename, "w") as f:
            json.dump(self.getParams(), f, cls=NumpyEncoder)

    @classmethod
    def load(cls, filename: str) -> Dist:
        """ Load parameters """
        with open(filename) as f:
            data = json.load(f)
        return cls.setParams(data)

    def getParams(self) -> dict[str, Any]:
        """ Get parameters """
        raise NotImplementedError

    @classmethod
    def setParams(cls, data: dict[str, Any]) -> Dist:
        """ Set parameters """
        raise NotImplementedError

    def plot(self, title: str = "") -> list[go.Figure]:
        """ Plot the model """
        raise NotImplementedError


class CNgroup(Dist):
    """
    Our method: CN_group

    Attributes:
      x_max:    For normalization
      base:     u1 (mean of CN=1 distribution)
      base_dev: SD1 (devation of CN=1 distribution)
      y0_dev:   SD0 (devation of CN=0 distribution)
      bin_num:  Bins number in [0, x_max]
      max_cn:   CN limits
    """
    def __init__(self):
        # const
        super().__init__()
        self.bin_num:  int   = 500
        self.max_cn:   int   = 6

        # parameters
        self.x_max:    float = 1
        self.base:     float = None
        self.base_dev: float = 0.08
        self.y0_dev:   float = 1.5

        # result (saved for plotting)
        self.data            = []
        self.likelihood      = []

    def getParams(self) -> dict[str, Any]:
        """ Save parameters """
        return {
            'method'  : "CNgroup",
            'x_max'   : self.x_max,
            'base'    : self.base,
            'base_dev': self.base_dev,
            'y0_dev'  : self.y0_dev,
            'bin_num' : self.bin_num,
            'max_cn'  : self.max_cn,
            'data'    : self.data,
            'likelihood' : self.likelihood,
        }

    @classmethod
    def setParams(cls, data: dict[str, Any]) -> CNgroup:
        """ Load parameters """
        assert data['method'] == "CNgroup"
        self = cls()
        self.base     = data['base']
        self.base_dev = data['base_dev']
        self.x_max    = data['x_max']
        self.y0_dev   = data['y0_dev']
        self.bin_num  = data['bin_num']
        self.max_cn   = data['max_cn']
        self.data     = data['data']
        self.likelihood = np.array(data['likelihood'])
        return self

    def fit(self, values: list[float]) -> None:
        """
        Find the maximum CN distributions to fit the values.

        1. Normalize by `x_max` (Maximum of depths)
        2. Find `base` (Mean of CN=1 distribution)
        """
        assert self.base is None
        # normalize
        max_depth = max(values) * 1.2
        self.base_dev *= max_depth
        self.x_max = max(max_depth, 1e-6)  # to avoid divided by 0
        self.data = values

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

    def assignCN(self, values: list[float]) -> list[int]:
        """ Assign CN group for each depths """
        assert self.base is not None
        cn_group = self.calcCNGroupProb(self.base)
        cn_max = cn_group.argmax(axis=0)
        space = self.x_max / self.bin_num
        return [cn_max[int(depth / space)] for depth in values]

    def calcCNGroupProb(self, base: float) -> np.ndarray:
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

    def plot(self, title: str = "") -> list[go.Figure]:
        assert self.base is not None

        # loss function of distribution
        fig_loss = px.line(x=self.likelihood[:, 0], y=self.likelihood[:, 1])
        fig_loss.add_vline(x=self.base, line_dash="dash", line_color="gray",
                           annotation_text=f"max value={self.base}")
        fig_loss.update_layout(
            xaxis_title="read-depth",
            yaxis_title="likelihood",
        )
        fig_loss.update_layout(title=title)

        # histogram of read depth
        fig_dist = make_subplots(specs=[[{"secondary_y": True}]])
        fig_dist.add_trace(go.Histogram(x=list(self.data),
                                        name="Observed",
                                        nbinsx=self.bin_num),
                           secondary_y=True)
        fig_dist.update_layout(
            title=title,
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
    """
    CN estimation via KDE

    Inspired by Sakaue, Saori, et al.
    "Decoding the diversity of killer immunoglobulin-like receptors "
    by deep sequencing and a high-resolution imputation method."
    Cell Genomics 2.3 (2022): 100101.

    Attributes:
      bandwidth: KDE's bandwidth
      points:    The number of bins in [0, 1]
      neighbor:  The thresholding points to cut local minimal
      x_max:     All values are normalize by this value into [0, 1]
      kde:       KernelDensity object
      local_min: list of positions of local minimal i.e. CN thesholds
    """

    def __init__(self):
        super().__init__()
        self.bandwidth: float       = 0.05
        self.points:    int         = 100
        self.neighbor:  int         = 5

        # data
        self.x_max:     float       = 0
        self.kde                    = None
        self.local_min: list[float] = []

        # for plot
        self.data: list[float]      = []  # normalized data

    def getParams(self) -> dict[str, Any]:
        """ Save parameters """
        return {
            'method'   : "KDEcut",
            'bandwidth': self.bandwidth,
            'points'   : self.points,
            'neighbor' : self.neighbor,
            'x_max'    : self.x_max,
            'kde'      : self.kde.get_params(),
            'local_min': self.local_min,
        }

    @classmethod
    def setParams(cls, data: dict[str, Any]) -> KDEcut:
        """ Load parameters """
        assert data['method'] == "KDEcut"
        self = cls()
        self.x_max     = data['x_max']
        self.bandwidth = data['bandwidth']
        self.neighbor  = data['neighbor']
        self.points    = data['points']
        self.kde       = KernelDensity().set_params(**data['kde'])
        self.local_min = data['local_min']
        return self

    def fit(self, values: list[float]) -> None:
        """
        1. Normalize to 0-1 by (the maximum depth `x_max`)
        2. Fit the data to KDE
        3. Find threshold for cutting local minimal
        4. Save the thresholds (`local_min`)
        """
        assert self.kde is None
        # normalize to 0 - 1
        self.x_max = np.max(values)
        data = np.array(values)[:, None] / self.x_max
        self.kde = KernelDensity(kernel='gaussian', bandwidth=self.bandwidth).fit(data)

        # cut
        x = np.linspace(0, 1.1, self.points)
        y = self.kde.score_samples(x[:, None])
        self.local_min = list(x[argrelextrema(y, np.less, order=self.neighbor)[0]])
        print("Threshold", self.local_min)

        # for plot
        self.data = values

    def assignCN(self, values: list[float]) -> list[int]:
        """ Assign CN for each depths """
        assert self.kde is not None
        x = np.array(values) / self.x_max
        cn = np.searchsorted(self.local_min, x)
        return list(cn)

    def plot(self, title: str = "") -> list[go.Figure]:
        x = np.linspace(0, 1.1, self.points)
        y = self.kde.score_samples(x[:, None])
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(go.Scatter(x=x, y=y, name="KDE"))
        fig.update_layout(yaxis_title="KDE score",
                          yaxis2_title="Fraction of samples",
                          xaxis_title="Normalized Depth")
        for cn, m in enumerate(self.local_min):
            fig.add_vline(x=m, line_width=2, line_dash="dash", annotation_text=f"cn={cn}")

        fig.add_trace(go.Histogram(
            x=np.array(self.data) / self.x_max,
            name="samples", nbinsx=100, histnorm="probability"), secondary_y=True)
        return [fig]


def loadCNModel(filename: str) -> Dist:
    """ Load model from json. """
    with open(filename) as f:
        data = json.load(f)
    if data['method'] == "KDEcut":
        return KDEcut.load(filename)
    if data['method'] == "CNgroup":
        return CNgroup.load(filename)
    raise NotImplementedError
