from typing import List
import numpy as np
from tqdm import tqdm


class EulerSolver:
    def __init__(self, alphas: List[float], tau: float, h: float, m: int, ns: List[int]):
        """
        采用隐式 Euler 方法求解多层复合材料的热传导方程的初边值问题。
        注意：空间被分为 n = sum(ns) + 1 个网格，其中第一个网格和最后一个网格为边界条件，仅有 n - 1 个网格待求解。
        Args:
            alphas (List[float]): 各层热扩散系数，单位 m^2/s
            tau (float): 时间步长
            h (float): 空间步长
            m (int): 时间步数
            ns (List[int]): 各层材料空间步数
        """
        for n in ns:
            assert n > 2
        self.alphas = alphas
        self.h = h
        self.tau = tau
        self.m = m
        self.ns = ns
        
    def create_disecret_matrix(self) -> np.array:
        """
        创建离散矩阵
        该方法生成一个三对角矩阵，用于离散化热方程。

        Returns:
            np.array: 离散矩阵
        """
        lowers = [- self.tau * alpha / self.h ** 2 for alpha in self.alphas]
        mids = [1 + 2 * self.tau * alpha / self.h ** 2 for alpha in self.alphas] 
        uppers = [- self.tau * alpha / self.h ** 2 for alpha in self.alphas]
        
        lower_diag = np.diag([lowers[i] for i in range(len(self.ns)) for _ in range(self.ns[i])][1: -1], -1)
        mid_diag = np.diag([mids[i] for i in range(len(self.ns)) for _ in range(self.ns[i])][: -1])
        upper_diag = np.diag([uppers[i] for i in range(len(self.ns)) for _ in range(self.ns[i])][1: -1], 1)
        
        tri_diag = mid_diag + lower_diag + upper_diag
        for i in range(1, len(self.ns)):
            row = sum(self.ns[:i]) - 1
            tri_diag[row, row - 1: row + 2] = [-1, 1 + self.alphas[i] / self.alphas[i - 1], -self.alphas[i] / self.alphas[i - 1]]
        return tri_diag
        

    def solve(self, ut0: np.array, ux0: np.array, ux1: np.array) -> np.array:
        """
        使用隐式 Euler 方法求解热传导方程。

        参数:
            ut0 (np.array): 时间区间起始温度分布。
            ux0 (np.array): 边界温度分布，为外界温度，恒定。
            ux1 (np.array): 边界温度分布，为体温，恒定。

        返回:
            np.array: 求解热传导方程后的温度分布。
        """
        assert ut0.shape[0] == sum(self.ns) + 1
        assert ux0.shape[0] == self.m + 1
        assert ux1.shape[0] == self.m + 1
        
        u = np.zeros((self.m + 1, sum(self.ns) + 1))  # 时间步长 * 空间步长
        u[0, :] = ut0
        u[:, 0] = ux0
        u[:, -1] = ux1
        u[1:, 1] = self.tau * self.alphas[0] / self.h ** 2 * ux0[1:]
        u[1:, -2] = self.tau * self.alphas[-1] / self.h ** 2 * ux1[1:]
        D = self.create_disecret_matrix()
        
        for j in tqdm(range(1, self.m + 1)):
            u[j, 1: -1] += u[j - 1, 1: -1]
            for i in range(1, len(self.ns)):
                u[j, sum(self.ns[:i])] = 0
            u[j, 1: -1] = np.linalg.solve(D, u[j, 1: -1])
            
        return u
    
    
def solve_on_composed_materials(u_outer: float, ls: List[float], lambdas: List[float], rhos: List[float], cs: List[float]):
    
    # 热扩散系数，单位 m^2/s
    alphas = [lambdas[i] / (rhos[i] * cs[i]) for i in range(len(lambdas))]
    
    h = 1e-4
    ns = [round(l / h) for l in ls]

    # 时间，单位 s
    T = 90 * 60
    m = T
    tau = T / m
    
    solver = EulerSolver(alphas, tau, h, m, ns)
    
    ut0 = np.array([u_outer] + [37] * sum(ns))
    ux0 = np.array([u_outer] * (m + 1))
    ux1 = np.array([37] * (m + 1))
    
    u = solver.solve(ut0, ux0, ux1)
    u_skin = u[:, sum(ns[:-1])]

    return u, u_skin


def binary_search_l2(a: float, 
                     b: float, 
                     tol: float, 
                     t1: int,
                     t2: int,
                     lim_t1: float, 
                     lim_t2: float, 
                     u_outer: float,
                     ls: List[float], 
                     lambdas: List[float], 
                     rhos: List[float], 
                     cs: List[float]):
    """
    二分法求解方程 f(x) = 0 的近似解。
    Args:
        f (function): 方程 f(x) = 0
        a (float): 区间左端点
        b (float): 区间右端点
        tol (float): 误差容限
    Returns:
        float: 方程 f(x) = 0 的近似解
    """
    def f(l2):
        ls[1] = l2
        _, u_skin = solve_on_composed_materials(u_outer, ls, lambdas, rhos, cs)
        return u_skin[t1 * 60] - lim_t1, u_skin[t2 * 60] - lim_t2
    
    fa55, fa60 = f(a)
    fb55, fb60 = f(b)
    a_s, b_s, fa55s, fb55s, fa60s, fb60s= [], [], [], [], [], []
    
    assert a < b
    assert (fa55 > 0 or fa60 > 0) and (fb55 < 0 and fb60 < 0)

    while b - a > tol:
        
        fa55s.append(fa55 + lim_t1)
        fa60s.append(fa60 + lim_t2)
        fb55s.append(fb55 + lim_t1)
        fb60s.append(fb60 + lim_t2)
        a_s.append(a * 1000)
        b_s.append(b * 1000)
        
        c = (a + b) / 2
        fc55, fc60 = f(c)
        print(c, fc55 + lim_t1, fc60 + lim_t2)
        
        if (fc55 == 0 and fc60 < 0) or (fc55 < 0 and fc60 == 0):
            fa55s.append(fc55 + lim_t1)
            fb55s.append(fc55 + lim_t1)
            fa60s.append(fc60 + lim_t2)
            fb60s.append(fc60 + lim_t2)
            break
        elif (fc55 < 0 and fc60 < 0):
            b = c
            fb55 = fc55
            fb60 = fc60
        else:
            a = c
            fa55 = fc55
            fa60 = fc60

    c = (a + b) / 2
    return c, a_s, b_s, fa55s, fb55s, fa60s, fb60s