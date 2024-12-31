import numpy as np
from matplotlib import pyplot as plt
from src.solve import solve_on_composed_materials, binary_search_l2
from src.plot import plot_temperature, plot_binary_search_l2


if __name__ == "__main__":
    # 热传导率，单位 W/(m*K)
    lambdas = [0.082, 0.37, 0.045, 0.028, 0.42]
    # 密度，单位 kg/m^3
    rhos = [300, 862, 74.2, 1.18, 1090]
    # 热容，单位 J/(kg*K)
    cs = [1377, 2100, 1726, 1005, 3350]
    # 厚度，单位 m
    l2 = 6e-3
    l4 = 5e-3
    ls = [0.6e-3, l2, 3.6e-3, l4, 2e-3]
    
    # 与测量值比较，验证算法正确性
    u_outer = 75
    u_composed, u_skin = solve_on_composed_materials(u_outer, ls, lambdas, rhos, cs)
    plot_temperature(u_skin, "appendices/CUMCM-2018-Problem-A-Chinese-Appendix.xlsx", "附件2")
    
    # 问题(1) 二分法求解 l2
    ## 画出皮肤温度随l2变化的曲线
    l2_values = np.linspace(0.6e-3, 25e-3, 20)
    u_skins_55 = []
    u_skins_60 = []
    for l2 in l2_values:
        ls[1] = l2
        u_outer = 65
        _, u_skin = solve_on_composed_materials(u_outer, ls, lambdas, rhos, cs)
        u_skins_55.append(u_skin[55 * 60])
        u_skins_60.append(u_skin[60 * 60])
    plt.plot(l2_values * 1000, u_skins_55, label='Skin Temp. at 55min')
    plt.plot(l2_values * 1000, u_skins_60, label='Skin Temp. at 60min')
    plt.legend()
    plt.title('Skin Temperature')
    plt.xlabel('l2/mm')
    plt.ylabel('Temperature/℃')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('figures/skin_temperature_vs_l2.png')
    plt.close()
    ## 二分法求解l2
    ls = [0.6e-3, None, 3.6e-3, 5.5e-3, 2e-3]
    l2, a_s, b_s, fa55s, fb55s, fa60s, fb60s = binary_search_l2(0.6e-3, 25e-3, 1e-6, 55, 60, 44.0, 47.0, 65, ls, lambdas, rhos, cs)
    print(f"l2 = {l2}")
    plot_binary_search_l2(a_s, b_s, fa55s, fb55s, 'figures/binary_search_l2_55min.png', 'Binary Search for l2 (65℃, 55min)')
    plot_binary_search_l2(a_s, b_s, fa60s, fb60s, 'figures/binary_search_l2_60min.png', 'Binary Search for l2 (65℃, 60min)')
    
    # 问题(2) l4最大时，二分法求解 l2
    ls = [0.6e-3, None, 3.6e-3, 6.4e-3, 2e-3]
    l2, a_s, b_s, fa55s, fb55s, fa60s, fb60s = binary_search_l2(6.4e-3, 25e-3, 1e-6, 25, 30, 44.0, 47.0, 80, ls, lambdas, rhos, cs)
    print(f"l2 = {l2}")
    plot_binary_search_l2(a_s, b_s, fa55s, fb55s, 'figures/binary_search_l2_25min.png', 'Binary Search for l2 (80℃, 25min)')
    plot_binary_search_l2(a_s, b_s, fa60s, fb60s, 'figures/binary_search_l2_30min.png', 'Binary Search for l2 (80℃, 30min)')