import pandas as pd
import matplotlib.pyplot as plt


def plot_temperature(u_estimated, file_path, sheet_name):
    # 加载Excel文件
    data = pd.read_excel(file_path, sheet_name=sheet_name)

    # 跳过第一行并使用第二行作为标题
    data = pd.read_excel(file_path, sheet_name=sheet_name, header=1)

    # 假设Excel文件有'时间 (s)'和'温度 (ºC)'列
    time = data['时间 (s)']
    u_gt = data['温度 (ºC)']

    # Plot the data
    # plt.figure(figsize=(10, 6))
    # plt.plot(time, temperature, marker='o', linestyle='-', color='b')
    plt.clf()
    plt.figure()
    plt.plot(time / 60, u_estimated, label='Estimated Temp.')
    plt.plot(time / 60, u_gt, label='Measured Temp.')
    plt.legend()
    plt.title('Skin Temperature')
    plt.xlabel('Time/min')
    plt.ylabel('Temperature/℃')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('figures/temperature_plot.png')
    plt.close()
    
    
def plot_binary_search_l2(a_s, b_s, fas, fbs, file_path, title):
    plt.clf()
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Skin Temperature/℃', color='tab:blue')
    ax1.plot(fas, marker='o', label='Skin Temp. >=', color='tab:blue')
    ax1.plot(fbs, marker='x', label='Skin Temp. <=', color='tab:cyan')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.legend(loc='lower right')
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.set_ylabel('l2/mm', color='tab:red')
    ax2.plot(a_s, marker='s', label='l2 >=', color='tab:red')
    ax2.plot(b_s, marker='d', label='l2 <=', color='tab:orange')
    ax2.tick_params(axis='y', labelcolor='tab:red')
    ax2.legend(loc='upper right')

    plt.title(title)
    fig.tight_layout()
    plt.savefig(file_path)
    plt.close()