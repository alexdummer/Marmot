import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("tsai_wu_plasticity_history.csv", delimiter=",")

# plot stress strain for each stress and associated strain for 6 compone
fig, axs = plt.subplots(2, 3, figsize=(15, 10))

sigma_labels = ["σ11", "σ22", "σ33", "σ12", "σ13", "σ23"]
eps_labels = ["ε11", "ε22", "ε33", "ε12", "ε13", "ε23"]


axs = axs.flatten()
for i in range(6):

    # 6 stress components ar in col 1,2,3,4,5,6
    # 6 strain components are in col 7,8,9,10,11,12

    # plot stress vs strain
    axs[i].plot(data[:, 7 + i], data[:, 1 + i], label=f"Stress Component {i + 1} vs Strain Component {i + 7}")
    axs[i].set_xlabel(eps_labels[i])
    axs[i].set_ylabel(sigma_labels[i])
    axs[i].grid()


plt.tight_layout()
plt.show()
