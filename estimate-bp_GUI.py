#!/usr/bin/env python3

import math
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, MultipleLocator

R = 8.314  # Gas constant (J/mol·K)
STD_PRESSURE = 760.0  # Atomospheric pressure (Torr)

# Plot style
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "font.size": 8,
        "axes.linewidth": 1.0,
    }
)


class Molecule:
    def __init__(self, name, Tb=273.0, cal_method="Trouton's rule"):
        self.name = name  # Name of the molecule
        self.known_bp = [Tb, STD_PRESSURE]  # 既知の圧力下での沸点 [K, Torr]
        self.is_set_from_known_bp = False
        self.Tb = Tb  # 大気圧下での沸点 (K)

        self.delta_Hs = {  # モル蒸発エンタルピーの関数。xは大気圧下での沸点Tb
            "Trouton's rule": lambda x: x * 85,
            "Methane": lambda x: x * 73,
            "Water": lambda x: x * 109,
            "T-H-E rule": lambda x: (4.4 + math.log(x)) * R * x,
        }
        self.cal_method = cal_method  # モル蒸発エンタルピーの計算方法

    def set_from_Tb(self, Tb):  # Tbが既知の場合のTbの設定
        self.Tb = Tb
        self.known_bp = [Tb, STD_PRESSURE]

    def set_from_known_bp(self, known_T, known_p):  # Tbが未知の場合のTbの設定
        self.is_set_from_known_bp = True
        if self.cal_method == "T-H-E rule":  # T-H-E則の場合は解析的に解けないのでTrouton's ruleを使用
            cal_method = "Trouton's rule"
        else:
            cal_method = self.cal_method
        self.Tb = known_T * (1 - R / self.delta_Hs[cal_method](1) * math.log(known_p / STD_PRESSURE))
        self.known_bp = [known_T, known_p]

    def recalc_Tb(self):  # ΔHの計算手法に変更がありかつTbが未知の場合にTbを再計算
        if self.is_set_from_known_bp:
            self.set_from_known_bp(self.known_bp[0], self.known_bp[1])

    def boiling_point_at(self, p):  # Clausius-Clapeyron式により圧力pにおける沸点を計算
        delta_H = self.delta_Hs[self.cal_method](self.Tb)

        return 1 / ((1 / self.Tb) - (math.log(p / STD_PRESSURE) * R / delta_H))


# Sidebar for molecule settings
def molecule_sidebar(mol: Molecule, idx: int):
    mol.name = st.sidebar.text_input("Name", value=mol.name, key=f"name_{idx}")

    use_other_pressure = st.sidebar.checkbox("Use other pressure", key=f"useP_{idx}")
    if use_other_pressure:
        P = st.sidebar.number_input("Pressure (Torr)", value=760.0, key=f"P_{idx}")
        T = st.sidebar.number_input("Temperature (℃)", value=mol.Tb - 273, key=f"T_{idx}") + 273
        mol.set_from_known_bp(T, P)
    else:
        Tb_input = st.sidebar.number_input("b.p. at 760 Torr (℃)", value=mol.Tb - 273, key=f"Tb_{idx}") + 273
        mol.set_from_Tb(Tb_input)

    mol.cal_method = st.sidebar.selectbox("Use ΔH of ...", list(mol.delta_Hs), key=f"method_{idx}")

    mol.recalc_Tb()


def main():
    # 初期設定. インスタンス化
    if "mol1" not in st.session_state:
        st.session_state.mol1 = Molecule("Molecule 1")
        st.session_state.mol2 = Molecule("Molecule 2")
        st.session_state.mol3 = Molecule("Molecule 3")

    # sidebar
    st.sidebar.title("Settings")
    molecule_sidebar(st.session_state.mol1, 1)
    st.sidebar.divider()
    add_second = st.sidebar.checkbox("Add molecule")
    if add_second:
        molecule_sidebar(st.session_state.mol2, 2)
        st.sidebar.divider()
    add_third = st.sidebar.checkbox("Add molecule", key="add3")
    if add_third:
        molecule_sidebar(st.session_state.mol3, 3)

    # main window
    st.title("Boiling Points Estimator")
    pressure_range = st.slider("Pressure range (Torr)", 0.01, 10.0, (0.1, 7.0))
    pressures = [p * 0.01 for p in range(int(pressure_range[0] * 100), int(pressure_range[1] * 100))]
    bp1 = [st.session_state.mol1.boiling_point_at(p) - 273 for p in pressures]

    # Print Tb
    st.divider()
    st.write(f"Boiling Points at {STD_PRESSURE} Torr")
    st.write(f"- {st.session_state.mol1.name}: {st.session_state.mol1.Tb - 273:.2f} ℃")
    if add_second:
        st.write(f"- {st.session_state.mol2.name}: {st.session_state.mol2.Tb - 273:.2f} ℃")
    if add_third:
        st.write(f"- {st.session_state.mol3.name}: {st.session_state.mol3.Tb - 273:.2f} ℃")

    # graph
    st.divider()
    fig, ax = plt.subplots()
    ax.set_title("Pressure vs Boiling Points")
    ax.plot(pressures, bp1, label=st.session_state.mol1.name)
    if add_second:
        bp2 = [st.session_state.mol2.boiling_point_at(p) - 273 for p in pressures]
        ax.plot(pressures, bp2, label=st.session_state.mol2.name)
    if add_third:
        bp3 = [st.session_state.mol3.boiling_point_at(p) - 273 for p in pressures]
        ax.plot(pressures, bp3, label=st.session_state.mol3.name)
    ax.set_xlim(pressure_range)
    ax.set_xlabel("Pressure (Torr)")
    ax.set_ylabel("Temperature (℃)")
    ax.xaxis.set_major_locator(LinearLocator(10))
    ax.xaxis.set_minor_locator(LinearLocator(50))
    ax.yaxis.set_major_locator(LinearLocator(10))
    ax.yaxis.set_minor_locator(LinearLocator(50))
    plt.gca().xaxis.set_major_formatter(plt.FormatStrFormatter("%.1f"))  # x軸小数点以下1桁表示
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter("%.1f"))  # y軸小数点以下1桁表示
    ax.legend(loc="lower right")
    ax.grid(which="minor")
    st.pyplot(fig)


if __name__ == "__main__":
    main()
