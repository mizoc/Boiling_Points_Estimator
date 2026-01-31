#!/usr/bin/env python3

import math
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, MultipleLocator

R = 8.314  # 気体定数
plt.rcParams["font.family"] = "sans-serif"  # 使用するフォント
plt.rcParams["xtick.direction"] = "in"  # x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams["ytick.direction"] = "in"  # y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams["xtick.major.width"] = 1.0  # x軸主目盛り線の線幅
plt.rcParams["ytick.major.width"] = 1.0  # y軸主目盛り線の線幅
plt.rcParams["font.size"] = 8  # フォントの大きさ
plt.rcParams["axes.linewidth"] = 1.0  # 軸の線幅edge linewidth。囲みの太さ


class Molecule:
    def __init__(self, name, Tb=273.0, cal_method="Trouton's rule"):
        self.name = name  # Name of the molecule
        self.cal_method = cal_method
        self.delta_Hs = {
            "Trouton's rule": lambda x: x * 85,
            "Methane": lambda x: x * 73,
            "Water": lambda x: x * 109,
            "T-H-E rule": lambda x: (4.4 + math.log(x)) * R*x,
        }
        self.known_bp = [Tb, 760]  # 既知の圧力 (Torr)
        self.Tb = Tb  # 沸点 (K)
        self.is_set_from_known_bp = False

    def set_from_Tb(self, Tb):
        self.Tb = Tb
        self.known_bp = [Tb, 760]

    def set_from_known_bp(self, known_T, known_p):
        self.is_set_from_known_bp = True
        if self.cal_method == "T-H-E rule":  # 解析的に解けないのでTrouton's ruleを使用
            cal_method = "Trouton's rule"
        else:
            cal_method = self.cal_method
        self.Tb = known_T * (1 - R / self.delta_Hs[cal_method](1) * math.log(known_p / 760))
        self.known_bp = [known_T, known_p]

    def re_cal_Tb(self):
        if self.is_set_from_known_bp:
            self.set_from_known_bp(self.known_bp[0], self.known_bp[1])

    def clausiusclapeyron(self, p):
        delta_H = self.delta_Hs[self.cal_method](self.Tb)  # エンタルピー (J/mol)

        return 1 / ((1 / self.Tb) - (math.log(p / 760) * R / delta_H))


def main():
    if "use_other_pressure" not in st.session_state:  # 初期設定
        st.session_state.use_other_pressure = [False, False]
        st.session_state.is_molecule2 = False
        st.session_state.molecule1 = Molecule(name="Molecule1")
        st.session_state.molecule2 = Molecule(name="Molecule2")

    # side bar
    st.sidebar.title("Setting")
    # molecular setting
    st.session_state.molecule1.name = st.sidebar.text_input("Name", value=st.session_state.molecule1.name)
    st.session_state.molecule1.set_from_Tb(
        st.sidebar.number_input(
            "b.p. at 760 Torr (℃)", key=0, value=st.session_state.molecule1.Tb - 273, disabled=st.session_state.use_other_pressure[0]
        )
        + 273
    )
    if use_other_pressure1 := st.sidebar.checkbox("Use other pressure", value=False, key=1):
        p = st.sidebar.number_input("Pressure (Torr)", key=2, value=760.0)
        T = st.sidebar.number_input("Temperature (℃)", key=3, value=st.session_state.molecule1.Tb - 273) + 273
        st.session_state.molecule1.set_from_known_bp(T, p)
    if use_other_pressure1 != st.session_state.use_other_pressure[0]:
        st.session_state.use_other_pressure[0] = use_other_pressure1
        if not st.session_state.is_molecule2:
            st.rerun()

    st.session_state.molecule1.cal_method = st.sidebar.selectbox("Use ΔH of ...", list(st.session_state.molecule1.delta_Hs), index=0, key=4)
    st.session_state.molecule1.re_cal_Tb()

    st.sidebar.divider()
    if is_molecule2 := st.sidebar.checkbox("Add molecule"):
        st.session_state.is_molecule2 = True
        st.session_state.molecule2.name = st.sidebar.text_input("Name", value=st.session_state.molecule2.name)
        st.session_state.molecule2.set_from_Tb(
            st.sidebar.number_input(
                "b.p. at 760 Torr (℃)", key=5, value=st.session_state.molecule2.Tb - 273, disabled=st.session_state.use_other_pressure[1]
            )
            + 273
        )
        if use_other_pressure2 := st.sidebar.checkbox("Use other pressure", value=False, key=6):
            p = st.sidebar.number_input("Pressure (Torr)", key=7, value=760.0)
            T = st.sidebar.number_input("Temperature (℃)", key=8, value=st.session_state.molecule2.Tb - 273) + 273
            st.session_state.molecule2.set_from_known_bp(T, p)
        if use_other_pressure2 != st.session_state.use_other_pressure[1]:
            st.session_state.use_other_pressure[1] = use_other_pressure2
            st.rerun()

        st.session_state.molecule2.cal_method = st.sidebar.selectbox("Use ΔH of ...", list(st.session_state.molecule2.delta_Hs), index=0, key=9)
        st.session_state.molecule2.re_cal_Tb()

    # main window
    st.title("b.p. estimater")
    pressure_range = st.slider("Pressure range (Torr)", 0.01, 10.0, (0.1, 7.0))
    pressure = [i * 0.01 for i in range(round(pressure_range[0] * 100), round(pressure_range[1] * 100))]
    st.session_state.molecule1_bp = [st.session_state.molecule1.clausiusclapeyron(p) - 273 for p in pressure]
    if is_molecule2:
        st.session_state.molecule2_bp = [st.session_state.molecule2.clausiusclapeyron(p) - 273 for p in pressure]

    # graph
    # plt.style.use("dark_background")
    st.write(f"{st.session_state.molecule1.name}: Tb = {round(st.session_state.molecule1.Tb-273, 3)} ℃")
    if is_molecule2:
        st.write(f"{st.session_state.molecule2.name}: Tb = {round(st.session_state.molecule2.Tb-273, 3)} ℃")
        # diff=np.where(np.abs(np.array(st.session_state.molecule1_bp)-np.array(st.session_state.molecule2_bp))>20)
        # print(diff[0])
        # print(pressure_range[diff[0]])

    fig, ax = plt.subplots()
    ax.plot(pressure, st.session_state.molecule1_bp, label=st.session_state.molecule1.name)
    ax.set_xlim(pressure_range)
    # ax.set_title("Title")
    ax.set_xlabel("Pressure (Torr)")
    ax.set_ylabel("b.p. (℃)")
    ax.xaxis.set_major_locator(LinearLocator(10))
    ax.xaxis.set_minor_locator(LinearLocator(50))
    ax.yaxis.set_major_locator(LinearLocator(10))
    ax.yaxis.set_minor_locator(LinearLocator(50))
    if is_molecule2:
        ax.plot(pressure, st.session_state.molecule2_bp, label=st.session_state.molecule2.name)
    plt.gca().xaxis.set_major_formatter(plt.FormatStrFormatter("%.1f"))  # y軸小数点以下1桁表示
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter("%.1f"))  # y軸小数点以下1桁表示
    ax.legend(loc=4)
    ax.grid(which="minor")
    st.pyplot(fig)


if __name__ == "__main__":
    main()
