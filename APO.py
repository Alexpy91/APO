import matplotlib.pyplot as plt
import math
import warnings

warnings.filterwarnings("ignore")

# kob = 1.6
# Tob1 = 0.225  # Данные объекта управления
# Tob2 = 6.37

kob = 0.15
Tob1 = 0.225
Tob2 = 1.05

count_s0 = 10 # количество итераций АПО
z = 0.5  # Задание регулятора
q1 = 0   # 2
q1_old = q1
q2 = 12.9  # Кэффициенты ШИМ регулятора # 12
q3 = 0.9  # 0.9
kim = 0.23  # значение коэф исполнительного механизма
dt = 0.02  # Шаг dt
L = 100  # конечное значение шага при моделировании
tau = 1.75  # Значение запаздывания для объекта
h = 0.018  # 0.018
T = 5  # Период ШИМ  # Корректируемые переменные

# ____________________
Result_I = {"Iter": 0, "I": 0, "q1": 0, "q2": 0, "q3": 0, "I_old": 0,
            "q1_old": q1, "q2_old": q2, "q3_old": q3}

flag1 = 0
I_min = 200
count_iter = 0
s0 = s1 = s2 = s3 = 0

data_y = []
data_tk = []
data_u = []
data_eps = []
data_eps_KT = []
data_s1, data_s2, data_s3 = [], [], []
data_dq1, data_dq2, data_dq3 = [], [], []
data_delta = []
data_y_delta = []
data_delta1 = []
data_y_delta1 = []
data_ksi1, data_ksi2, data_ksi3 = [], [], []
data_ksi1_KT = []
data_dutk = []
data_I = []
data_I_old = []
ureg_delta = 0
data_dtk_dq1, data_dtk_dq2, data_dtk_dq3 = [], [], []

dI = 0
dI_old = 0
SI = 0
Sdq1 = Sdq2 = Sdq3 = 0
nsi = ndi1 = ndi2 = ndi3 = 0
# ----------------

I = 0
I_old = 0
dutk = 0
dtk_dq1, dtk_dq2, dtk_dq3 = 0, 0, 0

ns = int(tau / dt)
ns_delta = int(tau / dt)  # Исходные переменные


# __________________

def delta_fun(x):  # Функция для вычисления весовой функции
    global ureg_delta, y1_delta, z1_delta, z2_delta, m_delta, y_delta, y_delta_old, delta
    #  u_delta - Единичное воздействие
    if x == 1:
        u_delta = 1
    else:
        u_delta = 0
    ureg_delta += kim * u_delta * dt  # исполнительный механизм (kim/P)
    y1_delta = z1_delta  # моделирование объекта по Рунге-Кутту
    k1_delta = dt * (z2_delta - Tob2 / Tob1 * y1_delta)
    m1_delta = dt * (kob / Tob1 * ureg_delta - y1_delta / Tob1)
    k2_delta = dt * (z2_delta + m1_delta / 2 - Tob2 / Tob1 * (y1_delta + k1_delta / 2))
    m2_delta = dt * (kob / Tob1 * ureg_delta - 1 / Tob1 * (y1_delta + k1_delta / 2))
    k3_delta = dt * (z2_delta + m2_delta / 2 - Tob2 / Tob1 * (y1_delta + k2_delta / 2))
    m3_delta = dt * (kob / Tob1 * ureg_delta - 1 / Tob1 * (y1_delta + k2_delta / 2))
    k4_delta = dt * (z2_delta + m3_delta - Tob2 / Tob1 * (y1_delta + k3_delta))
    m4_delta = dt * (kob / Tob1 * ureg_delta - 1 / Tob1 * (y1_delta + k3_delta))
    z1_delta = z1_delta + 1 / 6 * (k1_delta + 2 * k2_delta + 2 * k3_delta + k4_delta)
    z2_delta = z2_delta + 1 / 6 * (m1_delta + 2 * m2_delta + 2 * m3_delta + m4_delta)

    if m_delta >= ns_delta:  # моделирование запаздывания
        m_delta = 0
    y_delta = mas_delta[m_delta]
    mas_delta[m_delta] = y1_delta
    m_delta += 1

    delta = (y_delta - y_delta_old) / dt
    y_delta_old = y_delta

    # _________________


def model():  # Функция модели объекта с регулятором
    global ureg, y1, z1, z2, m, y
    ureg += kim * u * dt  # Исполнительный механизм kim/P u - выход регул
    y1 = z1  # Объект
    k1 = dt * (z2 - Tob2 / Tob1 * y1)
    m1 = dt * (kob / Tob1 * ureg - y1 / Tob1)
    k2 = dt * (z2 + m1 / 2 - Tob2 / Tob1 * (y1 + k1 / 2))
    m2 = dt * (kob / Tob1 * ureg - 1 / Tob1 * (y1 + k1 / 2))
    k3 = dt * (z2 + m2 / 2 - Tob2 / Tob1 * (y1 + k2 / 2))
    m3 = dt * (kob / Tob1 * ureg - 1 / Tob1 * (y1 + k2 / 2))
    k4 = dt * (z2 + m3 - Tob2 / Tob1 * (y1 + k3))
    m4 = dt * (kob / Tob1 * ureg - 1 / Tob1 * (y1 + k3))
    z1 = z1 + 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    z2 = z2 + 1 / 6 * (m1 + 2 * m2 + 2 * m3 + m4)

    if m >= ns:  # Запаздывание
        m = 0
    y = mas[m]
    mas[m] = y1
    m += 1


# ЗДЕСЬ НАЧИНАЕТСЯ АЛГОРИТМ ПРОГРАММЫ

while s0 < count_s0:  # (I > 82) or (s0 == 0):
    count_iter += 1
    I_old = I
    dI_old = dI

    #  ОБНУЛЕНИЕ ПЕРЕМЕННЫХ ПЕРЕД НОВОЙ ИТЕРАЦИЕЙ

    I = 0
    dI = 0

    k1, k2, k3, k4 = 0, 0, 0, 0
    m, m1, m2, m3, m4 = 0, 0, 0, 0, 0

    k1_delta, k2_delta, k3_delta, k4_delta = 0, 0, 0, 0
    m_delta, m1_delta, m2_delta, m3_delta, m4_delta = 0, 0, 0, 0, 0
    y, y_delta = 0, 0
    y1, y1_delta = 0, 0
    z1, z2 = 0, 0
    z1_delta, z2_delta = 0, 0
    ksi_1, ksi_2, ksi_3 = 0, 0, 0
    ksi_1_KT, ksi_2_KT, ksi_3_KT = 0, 0, 0
    dq1, dq2, dq3 = 0, 0, 0
    dq1_old, dq2_old, dq3_old = 0, 0, 0
    ureg = 0
    ureg_delta = 0
    dI = 0
    y_delta_old = 0
    u_old = 0
    mas = [0 for i in range(ns)]
    mas_delta = [0 for i in range(ns)]
    delta = 0
    eps_KT = 0
    eps = 0
    s1, s2, s3 = 0, 0, 0
    tk = 0
    u = 0
    u_old = 0

    #  ________________________

    while s1 <= L:

        delta_fun(0)  # запуск дельта
        data_delta.append(delta)
        data_y_delta.append(y_delta)
        if s2 < T:  # Цикл до T пока период не закончен

            if s3 < tk:  # Пока s3 не достигла tk  dt * T < s3 <= tk + dt * T:
                if eps_KT > 0:  # Определение значение сигнала регулятора 1, -1 или 0
                    u = 1
                    s3 += dt
                    data_u.append(u)
                else:
                    u = -1
                    s3 += dt
                    data_u.append(u)
            else:
                u = 0
                s3 += dt
                data_u.append(u)

            # определяем момент разрыва
            if (u != u_old) and (u == 0):  # если новое

                # здесь опрос всего для АПО в момент разрыва
                delta_fun(1)  # запуск дельта

                dutk = u - u_old  # определение скачка -1 или +1
                print(u)
                # dtk_dq1 = 1 - q2 * ksi_1_KT - 2 * q3 * eps_KT * ksi_1_KT  # частн пр по q1
                dtk_dq2 = eps_KT - q2 * ksi_2_KT - 2 * q3 * eps_KT * ksi_2_KT  # по q2
                dtk_dq3 = -q2 * ksi_3_KT + (eps_KT * eps_KT) - 2 * q3 * eps_KT * ksi_3_KT  # q3

                # dq1 += -2 * eps * ksi_1 * dt  # Направл градиента q1
                dq2 += -2 * eps * ksi_2 * dt  # Направл градиента q2
                dq3 += -2 * eps * ksi_3 * dt  # Направл градиента q3

                # ____________________________
                # print(f" Момент разрыва {tk} DuTK = {dutk} Значение U = {u} U предыдущее = {u_old}")

                data_delta1.append(delta)

                data_y_delta1.append(y_delta)
                data_dutk.append(dutk)
                data_dtk_dq1.append(dtk_dq1)
                data_dtk_dq2.append(dtk_dq2)
                data_dtk_dq3.append(dtk_dq3)
                data_dq1.append(dq1)
                data_dq2.append(dq2)
                data_dq3.append(dq3)

                # _______________________________________-

            u_old = u  # запоминаем предыдущее значение выхода регулятора
            s2 += dt


        else:   # Иначе если мы дошли до конца периода T
            eps_KT_TST = z - y  # Определ ошибку в точке KT
            # ksi_1_KT += - (dutk * dtk_dq1 * delta)  # функция чувствительности 1
            ksi_2_KT += - (dutk * dtk_dq2 * delta)  # функция чувствительности 2
            ksi_3_KT += - (dutk * dtk_dq3 * delta)  # функция чувствительности 3
              # Определ модуляционную характеристику для след T
            if abs(eps_KT_TST) < 0.005:
                eps_KT = 0
                q1 = 0
            else:
                eps_KT = eps_KT_TST
                q1 = q1_old
            tk = q2 * eps_KT + q3 * (eps_KT * eps_KT)
            s2 = 0
            s3 = 0

            # ____________________________

            data_ksi1_KT.append(ksi_1_KT)
            data_tk.append(tk)
            data_eps_KT.append(eps_KT)
            # print(f"eps = {eps_KT} EPS_TST = {eps_KT_TST}")

        model()  # модель с регулятором

        s1 += dt
        eps = z - y  # определение ошибки системы
        I += (eps * eps) * dt  # Рассчет интегрального критерия
        # ksi_1 += - (dutk * dtk_dq1 * delta)  # функция чувствительности 1
        ksi_2 += - (dutk * dtk_dq2 * delta)  # функция чувствительности 2
        ksi_3 += - (dutk * dtk_dq3 * delta)  # функция чувствительности 3

        # __________________________
        data_y.append(y)
        data_ksi1.append(ksi_1)
        data_ksi2.append(ksi_2)
        data_ksi3.append(ksi_3)
        data_I.append(I)
        data_eps.append(eps)

    s1 = 0

    # _______Для вывода_______

    if I < I_min:
        I_min = I
        Result_I["Iter"] = count_iter
        Result_I["I"] = I_min
        Result_I["q1"] = q1
        Result_I["q2"] = q2
        Result_I["q3"] = q3
        if flag1 == 0:
            Result_I['I_old'] = I
            flag1 = 1
    print(f"№ = {count_iter} I предыдущее = {I_old} I новое = {I}")
    print(f"eps = {eps_KT}")
    # _______КОРРЕКТОР___________

    # if I < I_old:
    #     h *= 1.5    #  Коррекция шага
    # else:
    #     h *= 0.5
    if dq1 != 0 or dq2 != 0 or dq3 != 0:
        # q1 = q1 + (h * (dq1 / (math.sqrt(dq1 ** 2 + dq2 ** 2 + dq3 ** 2))))
        q2 = q2 + (h * (dq2 / (math.sqrt(dq2 ** 2 + dq2 ** 2 + dq3 ** 2))))  # Коррекция q1,2,3
        q3 = q3 + (h * (dq3 / (math.sqrt(dq3 ** 2 + dq2 ** 2 + dq3 ** 2))))
    # ______________APO

    s0 += 1

# _______ВЫВОД РЕЗУЛЬТАТОВ___________

color_line = ['', 'blue', 'green', 'red', 'yellow']  # цвета для графиков

plt.plot(data_y, "b", color=color_line[2])  # TESTING
# plt.plot(data_y_delta1, "b", color=color_line[1])  # TESTING
plt.title('Импульсы ШИМ регулятора  ')
plt.ylabel('Amplitude')
plt.xlabel('Time(sec)')
plt.grid(True)
plt.show()

print(f"RESULTS_____________________________"
      f"\nНачальное значение I = {Result_I['I_old']} \nМинимальное значение I = {Result_I['I']} "
      f"\nНачальные значения коэффициентов: q1 = {Result_I['q1_old']} "
      f"q2 = {Result_I['q2_old']} q3 = {Result_I['q3_old']} "
      f"\nЗначения коффициентов: q1 = {Result_I['q1']} "
      f"q2 = {Result_I['q2']} q3 = {Result_I['q3']}  \nНомер итерации: {Result_I['Iter']} "
      f"\n_______________________"
      f"\nКоэффициенты после итоговой коррекции: q1 = {q1} q2 = {q2} q3 = {q3} "
      f" \nЗначение I итоговое = {I} \nИтоговое число итераций = {count_iter}")


