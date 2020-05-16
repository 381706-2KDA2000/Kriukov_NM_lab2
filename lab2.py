# -*- coding: utf-8 -*-
"""
Created on Wed May 13 13:04:59 2020

@author: dimentor
"""
import PySimpleGUI as sg

import matplotlib.pyplot as plt

from math import cos

pi = 3.1415926535

def func(x, l, f1, f2):
    return 1/l + f1 * cos((pi*x)/l) + f2 * cos(2*(pi*x)/l)

def bfunc(x, l, b0, b1, b2):
    return b0 + b1 * cos((pi*x)/l) + b2 * cos(2*(pi*x)/l)


def integrate(h, fu): #Фунуция численного интегрирования
    res = (h/3)*(fu[0] + fu[len(fu) - 1])
    for i in range(1, len(fu) - 1, 2):
        res += (h/3)*(4*fu[i] + 2*fu[i + 1])
    return res
        
def tridiagAlg(a, b, c, func, count):#Метод прогонки
    A = []
    B = []
    res = [0] * count
    A.append(-c[0]/b[0])
    B.append(func[0]/b[0])
    for i in range(1, count):
        A.append(-c[i] / (a[i] * A[i - 1] + b[i]))
        B.append((func[i] - a[i] * B[i - 1]) / (a[i] * A[i - 1] + b[i]))
    res[count-1] = B[count - 1]
    for i in range(count - 2, -1, -1):
        res[i] = (A[i] * res[i + 1] + B[i])
    return res

fig, ax = plt.subplots(figsize=(6, 4))
ax.grid()
ax.plot([0], [0])
fig.savefig('NM2plot.png')
sg.theme('DarkAmber')
sg.SetOptions(input_elements_background_color='#FFFFFF',
              input_text_color = '#000000')

layout = [[sg.Button('Plot'), sg.Button('Cancel'), sg.ProgressBar(1000,
           orientation='h', size=(20, 20), key='progressbar')],
            [sg.Button('Alternative solution')],
           [sg.Text('Длина стержня', size=(16, 1)), sg.Input('7', size=(8, 1)),
            sg.Text('', size=(4, 1)), sg.Text('b', size=(8, 1)),
            sg.Text('f', size=(8, 1))],
          [sg.Text('Время воздействия', size=(16, 1)),
           sg.Input('1', size=(8, 1)), sg.Text('0', size=(2, 1)),
           sg.Input('0', size=(8, 1)), sg.Input('1/len', size=(8, 1), key = 'f0',
                    disabled = True, background_color = 'DarkBlue')],
          [sg.Text('Шаг по времени', size=(16, 1)),
           sg.Input('0.01', size=(8, 1)), sg.Text('1', size=(2, 1)),
           sg.Input('-7', size=(8, 1)), sg.Input('0', size=(8, 1))],
          [sg.Text('Шаг по x', size=(16, 1)), sg.Input('0.2', size=(8, 1)),
           sg.Text('2', size=(2, 1)), sg.Input('0', size=(8, 1)),
           sg.Input('0', size=(8, 1))],
          [sg.Image(r'NM2plot.png', key='NM2plot')]]

window = sg.Window('NM lab2', layout)

func_val = []
resA = []
resB = []
x_val = []

flag = bool(0)

while True:
    event, values = window.Read(timeout = 100)
    
    if event in (None, 'Cancel'):
        break
    
    if event in (None, 'Alternative solution'):
        if(flag == 0):
            ax.plot(x_val, func_val, 'b')
            ax.plot(x_val, resA, 'g')
            ax.plot(x_val, resB, 'r')
            fig.savefig('NM2plot.png')
            window['NM2plot'].update(r'NM2plot.png')
            flag = 1
        else:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.grid()
            ax.plot(x_val, func_val, 'b')
            ax.plot(x_val, resB, 'r')
            fig.savefig('NM2plot.png')
            window['NM2plot'].update(r'NM2plot.png')
            flag = 0
        
    if event in ('Plot'):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.grid()
        fig.savefig('NM2plot.png')
        progress_bar = window['progressbar']
        window['NM2plot'].update(r'NM2plot.png')
        try:
            _len = float(values[0])
            window['f0'].update(round(1/_len, 4))
            time = float(values[1])
            b0 = float(values[2])
            delta_t = float(values[3])
            b1 = float(values[4])
            f1 = float(values[5])
            delta_x = float(values[6])
            b2 = float(values[7])
            f2 = float(values[8])
        except:
            continue
        func_val = []
        bfunc_val = []

        slices1 = [[]]
        slices2 = [[]]
        count_N = int(_len/delta_x)
        count_T = int(time/delta_t)
        prbar_step = 1000/count_T
        
        # Вычисление значений функции и заполнение первого слоя сетки
        for i in range(0, count_N):
            func_val.append(func(i*delta_x, _len, f1, f2))
            bfunc_val.append(bfunc(i*delta_x, _len, b0, b1, b2))
            slices1[0].append(func_val[i])
            slices2[0].append(func_val[i])

        # Заполнение матрицы коэффициентов для метода прогонки
        coeff_a = [0.0]
        coeff_b = [1.0]
        coeff_c = [-1.0]
        for i in range(1, count_N - 1):
            coeff_a.append(delta_t / (delta_x * delta_x))
            coeff_b.append(-1 - 2*delta_t / (delta_x * delta_x))
            coeff_c.append(delta_t / (delta_x * delta_x))
        coeff_a.append(-1.0)
        coeff_b.append(1.0)
        coeff_c.append(0.0)
        
        #Вычисление последующих слоев сетки
        for i in range(1, count_T):
            I = integrate(delta_x, bfunc_val)
            fu = [0]
            fu2 = [0]
            slices1.append([])
            slices2.append([])
            
            #Вычисляем правую часть системы для прогонки
            for j in range(1, count_N - 1):
                fu.append(-slices1[i - 1][j] * ((bfunc_val[j] - I) * delta_t * delta_t  + 1.0))
                fu2.append(-slices2[i - 1][j] * (bfunc_val[j] * delta_t * delta_t + 1.0))
            fu.append(0)
            fu2.append(0)
        
            #Метод прогонки для системы из B
            res = tridiagAlg(coeff_a, coeff_b, coeff_c, fu, count_N)
            for j in range(0, count_N):
                slices1[i].append(res[j])
            
            #Метод прогонки для системы из A
            res2 = tridiagAlg(coeff_a, coeff_b, coeff_c, fu2, count_N)
            for j in range(0, count_N):
                slices2[i].append(res2[j])
            progress_bar.UpdateBar(i * prbar_step)
                
        I = integrate(delta_x, slices2[count_T - 1])
        
        resB = []
        resA = []
        for j in range(0, count_N):
            resA.append(slices2[count_T - 1][j] / I)
            resB.append(slices1[count_T - 1][j])
        
        x_val = []
        
        for i in range(0, count_N):
            x_val.append(i * delta_x)
        
        ax.plot(x_val, func_val, 'b')
        ax.plot(x_val, slices1[count_T - 1], 'r')
        fig.savefig('NM2plot.png')
        window['NM2plot'].update(r'NM2plot.png')
        progress_bar.UpdateBar(1000)
        flag = 0
        
window.close()
