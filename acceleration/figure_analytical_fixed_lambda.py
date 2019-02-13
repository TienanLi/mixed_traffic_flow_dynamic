import matplotlib.pyplot as plt
from math import floor
import numpy as np
from expectation_function import multi_value_simulation

def get_lamb(q, m): # unit in h_0
    return (1 / q - 1) * m

def get_o(u, v_0, a):
    return (u-v_0) ** 2 / (2 * u *a)

def get_term1(u, v_0, a_s, a_l, p):
    return ((u - v_0) ** 2) / 2 / u * (1 / a_s - p * (a_l - a_s) / a_l / a_s)

def get_term2(h_0, a_s, a_l, p, o_l, lamb, m, K_star_2, type, K_star_1):
    R_ls=(a_l / a_s - 1)
    if type == 'homo':
        coef = 1
    if type == 'heter':
        coef = m

    p=p**coef
    term_1 = o_l * (1 - p ** (K_star_1))
    term_2 = (o_l + h_0) * (p ** (K_star_1) - p ** (K_star_2))
    term_3 = lamb * (( p ** (K_star_1) - p ** (K_star_2)) / (1 - p ) + K_star_1 * p ** (K_star_1) - K_star_2 * p ** (K_star_2))
    term_3_numerical=0
    for k in range(K_star_1+1,K_star_2+1):
        term_3_numerical += k*(p**(k-1))*(1-p)
    print(K_star_1,K_star_2,p,term_3,term_3_numerical)

    return R_ls * p * (term_1 + term_2 - term_3_numerical)

def analytical_result(u, v_0, a_s, a_l, q, m, p, h_0, type):
    lamb = get_lamb(q, m) * h_0
    o_l = get_o(u, v_0, a_l)
    K_star_2 = floor((o_l  +h_0) / lamb)
    K_star_1 = floor(h_0 / lamb)

    term_1 = get_term1(u, v_0, a_s, a_l, p)
    term_2 = get_term2(h_0, a_s, a_l, p, o_l, lamb, m, K_star_2, type, K_star_1)
    return [term_1/h_0, term_2/h_0, (term_1 + term_2) / h_0]

def drawing_analytical():
    p_range=np.arange(0,1,0.1)
    homo_m_4=[analytical_result(30, 10, 1., 3., .8, 4, p, 3600 / 2200, 'homo') for p in p_range]
    heter_m_4=[analytical_result(30, 10, 1., 3., .8, 4, p, 3600 / 2200, 'heter') for p in p_range]
    y_list = [[a[0] for a in homo_m_4],
              [a[0] for a in heter_m_4],
              [a[1] for a in homo_m_4],
              [a[1] for a in heter_m_4]]
    name_list=['E[$o_{in}$], homo','E[$o_{in}$],heter',r'E[$\tilde{\omega}$], homo',r'E[$\tilde{\omega}$], heter']
    style_list=['r-*','g-*','r:v','g:v']
    x_label='$m$ ($veh$)'
    y_label=r'E[$o_{in}$] and E[$\tilde{\omega}$] ($h_0$)'
    customized_x_ticks=[0,0.3,0.6,0.9]
    customized_axis=[0, 0.9, 0, 4]
    draw_fig('p_range.png', p_range, y_list, name_list, style_list, x_label, y_label, customized_x_ticks, customized_axis)

    m_range = range(1, 12)

    homo_2=[analytical_result(30, 10, 1., 3., .8, m, 0.2, 3600 / 2200, 'homo') for m in m_range]
    homo_5=[analytical_result(30, 10, 1., 3., .8, m, 0.5, 3600 / 2200, 'homo') for m in m_range]
    homo_8=[analytical_result(30, 10, 1., 3., .8, m, 0.8, 3600 / 2200, 'homo') for m in m_range]
    y_list=[[a[0] for a in homo_2],[a[0] for a in homo_5],[a[0] for a in homo_8],[a[1] for a in homo_2],[a[1] for a in homo_5],[a[1] for a in homo_8]]
    name_list=['E[$o_{in}$], p=0.2','E[$o_{in}$], p=0.5','E[$o_{in}$], p=0.8',r'E[$\tilde{\omega}$], p=0.2',r'E[$\tilde{\omega}$], p=0.5',r'E[$\tilde{\omega}$], p=0.8']
    style_list=['r-*','g-*','b-*','r:v','g:v','b:v']
    customized_x_ticks=[0, 2, 4, 6, 8, 10]
    customized_axis=[1, 11, 0, 4]
    draw_fig('homo.png', m_range, y_list, name_list, style_list, x_label, y_label, customized_x_ticks, customized_axis)

    heter_2 = [analytical_result(30, 10, 1., 3., .8, m, 0.2, 3600 / 2200, 'heter') for m in m_range]
    heter_5 = [analytical_result(30, 10, 1., 3., .8, m, 0.5, 3600 / 2200, 'heter') for m in m_range]
    heter_8 = [analytical_result(30, 10, 1., 3., .8, m, 0.8, 3600 / 2200, 'heter') for m in m_range]
    y_list = [[a[0] for a in heter_2], [a[0] for a in heter_5], [a[0] for a in heter_8], [a[1] for a in heter_2],
              [a[1] for a in heter_5], [a[1] for a in heter_8]]
    draw_fig('heter.png', m_range, y_list, name_list, style_list, x_label, y_label, customized_x_ticks, customized_axis)

    y_list = [[a[2] for a in homo_2],
              [a[2] for a in homo_5],
              [a[2] for a in homo_8],
              [a[2] for a in heter_2],
              [a[2] for a in heter_5],
              [a[2] for a in heter_8]]
    name_list=['homo, p=0.2',
               'homo, p=0.5',
               'homo, p=0.8',
               'mixed, p=0.2',
               'mixed, p=0.5',
               'mixed, p=0.8']
    customized_axis=[1, 11, 1.8, 4.2]
    y_label=r'E[$o_{cum}$] ($h_0$)'
    style_list=['r-*','g-*','b-*','r:v','g:v','b:v']
    customized_x_ticks=[0, 2, 4, 6, 8, 10]
    draw_fig('results.png', m_range, y_list, name_list, style_list, x_label, y_label, customized_x_ticks, customized_axis)

def drawing_simulaion():
    m_range = range(1, 3)

    homo_2 = [multi_value_simulation((30-10)**2/2/30/3/3600*2200+1, .25*m, [1,3], [0.8, 0.2], m, 'homo') for m in m_range]
    homo_5 = [multi_value_simulation((30-10)**2/2/30/3/3600*2200+1, .25*m, [1,3], [0.5, 0.5], m, 'homo') for m in m_range]
    homo_8 = [multi_value_simulation((30-10)**2/2/30/3/3600*2200+1, .25*m, [1,3], [0.2, 0.8], m, 'homo') for m in m_range]
    heter_2 = [multi_value_simulation((30-10)**2/2/30/3/3600*2200+1, .25*m, [1,3], [0.8, 0.2], m, 'heter') for m in m_range]
    heter_5 = [multi_value_simulation((30-10)**2/2/30/3/3600*2200+1, .25*m, [1,3], [0.5, 0.5], m, 'heter')  for m in m_range]
    heter_8 = [multi_value_simulation((30-10)**2/2/30/3/3600*2200+1, .25*m, [1,3], [0.2, 0.8], m, 'heter')  for m in m_range]

    probability_range=np.arange(0,1.01,0.1)
    homo_p_4=[multi_value_simulation((30-10)**2/2/30/3/3600*2200+1, .25*4, [1,3], [1-p, p], 4, 'homo') for p in probability_range]
    heter_p_4=[multi_value_simulation((30-10)**2/2/30/3/3600*2200+1, .25*4, [1,3], [1-p, p], 4, 'heter') for p in probability_range]

    fig = plt.figure(figsize=(20, 12), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(probability_range, [a[0] for a in homo_2], 'r-*',
            lw=5, markersize=18, label='E[$o_{in}$], p=0.2')
    ax.plot(probability_range, [a[0] for a in homo_5], 'g-*',
            lw=5, markersize=18, label='E[$o_{in}$], p=0.5')

    ax.plot(probability_range, [a[1] for a in homo_2], 'r:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], p=0.2')
    ax.plot(probability_range, [a[1] for a in homo_5], 'g:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], p=0.5')

    # ax.legend(fontsize=32, loc=9, ncol=3, columnspacing=.5)
    plt.xticks([0, 2, 4, 6, 8, 10], [0, 2, 4, 6, 8, 10])
    ax.axis([1, 11, 0, 4])
    ax.tick_params(labelsize=38)
    plt.xlabel('$m$ ($veh$)', fontsize=44)
    plt.ylabel(r'E[$o_{in}$] and E[$\tilde{\omega}$] ($h_0$)', fontsize=44)
    # plt.ylabel('$\lambda$ ($h_0$)')
    ax.set_position([0.1, 0.1, 0.6, 0.85])
    ax.legend(fontsize=36, loc='center left', bbox_to_anchor=(1, 0, 1.2, 1))
    plt.savefig('probability.png')
    plt.close()
    plt.clf()



    fig = plt.figure(figsize=(20, 12), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(m_range, [a[0] for a in homo_2], 'r-*',
            lw=5, markersize=18, label='E[$o_{in}$], p=0.2')
    ax.plot(m_range, [a[0] for a in homo_5], 'g-*',
            lw=5, markersize=18, label='E[$o_{in}$], p=0.5')
    ax.plot(m_range, [a[0] for a in homo_8], 'b-*',
            lw=5, markersize=18, label='E[$o_{in}$], p=0.8')
    ax.plot(m_range, [a[1] for a in homo_2], 'r:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], p=0.2')
    ax.plot(m_range, [a[1] for a in homo_5], 'g:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], p=0.5')
    ax.plot(m_range, [a[1] for a in homo_8], 'b:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], p=0.8')

    # ax.legend(fontsize=32, loc=9, ncol=3, columnspacing=.5)
    plt.xticks([0, 2, 4, 6, 8, 10], [0, 2, 4, 6, 8, 10])
    ax.axis([1, 11, 0, 4])
    ax.tick_params(labelsize=38)
    plt.xlabel('$m$ ($veh$)', fontsize=44)
    plt.ylabel(r'E[$o_{in}$] and E[$\tilde{\omega}$] ($h_0$)', fontsize=44)
    # plt.ylabel('$\lambda$ ($h_0$)')
    ax.set_position([0.1, 0.1, 0.6, 0.85])
    ax.legend(fontsize=36, loc='center left', bbox_to_anchor=(1, 0, 1.2, 1))
    plt.savefig('homo.png')
    plt.close()
    plt.clf()

    y_list=[[a[0] for a in homo_2],[a[0] for a in homo_5],[a[0] for a in homo_8],[a[1] for a in homo_2],[a[1] for a in homo_5],[a[1] for a in homo_8]]
    name_list=['E[$o_{in}$], p=0.2','E[$o_{in}$], p=0.5','E[$o_{in}$], p=0.8',r'E[$\tilde{\omega}$], p=0.2',r'E[$\tilde{\omega}$], p=0.5',r'E[$\tilde{\omega}$], p=0.8']
    style_list=['r-*','g-*','b-*','r:v','g:v','b:v']
    x_label='$m$ ($veh$)'
    y_label=r'E[$o_{in}$] and E[$\tilde{\omega}$] ($h_0$)'
    customized_x_ticks=[0, 2, 4, 6, 8, 10]
    customized_axis=[1, 11, 0, 4]
    draw_fig('homo_test.png', m_range, y_list, name_list, style_list, x_label, y_label, customized_x_ticks, customized_axis)

    fig = plt.figure(figsize=(20, 12), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(m_range, [a[0] for a in heter_2], 'r-*',
            lw=5, markersize=18, label='E[$o_{in}$], p=0.2')
    ax.plot(m_range, [a[0] for a in heter_5], 'g-*',
            lw=5, markersize=18, label='E[$o_{in}$], p=0.5')
    ax.plot(m_range, [a[0] for a in heter_8], 'b-*',
            lw=5, markersize=18, label='E[$o_{in}$], p=0.8')
    ax.plot(m_range, [a[1] for a in heter_2], 'r:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], p=0.2')
    ax.plot(m_range, [a[1] for a in heter_5], 'g:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], p=0.5')
    ax.plot(m_range, [a[1] for a in heter_8], 'b:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], p=0.8')
    plt.xticks([0, 2, 4, 6, 8, 10], [0, 2, 4, 6, 8, 10])
    ax.axis([1, 11, 0, 4])
    ax.tick_params(labelsize=38)
    plt.xlabel('$m$ ($veh$)', fontsize=44)
    plt.ylabel(r'E[$o_{in}$] and E[$\tilde{\omega}$] ($h_0$)', fontsize=44)
    ax.set_position([0.1, 0.1, 0.6, 0.85])
    ax.legend(fontsize=36, loc='center left', bbox_to_anchor=(1, 0, 1.2, 1))
    plt.savefig('heter.png')
    plt.close()
    plt.clf()

    fig = plt.figure(figsize=(20, 12), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(m_range, [a[2] for a in homo_2], 'r-*',
            lw=5, markersize=18, label='homo, p=0.2')
    ax.plot(m_range, [a[2] for a in homo_5], 'g-*',
            lw=5, markersize=18, label='homo, p=0.5')
    ax.plot(m_range, [a[2] for a in homo_8], 'b-*',
            lw=5, markersize=18, label='homo, p=0.8')
    ax.plot(m_range, [a[2] for a in heter_2], 'r:v',
            lw=5, markersize=18, label='mixed, p=0.2')
    ax.plot(m_range, [a[2] for a in heter_5], 'g:v',
            lw=5, markersize=18, label='mixed, p=0.5')
    ax.plot(m_range, [a[2] for a in heter_8], 'b:v',
            lw=5, markersize=18, label='mixed, p=0.8')
    plt.xticks([0, 2, 4, 6, 8, 10], [0, 2, 4, 6, 8, 10])
    ax.axis([1, 11, 1.8, 4.2])
    ax.tick_params(labelsize=38)
    plt.xlabel('$m$ ($veh$)', fontsize=44)
    plt.ylabel('E[$o_{cum}$] ($h_0$)', fontsize=44)
    ax.set_position([0.1, 0.1, 0.6, 0.85])
    ax.legend(fontsize=36, loc='center left', bbox_to_anchor=(1, 0, 1.2, 1))
    # plt.ylabel('$\lambda$ ($h_0$)')
    plt.savefig('results.png')
    plt.close()
    plt.clf()

def draw_fig(fig_name,x,y_list,name_list,style_list,x_label,y_label,customized_x_ticks,customized_axis):
    fig = plt.figure(figsize=(20, 12), dpi=100)
    ax = fig.add_subplot(111)
    for i in range(len(y_list)):
        ax.plot(x, y_list[i], style_list[i],
            lw=5, markersize=18, label=name_list[i])
    plt.xticks(customized_x_ticks, customized_x_ticks)
    ax.axis(customized_axis)
    ax.tick_params(labelsize=38)
    plt.xlabel(x_label, fontsize=44)
    plt.ylabel(y_label, fontsize=44)
    ax.set_position([0.1, 0.1, 0.6, 0.85])
    ax.legend(fontsize=36, loc='center left', bbox_to_anchor=(1, 0, 1.2, 1))
    plt.savefig(fig_name)
    plt.close()
    plt.clf()

def comparing_n_value():
    m_range = range(1, 12)

    homo_5 = [
        multi_value_simulation((30 - 10) ** 2 / 2 / 30 / max([0.5, 1, 2, 3, 3.5]) / 3600 * 2200 + 1, .25 * m, [0.5, 1, 2, 3, 3.5],
                               [0.2, 0.2, 0.2, 0.2, 0.2], m, 'homo')
        for m in m_range]
    homo_2 = [
        multi_value_simulation((30 - 10) ** 2 / 2 / 30 / 3 / 3600 * 2200 + 1, .25 * m, [1,  3], [0.5, 0.5], m, 'homo')
        for m in m_range]


    heter_2 = [multi_value_simulation((30 - 10) ** 2 / 2 / 30 / 3 / 3600 * 2200 + 1, .25 * m, [1, 3], [0.5, 0.5], m,
                                      'heter') for m in m_range]
    heter_5 = [multi_value_simulation((30 - 10) ** 2 / 2 / 30 / max([0.5, 1, 2, 3, 3.5])/ 3600 * 2200 + 1, .25 * m, [0.5, 1, 2, 3, 3.5], [0.2, 0.2, 0.2, 0.2, 0.2], m,
                                      'heter') for m in m_range]

    fig = plt.figure(figsize=(20, 12), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(m_range, [a[0] for a in homo_2], 'r-*',
            lw=5, markersize=18, label='E[$o_{in}$], 2-value')
    ax.plot(m_range, [a[0] for a in homo_5], 'g-*',
            lw=5, markersize=18, label='E[$o_{in}$], 5-value')
    ax.plot(m_range, [a[1] for a in homo_2], 'r:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], 2-value')
    ax.plot(m_range, [a[1] for a in homo_5], 'g:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], 5-value')

    # ax.legend(fontsize=32, loc=9, ncol=3, columnspacing=.5)
    plt.xticks([0, 2, 4, 6, 8, 10], [0, 2, 4, 6, 8, 10])
    ax.axis([1, 11, 0, 4])
    ax.tick_params(labelsize=38)
    plt.xlabel('$m$ ($veh$)', fontsize=44)
    plt.ylabel(r'E[$o_{in}$] and E[$\tilde{\omega}$] ($h_0$)', fontsize=44)
    # plt.ylabel('$\lambda$ ($h_0$)')
    ax.set_position([0.1, 0.1, 0.6, 0.85])
    ax.legend(fontsize=36, loc='center left', bbox_to_anchor=(1, 0, 1.2, 1))
    plt.savefig('homo.png')
    plt.close()
    plt.clf()

    fig = plt.figure(figsize=(20, 12), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(m_range, [a[0] for a in heter_2], 'r-*',
            lw=5, markersize=18, label='E[$o_{in}$], 2-value')
    ax.plot(m_range, [a[0] for a in heter_5], 'g-*',
            lw=5, markersize=18, label='E[$o_{in}$], 5-value')
    ax.plot(m_range, [a[1] for a in heter_2], 'r:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], 2-value')
    ax.plot(m_range, [a[1] for a in heter_5], 'g:v',
            lw=5, markersize=18, label=r'E[$\tilde{\omega}$], 5-value')
    plt.xticks([0, 2, 4, 6, 8, 10], [0, 2, 4, 6, 8, 10])
    ax.axis([1, 11, 0, 4])
    ax.tick_params(labelsize=38)
    plt.xlabel('$m$ ($veh$)', fontsize=44)
    plt.ylabel(r'E[$o_{in}$] and E[$\tilde{\omega}$] ($h_0$)', fontsize=44)
    ax.set_position([0.1, 0.1, 0.6, 0.85])
    ax.legend(fontsize=36, loc='center left', bbox_to_anchor=(1, 0, 1.2, 1))
    plt.savefig('heter.png')
    plt.close()
    plt.clf()

    fig = plt.figure(figsize=(20, 12), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(m_range, [a[2] for a in homo_2], 'r-*',
            lw=5, markersize=18, label='homo, 2-value')
    ax.plot(m_range, [a[2] for a in homo_5], 'g-*',
            lw=5, markersize=18, label='homo, 5-value')
    ax.plot(m_range, [a[2] for a in heter_2], 'r:v',
            lw=5, markersize=18, label='mixed, 2-value')
    ax.plot(m_range, [a[2] for a in heter_5], 'g:v',
            lw=5, markersize=18, label='mixed, 5-value')
    plt.xticks([0, 2, 4, 6, 8, 10], [0, 2, 4, 6, 8, 10])
    ax.axis([1, 11, 2.5, 6.5])
    ax.tick_params(labelsize=38)
    plt.xlabel('$m$ ($veh$)', fontsize=44)
    plt.ylabel('E[$o_{cum}$] ($h_0$)', fontsize=44)
    ax.set_position([0.1, 0.1, 0.6, 0.85])
    ax.legend(fontsize=36, loc='center left', bbox_to_anchor=(1, 0, 1.2, 1))
    # plt.ylabel('$\lambda$ ($h_0$)')
    plt.savefig('results.png')
    plt.close()
    plt.clf()

# comparing_n_value()
drawing_analytical()