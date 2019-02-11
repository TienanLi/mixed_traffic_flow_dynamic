from math import e, factorial, sqrt,pi,log
from scipy.special import gammainc, gamma, hyp1f1
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
from resolving import Approximate_solution, ran_type, ran_type_multiple
import random
import statistics
import matplotlib

font = {'family' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)

def approximate_dist(k_theta_list):
    numerator = 0
    denominator = 0
    for dist in k_theta_list:
        numerator += dist[0] * dist[1]
        denominator += dist[0] * dist[1] * dist[1]
    k_sum = numerator * numerator / denominator
    theta_sum = numerator / k_sum
    return k_sum, theta_sum

def term(o,t,j):
    x=o/t
    return (x**j)/gamma(j+1)*(e**(-x))

def general_term(o, t, s, k):
    # print ('\notsk:',o,t,s,k)
    term_1 = (e ** (-o / s))
    term_2 = ((1 - t / s) ** (-k))
    term_3 = gamma(k)
    term_4 = gammainc(k, o * (1 / t - 1 / s)) #=lower part/full gamma
    # print('term1234',term_1, term_2, term_3, term_4)
    return  term_1 * term_2 * term_4

def F_MINUS_F_term(o, lamb, k, j, amp_index):
    F_jo = integrate.quad(sum_pdf, 0, amp_index * o, args=(amp_index, k, j, lamb))[0]
    F_jp1o = integrate.quad(sum_pdf, 0, amp_index * o, args=(amp_index, k, j + 1, lamb))[0]
    term = F_jo - F_jp1o
    if term < -1e-08:
       print(k,j,F_jo, F_jp1o)
    return term

def sum_pdf(x, amp_index, k , j, lamb):
    if k < j:
        pdf=(amp_index ** (k)) * (e ** (-x / lamb)) * (x ** (j - 1)) * ((lamb) ** (-j)) * hyp1f1(k, j, - (amp_index - 1) * x / (lamb)) / gamma(j)
    else:
        pdf = (e ** (- x / (lamb / amp_index))) * (x ** (j - 1)) * ((lamb / amp_index) ** (-j)) / gamma(j)
    return pdf

def general_approximate_term(o,lamb,k,j,amp_index):
    # hyper_term = (amp_index**(k))*(e**(-o/lamb))*((o/lamb)**(j))*hyp1f1(k, j + 1, - (amp_index - 1) * o / (lamb)) / gamma(j + 1)
    hyper_term = integrate.quad(integrated_function, 0, o, args=(amp_index, k, j, o, lamb))
    return hyper_term

def integrated_function(x, amp_index, k, j, z, lamb):
    return (amp_index ** (k)) * (e ** (-z / lamb)) * ((z / lamb) ** (j)) * hyp1f1(k, j, - (amp_index - 1) * z / (lamb)) / gamma(j + 1) * (e ** (- (x - z) / lamb))

def hypergeometric1F1Regularized(a,b,x):
    term = 1 / gamma(b)
    sum = term
    for i in range(75):
        term = term * gamma(b + i) / gamma(b + i + 1) * (a + i) * x
        sum += term
    return sum

def sum_to_j(o,lamb,p,j,amp_index):
    sum=0
    for k in range(1, j + 1):

        # sum_k, sum_lamb = approximate_dist([(k, lamb / amp_index),(j - k, lamb)])
        # sum_k_1, sum_lamb_1 = approximate_dist([(k, lamb / amp_index), (j - k + 1, lamb)])

        ancient_term_k=(p ** (k - 1)) * (1 - p) * term(o * amp_index, lamb, j + (amp_index - 1) * k)
        revised_ancient_term_k=(p ** (k - 1)) * (1 - p) * (term(o * amp_index, lamb, j + (amp_index - 1) * k) * (1 - gammainc(k, 1 / lamb)) + term(o * amp_index, lamb, j) * gammainc(k, 1 / lamb))
        # h_0=1 (standardized)
        # if ancient_term_k>0:
        #     print(ancient_term_k,revised_ancient_term_k,k,j)
        #
        # o_term_k=(p ** (k - 1)) * (1 - p) * term(o * amp_index, sum_lamb, sum_k)
        # o_term_k_2=(p ** (k - 1)) * (1 - p) * term(o * amp_index, sum_lamb_1, sum_k_1)
        #
        # term_k = (p ** (k - 1)) * (1 - p) * general_term(o * amp_index, sum_lamb, lamb, sum_k)
        #
        # n_term_k = (p ** (k - 1)) * (1 - p) * (general_approximate_term(o , lamb, k, j, amp_index)[0])

        # F_term_k = (p ** (k - 1)) * (1 - p) * F_MINUS_F_term(o , lamb, k, j, amp_index)
        # if term_k > 1e-04:
            # print('j,k,k_sum,lamb_sum,term_k,o_term_k,ancient_term_k,n_term_k:',
            #       j, k, sum_k, sum_lamb, term_k, o_term_k, ancient_term_k, n_term_k)

        sum+=revised_ancient_term_k
        # if term_k<1e-08:
        #     break
    return sum

def simulation_homo(o,lamb,p,amp_index,platoon_length,stochasticity_gap,stochasticity_m, CF_type):
    print('o',o,'p_al', p, 'homo')

    amplification_coef=amp_index

    initial_void=o
    a_l_penetration_rate=p

    L_c_range=[platoon_length]
    beta_range=[a_l_penetration_rate]

    X,Y = np.meshgrid(L_c_range,beta_range)

    record=[]

    platoon_number = []
    vehicle_number = []
    used_gap = []
    o_cum = []
    o_in_g=[]
    omega_g=[]

    lambda_expectation = lamb
    L_c=L_c_range[0]
    beta=beta_range[0]

    for i in range(10000):
        o_CF = 0

        o_l = initial_void
        num=0. #used gap
        j=0 #affected platoon
        num_veh=0. #affected vehicle number
        un_amp=True

        type = ran_type(beta)  # gap amplification
        if type & un_amp:
            o_l = 1 + (o_l - 1) * amplification_coef
            un_amp = False

        o_in = o_l - 1
        omega = 0

        while o_l>0:
            if stochasticity_gap == 'd':
                gap = lambda_expectation
            else:
                gap=np.random.exponential(lambda_expectation) #gap size

            # if (num==0) & (gap < h_0):
            #     continue

            if o_l>gap:#counted the used gap
                num = num + 1.
            else:
                num=num+o_l/gap

            o_l = o_l - gap #gap resolve

            if o_l>0: #if o_l still larger than 0, means it will influence the following platoon, count the platoon info
                j+=1
                if stochasticity_m == 'd':
                    this_platoon_length = L_c
                else:
                    this_platoon_length = np.random.random_integers(1, 2 * L_c - 1)
                num_veh += this_platoon_length

                type = ran_type(beta) #gap amplification
                if type&un_amp:
                    amp_part = min(o_l, o_in)
                    omega = amp_part * (amplification_coef - 1)
                    o_l = o_l + omega
                    un_amp=False
                else:
                    if CF_type == 'exp':
                        CF_coef_a_s=1.1**this_platoon_length
                        CF_coef_a_l=1.05**this_platoon_length
                        if type:  # CF effects
                            o_CF += o_l * (CF_coef_a_s - 1)
                            o_l = o_l * CF_coef_a_s

                        if not type:
                            o_CF += o_l * (CF_coef_a_l - 1)
                            o_l = o_l * CF_coef_a_l
                        if o_l>100:
                            break
                    if CF_type == 'linear':
                        CF_coef_a_s=0.1*this_platoon_length
                        CF_coef_a_l=.05*this_platoon_length
                        if type:  # CF effects
                            o_CF += CF_coef_a_s
                            o_l = o_l + CF_coef_a_s
                        if not type:
                            o_CF += CF_coef_a_l
                            o_l = o_l + CF_coef_a_l
                        if o_l>100:
                            break

        o_in_g.append(o_in)
        omega_g.append(omega)
        o_cum.append(o_in+omega+o_CF)
        used_gap.append(num)
        platoon_number.append(j)
        vehicle_number.append(num_veh)
    return (statistics.mean(o_in_g),statistics.mean(omega_g),statistics.mean(o_cum),statistics.mean(platoon_number),statistics.mean(vehicle_number))

def simulation_heter(o,lamb,p,amp_index,platoon_length,stochasticity_gap,stochasticity_m,CF_type):
    print('o',o,'p_al', p, 'homo')

    amplification_coef=amp_index

    initial_void=o
    a_l_penetration_rate=p

    L_c_range=[platoon_length]
    beta_range=[a_l_penetration_rate]

    X,Y = np.meshgrid(L_c_range,beta_range)

    record=[]

    platoon_number = []
    vehicle_number = []
    used_gap = []
    o_cum = []
    o_in_g=[]
    omega_g=[]

    lambda_expectation = lamb
    L_c=L_c_range[0]
    beta=beta_range[0]

    for i in range(10000):
        o_CF = 0
        o_l = initial_void
        num=0. #used gap
        j=0 #affected platoon
        num_veh=0. #affected vehicle number
        un_amp=True

        type = ran_type(beta)  # gap amplification
        if type & un_amp:
            o_l = 1 + (o_l - 1) * amplification_coef
            un_amp = False

        o_in = o_l - 1
        omega = 0

        while o_l>0:
            if stochasticity_gap == 'd':
                gap = lambda_expectation
            else:
                gap=np.random.exponential(lambda_expectation) #gap size

            # if (num==0) & (gap < h_0):
            #     continue

            if o_l>gap:#counted the used gap
                num = num + 1.
            else:
                num=num+o_l/gap

            o_l = o_l - gap #gap resolve

            if o_l>0: #if o_l still larger than 0, means it will influence the following platoon, count the platoon info
                j+=1

                if stochasticity_m == 'd':
                    this_platoon_length = L_c
                else:
                    this_platoon_length = np.random.random_integers(1, 2 * L_c - 1)
                num_veh += this_platoon_length

                for i in range(this_platoon_length):
                    type = ran_type(beta) #gap amplification
                    if type&un_amp:
                        amp_part = min(o_l, o_in)
                        omega = amp_part * (amplification_coef - 1)
                        o_l = o_l + omega
                        un_amp=False
                    else:
                        if CF_type == 'exp':
                            CF_coef_a_s = 1.1
                            CF_coef_a_l = 1.05
                            if type:  # CF effects
                                o_CF += o_l * (CF_coef_a_s - 1)
                                o_l = o_l * CF_coef_a_s
                            if not type:
                                o_CF += o_l * (CF_coef_a_l - 1)
                                o_l = o_l * CF_coef_a_l

                            if o_l > 100:
                                break
                        if CF_type == 'linear':
                            CF_coef_a_s = 0.1
                            CF_coef_a_l = .05
                            if type:  # CF effects
                                o_CF += CF_coef_a_s
                                o_l = o_l + CF_coef_a_s
                            if not type:
                                o_CF += CF_coef_a_l
                                o_l = o_l + CF_coef_a_l
                            if o_l > 100:
                                break

        o_in_g.append(o_in)
        omega_g.append(omega)
        o_cum.append(o_in + omega + o_CF)
        used_gap.append(num)
        platoon_number.append(j)
        vehicle_number.append(num_veh)

    return (statistics.mean(o_in_g),statistics.mean(omega_g),statistics.mean(o_cum),statistics.mean(platoon_number),statistics.mean(vehicle_number))

def ran_type_multiple(a_group,p_group):
    cdf=[sum(p_group[:i+1]) for i in range(len(p_group))]

    type_random = random.random()
    for i in range(len(cdf)):
        if type_random <= cdf[i]:
            return a_group[i]

def multi_value_simulation(initial_void,lamb,a_group,p_group,platoon_length,scenario):
    print(a_group,p_group,platoon_length, scenario)

    L_c_range=[platoon_length]
    platoon_number = []
    vehicle_number = []
    used_gap = []
    o_cum = []
    o_in_g=[]
    omega_g=[]

    lambda_expectation = lamb
    L_c=L_c_range[0]

    for i in range(10000):
        current_a = max(a_group)
        o_l = initial_void
        num=0. #used gap
        j=0 #affected platoon
        num_veh=0. #affected vehicle number

        follower_a = ran_type_multiple(a_group, p_group)  # initial gap amplification (LC vehicle)
        if follower_a < current_a:
            o_l = 1 + (o_l - 1) * (current_a / follower_a) #h_0 + void
            current_a = follower_a

        o_in = o_l - 1 # without the h_0 part, only this part will be amplified
        omega = 0

        while o_l>0:
            gap=np.random.exponential(lambda_expectation) #gap size

            if o_l>gap:#counted the used gap
                num = num + 1.
            else:
                num=num+o_l/gap

            o_l = o_l - gap #gap resolve

            if o_l>0: #if o_l still larger than 0, means it will influence the following platoon, count the platoon info
                j+=1

                this_platoon_length=np.random.random_integers(1,2*L_c-1)
                num_veh+=this_platoon_length

                if scenario == 'heter':
                    this_platoon_type=this_platoon_length
                if scenario == 'homo':
                    this_platoon_type=1

                for i in range(this_platoon_type):
                    follower_a = ran_type_multiple(a_group, p_group)  # gap amplification
                    if follower_a < current_a:
                        omega = min(o_l, o_in) * (current_a / follower_a - 1) #only the void part is amplified, if h_0 is not fully resolved yet
                        o_l = o_l + omega
                        o_in = o_in * (current_a / follower_a) # basic magnitude of the void part
                        current_a = follower_a

        o_in_g.append(o_in)
        omega_g.append(omega)
        o_cum.append(o_in+omega)
        used_gap.append(num)
        platoon_number.append(j)
        vehicle_number.append(num_veh)

    return (statistics.mean(o_in_g),statistics.mean(omega_g),statistics.mean(o_cum))

def not_amplified(o,lamb):
    x=1
    j=1
    sum=0
    sum_p=0
    while j<75:
        prob=term(o,lamb,j)

        # print('j:',j,'probability:',round(prob,3))
        x=prob*j
        j+=1
        sum+=x
        sum_p+=prob

    # print('j count to:', j, 'probability count to:', sum_p)
    return sum

def amplified(o,lamb,p,amp_index):
    x=0
    j=1
    sum=0
    sum_p=0
    while j<75:
        prob=term(o,lamb,j)*(p**j)+sum_to_j(o,lamb,p,j,amp_index)
        # print('j:',j,'probability:',round(prob,3))
        x=prob*j
        # print(j , prob)
        j+=1
        sum+=x
        sum_p+=prob

    # print('j count to:', j, 'probability count to:', sum_p)
    return sum

def tri_amplified(o,lamb,p,r,a_l,a_m,a_s):
    x = 0
    j = 1
    sum = 0
    sum_p = 0
    while j < 75:
        prob = term(o, lamb, j) * (p ** j)
        for k in range(1, j + 1):
            prob += (p ** (k - 1)) * (r ** (j-k+1)) * term(o * a_l / a_m, lamb, j + (a_l / a_m - 1) * k)
            prob += (p ** (k - 1)) * (1 - p - r) * term(o * a_l / a_s, lamb, j + (a_l / a_s - 1) * k)
        for k in range(1, j):
            for l in range(1 , k + 1):
                prob += (p ** (k - l)) * (r ** l) * (1 - p - r) * term(o * a_l / a_s, lamb, j + (a_l / a_s - 1) * k + (a_m / a_s - 1))
        # print(j , prob)
        x = prob * j

        j += 1
        sum += x
        sum_p += prob

    # print('j count to:', j, 'probability count to:', sum_p)
    return sum

def drawing_sensitivity():

    color=['r','g','b','y','c','k']
    frequency=.5

    line=0
    o=(30-10)**2/2/30/3/3600*2200+1
    p_l = 4

    homo_chain = []
    heter_chain = []
    a_l = 3.
    a_s = 1.
    tested_probability = np.arange(0.2, 1.01, 0.3)
    for a_l_penetration_rate in tested_probability:
        homo_chain.append([])
        heter_chain.append([])

        for p_l in range(1,11):
            # expectation_chain.append(not_amplified(o,.25*p_l))
            homo_chain[line].append(simulation_homo(o, .25*p_l, p=a_l_penetration_rate, amp_index=a_l / a_s,
                                                    platoon_length = p_l,stochasticity_gap='r',stochasticity_m='d', CF_type='no'))
            heter_chain[line].append(simulation_heter(o, .25*p_l, p=a_l_penetration_rate, amp_index=a_l / a_s,
                                                      platoon_length = p_l,stochasticity_gap='r',stochasticity_m='d', CF_type='no'))
        line += 1

    label_1='homo'
    label_2='heter'
    a_l_penetration_rate_group = tested_probability
    fig = plt.figure(figsize=(6, 6), dpi=100, tight_layout=True)
    bx = fig.add_subplot(111)
    # bx.plot(np.arange(1, 4.1, frequency),expectation_chain,c=color[line],label='analytical ' + str(a_l_penetration_rate))
    for i in range(line):
        bx.plot(range(1,11),[s[3] for s in homo_chain[i]],'-o', c=color[i], label=label_1 + ' ' + str(a_l_penetration_rate_group[i]))
        bx.plot(range(1,11),[s[3] for s in heter_chain[i]],'--*', c=color[i], label=label_2 + ' ' + str(a_l_penetration_rate_group[i]))

    # bx.axis([1, 4, 1, 4])
    plt.title('influenced platoon ($j$)')
    plt.legend()
    plt.xlabel('m')
    plt.savefig('j.png')
    # plt.show()


    fig = plt.figure(figsize=(6, 6), dpi=100, tight_layout=True)
    bx = fig.add_subplot(111)
    # bx.plot(np.arange(1, 4.1, frequency),expectation_chain,c=color[line],label='analytical ' + str(a_l_penetration_rate))
    for i in range(line):
        bx.plot(range(1, 11), [s[2] for s in homo_chain[i]], '-o', c=color[i],
                label=label_1 + ' ' + str(a_l_penetration_rate_group[i]))
        bx.plot(range(1, 11), [s[2] for s in heter_chain[i]], '--*', c=color[i],
                label=label_2 + ' ' + str(a_l_penetration_rate_group[i]))

    # bx.axis([1, 4, 1, 4])
    plt.title('cumulative void ($o_{cum}$)')
    plt.legend()
    plt.xlabel('m')
    plt.savefig('o_cum.png')
    # plt.show()

    fig = plt.figure(figsize=(6, 6), dpi=100, tight_layout=True)
    bx = fig.add_subplot(111)
    # bx.plot(np.arange(1, 4.1, frequency),expectation_chain,c=color[line],label='analytical ' + str(a_l_penetration_rate))
    for i in range(line):
        bx.plot(range(1, 11), [s[4] for s in homo_chain[i]], '-o', c=color[i],
                label=label_1 + ' ' + str(a_l_penetration_rate_group[i]))
        bx.plot(range(1, 11), [s[4] for s in heter_chain[i]], '--*', c=color[i],
                label=label_2 + ' ' + str(a_l_penetration_rate_group[i]))

    # bx.axis([1, 4, 1, 4])
    plt.title('influenced vehicle ($\pi$)')
    plt.legend()
    plt.xlabel('m')
    plt.savefig('pi.png')
    # plt.show()


# for j in range(75):
#     for k in range(10):
#         print(j, k)
#         print(approximate_dist([(k, 2.5 / 2), (j , 2.5)]))
#         print(approximate_dist([(k, 2.5 / 2), (j + 1 , 2.5)]))


# print('non-amplified: ',not_amplified(40/9+1,2.5))
# print(amplified(40/9+1, 2.5, p=0.5, amp_index=2))
# print(tri_amplified(40/9+1, 2.5, p=0.33, r=0.33, a_l=2, a_m=1.5, a_s=1))
#
drawing_sensitivity()

# print(term(40 / 9 + 1, 2.5, 10))
# print(general_term(40 / 9 + 1, 2.45 , 2.5, 10))
