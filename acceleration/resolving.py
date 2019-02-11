import random
import statistics
import numpy as np

def ran_type(a_l_penetration):
    type_random = random.random()
    if type_random >= a_l_penetration:
        type = True  # 1 indicates a_s
    else:
        type = False
    return type

def ran_type_multiple(a_l_penetration, m):
    vehicle_types=[]
    for i in range(m):
        vehicle_types.append(ran_type(a_l_penetration))
    return vehicle_types

def ran_n_acc_type(rate_group,m):
        s = sum(rate_group)
        rate_group = [r / s for r in rate_group]
        rate_group_order = [sum(rate_group[:i+1]) for i in range(len(rate_group))]
        rate_group_order[-1] = 1
        type_random = random.random()
        for i in range(len(rate_group_order)):
            if rate_group_order[i]>=type_random:
                return i

def Approximate_solution(amplification_coef,initial_void,q,h_0,a_s_penetration_rate,N_p_c,beta):
    if (a_s_penetration_rate==0)|(a_s_penetration_rate==1):
        affected_vehicle=initial_void*q/(1-q*h_0)-1/2*N_p_c
    else:
        affected_vehicle=amplification_coef*initial_void*q/(1-q*h_0)-1/a_s_penetration_rate*N_p_c\
                                                                 *((amplification_coef-3./2.)*beta-amplification_coef+1)
    return affected_vehicle

def main():
    amplification_coef=2#a_l/a_s
    # amplification_coef=np.arange(1,2.1,0.1)

    CF_coef_a_s=1.00
    CF_coef_a_l=1.00
    h_0=1

    initial_void=(40/9.)+h_0
    # initial_void=0.25



    q_range=[4/5] #ratio
    # q=[i/3600. for i in q]#transfer to vehicle/sec, corresponding to sec scale

    #independent
    a_l_penetration_rate=0.1


    L_c_range=[10]
    beta_range=[a_l_penetration_rate]#a_l_platoon_ratio
    # N_p_c_range=np.arange(4,8,4)
    # beta_range=np.arange(0,0.05,0.1)

    X,Y = np.meshgrid(L_c_range,beta_range)

    record=[]

    for index_p in range(len(L_c_range)):
        for index_b in range(len(beta_range)):
            for q in q_range:
                platoon_number = []
                vehicle_number = []
                used_gap = []

                L_c=L_c_range[index_p]
                beta=beta_range[index_b]
                #dependent variable
                if (a_l_penetration_rate==0)|(a_l_penetration_rate==1):
                    lambda_expectation=round((1-h_0*q)/(q)*L_c,1)
                    L_r=L_c
                else:
                    lambda_expectation=(1-h_0*q)/(q*a_l_penetration_rate)*beta*L_c
                    L_r=(1-a_l_penetration_rate)/a_l_penetration_rate*beta/(1-beta)*L_c

                for i in range(100000):
                    o_l = initial_void
                    num=0. #used gap
                    j=0 #affected platoon
                    num_veh=0. #affected vehicle number
                    un_amp=True

                    # type=ran_type(beta) #general term of generated void, no need to define O_r/O_c
                    # if type & un_amp:
                    #     o_l = o_l * amplification_coef
                    #     un_amp = False

                    while o_l>0:
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

                            type = ran_type(beta) #gap amplification
                            if type&un_amp:
                                o_l=o_l*amplification_coef
                                un_amp=False

                            if type:
                                num_veh+=random.randint(1,2*L_c-1)
                            else:
                                num_veh+=random.randint(1,2*L_r-1)

                            if type&(not un_amp): #CF effects
                                o_l=o_l*CF_coef_a_s
                            if not type:
                                o_l=o_l*CF_coef_a_l

                    used_gap.append(num)
                    platoon_number.append(j)
                    vehicle_number.append(num_veh)

                # print('void size:',initial_void)
                # print('platoon len:',N_p_c)
                # print('beta:', beta)
                # print('used gap:',statistics.mean(used_gap))
                # print('lambda expectation:', lambda_expectation)
                print('affected platoon:',statistics.mean(platoon_number))
                # print('affected vehicle:',statistics.mean(vehicle_number))
                # print('estimated affected vehicle:', Approximate_solution(amplification_coef,initial_void,q,h_0,a_s_penetration_rate,N_p_c,beta))
                # print('vehicle number variance:',statistics.variance(vehicle_number))
                record.append((beta,L_c,q*3600,statistics.mean(vehicle_number),lambda_expectation))


    print('p_slow','exp_length','flow','sim','lam_exp')
    for rec in record:
        print(rec[0],rec[1],rec[2],rec[3],rec[4])

if __name__ == '__main__':
    main()