from random import randint
import statistics
import matplotlib.pyplot as plt

def ran_initial_vehicle(n_l,n_s):
    veh=randint(1,n_l+n_s)
    return veh

def veh_a(n_l,n_s,num,a_l,a_s):
    num=(num-1)%(n_l+n_s)
    if num<n_l:
        return a_l
    return a_s

def simple_simulation(initial_void_l,lambda_expectation,a_l,a_s,n_l,n_s,platoon_length,scenario):
    print('p_l=',n_l,'/',n_s+n_l)

    platoon_number = []
    vehicle_number = []
    used_gap = []
    o_cum = []
    o_initial_g=[]
    omega_g=[]

    for i in range(10000):
        current_a = a_l
        o_l = initial_void_l #void in unit h_0
        num=0. #used gap
        j=0 #affected platoon
        num_veh=0. #affected vehicle number

        LC_veh=ran_initial_vehicle(n_l,n_s)
        initial_a=veh_a(n_l,n_s,LC_veh,a_l,a_s)

        if initial_a < current_a:
            o_l = 1 + (o_l - 1) * (current_a / initial_a) #h_0 + void
            current_a = initial_a

        o_initial=o_l-1
        o_in = o_l - 1 # without the h_0 part, only this part will be amplified
        omega = 0

        current_veh=ran_initial_vehicle(n_l,n_s)

        while o_l>0:
            gap=lambda_expectation

            if o_l>gap:#counted the used gap
                num = num + 1.
            else:
                num=num+o_l/gap

            o_l = o_l - gap #gap resolve

            if o_l>0: #if o_l still larger than 0, means it will influence the following platoon, count the platoon info
                j+=1

                this_platoon_length=platoon_length
                num_veh+=this_platoon_length

                if scenario == 'heter':
                    this_platoon_type=this_platoon_length
                if scenario == 'homo':
                    this_platoon_type=1

                for i in range(this_platoon_type):
                    follower_a=veh_a(n_l,n_s,current_veh,a_l,a_s)

                    if follower_a < current_a:
                        omega = min(o_l, o_in) * (current_a / follower_a - 1) #only the void part is amplified, if h_0 is not fully resolved yet
                        o_l = o_l + omega
                        o_in = o_in * (current_a / follower_a) # basic magnitude of the void part
                        current_a = follower_a
                current_veh+=1

        o_initial_g.append(o_initial)
        omega_g.append(omega)
        o_cum.append(o_initial+omega)
        used_gap.append(num)
        platoon_number.append(j)
        vehicle_number.append(num_veh)
    print(statistics.mean(o_initial_g),statistics.mean(omega_g),statistics.mean(o_cum))
    return statistics.mean(o_cum)

def line(n_s_set,n_l_set):
    r_ls=(a_l/a_s-1)
    h_0=3600/capacity
    lamb=(1/s-1)*platoon_length*h_0
    o_in_s=(30-10)**2/2/30/a_s
    o_in_l=(30-10)**2/2/30/a_l
    print(o_in_s,o_in_l)
    x=[]
    y=[]
    for i in range(len(n_s_set)):
        veh_period=n_s_set[i]+n_l_set[i]
        p_l=n_l_set[i]/veh_period
        p_s=1-p_l
        x.append(p_l)
        y.append((o_in_s+p_l*(h_0*r_ls-lamb*r_ls*(p_s+p_l/2+3/2/veh_period)))/h_0)
    print(x,y)
    return x,y

def draw_figure(n_s_set,n_l_set,o_cum):
    fig = plt.figure(figsize=(6, 6), dpi=100, tight_layout=True)
    bx = fig.add_subplot(111)
    for i in range(len(n_s_set)):
        bx.scatter(n_l_set[i]/(n_l_set[i]+n_s_set[i]), o_cum[i])
    x,y=line(n_s_set,n_l_set)
    bx.plot(x,y)
    plt.xlabel('p_l')
    plt.savefig('o_cum.png')
    # plt.show()

def main():
    global a_l,a_s,s,platoon_length,capacity
    capacity=2200
    a_l=1
    a_s=0.5
    s=0.8
    platoon_length=8

    o_cum=[]
    n_s_set=[]
    n_l_set=[]
    n_l=1
    for n_s in range(9,0,-1):
        o_cum.append(simple_simulation((30-10)**2/2/30/a_l/3600*capacity+1,(1/s-1)*platoon_length,a_l,a_s,n_l,n_s,4,'homo'))
        n_s_set.append(n_s)
        n_l_set.append(n_l)
    n_s = 1
    for n_l in range(2, 10, 1):
        o_cum.append(
            simple_simulation((30 - 10) ** 2 / 2 / 30 / a_l / 3600 * capacity + 1, (1/s-1) * platoon_length, a_l, a_s, n_l,
                              n_s, 4, 'homo'))
        n_s_set.append(n_s)
        n_l_set.append(n_l)
    draw_figure(n_s_set,n_l_set,o_cum)


main()