# PNJL 常数

# 热力学常数

hc=197.33

rho0=0.16

T0 = 210 / hc 

# PNJL 模型常数
Nc=3.0
Nf=3.0

Lambda=602.3    # !(MeV)

G_Lam2=1.835

K_Lam5=12.36

m0_q=5.5       # !(MeV)
m0_s=140.7    # !(MeV)



Lambda_f=Lambda/hc  #    !(fm**(-1))

G_f=G_Lam2/Lambda_f^2   # !(fm**(2))


K_f=K_Lam5/Lambda_f^5   #  !(fm**(5))

m0_q_f=m0_q/hc       #    !(fm**(-1))
m0_s_f=m0_s/hc

m0 = [m0_q_f, m0_q_f, m0_s_f] # 3×1 夸克质量数组


# Polyakov-loop 势能常数

a0 = 3.51
a1 = -2.47
a2 = 15.2

b3 = -1.75
b4=7.555
