# PNJL 常数

# 热力学常数
pi=3.141592653
hc=197.33

rho0=0.16

T0 = 210 / hc 


mixmix=1 

# PNJL 模型常数
Nc=3.0
Nf=3.0

Lambda=602.3    # !(MeV)

Lambda_QCD = 300/hc

G_Lam2=1.835

K_Lam5=12.36

m0_q=5.5       # !(MeV)
m0_s=140.7



Lambda_f=Lambda/hc  #    !(fm**(-1))

G_f=G_Lam2/Lambda_f^2   # !(fm**(2))


K_f=K_Lam5/Lambda_f^5   #  !(fm**(5))

m0_q_f=m0_q/hc       #    !(fm**(-1))
m0_s_f=m0_s/hc

m0 = [m0_q_f, m0_q_f, m0_s_f] # 真实流夸克质量
alpha_D = [1e10, 1e10, 1e10]  # Dirichlet boundary condition
alpha_N = [1e-3, 1e-3, 1e-3] # Neumann boundary condition

alpha_Dr = [-1e10, -1e10, -1e10]  #rev Dirichlet boundary condition

qf = [2/3, 1/3, 1/3] #charge of quarks
# emm常数

a0 = 3.51
a1 = -2.47
a2 = 15.2
b3 = -1.75



a = 0.0108805
b = -1.0133e-4
c = 0.02228
d = 1.84558e-4


