import matplotlib.pyplot as plt
from modsim import *

def make_system( alpha, beta, gamma, delta, epsilon, theta, zeta, eta, mu, nu, tau, l, kappa, xi, rho,
sigma):
    init = State(S=385600, I=2410, D=0, A=0, R=0, T=0, H=0, E=0)
    init /= sum(init)
    t0 = 0
    t_end = 7 * 36
    return System(init=init, t0=t0, t_end=t_end,alpha = alpha,beta = beta,gamma = gamma,delta = delta,epsilon = epsilon,theta = theta,
                  zeta = zeta,eta = eta,mu = mu,nu = nu,tau = tau,l = l,kappa = kappa,xi = xi,rho = rho,sigma = sigma)

tc = 9
tr = 10
alpha, beta, gamma, delta = .4, .3, .2, .1
epsilon, theta = .2,.7
zeta, eta = 1/3, 1/7
mu, nu = 1/14,1/7
tau = .02
l, kappa, xi, rho, sigma = .8, .5, .5, .4, .3
system = make_system(alpha, beta, gamma, delta, epsilon, theta, zeta, eta, mu, nu, tau, l, kappa, xi, rho, sigma)

def update_func(state, T, system):
    constants  = ["alpha", "beta", "gamma", "delta", "epsilon", "theta", "zeta", "eta", "mu", "nu", "tau", "l", "kappa", "xi", "rho", "sigma"]
    (alpha, beta, gamma, delta, epsilon, theta, zeta, eta, mu, nu, tau, l, kappa, xi, rho, sigma) = [system[key] for key in constants]
    S,I,D,A,R,T,H,E = state
    infected = S*(alpha*I+beta*D+gamma*A*delta*R) - (epsilon + zeta+l )*I
    diagnosed = epsilon*I-(eta+rho)*D
    ailing = zeta*I-(theta+mu+kappa)*A
    recognized = eta*D+theta*A-(nu+xi)*R
    threatened = mu*A*nu*R-(sigma+tau)*T
    healed = l*I+rho*D+kappa*A+xi*R+sigma*T
    extinct = tau*T
    S -= infected
    I += infected
    D += diagnosed
    A += ailing 
    R += recognized 
    T += threatened 
    H += healed
    E += extinct
    return State(S=S,I=I,D=D,A=A,R=R,T=T,H=H,E=E)

init = State(S=385600, I=2410, D=0, A=0, R=0, T=0, H=0, E=0)
init /= sum(init)
state = update_func(init, 0, system)

def run_simulation(system, update_func):
    frame = TimeFrame(columns=system.init.index)
    frame.row[system.t0] = system.init
    state = system.init
    t0 = system.t0
    for t in linrange(system.t0, system.t_end): frame.row[t+1] = update_func(frame.row[t], t, system)
    return frame

tc = 9
tr = 10
alpha = .03*3
beta = .005*2
gamma = .005*3
delta = .005*3
epsilon = .02/7
theta = .1*1/7
zeta = .1/7 
eta = .1/7
mu = 1/7*.02
nu = .05/7
tau = .02
l,kappa,xi,rho,sigma = np.array([.3,.25,.2,.2,.05])*.05*1/14
system = make_system(alpha, beta, gamma, delta, epsilon, theta, zeta, eta, mu, nu, tau, l, kappa, xi, rho, sigma)
results = run_simulation(system, update_func)
results.head()

def plot_results(S, I, D, A, R, T, H,E):
    plot(S, 'yellow', label='Susceptible')
    plot(I,'red', label = "Infected")
    plot(D,'magenta', label = "Diagnosed")
    plot(A,'purple', label = "Ailing")
    plot(R,'cyan', label = "Recognized")
    plot(T,'orange', label = "Threatened")
    plot(H,'green', label = "Healed")
    plot(E,'black', label = "Extinct")
    decorate(xlabel='Time (days)',
             ylabel='Population (1.0 = 100%)')
    
plot_results(results.S, results.I, results.D, results.A, results.R, results.T, results.H,results.E)