
import sympy as sp 
"""Define algebraic symbols"""
t = sp.symbols("t")
tau = sp.symbols('tau')
r= sp.symbols("r", positive = True)
theta = sp.symbols("theta") 
phi = sp.symbols("phi")
x = sp.symbols("x")
y = sp.symbols("y")
z = sp.symbols("z")
u = sp.symbols("u")
v = sp.symbols("v")
U = sp.symbols("U")
V = sp.symbols("V")


"""Define function version of symbols for use in calculation of Geodesics"""
functions = {
t : sp.Function('t')(tau),
x : sp.Function('x')(tau),
y : sp.Function('y')(tau),
z : sp.Function('z')(tau),
r : sp.Function('r')(tau),
theta : sp.Function('theta')(tau),
phi : sp.Function('phi')(tau),
u : sp.Function('u')(tau),
v : sp.Function('v')(tau),
U : sp.Function('U')(tau),
V : sp.Function('V')(tau)
}

"""Define indices for use in tensor (and non tensor) calculations"""
mu = sp.symbols("mu")
nu = sp.symbols("nu")
sigma = sp.symbols("sigma")
rho = sp.symbols("rho")
epsilon = sp.symbols("epsilon")
alpha = sp.symbols("alpha")

"""Define Physical Qauntities """
c = sp.symbols('c')
G = sp.symbols("G")
M = sp.symbols("M", positive = True)
J = sp.symbols("J")
A = sp.symbols("A")
Lambda = sp.symbols('Lambda')
a = sp.Function('a')(t)
A = J/M
Sigma = sp.symbols('Sigma')
Delta = sp.symbols('Delta')


f = 1 - (Lambda*r**2)/3

"""Define Local symbols for use in inputs, so string characters can be understood as sympy symbols"""
local_symbols = {
    't' : t,
    'r': r,
    'theta' : theta,
    'phi': phi,
    'x' : x,
    'y' : y,
    'z' : z,
    'tau':'tau',
    'sin' : sp.sin,
    'cos' : sp.cos,
    'G' : G,
    'M' : M,
    'exp' : sp.exp
}

"""Dict that defines different scale factors for use in the FLRW metric and allows for manual input"""
scale_factor = {
    'Manual Entry' : ' ',
    'Ambiguous Scale Function a(t)': a,
    'Matter Dominated' : t**sp.Rational(2,3),
    'Radiation Dominated' : t**sp.Rational(1,2),
    'Curvature Dominated' : t
}

"""Define different kinds of coordinates for use in metric calculations, if a manual metric is chosen, a user can choose, or manually enter coordinates"""
coordinates_dict = {
    'Manual' : '[ , , , ]',
    'Cartesian' : [t, x, y, z],
    'Spherical' : [t, r, theta, phi],
    'Cylindrical' : [t, r, phi, z],
    'Null (Cartesian)' : [u, v, y, z],
    'Null (Kruskal)' : [U, V, theta, phi],
}

Metrics = {
    'Manual Entry' : {'Metric' : sp.zeros(4),
    'coordinates':[]
    },
    'Minkowksi Cartesian' : {'Metric' : sp.Matrix([
        [-1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1]
    ]),
    'coordinates' : coordinates_dict['Cartesian']
    },
    'Minkowski Spherical' : { 'Metric' : sp.Matrix([
        [-1,0,0,0],
        [0,1,0,0],
        [0,0,r**2,0],
        [0,0,0,r**2 * sp.sin(theta)**2]
        
    ]),
    'coordinates': coordinates_dict['Spherical']
    },
    'Schwarzchild':{'Metric': sp.Matrix([
        [-(1-2*G*M/r),0,0,0],
        [0,(1-2*G*M/r)**-1,0,0],
        [0,0,r**2,0],
        [0,0,0,r**2*sp.sin(theta)**2]
    ]),
    'coordinates': coordinates_dict['Spherical']
    },
    'Kerr':{'Metric':sp.Matrix([
        [-(1-(2*M*r)/Sigma),0,0,-(2*M*(A)*r*sp.sin(theta)**2)/Sigma],
        [0,Sigma/Delta,0,0],
        [0,0,Sigma,0],
        [-(2*M*(A)*r*sp.sin(theta)**2)/Sigma,0,0,sp.sin(theta)**2*(r**2+(A)**2)+(2*M*(A)**2*r*sp.sin(theta)**2)/Sigma]]),
    'coordinates': coordinates_dict['Spherical']
    },
    'FLRW':{'Metric': sp.Matrix([
        [-1,0,0,0],
        [0,a**2,0,0],
        [0,0,a**2,0],
        [0,0,0,a**2]
    ]),
    'coordinates': coordinates_dict['Cartesian']
    },

    'de Sitter Spherical': {'Metric': sp.Matrix([
        [-(f),0,0,0],
        [0,f**-1,0,0],
        [0,0,r**2,0],
        [0,0,0,r**2 * sp.sin(theta)**2]
    ]),
    'coordinates':coordinates_dict['Spherical']
    }
}

"""Lists metrics from dictionary for user, allows them to choose a number that picks a metric and coordinate system to work with"""
for i, name in enumerate(Metrics, start = 0):
     print(f'{i}. {name}')
choice = int(input('Enter number to choose Metric: '))
metric_name = list(Metrics.keys())[choice]
metric = Metrics[metric_name]['Metric']
coordinates = Metrics[metric_name]['coordinates']
"""Starts Branch for manual entry of metric tensor, lists choices of coordinates"""
if metric_name == 'Manual Entry':
    coordinates = []
    for i, (name,value) in enumerate(coordinates_dict.items(), start = 0):
       print(f'{i}. {name} : {value}') 
    coordinate_choice = int(input('Choose coordinates (press "0" to manually define): '))
    """Manual entry of coordinates"""
    if coordinate_choice == 0:
        coordinates = (input('Enter coordinates seperated by spaces: ').split())
        coordinates = [sp.sympify(entry, locals = local_symbols) 
                       for entry in coordinates]
    else:
        coordinates = list(coordinates_dict.values())[coordinate_choice]
    
    """Prompts user to make entering metric easier, if manual metric is diagonal, only 4 entries are needed"""
    diagonal_list = []
    new_diag_list = []
    yes_or_no = input("Diagonal metric, y/n?")
    if yes_or_no == "y":
        diagonal_list = input('Enter g_00, g_11, g_22, g_33 seperated by spaces: ').split()
        if len(diagonal_list) != 4:
            raise ValueError('Error: you must enter exactly 4 components!')
        for entry in diagonal_list:
            new_diag_list.append(sp.sympify(entry, locals = local_symbols))
        for mu in range (4):
            metric[mu,mu] = new_diag_list[mu] 

    """allows User to manually define their own metric, one with 10 independnet 
    components, uses metric symmetry to make it easier for the user to enter componennts."""
    if yes_or_no == "n":
        print('Enter upper right triangle of metric:' )
        Upper_row_list = input('Enter g_00, g_01, g_02, g_03 seperated by spaces: ' ).split()
        if len(Upper_row_list) != 4:
            raise ValueError('You must enter 4 values!')
        Second_row_list = input('Enter g_11, g_12, g_13 seperated by spaces: ').split()
        if len(Second_row_list) != 3:
            raise ValueError('You must enter 3 values!')
        Third_row_list = input('Enter g_22, 23 seperated by spaces: ').split()
        if len(Third_row_list) !=2:
            raise ValueError('You must enter 2 values!')
        Final_row_list = input('Enter g_33: ')
        if len(Final_row_list) != 1:
            raise ValueError('You must enter 1 value!')
        Upper_row_list = [
            sp.sympify(entry, locals = local_symbols)
            for entry in Upper_row_list]
        Second_row_list = [sp.sympify(entry, locals = local_symbols) 
                           for entry in Second_row_list]
        Third_row_list = [sp.sympify(entry, locals = local_symbols) 
                           for entry in Third_row_list]
        Final_row_list = [sp.sympify(entry, locals = local_symbols)
                          for entry in Final_row_list]
        for mu in range (4):
            metric[0,mu] = Upper_row_list[mu]
        for mu in range(3):
            metric[1,mu+1] = Second_row_list[mu]
        for mu in range (2):
            metric[2,mu+2] = Third_row_list[mu]
        metric[3,3] = Final_row_list
        for mu in range(3):
            metric[mu+1,0] = metric[0,mu+1]
        for mu in range(2):
            metric[mu+2,1] = metric[1,mu+2]
        metric[3,2] = metric[2,3]

"""Special case for FLRW metric, allows choice of scale factor from scale factor dictionary"""
if metric_name == 'FLRW':
    for i, (name,value) in enumerate(scale_factor.items(), start = 0):
        print(f'{i}. {name} : {value}') 
    scale_choice = int(input('Choose scale factor (press "0" for manual entry): '))

    if scale_choice == 0:
        manual_entry = input('Input scale factor of choice: ')
        s_s_f = sp.sympify(manual_entry)
    else:
        s_s_f = list(scale_factor.values())[scale_choice]
    metric = metric.subs(a, s_s_f)

if metric_name == 'Kerr':
    input('Press "f" to pay respects for your computer')

"""Start of calculation"""
g = metric
g_inv = g.inv()
sp.pprint(g)
gamma = {}
"""Find Connection Coefecients"""

"""FIX CONNECTION COEFFECIENTS USING SYMBOLS AND NOT FUNCTIONS WHEN IN USE BY GEODESIC FUNC"""
for sigma in range(4):
    for mu in range(4):
        for nu in range (4):
            gamma_x_two = 0
            for rho in range (4):
                gamma_x_two += g_inv[sigma, rho] * (
                sp.diff(g[rho, nu], coordinates[mu]) 
                + sp.diff(g[rho, mu], coordinates[nu]) 
                - sp.diff(g[mu, nu], coordinates[rho]))
            gamma[sigma, mu, nu] = sp.simplify(gamma_x_two/2)
for indices, value in gamma.items():
    if value != 0:
        sigma, mu, nu = indices
        print(f"Gamma^{sigma}_{mu}{nu}=", value)

"""Calculate Geodesics calls on dictionary coordinate func, which maps sympy symbols to sympy functions of the same name for use in differentiaon with repsect to proper time"""
coordinates_func = [functions[nu] for nu in coordinates] 

d2xdtau2 = {}
for mu in range(4):
    geodesic = 0
    for sigma in range (4):
        for rho in range(4):
            geodesic += gamma[mu,sigma,rho] * sp.diff(coordinates_func[sigma], tau) * sp.diff(coordinates_func[rho], tau)
        d2xdtau2[mu] = -(sp.simplify(geodesic))
for index, value in d2xdtau2.items():
    if value !=0:
        mu = index
        print(f'd2x^{mu}/dtau2 =', value)           


"""Find Reimann Tesnsor Using Connection Coeffecients, calculated above"""
R = {}
for mu in range(4):
    for nu in range(4):
        for sigma in range(4):
            for rho in range(4):
                
                Reimann = (sp.diff(gamma[rho,mu,sigma],coordinates[nu]) 
                - sp.diff(gamma[rho,mu,nu],coordinates[sigma]))
                
                for alpha in range(4):
                    Reimann += (gamma[alpha,mu,sigma] * gamma[rho,nu,alpha]
                    - gamma[alpha,mu,nu] * gamma[rho,sigma,alpha])
                R[rho,mu,nu,sigma] = sp.simplify(Reimann)
for index, value in R.items():
    if value !=0:
        rho,mu,nu,sigma = index
        print(f'R^{rho}_{mu}{nu}{sigma} =', value)

"""Find Ricci Tensor Using Connection Coeffecients, Same as Reimann tesnor except first and third indicies are contracted"""
Ricci = {}
for mu in range(4):
    for nu in range(4):
        Ric = 0
        for sigma in range(4):
            Ric += (sp.diff(gamma[sigma,mu,nu],coordinates[sigma]) 
            - sp.diff(gamma[sigma,mu,sigma],coordinates[nu]))
            for rho in range(4):  
                Ric += ((gamma[rho,mu,nu] * gamma[sigma,sigma,rho])
                -(gamma[rho,mu,sigma] * gamma[sigma,nu,rho]))
        Ricci[mu,nu] = sp.simplify(Ric)
for index, value in Ricci.items():
    if value != 0:
        mu,nu = index
        print(f'R_{mu}_{nu} =', value)

"""Calculate Ricci Scalar"""
Ricci_Scalar = 0
for mu in range (4):
    for nu in range(4):
        Ricci_Scalar += g_inv[mu,nu] * Ricci[mu,nu]
print(f'Ricci Scalar =', sp.simplify(Ricci_Scalar))

"""Calulate Einstein Tensor"""
E={}
for mu in range(4):
    for nu in range(4):
        E[mu,nu] = Ricci[mu,nu] - sp.Rational(1,2) * Ricci_Scalar * g[mu,nu]
for index, value in E.items():
    if value !=0:
        mu,nu = index
        print(f'G_{mu}_{nu} = ', value)

"""Calculate Stess Energy Tensor"""
T = {}
for mu in range(4):
    for nu in range(4):
        T[mu,nu] = c**4/(8*sp.pi*G) * E[mu,nu]
for index, value in T.items():
    if value !=0:
        mu,nu = index
        print(f'T_{mu}_{nu} = ', value)