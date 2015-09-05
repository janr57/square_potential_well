# square_potential_well.sage
# Bound states of a particle in a square potential well
# (C) 2015 José A. Navarro Ramón <josea.navarro@murciaeduca.es>

# ############################################################################
# ###################### VARIABLES GLOBALES ##################################
# ############################################################################
NODEBUG = false
DEBUG = true

# ############################################################################
# CALCULATIONS PRECISION
# 53-bit is a standard precision though we can change it.
# ============================================================================
precision_bits = 53
# ============================================================================
# RealNumber field precision bits declaration
R = RealField(precision_bits)
RealNumber = R

# ############################################################################
# A particle in a square potential well:
# mr = R(1.0) # Reduced mass of the particle;
#               Relative mass of the particle with respect to that of an
#               electron (dimensionless); mr = particle mass/electron mass.
# u0ev = R(10) # Energía potencial del pozo (eV)
# lnm = R(1) # Anchura del pozo (nm)
# ============================================================================
mr = R(1.0) # Reduced mass of the particle (dimensionless)
u0ev = R(10) # Energía potencial del pozo (eV)
lnm = R(1) # Anchura del pozo (nm)
# ============================================================================

# ############################################################################
# Some independent physical constants
# h = R(6.626069311e-34) # Planck's constant in SI units (J·s)
# qe = R(1.6021765314e-19) # Charge of an electron in SI units (C)
# me = R(9.109382616e-31) # Mass of an electron in SI units (kg)
# ============================================================================
h = R(6.626069311e-34) # Plancks constant in SI units (J·s)
qe = R(1.6021765314e-19) # Absolute charge of an electron in SI units (C)
me = R(9.109382616e-31) # Mass of an electron in SI units (kg)
# ============================================================================

# Other physical constants
hbar = R(h/(2*pi))

# Quantities in SI units
m = mr * me
mSI = m
u0 = R(u0ev) * qe # Potential energy of the well (J)
u0SI = u0
l = R(lnm * 1e-9) # Well's width (m)
lSI = l

# Other reduced quantities apart from mr (reduced mass)
B = R(hbar / (sqrt(m*u0)))
C = R(u0)
lr = R(l/B)
u0r = R(1.0)

# ###########################################################################
# Variable declarations
# ###########################################################################
# xr -> reduced position: -oo < xr < oo
# epsilonr -> reduced energy: 0 < epsilonr < 1
xr,epsilonr = var('xr,epsilonr')

# ############################################################################
# ###################### FUNCTIONS ###########################################
# ############################################################################

# ############################################################################
# GLOBAL FUNCTIONS
# ############################################################################

# Ecuación de valores propios
# Primera parte de la ecuación
eigenvalfunc1(epsilonr) = tan(lr*sqrt(2*epsilonr))
# Segunda parte de la ecuación
eigenvalfunc2(epsilonr) = 2*epsilonr*sqrt(1/epsilonr-1)/(2*epsilonr-1)
# Función de valores propios que debe ser igual a cero
eigenvalfunc(epsilonr) = eigenvalfunc1(epsilonr) - eigenvalfunc2(epsilonr)

# Encuentra todas las raíces de 'f==0' en el intervalo (a,b)
def find_all_roots(f, a, b, roots, eps=1e-12):
	intervals_to_check = [(a,b)]
	while intervals_to_check:
		start, end = intervals_to_check.pop()
		try:
			root = find_root(f, start, end)
		except RuntimeError:
			continue
		if root in roots:
			continue
		if abs(f(root)) < 1:
			roots.append(RealNumber(root))
		intervals_to_check.extend([(start,root-eps), (root+eps,end)])
	roots.sort()
	return roots

# ############################################################################

def initial_quantities():
	print('-'*79)
	print('CONFIGURABLE QUANTITIES (lines 10, 24, 25 and 26):')
	print('Calculation precision: {0} bits'.format(precision_bits))
	print("Reduced particle mass -> mr= {0}".format(mr.n(digits=8)))
	print('Potential energy of the well -> u0= {0} eV'.format(u0ev.n(digits=8)))
	print('Width of the potential well -> l= {0} nm'.format(lnm.n(digits=8)))
	print('-'*79)
	print('REDUCED QUANTITIES:')
	print("Mass of the particle -> mr= {0}".format(mr.n(digits=8)))
	print('Potential energy of the well -> u0r= {0}'.format(u0r.n(digits=8)))
	print('Width of the potential well -> lr= {0} nm'.format(lr.n(digits=8)))
	print("Run:")
	print(" 'start()/start(DEBUG)' to calculate energy levels and wavefunctions")
	print(" 'change_initial_quantities()' to change the configurable quantities")
	print(" 'initial_quantities()' to see this values again")
	print('-'*79)
#print('Run: helpme() for options')


# Initial output:
initial_quantities()

#..............................................................................
# Gráfica de la función tangente en la ecuación de valores propios
pf1 = plot(eigenvalfunc1(epsilonr), (epsilonr,0,1), detect_poles='show', 
	color="blue")
# Gráfica de la función que no es tangente en la ecuación de valores propios
pf2 = plot(eigenvalfunc2(epsilonr), (epsilonr,0,1), detect_poles='show',
	color="red")
# Donde se crucen se encuentran los valores propios
#(pf1+pf2).show(ymin=-4*pi, ymax=4*pi)

#------------------------------------------------------------------------------

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#print('-------------------------------------------------------------')
#print('Estados ligados detectados -> len(ers)= {0}'.format(len(ers)))
#print("CORRECTO")
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
def start(mode=NODEBUG):
	print("CORRECT!")
# Valores propios
#ers = []
#find_all_roots(eigenvalfunc, R(1e-8), R(1-1e-8), ers)


#------------------------------------------------------------------------------
# Funciones propias
norms = []
#subnorms = []
funcs = []
for er in ers:
	#novale = []	
	# Normalization
	I1 = R(1/(2*sqrt(2*(1-er))))
	I3 = I1
	I21 = R(lr)
	I22 = R((1-2*er)*(2*lr*sqrt(2*er)-sin(2*lr*sqrt(2*er)))/(4*er*sqrt(2*er)))
	I23 = R(sqrt((1-er)/er)*(1-cos(2*lr*sqrt(2*er)))/(2*sqrt(2*er)))
	I2 = I21 + I22 + I23
	#novale.append(I1)
	#novale.append(I2)
	#novale.append(I3)
	#subnorms.append(novale)
	A = R(1/sqrt(I1 + I2 + I3))
	norms.append(A)
	zones = []
	f1(xr) = A * exp(sqrt(2*(1-er))*xr)
	f2(xr) = A * (cos(sqrt(2*er)*xr) + sqrt((1-er)/er)*sin(sqrt(2*er)*xr))
	c3 = R(A * (cos(lr*sqrt(2*er))+sqrt((1-er)/er)*sin(lr*sqrt(2*er))))
	f3(xr) = c3 * exp(-sqrt(2*(1-er))*(xr-lr))
	zones.append(f1(xr))
	zones.append(f2(xr))
	zones.append(f3(xr))
	funcs.append(zones)


# ########################## GRÁFICAS #########################################
# ...........................................................................
# Gráfica del pozo de potencial teniendo en cuenta los valores u0r y lr
ppot = line([(-2*lr,1),(0,1),(0,0),(lr,0),(lr,1),(3*lr,1)], color="black", thickness=2)

# ...........................................................................
# Gráfica de las energías propias junto con pozo de potencial
plen = ppot
for energy in ers:
	plen += line([(0,energy),(lr,energy)],color="red")

# ...........................................................................
# Gráfica de las funciones propias, energías propias y el pozo de potencial
plenfunc = plen
scal=0.10
for i in range(len(ers)):
	plenfunc += plot(scal*funcs[i][0]+ers[i], (x,-2*lr,0))
	plenfunc += plot(scal*funcs[i][1]+ers[i], (x,0.0,lr))
	plenfunc += plot(scal*funcs[i][2]+ers[i], (x,lr,3.0*lr))

# ...........................................................................
# Gráfica de los cuadrados de las funciones propias, energías propias
# y el pozo de potencial
plenfunc2 = plen
scal=0.50
for i in range(len(ers)):
	plenfunc2 += plot(scal*funcs[i][0]^2+ers[i], (x,-2*lr,0))
	plenfunc2 += plot(scal*funcs[i][1]^2+ers[i], (x,0.0,lr))
	plenfunc2 += plot(scal*funcs[i][2]^2+ers[i], (x,lr,3.0*lr))

# ...........................................................................
plf = []
for f in funcs:
	plno = ppot
	plno += plot(f[0], (x,-2*lr,0))
	plno += plot(f[1], (x,0.0,lr))
	plno += plot(f[2], (x,lr,3.0*lr))
	plf.append(plno)

# ...........................................................................
plf2 = []
for f in funcs:
	plno = ppot
	plno += plot(f[0]^2, (x,-2*lr,0))
	plno += plot(f[1]^2, (x,0.0,lr))
	plno += plot(f[2]^2, (x,lr,3.0*lr))
	plf2.append(plno)


# ...........................................................................
pls2 = []
for f in funcs:
	novale = []
	p1 = plot(f[0]^2, (x,-2*lr,0))
	p2 = plot(f[1]^2, (x,0.0,lr))
	p3 = plot(f[2]^2, (x,lr,3*lr))
	novale.append(p1)
	novale.append(p2)
	novale.append(p3)
	pls2.append(novale)

#in01 = integral(funcs[0][0]^2,(xr,-oo,0))
#in02 = integral(funcs[0][1]^2,(xr,0,lr))
#in03 = integral(funcs[0][2]^2,(xr,lr,oo))
#
#in11 = integral(funcs[1][0]^2,(xr,-oo,0))
#in12 = integral(funcs[1][1]^2,(xr,0,lr))
#in13 = integral(funcs[1][2]^2,(xr,lr,oo))
#
#in21 = integral(funcs[2][0]^2,(xr,-oo,0))
#in22 = integral(funcs[2][1]^2,(xr,0,lr))
#in23 = integral(funcs[2][2]^2,(xr,lr,oo))
#
#in31 = integral(funcs[3][0]^2,(xr,-oo,0))
#in32 = integral(funcs[3][1]^2,(xr,0,lr))
#in33 = integral(funcs[3][2]^2,(xr,lr,oo))
#
#in41 = integral(funcs[4][0]^2,(xr,-oo,0))
#in42 = integral(funcs[4][1]^2,(xr,0,lr))
#in43 = integral(funcs[4][2]^2,(xr,lr,oo))


# ########################## FUNCIONES ########################################
# FUNCIONES DE USUARIO
# ...........................................................................
def ayuda():
	print("FUNCIONES DISPONIBLES:")
	print("configini() -> Valores configurables iniciales")
	print("configred() -> Valores configurables reducidos")
	print("plot_ecvalprop() -> Gráfica de la ecuación de valores propios")
	print("valprop() ->  Energías reducidas de los estados ligados")
	print("valpropSI() ->  Energías de los estados ligados en el SI")
	print("valpropeV() ->  Energías de los estados ligados en el SI")
	print("funcpropshow() ->  Presenta variable que contiene a las funciones propias")
	print("plot_ppot() -> Gráfica del pozo de potencial")
	print("plot_plen() -> Gráfica de valores propios en pozo de potencial")
	print("plot_plenfunc() -> Gráfica de valores propios y funciones propias")
	print("plot_plenfunc2() -> Gráfica de valores propios y cuadrados de funciones propias")
	print("plot_func(i) -> Gráfica de función propia i y pozo de potencial")
	print("plot_func2(i) -> Gráfica del cuadrado de función propia i y pozo de potencial")

# ...........................................................................
def configini():
	print('VALORES INICIALES CONFIGURABLES:')
	print('Precisión de los cálculos -> precision_bits= {0} bits'.format(precision_bits))
	print('Masa de la partícula -> mr= {0}'.format(1.n(digits=8)))
	print('Profundidad del pozo de potencial -> Ur= {0}'.format(1.n(digits=8)))
	print('Anchura del pozo de potencial -> lr= {0}'.format(lr.n(digits=8)))

# ...........................................................................
def configred():
	print('Masa de la partícula -> mr= {0}'.format(1.n(digits=8)))
	print('Profundidad del pozo de potencial -> Ur= {0}'.format(1.n(digits=8)))
	print('Anchura del pozo de potencial -> lr= {0}'.format(lr.n(digits=8)))

# ...........................................................................
def plot_ecvalprop():
	(pf1+pf2).show(ymin=-4*pi, ymax=4*pi)

# ...........................................................................
def valprop():
	print("VALORES PROPIOS REDUCIDOS -> ers")
	i = 0
	for en in ers:
		print('ers[{0}]= {1}'.format(i, ers[i].n(digits=5)))
		i = i + 1

# ...........................................................................
def valpropSI():
	print("VALORES PROPIOS EN EL SI: u0 * ers[i]")
	i = 0
	for en in ers:
		print('E{0}= {1} J'.format(i, u0 * ers[i].n(digits=5)))
		i = i + 1

# ...........................................................................
def valpropeV():
	print("VALORES PROPIOS EN eV: u0ev * ers[i]")
	i = 0
	for en in ers:
		print('E{0}= {1} eV'.format(i, u0ev * ers[i].n(digits=5)))
		i = i + 1

# ...........................................................................
def funcpropshow():
	print("FUNCIONES PROPIAS EN CADA ZONA, EN LA VARIABLE 'funcs' -> funcs[i][j]")
	print("i= 0,1,... -> función propia 0, función propia 1, etc.")
	print("j= 0,1,2 -> zona 1, zona 2 y zona 3.")

# ...........................................................................
def plot_ppot():
	ppot.show()

# ...........................................................................
def plot_plen():
	plen.show()

# ...........................................................................
def plot_plenfunc():
	plenfunc.show()

# ...........................................................................
def plot_plenfunc2():
	plenfunc2.show()

# ...........................................................................
def plot_func(i):
	if i < nmax:
		plf[i].show()
	else:
		print("Error: máximo valor de 'i' es {0}".format(nmax))

# ...........................................................................
def plot_func2(i):
	if i < nmax:
		plf2[i].show()
	else:
		print("Error: máximo valor de 'i' es {0}".format(nmax))

# Cálculo del número de estados ligados de energía
#nmax = floor(R(lr*sqrt(2)/pi + 0.5-1e-8))
