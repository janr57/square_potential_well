# pozo_potencial_rectangular.jl
# Estados ligados de una partícula en un pozo rectangular de potencial.
# (C) 2015 José A. Navarro <josea.navarro@murciaeduca.es>

# Código unicode útil
# -------------------
# pi \u03c0
# alpha \u03b1
# delta \u03b4
# epsilon \u03b5
# hbar \u0127

# ############################################################################
# ###################### VARIABLES GLOBALES ##################################
# ############################################################################
NODEBUG= 0
DEBUG= 1

# Partícula en un pozo de potencial rectangular
# mr = 1.0   # Masa reducida de la partícula;
#            # Masa relativa de la partícula con respecto de la del electrón
#            # (adimensional);
#            #  mr = masa de la partícula / masa del electrón
# u0ev = 10.0 # Energía potencial del pozo (eV)
# lnm = 1.0 # Anchura del pozo (nm)
# ============================================================================
mr = 1.0; # Masa reducida de la partícula (adimensional)
u0ev = 10.0; # Energía potencial del pozo (eV)
lnm = 1.0; # Anchura del pozo (nm)
# ============================================================================

# ############################################################################
# Algunas constantes físicas independientes
# h = 6.626069311e-34 # Constante de Planck en el SI (J·s)
# qe = 1.6021765314e-19 # Carga del electrón en el SI (C)
# me = 9.109382616e-31 # Masa de la partícula (kg)
# ============================================================================
h = 6.626069311e-34; # Constante de Planck en el SI (J·s)
qe = 1.6021765314e-19; # Carga del electrón en el SI (C)
me = 9.109382616e-31; # Masa de la partícula (kg)
# ============================================================================

# Otras constantes físicas
ħ = h/(2 * π);

# Magnitudes en el SI
m = mr * me;
mSI = m;
u0 = u0ev * qe # Energía potencial del pozo (J)
u0SI = u0;
l = lnm * 1e-9; # Anchura del pozo (m)
lSI = l;

# Otras variables reducidas, aparte de la masa reducida (mr)
B = ħ / (sqrt(m*u0));
C = u0;
lr = l / B;
u0r = 1.0;

# ############################################################################
# ###################### FUNCIONES ###########################################
# ############################################################################

# ############################################################################
# FUNCIONES GLOBALES
# ############################################################################

# Ecuación de valores propios
# ---------------------------
# Primera parte de la ecuación
eigenvalfunc1(εr::Float64) = tan(lr*sqrt(2εr))
## Segunda parte de la ecuación
eigenvalfunc2(εr::Float64) = 2εr*sqrt(1/εr-1)/(2εr-1)
## Función de valores propios que debe ser igual a cero
eigenvalfunc(εr::Float64) = eigenvalfunc1(εr) - eigenvalfunc2(εr)
#eigenvalfunc(arr::Array{Float64,1}} = 

# ############################################################################
# LIBRERÍAS
# ############################################################################
using Roots


##..............................................................................
## Gráfica de la función tangente en la ecuación de valores propios
#pf1 = plot(eigenvalfunc1(εr), (εr,0,1), detect_poles='show', 
#	color="blue")
## Gráfica de la función que no es tangente en la ecuación de valores propios
#pf2 = plot(eigenvalfunc2(εr), (εr,0,1), detect_poles='show',
#	color="red")
## Donde se crucen se encuentran los valores propios
##(pf1+pf2).show(ymin=-4*pi, ymax=4*pi)
#
##------------------------------------------------------------------------------
## Valores propios
ers = Float64[]
nmax = floor(Int64, lr*sqrt(2)/π + 0.5-1e-8);
if v"0.3" <= VERSION < v"0.4-"
	sizehint(ers, nmax+1)
elseif VERSION >= v"0.4-"
	sizehint!(ers, nmax+1)
end

# ############################################################################
# ####################### FUNCIONES ##########################################
# ############################################################################

# ############################################################################
# FUNCIONES INTERNAS
# Produce un array de asíntotas de las dos funciones que forman la
# ecuación de valores propios: la función tangente y la otra
function find_all_asymptotes!(asymptotes, lr)
	# Cálculo del número total de nodos de la función tangente
	#tangentFuncAsymptNum = int64(floor(lr*sqrt(2)/π + 0.5-1e-8));
	tangentFuncAsymptNum = floor(Int64, lr*sqrt(2)/π + 0.5-1e-8);
	#tangentFuncAsymptNum = iround(lr*sqrt(2)/π + 0.5-1e-8);
	# Construye array de asíntotas de la función tangente
	for i in 1:tangentFuncAsymptNum
		push!(asymptotes, 0.5*((2i-1)π/(2lr))^2)
	end
	# Se añade asíntota de la segunda función	
	push!(asymptotes,0.5)
	# Se añade el máximo valor posible
	push!(asymptotes,1.0)
	sort!(asymptotes)
end

# Se han encontrado una o más posibles raíces.
# Hay que elegir la verdadera raíz.
function chooseTheRoot!(roots, probRoots, f, almostZero, left, right, δ, mode)
	errMsg = "interval ($left, $right); delta= δ; without a root"

	values =filter(x -> x < almostZero, abs(map(f,probRoots)))
	if mode == DEBUG
		println("Values:")
		println(values)
	end

	if length(values) > 0
		leastValue = minimum(values)
		iVal = findfirst(values, leastValue)
		# Se ha encontrado la raíz correcta
		root = probRoots[iVal]
		if mode == DEBUG
			println("Raíz encontrada: $root")
		end
		push!(roots, root)
	else
		if mode == DEBUG
			println("[SOME] $msgNoRoot")
		end
		error("\n$errMsg")
	end
end

# No se ha encontrado ninguna raíz.
# Busca de nuevo, utilizando ahora el método Roots.fzero
function searchRootAgain!(roots, f, f1, almostZero, left, right, δ, mode)
	errMsg= """
	en el intervalo ($left, $right); con delta= $δ; no hay ninguna raíz"""
	noRootMsg= "No hay raíz en este intervalo, pero no es un error"
	try
		root = Roots.fzero(f, left+δ, right-δ)
		if abs(f(root)) < almostZero
			mode == DEBUG && println("Raíz encontrada: $root")
			push!(roots,root)
		else
			mode == DEBUG && println("[ZERO1] $msgNoRoot")
			error("\n$errMsg")
		end
	catch e
		if isa(e, ErrorException)
			#println(e)
			if right == 1.0 && f1(1.0) > 0
				# Debe haber una raíz, pero esta
				# excepción se lanzó porque no se encuentra
				mode == DEBUG && println("[ZERO2] $msgNoRoot")
				error("\n$errMsg")
			else
				mode == DEBUG && println("$noRootMsg")
			end
		else
			throw(e)
		end
	end
end

# Ejecuta distintos algoritmos según se haya encontrado alguna posible raíz
# o ninguna, en un intervalo dado por: (left, right)
function analizeProbRoots!(roots, probRoots, f, f1, left, right, δ, mode)
	# Aproximmación: f(root) == 0 ----> f(root) < almostZero
	almostZero= 0.1

	if length(probRoots) > 0
		# Hay una o más raíces posibles.
		# Se debe detectar la correcta (si existe).
		# Se descartan aquellas que no cumplan abs(f()) < 0.1
		chooseTheRoot!(roots,probRoots,f,almostZero,left,right,δ,mode)
	else
		# Prueba con la función fzero si no estamos en el
		# último intervalo
		searchRootAgain!(roots, f, f1, almostZero, left, right, δ, mode)
	end
end

# Encuentra todas las raíces de 'f==0'
# TODO: Comprobar que en cada intervalo entre asíntotas de la tangente
# hay una raíz (quizás menos en el último y/o en uno con epsilonr cerca de 0.5)
# En cada intervalo debe haber una raíz (excepto quizás el ultimo ¿y/o alguno
# cerca de la asíntota 0.5?)
function findAllRoots!(roots, f, f1, lr, mode)
	# Intervalos sin raíz
	asymptotes = Float64[]
	find_all_asymptotes!(asymptotes, lr)
	if mode == DEBUG
		println("asymptotes:")
		println("$asymptotes")
	end
	# Éste será el elemento de la izquierda del primer intervalo
	right = shift!(asymptotes)
	while length(asymptotes) != 0
		# left y right son los extremos del intervalo
		left = right
		right = shift!(asymptotes)
		long = right - left
		δ = long * 1e-7
		if mode == DEBUG
			println(repeat("-", 79))
			println("Intervalo: ($left, $right)")
			println("Longitud intervalo= $long")
			println("Delta= $δ")
		end
		probRoots = Roots.fzeros(f, left+δ, right-δ)
		if mode == DEBUG
			println("Raíces probables:")
			println("$probRoots")
		end
		# Analiza el array 'probRoots'
		analizeProbRoots!(roots, probRoots, f, f1, left, right, δ, mode)
	end
	if mode == DEBUG
		println(repeat("-", 79))
	end
	sort!(roots)
	return
end

function start(mode=NODEBUG)
	# Vacía el array de valores propios 'ers'
	while length(ers) > 0
		shift!(ers)
	end
	# Busca valores propios
	findAllRoots!(ers, eigenvalfunc, eigenvalfunc1, lr, mode)
	# Personaliza mensaje resumen
	msg1 = "Encontrado 1 valor propio"
	msgs = "Encontrados $(length(ers)) valores propios"
	length(ers) == 1 ? msg= msg1 : msg= msgs
	println(msg)

	nmax = floor(Int64, lr*sqrt(2)/π + 0.5-1e-8);
	println("nmax= $nmax")
	println("CORRECTO")
	println(repeat("*",79))
	println("Ejecute: ayuda() para ver opciones disponibles")
	println(repeat("*",79))
end

function valPropios()
	headMsg= "Valores propios: $(length(ers))"
	println(repeat("-", 79))
	println(headMsg)
	println(repeat("-", 79))
	for i in 1:length(ers)
		println("ers[$i]= $(ers[i])")
	end
end

function ayuda()
	println("valPropios() -> Listado de valores propios")
end

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Salida inicial de información
println(repeat("*",79))
println("CONFIGURABLE QUANTITIES (lines 10, 24, 25 and 26):")
println("Reduced particle mass -> mr= $mr")
println("Potential energy of the well -> u0= $u0ev eV")
println("Width of the potential well -> l= $lnm nm")
println(repeat("*",79))
println("VALORES REDUCIDOS INICIALES:")
println("Masa de la partícula-> mr= $mr")
println("Profundidad del pozo de potencial -> U0r= $u0r")
println("Longitud del pozo de potencial -> lr= $lr")
println(repeat("*",79))
println("Run 'start()/start(DEBUG)' to calculate energy levels and wavefunctions")
println(repeat("*",79))
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||











