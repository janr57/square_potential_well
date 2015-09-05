# editar_todo.sh
# Uso: $ sh editar_todo.sh <editor>
# 

files="Makefile\
	pozo_potencial_rectangular.tex\
	texto/ecuacion_schrodinger.tex\
	texto/estados_ligados.tex\
	texto/estados_no_ligados.tex\
	texto/solucion_numerica.tex"

if [ $# -ne 1 ]
then
	echo "USAGE: $0 <editor>"
	exit 1
elif [ $1 = 'gvim' ]
then
	$1 -geom 80x61 $files
elif [ $1 = 'vim' ]
then
	$1 $files
else
	$1 $files &
fi
