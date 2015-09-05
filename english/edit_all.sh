# edit_all.sh
# Usage: $ sh edit_all <editor>
# 

files="Makefile\
	square_potential_well.tex\
	text/schrodinger_equation.tex"

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
