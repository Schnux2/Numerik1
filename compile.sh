OPTIONS="-std=c11 -Wall -Wimplicit-fallthrough -Wpedantic"
LIBS="-lm -lncurses"

INFILE=test.c
OUTFILE="${INFILE/.c/.exe}"

rm -f $OUTFILE
echo ==================================
echo attempting to compile
echo gcc $OPTIONS -o $OUTFILE $INFILE matrix.c LR_Zerlegung.c QR_Zerlegung.c Bisektion_Newton.c Diskretisierung.c Interpolation.c Integration.c SOR.c $LIBS
echo ..................................
echo ""

gcc $OPTIONS -o $OUTFILE $INFILE matrix.c LR_Zerlegung.c QR_Zerlegung.c Bisektion_Newton.c Diskretisierung.c Interpolation.c Integration.c SOR.c $LIBS

if [[ -f $OUTFILE ]];then
	echo ""
	echo ..................................
	echo Compilation successful.
	echo ..................................
	echo ""
	./$OUTFILE
else
	echo "Some compilation error occurred."
fi