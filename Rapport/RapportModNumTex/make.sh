./clean.sh; 
pdflatex -shell-escape rapport.tex;
pdflatex -shell-escape rapport.tex;

./clean.sh; 
echo
echo
echo "================================================="
echo "Taper evince rapport.pdf" pour voir votre rapport
echo "================================================="
echo
echo

