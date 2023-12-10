#! /bin/bash

# cp ../alpha_of_mu_analytic.pdf ./   #it looks good already
# cp ../method_I_stability.pdf   ./   #it looks good already
# cp ../method_II_stability.pdf  ./   # it looks good already

OBJ="manuscript_s1328668"

if [ -f ${OBJ}.pdf ] ; then
  mv ${OBJ}.pdf .${OBJ}.pdf
fi
if [ -f ${OBJ}.tex ] ; then
  pdflatex ${OBJ}.tex
fi
