#! /bin/bash

#cp ../mu_over_lambda_analytic.pdf ./

#cp ../alpha_of_mu_analytic.pdf ./

OBJ="short_proof"

if [ -f ${OBJ}.pdf ] ; then
  mv ${OBJ}.pdf .${OBJ}.pdf
fi
if [ -f ${OBJ}.tex ] ; then
  pdflatex ${OBJ}.tex
fi
