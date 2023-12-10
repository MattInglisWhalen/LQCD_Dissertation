#! /bin/bash

FILE="for_roger_aug_01.pdf"

if [ -f ${FILE} ] ; then
  rm ${FILE}
fi

#pdftk methodI_lat_quenched_direct.pdf methodI_lat_quenched_recip.pdf methodIquenched.pdf methodIIquenched.pdf methodIIIquenched.pdf methodIIPquenched.pdf methodIIIPquenched.pdf Plaquette_Chiral_Extrap.pdf r0_Over_a_Chiral_Extrap.pdf methodI_nf=2.pdf methodIIP_nf=2.pdf methodIIIP_nf=2.pdf cat output ${FILE}#pdftk methodI_lat_quenched_direct.pdf methodI_lat_quenched_recip.pdf methodIquenched.pdf methodIIquenched.pdf methodIIIquenched.pdf methodIIPquenched.pdf methodIIIPquenched.pdf Plaquette_Chiral_Extrap.pdf r0_Over_a_Chiral_Extrap.pdf methodI_nf=2.pdf methodIIP_nf=2.pdf methodIIIP_nf=2.pdf cat output ${FILE}#pdftk methodI_lat_quenched_direct.pdf methodI_lat_quenched_recip.pdf methodIquenched.pdf methodIIquenched.pdf methodIIIquenched.pdf methodIIPquenched.pdf methodIIIPquenched.pdf Plaquette_Chiral_Extrap.pdf r0_Over_a_Chiral_Extrap.pdf methodI_nf=2.pdf methodIIP_nf=2.pdf methodIIIP_nf=2.pdf cat output ${FILE}
#pdftk methodIII_nf=2.pdf methodIIIP_nf=2.pdf rat32_v_rat20.pdf rat32_v_rat20_small.pdf cat output ${FILE}
pdftk Plaquette_Chiral_Extrap_new.pdf r0_Over_a_Chiral_Extrap_new.pdf methodI_nf=2.pdf methodII_nf=2.pdf methodIIP_nf=2.pdf methodIII_nf=2.pdf methodIIIP_nf=2.pdf cat output ${FILE}
