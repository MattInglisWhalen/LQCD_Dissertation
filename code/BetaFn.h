#ifndef BETAFN_HH
#define BETAFN_HH

#include <vector>

#include "TMath.h"
#include "TComplex.h"

#define RZ3 1.2020569031595942853997  // Value of the Riemann zeta function at z=3.0

//class for the QCD beta function, defined as \partial \alpha / \partial \ln \mu^2 = -\beta_0\alpha^2 - beta_1\alpha^3 - ... = \beta(\alpha)

class BetaFn{

 public:

  // constructors and destructors
  BetaFn();                    // default constructor
  BetaFn(Int_t,Int_t,Int_t=1); // _nloops, _nflav, _variable
  ~BetaFn();                   // default destructor

  // get functions
  Int_t nloops() const;
  Int_t nflav() const;
  Int_t variable() const;
  std::vector<Double_t> coeffs() const;

  // set functions
  void set_nloops(Int_t);
  void set_nflav(Int_t);
  void set_variable(Int_t);
  void set_coeff(Int_t,Double_t);

  //operators
  BetaFn& operator=(const BetaFn&);       //assignment operator
  Double_t& operator[](const Int_t);      //set index
  Double_t operator[](const Int_t) const; //get index

  //useful functions
  Double_t eval(Double_t);  //gets the value of the beta function at a certain alpha
  std::vector< TComplex > roots(); //get the (complex) roots of the beta function
  std::vector< TComplex > nonzeroroots(); //get the (complex) nonzero roots of the beta function

 private:
  Int_t _nloops;    // up to 4
  Int_t _nflav;     // between 0 and 6, inclusive
  Int_t _variable;  // 2==g, 1==alpha (=g^2/4pi) , 0==a (=alpha/4pi)
  std::vector<Double_t> _coeffs; //coefficients of the beta function

};

#endif

