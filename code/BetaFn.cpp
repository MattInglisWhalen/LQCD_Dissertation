#include "BetaFn.h"

#include <iostream>

// ---- constructors and destructors ----

BetaFn::BetaFn(){
  _nloops=1;
  _nflav=0;
  _variable=2;
  this->set_coeff(0,1.);
};

BetaFn::BetaFn(Int_t nl, Int_t nf, Int_t var) : _nloops(nl), _nflav(nf), _variable(var){

  // ----------- universal first two coefficients ----------------

  // index 0
  if( nl>0 ) {
    _coeffs.push_back( ( 11.- nf*(2./3.) )*TMath::Power( 4.*TMath::Pi() , -1.*((Double_t)var) ) );
  }
  else{
    _coeffs.push_back( 0. );
  }
  //index 1
  if( nl>1 )   {
    _coeffs.push_back( ( 102.- nf*(38./3.)  )*TMath::Power( 4.*TMath::Pi() , -2.*((Double_t)var) ) );
  }
  else{
    _coeffs.push_back( 0. );
  }

  // ----------- scheme dependent remaining coefficients ------------

  if( var>=0 ){ 
    //index 2
    if ( nl>2 )  {
      _coeffs.push_back( ( 2857./2. -nf*(5033./18.) + nf*nf*(325./54.)  )*TMath::Power( 4.*TMath::Pi() , -3.*((Double_t)var) ) );
    }
    else{
      _coeffs.push_back( 0. );
    }
    //index 3
    if( nl>3 )   {
      _coeffs.push_back( ( (149753./6. + 3564.*RZ3) - nf*(1078361./162.+6508.*RZ3/27.) 
          + nf*nf*(50065./162.+6472.*RZ3/81.) + nf*nf*nf*(1093./729.) )*TMath::Power( 4.*TMath::Pi() , -4.*((Double_t)var) ) );
    }
    else{
      _coeffs.push_back( 0. );
    }
  } 
};

BetaFn::~BetaFn(){
  // do nothing
};

// ---- get and set functions ----
Int_t BetaFn::nloops() const{
  return _nloops;
};
Int_t BetaFn::nflav() const{
  return _nflav;
};
Int_t BetaFn::variable() const{
  return _variable;
};
std::vector<Double_t> BetaFn::coeffs() const{
  return _coeffs;
};

void BetaFn::set_nloops(Int_t nl){
  if( (nl<1) || (nl>4) ){
    std::cout << "BetaFn: invalid number of loops\n";
    exit(-1);
  }
  _nloops=nl;
};
void BetaFn::set_nflav(Int_t nf){
  if( (nf<0) || (nf>6) ){
    std::cout << "BetaFn: invalid number of flavours\n";
    exit(-1);
  }
  _nflav=nf;
};
void BetaFn::set_variable(Int_t var){
  if( (var<0) || (var>3) ){
    std::cout << "BetaFn: invalid variable, g==2, alpha==1, a==0\n";
    exit(-1);
  }
  _variable=var;
};
//this function should not be used
void BetaFn::set_coeff(Int_t ind,Double_t val){
  if(ind<2){
    std::cout<<"Note: the beta function coefficient " << ind << " is universal -- be sure you know what you're doing!\n";
  }
  while(_coeffs.size()<(unsigned)ind ){
    _coeffs.push_back(0.);
  }
  _coeffs[ind]=val;
};

// ---- operators ----
//assignment operator
BetaFn& BetaFn::operator= (const BetaFn& bf){
  _nloops=bf.nloops();
  _nflav=bf.nflav();
  _variable=bf.variable();
  _coeffs=bf.coeffs();
  return *this;
};
//set index
Double_t& BetaFn::operator[] (const Int_t ind){
  if( (ind<0)||(ind>=4) ){
    std::cout<<"BetaFn: invalid index for [] get operator.\n";
    exit(-1);
  }
  return _coeffs[ind];
};

//get index
Double_t BetaFn::operator[] (const Int_t ind) const{
  if( (ind<0)||(ind>=4) ){
    std::cout<<"BetaFn: invalid index for [] get operator.\n";
    exit(-1);
  }
  return _coeffs[ind];
};

// ---- useful functions ---
Double_t BetaFn::eval(Double_t var){
  Double_t prefac;
  if(_variable==2){
    prefac=-var*var*var;
    var*=var;
  }
  else{
    prefac=-var*var;
  }
  return prefac*( (*this)[0] + (*this)[1]*var +(*this)[2]*var*var + (*this)[3]*var*var*var );
};
std::vector<TComplex> BetaFn::roots(){

  Double_t a,b,c,d;
  a=b=c=d=0.;

  std::vector<TComplex> v;
  v.push_back( 0. ); //there are always two roots at x=0
  if(_variable==2){
    v.push_back( 0. );  //it has one more root at x=0 if written in terms of g
  }
  if( _nloops==1 ){
    a=(*this)[0];
  }
  else if( _nloops==2 ){
    a=(*this)[1]; b=(*this)[0];
    v.push_back( -b/a );
  }
  else if( _nloops==3 ){
    a=(*this)[2]; b=(*this)[1]; c=(*this)[0];
    v.push_back( -b/(2.*a) + TComplex::Sqrt(b*b-4.*a*c)/(2.*a) );
    v.push_back( -b/(2.*a) - TComplex::Sqrt(b*b-4.*a*c)/(2.*a) );
  }
  else if( _nloops==4 ){
    a=(*this)[3]; b=(*this)[2]; c=(*this)[1]; d=(*this)[0];
    Double_t Delta_0 = b*b-3.*a*c;
    Double_t Delta_1 = 2.*b*b*b - 9.*a*b*c +27.*a*a*d;
    TComplex C=TComplex::Power( 0.5*Delta_1 + 0.5*TComplex::Power( TComplex(Delta_1*Delta_1 -4.*Delta_0*Delta_0*Delta_0), 1./2. ) , 1./3. );
    TComplex u_1(1.,0.);
    TComplex u_2(-0.5,+0.5*TMath::Sqrt(3));
    TComplex u_3(-0.5,-0.5*TMath::Sqrt(3)); 
    v.push_back( (-1./(3.*a))*(b+u_1*C+Delta_0/(u_1*C) ) );
    v.push_back( (-1./(3.*a))*(b+u_2*C+Delta_0/(u_2*C) ) );
    v.push_back( (-1./(3.*a))*(b+u_3*C+Delta_0/(u_3*C) ) );
  }
  else{
    v[0]=-1.;
  }
  return v;
};

std::vector<TComplex> BetaFn::nonzeroroots(){
  std::vector<TComplex> v,temp=this->roots();
  for(Int_t i=0;i<(signed)temp.size();++i){
    if(temp[i].Rho()>0.000000001){
      v.push_back(temp[i]);
    }
  }
  return v;
};


