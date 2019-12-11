#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <cmath> 
#include <math.h>
#include <iostream> 
#include <complex> 
#include <stdlib.h> 
#include <vector> 
#include <algorithm>
#include <unordered_map> 
#include <tuple>
#include <utility>
#include "TVectorD.h"
#include "TGraphErrors.h"

std::vector<std::vector<int>> getCombinations(int n) {
   std::vector<int> set;
   std::vector<double> millers;
 
  for (int i = 0; i < n; i++) {
    set.emplace_back(i);
  }
  
  std::vector<std::vector<int>> permutaions; 
  for (int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      for (int k = j; k < n; k++) {
        int millerH = set[i];
        int millerK = set[j];
        int millerL = set[k];
  
        if (millerH == 0 && millerK == 0 && millerL == 0) continue; // No (000) Miller Plane Exists
        double quantity = pow(millerH, 2) + pow(millerK, 2) + pow(millerL, 2); 
        bool alreadySeen = std::find(millers.begin(), millers.end(), quantity) != millers.end(); 
        if (alreadySeen) continue; 
        std::vector<int> perm = {set[i], set[j], set[k]}; 
        permutaions.emplace_back(perm);
        millers.emplace_back(quantity);  
      }
    }
  }
  return permutaions;
}


double pickNeutronWavelength(double T) {
  // Returns a neutron wavelength in Angstroms, temperature in K
  double h = 6.6e-34; // SI
  double c = 3e8; // ms-1
  double kb = 1.38e-23; // JK-1
  double m_n = 1.67e-27; // kg
 
  //TF1* maxwellBoltzmann = new TF1("maxwell", "([0]/(2*[1]*[2]*[3]))^(1.5) * exp(-x/(2*[2]*[3]))", 0.0, 1e-19);
 
  TF1* maxwellBoltzmann = new TF1("maxwell", "([0]/(2*[1]*[2]*[3]))^(0.5) * exp(-([0]*x*x)/(2*[2]*[3]))", 0.0, 5000.0);
  maxwellBoltzmann->SetParameters(m_n, M_PI, kb, T);
  double velocity = maxwellBoltzmann->GetRandom();  
  double energy = 0.5 * m_n * pow(velocity, 2); 
  double wavelength = h/sqrt(2*m_n*energy);
  return wavelength * 1e10; 
}

  //double neutronTemp = 290;
  //double wavelength = pickNeutronWavelength(neutronTemp); // drawn from blackbody distribution for thermal neutrons. 
  //double wavelength = 1.0;
  //double k = 2*M_PI/wavelength; 
  
  // Parameters
  int Z = 54;
  // double wavelength = 1.0; 
  double lattice_spacing = 100.0;
  double mu = 200.0; // eV
  double q_0 = 6.86; // A-1
  
  double b_F = -1.468e-3; // Foldy scattering length, fm
  double r_N_square = -0.1161; // PRL, fm^2
  double a_0 = 5.2917721067e4; // Scipy value, fm
  double b_I = (1/3)*(1.6726219e-27/9.10938356e-31)*(r_N_square/a_0); // fm Intrinsic scattering length
  double b_c = 4.92; // fm, Xe coherent scattering length from NIST database

  // Dimensionless constants from scattering lengths
  double x_em = -Z*(b_F + b_I)/b_c; // Need - sign to make a sensible em x-section, also Sears Eq 3.19 has - out front
 
  
  // New interaction stuff 
  double h_bar = 6.58211956e-16; // eV s
  double c_angstrom_per_sec = 2.998e+18; // speed of light in angstrom per second
  

  double mu_gev = mu*1e-9; // GeV
  double g_squared = 1e-16; // GeV^-2
  double xe_nuclear = 122; // GeV c^-2
  double neutron_nuclear = 0.93957; // GeV c^-2
  double b_c_ev = 25; // GeV^-1 
  double m_n = 0.93957; // ?? weird m_n term in GeV, neutron mass?? Discussed with Alex he thinks it's neutron mass also which 
                       // is weird because then why have it in Q_1 as well? 
  double x_new = m_n*g_squared*xe_nuclear*neutron_nuclear/(2*M_PI*b_c_ev*pow(mu_gev,2));


double f(double q) {
  return 1/pow(1+3*(pow(q/q_0, 2)), 0.5); 
}

double vanillaCrossSection(Double_t *angle, Double_t *par) {
  double scatAngle = angle[0]/2.0;
  double k = 2*M_PI/angle[1]; 
  double constTerm = 1;
  double emTerm = 2*x_em*(1-f(2*k*sin(scatAngle)));
  //return scatAngle; 
 return par[0]*2*M_PI*pow(b_c,2)*(constTerm + emTerm); 
}

double newCrossSection(Double_t *angle, Double_t *par) {
  double scatAngle = angle[0]/2.0;
  double k = 2*M_PI/angle[1]; 
  double q = 2*k*sin(scatAngle); 
  double constTerm = 1;
  double emTerm = 2*x_em*(1-f(q));  
  double mu_sq = pow(mu, 2); 
  double q_transfer_sq = pow(q*h_bar*c_angstrom_per_sec, 2);
  double newTerm = 2*x_new * mu_sq / (mu_sq + q_transfer_sq); 
  std::cout << "NEW TERM " << newTerm << std::endl;
  return par[0]*2*M_PI*pow(b_c,2)*(constTerm + emTerm + newTerm); 
 
}

std::complex<double> complexExp(double x) {
  const std::complex<double> i(0, 1); 
  return std::exp(-2*M_PI*i*x); 
}

std::complex<double> F_orthorhombic_cubic(int h, int k, int l, double f, double a, double b, double c) {
  std::complex<double> p0 = complexExp(0.0); 
  std::complex<double> p1 = complexExp(a*h); 
  std::complex<double> p2 = complexExp(b*k); 
  std::complex<double> p3 = complexExp(c*l);
  return f*(p0 + p1 + p2 + p3);
}

std::complex<double> F_primitive_cubic(int h, int k, int l, double f, double a, double b, double c) {
  return f; 
}

double square(double x) {
  return pow(x, 2); 
}

double hklThetaOrthorhombic(int h, int k, int l, double wavelength, double a, double b, double c) {
  double sinx = (wavelength/2.0) * pow(square(h/a) + square(k/b) + square(l/c), 0.5);
  double result = std::asin(sinx);
  if (!std::isnan(result)) {
    return result; 
  } 
  //std::cout << "Incorrect Argument for arcsin in hklTheta, putting at theta = 0" << std::endl;
  return -1; 
}

double I_ortho(int h, int k, int l, double f, double a, double b, double c) {
  return std::norm(F_primitive_cubic(h, k, l, f, a, b, c)); 
}

double angleToWavelength(double theta, int h, int k, int l, double a, double b, double c) {
  double d_hkl_reciprocal = pow(square(h/a) + square(k/b) + square(l/c), 0.5);  

  return 2*sin(theta)*1/d_hkl_reciprocal; 
} 

void neutrons() {
  
  gStyle->SetOptStat(0);  
  int nentries = int(1e6); 
  
  //TF1* vanilla_xsection = new TF1("vanilla", "6.283*[0]**2 * (1 + 2 * [1] * (1 - (1 + 3 * (2*[2]*sin(x/2)/[3])**2)**(-0.5)))", 0.0, 0.5);
  TF1* vanilla_xsection = new TF1("vanilla", vanillaCrossSection, 0.0, 0.5, 1); 
  vanilla_xsection->SetParameter(0, 1.0);
  vanilla_xsection->SetLineColor(kBlue);
 // TF1* new_xsection = new TF1("new", "6.283*[0]**2 * (1 + 2 * [1] * (1 - (1 + 3 * (2*[2]*sin(x/2)/[3])**2)**(-0.5)) + 2* [4] *[5]**2/([5]**2 + (2*[2]*[6]*[7]*sin(x/2))**2))", 0.0, 0.5);
  TF1* new_xsection = new TF1("new", newCrossSection, 0.0, 0.5, 1); 
  new_xsection->SetParameter(0, 1.0); 
  //new_xsection->SetParameters(b_c, x_em, k, q_0, x_new, mu, h_bar, c_angstrom_per_sec);
  //vanilla_xsection->SetParameters(b_c, x_em, k, q_0);
  
  //double integralVanilla = vanilla_xsection->Integral(0.0, 0.5); 
  //TF1* vanillaNormed = new TF1("vanilla_normed", "[0]*vanilla", 0.01, 0.5);
  /*
  TH1F* plain_xsect = new TH1F("h", "Neutron Nucleus Cross Section", 100, 0.0, 0.5);
  TH1F* modified_xsect = new TH1F("new", "Cross Section", 100, 0.0, 0.5);
  modified_xsect->SetLineColor(kRed);
  plain_xsect->SetLineColor(kBlue);

  for (int i = 0; i < nentries; i++) {
    modified_xsect->Fill(new_xsection->GetRandom());
    plain_xsect->Fill(vanilla_xsection->GetRandom());
  }
  TCanvas* canv = new TCanvas("c", "c", 1000, 600); 
  plain_xsect->Draw("E1"); 
  modified_xsect->Draw("E1 same"); 
  plain_xsect->Fit("vanilla"); 
  modified_xsect->Fit("new"); 
  std::cout << "Vanilla x^2/n = " << vanilla_xsection->GetChisquare()/vanilla_xsection->GetNDF() << std::endl;
  std::cout << "New x^2/n = " << new_xsection->GetChisquare()/new_xsection->GetNDF() << std::endl;
  */

  // -------------------- LATTICE -------------------- 
  
  TH1F* braggAnglesOriginal = new TH1F("bragg", "bragg", 1000, 0.0, 1.57); 
  braggAnglesOriginal->SetLineColor(kBlack);
  braggAnglesOriginal->GetXaxis()->SetTitle("Angle / radians"); 
  braggAnglesOriginal->GetYaxis()->SetTitle("Number of neutrons"); 
  TH1F* braggAnglesNew = new TH1F("braggnew", "braggnew", 100, 0.0, 1.57); 
  braggAnglesNew->SetLineColor(kRed);
  braggAnglesNew->GetXaxis()->SetTitle("Angle / radians"); 
  braggAnglesNew->GetYaxis()->SetTitle("Number of neutrons"); 
  TH1F* vanillaCrystal = new TH1F("vanillaCrystal", "Plain Crystal Distribution", 50, 0.0, 1.57);   
  
  TH1F* wavelengths = new TH1F("wavelengths", "Wavelength Distribution", 100, 0.0, 5.0);   
  double a = 10.0; // Angstrom
  double b = 10.0; // Angstrom
  double c = 10.0; // Angstrom
  
  std::vector<double> angles;
  std::vector<double> vanillaDistribution; 
  std::vector<double> newDistribution; 
  std::vector<double> millerPlanes;
 
  int miller = 10; 
  TRandom* rand = new TRandom(); 
  
  // Not normalising these because then lose units (discussion with Alex from 25/11/19) 
 // double momRatio = square(k/q_0); 
 // TF1* normalisedVanilla = new TF1("normalVanilla", vanillaCrossSection, 0.0, 0.5, 2);
 // double vanillaIntegral = 2*M_PI*square(b_c) * (2 * (1 + 2*x_em) - 2*x_em*((2*pow(12*momRatio+1, 0.5) - 1)/(6*momRatio)));
 // normalisedVanilla->SetParameter(0, 1);//1/vanillaIntegral);
  
  //TF1* newXSection = new TF1("newXSection", newCrossSection, 0.0, 0.5, 2); 
 // newXSection->SetParameter(0, 1); 

  typedef std::tuple<int, int, int> miller_plane_t; 
  
  std::unordered_map<double, miller_plane_t> angleToIndices; 
  std::vector<std::vector<int>> combos = getCombinations(miller); 
 
  for (int i = 0; i < nentries; i++) {

    double neutronTemp = 290; // K 
    //double wavelength = pickNeutronWavelength(neutronTemp);
    double wavelength = 1.0;
    wavelengths->Fill(wavelength); 
    int index = rand->Integer(combos.size());
    std::vector<int> combo = combos[index]; 
          
    int h = combo[0]; 
    int k = combo[1]; 
    int l = combo[2];

 
    double angle = hklThetaOrthorhombic(h, k, l, wavelength, a, b, c); 
    braggAnglesOriginal->Fill(angle);
    miller_plane_t plane = std::make_tuple(h, k, l);
    bool inAngles = std::find(angles.begin(), angles.end(), angle) != angles.end(); 
    if (!inAngles) angles.emplace_back(angle); 
    
    TAxis* xaxis = braggAnglesOriginal->GetXaxis(); 
    Int_t binx = xaxis->FindBin(angle); 
    double centre = braggAnglesOriginal->GetBinCenter(binx); 
    angleToIndices[centre] = plane;
     
   }

   int nbins = braggAnglesOriginal->GetNbinsX();
   TVectorD grAngles(nbins); 
   TVectorD xsectionAngles(nbins); 
   TVectorD errsX(nbins); 
   TVectorD errsY(nbins); 

   TVectorD newXsectionAngles(nbins); 
   TVectorD errsNewX(nbins); 
   TVectorD errsNewY(nbins); 

   // Normalisation Constant Calculation (only for primitive cubic, can change it for general F_hkl later) 

   TF1* normalisedVanilla = new TF1("normalVanilla", vanillaCrossSection, 0.0, 1.57, 1);
   TF1* newXSection = new TF1("newXSection", newCrossSection, 0.0, 1.57, 1); 
   double normalisationVanilla = 0.0; 
   double normalisationNew = 0.0; 
   
  for (int i = 0; i < nbins; i++) {
     double centre = braggAnglesOriginal->GetBinCenter(i); 
     double content = braggAnglesOriginal->GetBinContent(i); 
     if (content == 0) continue; 
    
       miller_plane_t plane = angleToIndices[centre];
       int h = std::get<0>(plane); 
       int k = std::get<1>(plane); 
       int l = std::get<2>(plane); 
              
       double wavelength = angleToWavelength(centre, h, k, l, a, b, c); 
       double args[2] = {centre, wavelength}; 
       normalisedVanilla->SetParameter(0, 1.0); 
       double contribution = normalisedVanilla->EvalPar(args); 

       newXSection->SetParameter(0, 1.0);
       double newXSectionContribution = newXSection->EvalPar(args); 
       double newContribution = newXSection->EvalPar(args);

       normalisationVanilla += I_ortho(h, k, l, newContribution, a, b, c); 
       normalisationNew += I_ortho(h, k, l, newContribution, a, b, c);
   }
   
   double height = -120;
   
   for (int i = 0; i < nbins; i++) {
       double centre = braggAnglesOriginal->GetBinCenter(i); 
       double error = braggAnglesOriginal->GetBinError(i);
       if (error == 0.0) continue;  
       grAngles[i] = centre;  
       
       miller_plane_t plane = angleToIndices[centre];
       int h = std::get<0>(plane); 
       int k = std::get<1>(plane); 
       int l = std::get<2>(plane); 
              
       double wavelength = angleToWavelength(centre, h, k, l, a, b, c); 
       //std::cout << "Wavelength " << wavelength << std::endl;
       //double wavelength = 1.0; 

       normalisedVanilla->SetParameter(0, 1.0); 
       newXSection->SetParameter(0, 1.0); 

       double args[2] = {centre, wavelength}; 
       double xsectionF = normalisedVanilla->EvalPar(args);
       double xsectionNewF = newXSection->EvalPar(args); 
       double finalVal = I_ortho(h, k, l, xsectionF, a, b, c)/normalisationVanilla;
       double finalValNew = I_ortho(h, k, l, xsectionNewF, a, b, c)/normalisationNew;

       xsectionAngles[i] = finalVal * nentries;
       newXsectionAngles[i] =  finalValNew * nentries; 
       errsY[i] = finalVal * error;
       errsNewY[i] = finalValNew * error; 
       errsX[i] = 0.0;  
       errsNewX[i] = 0.0;
  
   } 
   
  TGraphErrors* vanillaGr = new TGraphErrors(grAngles, xsectionAngles, errsX, errsY); 
  TGraphErrors* newGr = new TGraphErrors(grAngles, newXsectionAngles, errsNewX, errsNewY); 

  newGr->SetLineColor(kRed); 
  newGr->SetMarkerColor(kRed); 
  //newGr->SetLineWidth(4); 

  vanillaGr->SetMarkerStyle(kPlus); 
  //vanillaGr->SetLineWidth(4); 
  newGr->SetMarkerStyle(kMultiply); 
   
  //vanillaGr->SetMarkerSize(4); 
  //newGr->SetMarkerSize(4); 

  
 //newGr->GetYaxis()->SetRange(long(3700e9), long(3750e9));
 // vanillaGr->GetYaxis()->SetRange(long(3700e9), long(3750e9)); 
  TCanvas* lambda = new TCanvas("l", "l", 1000, 600);
  wavelengths->Draw("E1");
  TCanvas* br = new TCanvas("br", "br", 1000, 600); 
  braggAnglesOriginal->Draw("E1");
  TCanvas* vanilla =  new TCanvas("van", "van", 1000, 600); 
   //vanillaGr->Draw("A*");
  //newGr->Draw("A* same");  
  vanillaGr->Draw("AP"); 
  vanillaGr->GetXaxis()->SetTitle("Angles / Radians"); 
  vanillaGr->GetYaxis()->SetTitle("Number of Neutrons"); 
  newGr->Draw("P"); 
  
  double param = 0.0;  
  for (int i = 0; i < newXsectionAngles.GetNoElements(); i++) {
    if (newXsectionAngles[i] == 0) continue;
    param = newXsectionAngles[i]; 
    break;  
 }


  std::cout << "Lambda Yukawa: " << mu << " eV, Height: " << param << std::endl;
}
