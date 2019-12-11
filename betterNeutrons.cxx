#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom3.h>
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


/* Global Parameters for both vanilla and new cross section. */
  // Parameters
  int Z = 54;
  // double wavelength = 1.0; 
  double lattice_spacing = 100.0;
  //double mu = 200.0e-15; // eV
  double q_0 = 6.86; // A-1
  
  double b_F = -1.468e-3; // Foldy scattering length, fm
  double r_N_square = -0.1161; // PRL, fm^2
  double a_0 = 5.2917721067e4; // Scipy value, fm
  double b_I = (1/3)*(1.6726219e-27/9.10938356e-31)*(r_N_square/a_0); // fm Intrinsic scattering length
  double b_c = 4.92; // fm, Xe coherent scattering length from NIST database

  // Dimensionless constants from scattering lengths
  double x_em = -Z*(b_F + b_I)/b_c; // Need - sign to make a sensible em x-section, also Sears Eq 3.19 has - out front


/* Utility Methods */ 

double square(double x) {
  return pow(x, 2); 
}

std::complex<double> complexExp(double x) {
  const std::complex<double> i(0, 1); 
  return std::exp(-2*M_PI*i*x); 
}

// Need to define special contain method because of floating point errors
bool contains(long double x, std::vector<long double> vec, double tolerance=1e-3) {
  for (double num : vec) {
    if (fabs(num - x) < tolerance) return true; 
  }
  return false; 
}

// Miller Permutations
std::vector<std::vector<int>> getCombinations(int n) {
   std::vector<int> set;
   std::vector<long double> millers;
 
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
        long double quantity = pow(millerH, 2) + pow(millerK, 2) + pow(millerL, 2); 
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

std::vector<std::vector<int>> getMillerCombinations(int n, double a, double b, double c) {
  std::vector<std::vector<int>> perms; 
  std::vector<long double> millers;
  std::vector<int> set;
  
  for (int i = 0; i < n; i++) {
    set.emplace_back(i);
  }
  
  for (int h = 0; h < n; h++) {
    for (int k = 0; k < n; k++) {
      for (int l = 0; l < n; l++) {
        if (h == 0 && k == 0 && l == 0) continue;
        int millerH = set[h]; 
        int millerK = set[k];
        int millerL = set[l];

        long double quantity = square(millerH/a) + square(millerK/b) + square(millerL/c);
        bool alreadySeen = contains(quantity, millers); 
        //bool alreadySeen = std::find(millers.begin(), millers.end(), quantity) != millers.end();
        if (alreadySeen) continue; 
        std::vector<int> perm = {millerH, millerK, millerL};
        perms.emplace_back(perm); 
        millers.emplace_back(quantity);
      }
    }
  }
  return perms; 
}

// Neutron Wavelength Picker
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

// Form Factor function from PRL
double f(double q) {
  return 1/pow(1+3*(pow(q/q_0, 2)), 0.5); 
}

// Vanilla Cross Section, needs theta and wavelength
double vanillaCrossSection(Double_t *args, Double_t *par) {
  double scatAngle = args[0]/2.0;
  double k = 2*M_PI/args[1]; 
  double constTerm = 1;
  double emTerm = 2*x_em*(1-f(2*k*sin(scatAngle)));
  //return scatAngle; 
 return par[0]*2*M_PI*pow(b_c,2)*(constTerm + emTerm); 
}

// New cross section function, needs theta, wavelength, mu in eV, and g in GeV^-2
double newCrossSection(Double_t *args, Double_t *par) {
  double mu = par[2]; // eV
  double g_squared = par[1]; // GeV^-2 
  //std::cout << "mu: " << mu << " g^2: " << g_squared << std::endl;
  double mu_gev = mu*1e-9; // GeV
  double xe_nuclear = 122; // GeV c^-2
  double neutron_nuclear = 0.93957; // GeV c^-2
  double b_c_ev = 25; // GeV^-1 
  double m_n = 0.93957; // ?? weird m_n term in GeV, neutron mass?? Discussed with Alex he thinks it's neutron mass also which 
                       // is weird because then why have it in Q_1 as well? 
  double x_new = m_n*g_squared*xe_nuclear*neutron_nuclear/(2*M_PI*b_c_ev*pow(mu_gev,2)); 
  
  double h_bar = 6.58211956e-16; // eV s
  double c_angstrom_per_sec = 2.998e+18; // speed of light in angstrom per second
  double scatAngle = args[0]/2.0;
  
  double k = 2*M_PI/par[3]; 
  double q = 2*k*sin(scatAngle); 
  double constTerm = 1;
  double emTerm = 2*x_em*(1-f(q));  
  double mu_sq = pow(mu, 2); 
  double q_transfer_sq = pow(q*h_bar*c_angstrom_per_sec, 2);
  double newTerm = 2*x_new * mu_sq / (mu_sq + q_transfer_sq); 
  //std::cout << "NEW TERM: " << newTerm << std::endl;
  return par[0]*2*M_PI*pow(b_c,2)*(constTerm + emTerm + newTerm); 
 
}

double squaredXnewSection(Double_t* args, Double_t *par) {
  return square(newCrossSection(args, par));
}

/* Lattice Physics  */

// Primitive cubic structure factor
std::complex<double> F_primitive_cubic(int h, int k, int l, double f, double a, double b, double c) {
  return f; 
}

// Orthorhombic cubic structure factor
std::complex<double> F_tetragonal_cubic(int h, int k, int l, double f, double a, double b, double c) {
  std::complex<double> p0 = complexExp(0.0); 
  std::complex<double> p1 = complexExp(a*h); 
  std::complex<double> p2 = complexExp(b*k); 
  std::complex<double> p3 = complexExp(c*l);
  //return f*(p0 + p1 + p2 + p3);
  return f*p0; // Tetragonal structure factor same as primitive cubic (from Reddit: https://www.reddit.com/r/AskPhysics/comments/7cnbau/how_to_calculate_the_structure_factor_of_a/dpthn8g/) 
}

// Convert miller indices to angle using Bragg relation 
double hklThetaOrthorhombic(int h, int k, int l, double wavelength, double a, double b, double c) {
  double sinx = (wavelength/2.0) * pow(square(h/a) + square(k/b) + square(l/c), 0.5);
  double result = std::asin(sinx);
  if (!std::isnan(result)) {
    return result; 
  } 
  return -1; 
}

// Intensity from structure factor I_hkl = |F_hkl|^2
double I_ortho(int h, int k, int l, double f, double a, double b, double c) {
  return std::norm(F_primitive_cubic(h, k, l, f, a, b, c)); 
}

// Convert lattice output angle to wavelength using Bragg relation
double angleToWavelength(double theta, int h, int k, int l, double a, double b, double c) {
  double d_hkl_reciprocal = pow(square(h/a) + square(k/b) + square(l/c), 0.5);  
  return 2*sin(theta)*1/d_hkl_reciprocal; 
} 


double doNeutronSimulation(int nentries, int miller, double mu, double g_squared, double a, double b, double c) {
  
    gStyle->SetOptStat(0);  
  //int nentries = int(1e5); 
  //double mu = 200.0; // eV
  //double g_squared = 1e-16; //GeV^-2
    // -------------------- LATTICE -------------------- 
  
  TH1F* braggAnglesOriginal = new TH1F("bragg", "bragg", 100, 0.0, 0.0); 
  braggAnglesOriginal->SetLineColor(kBlack);
  braggAnglesOriginal->GetXaxis()->SetTitle("Angle / radians"); 
  braggAnglesOriginal->GetYaxis()->SetTitle("Number of neutrons"); 
  TH1F* braggAnglesNew = new TH1F("braggnew", "braggnew", 100, 0.0, 1.57); 
  braggAnglesNew->SetLineColor(kRed);
  braggAnglesNew->GetXaxis()->SetTitle("Angle / radians"); 
  braggAnglesNew->GetYaxis()->SetTitle("Number of neutrons"); 
  TH1F* vanillaCrystal = new TH1F("vanillaCrystal", "Plain Crystal Distribution", 50, 0.0, 1.57);   
  
  TH1F* wavelengths = new TH1F("wavelengths", "Wavelength Distribution", 100, 0.0, 5.0);   
  //double a = 10.0; // Angstrom
  //double b = 10.0; // Angstrom
  //double c = 50.0; // Angstrom
  
  std::vector<double> angles;
  std::vector<double> vanillaDistribution; 
  std::vector<double> newDistribution; 
  std::vector<double> millerPlanes;
  
  TRandom3* rand = new TRandom3(0); 
  rand->SetSeed(0);  
  // Not normalising these because then lose units (discussion with Alex from 25/11/19) 
 // double momRatio = square(k/q_0); 
 // TF1* normalisedVanilla = new TF1("normalVanilla", vanillaCrossSection, 0.0, 0.5, 2);
 // double vanillaIntegral = 2*M_PI*square(b_c) * (2 * (1 + 2*x_em) - 2*x_em*((2*pow(12*momRatio+1, 0.5) - 1)/(6*momRatio)));
 // normalisedVanilla->SetParameter(0, 1);//1/vanillaIntegral);
  
  //TF1* newXSection = new TF1("newXSection", newCrossSection, 0.0, 0.5, 2); 
 // newXSection->SetParameter(0, 1); 

  typedef std::tuple<int, int, int> miller_plane_t; 
  
  std::unordered_map<double, miller_plane_t> angleToIndices; 
  std::vector<std::vector<int>> combos = getMillerCombinations(miller, a, b, c); 
  //std::vector<std::vector<int>> combos = getCombinations(miller); 
  for (size_t i = 0; i < nentries; i++) {

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
    if (angle <= 0) continue; // Filter out incompatible events
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
   
   // Get Number of bins with nonzero content for vectors 
   // Normalisation Constant Calculation (only for primitive cubic, can change it for general F_hkl later) 

   TF1* normalisedVanilla = new TF1("normalVanilla", vanillaCrossSection, 0.0, 1.57, 1);
   TF1* newXSection = new TF1("newXSection", newCrossSection, 0.0, 1.57, 4);
   newXSection->SetParameter(1, mu); 
   newXSection->SetParameter(2, g_squared); 
   
   double normalisationVanilla = 0.0; 
   double normalisationNew = 0.0; 
   int numAngles = 0;  
  for (size_t i = 0; i < nbins; i++) {
     double centre = braggAnglesOriginal->GetBinCenter(i); 
     double content = braggAnglesOriginal->GetBinContent(i); 
     if (content == 0) continue; 
    
       miller_plane_t plane = angleToIndices[centre];
       int h = std::get<0>(plane); 
       int k = std::get<1>(plane); 
       int l = std::get<2>(plane); 
              
       double wavelength = angleToWavelength(centre, h, k, l, a, b, c); 
       double args[2] = {centre, wavelength};
       double argsNew[3] = {centre}; 
       normalisedVanilla->SetParameter(0, 1.0); 
       double contribution = normalisedVanilla->EvalPar(args); 

       newXSection->SetParameter(0, 1.0);
       newXSection->SetParameter(1, mu); 
       newXSection->SetParameter(2, g_squared);
       double pars[4] = {1.0, g_squared, mu, wavelength};
       double newContribution = newXSection->EvalPar(argsNew, pars);

       normalisationVanilla += I_ortho(h, k, l, contribution, a, b, c) * content; 
       normalisationNew += I_ortho(h, k, l, newContribution, a, b, c) * content;
       numAngles++; 
   }
    
   TVectorD grAngles(numAngles); 
   TVectorD xsectionAngles(numAngles); 
   TVectorD errsX(numAngles); 
   TVectorD errsY(numAngles); 

   TVectorD newXsectionAngles(numAngles); 
   TVectorD errsNewX(numAngles); 
   TVectorD errsNewY(numAngles);
   
   double height = -120;
   int vectorIndex = 0; 
   for (size_t i = 0; i < nbins; i++) {
       double centre = braggAnglesOriginal->GetBinCenter(i); 
       double error = braggAnglesOriginal->GetBinError(i);
       int content = braggAnglesOriginal->GetBinContent(i); 

       if (content == 0.0) continue;  
       grAngles[vectorIndex] = centre;  
       
       miller_plane_t plane = angleToIndices[centre];
       int h = std::get<0>(plane); 
       int k = std::get<1>(plane); 
       int l = std::get<2>(plane); 
              
       double wavelength = angleToWavelength(centre, h, k, l, a, b, c); 
       //std::cout << "Wavelength " << wavelength << std::endl;
       //double wavelength = 1.0; 

       double args[2] = {centre, wavelength};
       double argsNew[3] = {centre}; 
       double pars[4] = {1.0, g_squared, mu, wavelength}; 
       double xsectionF = normalisedVanilla->EvalPar(args);
       double xsectionNewF = newXSection->EvalPar(argsNew, pars);
      
       double finalVal = I_ortho(h, k, l, xsectionF, a, b, c) * (nentries/normalisationVanilla);
       double finalValNew = I_ortho(h, k, l, xsectionNewF, a, b, c) * (nentries/normalisationNew);
       // here nentries/normalisation since A * sum_i (F_j^2 * N_i) = N so A = N/(sum_i ...) 

       if (finalVal == 0 || finalValNew == 0) continue; 
       xsectionAngles[vectorIndex] = finalVal * content;
       newXsectionAngles[vectorIndex] =  finalValNew * content; 
       errsY[vectorIndex] = finalVal * error;
       errsNewY[vectorIndex] = finalValNew * error; 
       errsX[vectorIndex] = 0.0;  
       errsNewX[vectorIndex] = 0.0;
       vectorIndex++; 
  
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

  /* Suppress drawing for now. */  

 //newGr->GetYaxis()->SetRange(long(3700e9), long(3750e9));
 // vanillaGr->GetYaxis()->SetRange(long(3700e9), long(3750e9));
  
  TCanvas* lambda = new TCanvas("l", "l", 1000, 600);
  wavelengths->Draw("E1");
  TCanvas* br = new TCanvas("br", "br", 1000, 600); 
  braggAnglesOriginal->Draw("E1");
  TCanvas* vanilla =  new TCanvas("van", "van", 1000, 600); 
  //vanillaGr->Draw("A*");
  //newGr->Draw("A*");  
  vanillaGr->Draw("AP"); 
  vanillaGr->GetXaxis()->SetTitle("Angles / Radians"); 
  vanillaGr->GetYaxis()->SetTitle("Number of Neutrons"); 
  newGr->Draw("P"); 
 

  std::cout << "Angles: ";
  for (int i = 0; i < newGr->GetN(); i++) { 
    std::cout << newGr->GetX()[i] << ", "; 
  } 
  std::cout << "\n Intensities "; 
  for (int i = 0; i < newGr->GetN(); i++) { 
    std::cout << newGr->GetY()[i] << ", "; 
  } 
  std::cout << "\n Errors "; 
  for (int i = 0; i < newGr->GetN(); i++) { 
    std::cout << newGr->GetEY()[i] << ", "; 
  } 
  std::cout << "\n"; 
//  TCanvas* test = new TCanvas("test", "test", 1000, 600); 
 
//  vanillaGr->Draw("AP"); 
  // Fitting of graph 
 
  TF1* newXsectionFit = new TF1("newXsectFit", squaredXnewSection, 0.0, 1.0, 4);
  //newXsectionFit->FixParameter(1, 1e-16);  
  newXsectionFit->FixParameter(2, mu); // mu
  newXsectionFit->SetParameter(0, nentries); // normalisation
  newXsectionFit->FixParameter(3, 1.0); // Wavelength
  newXsectionFit->SetParameter(1, g_squared); // gsquare
  //newXsectionFit->SetParLimits(1, 5e-17, 5e-16);
  newGr->Fit("newXsectFit", "EM"); 
  TF1* fitFunc = newGr->GetFunction("newXsectFit");
  //TCanvas* fit = new TCanvas("fit", "fit", 1000, 600); 
  //fitFunc->Draw("C");  
  double g_fit = fitFunc->GetParameter(1);
  std::cout << "Normalisation " << fitFunc->GetParameter(0) << std::endl;  
  //delete braggAnglesOriginal;
  delete braggAnglesNew;
  delete vanillaCrystal;
  delete wavelengths; 
  delete rand;
  return g_fit;  

}
void betterNeutrons() {
  size_t nIter = 100000000; 
  int miller = 5; 
  double mu = 200.0;
  double g_squared = 7.5e-16; 
  double sameDim = 10.0; 
  double prismDim = 10.0;
  doNeutronSimulation(nIter, miller, mu, g_squared, sameDim, sameDim, prismDim); 
 
  /*
  std::vector<double> couplings;
  couplings.emplace_back(0.0);
  
  THStack* stack = new THStack("stk", "stk"); 
  
  for (int j = 0; j < couplings.size(); j++) {
    double coupling = couplings[j];  
    TH1F* g_fit_histo = new TH1F("gFit", "gFit", 1000, 0.0e-12, 0.0e-12); 
    g_fit_histo->SetLineColor(j+1); 
    for (int i = 0; i < 1000; i++) {
      double g_squareFit = doNeutronSimulation(nIter, miller, mu, coupling); // n_events, miller index, eV, GeV^-2
      std::cout << "Done iteration " << i << std::endl;
      g_fit_histo->Fill(g_squareFit);
    }
    stack->Add(g_fit_histo);
  }
  TCanvas* yeet = new TCanvas("feldman", "feldman", 1000, 600); 
  stack->Draw("E1");*/ 
}








