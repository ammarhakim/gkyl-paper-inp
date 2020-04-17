/** Compute exact solution to Euler Riemann problem. */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <limits>

class NameValuePair {
  public:
/** 
 * @param fname File containing data.
 */
    NameValuePair(const std::string& fname) {
      std::string line;
// open file for reading
      std::ifstream file(fname.c_str());
      if (!file)
        throw std::runtime_error("NameValuePair::NameValuePair: Unable to open input file");
// loop over each line in file
      while (std::getline(file, line)) {
        unsigned commentLoc = line.find('#');
        std::string myLine = line.substr(0, commentLoc);
        if (myLine != "") {
// only work on non-blank line: horribly ugly, but works
          unsigned eqLoc = myLine.find('=');
          std::string lhs, lhsStr, valueStr;
          lhsStr = myLine.substr(0, eqLoc);
          valueStr = myLine.substr(eqLoc+1, std::string::npos);
          std::istringstream lhsToken(lhsStr);
          lhsToken >> lhs;

          double value;
          std::istringstream tokens(valueStr);
          tokens >> value;
// insert it into map
          nvmap[lhs] = value;
        }
      }
    }

/** 
 * @param name Name to lookup
 * @return value associated with name
 */
    double getValue(const std::string& name) const {
      std::map<std::string, double>::const_iterator itr =
        nvmap.find(name);
      if (itr != nvmap.end())
        return itr->second;
      throw std::runtime_error("NameValuePair::getValue: not found");
    }

/**
 * @param name Name to lookup
 * @return true if it exists, false otherwise
 */
    bool hasValue(const std::string& name) const {
      std::map<std::string, double>::const_iterator itr =
        nvmap.find(name);
      if (itr != nvmap.end())
        return true;
      return false;
    }

  private:
    std::map<std::string, double> nvmap;
};

/** Problem state */
class ProblemState {
  public:
/** Density, velocity, pressure, sound-speed on left */
    double dl, ul, pl, cl;
/** Density, velocity, pressure, sound-speed on right */
    double dr, ur, pr, cr;
/** Normalising constant */
    double pscale;
/** Lower and upper bounds */
    double lower, upper;
/** Domain length, time at which solution is needed */
    double domlen, tEnd;
/** Ratio of specific heat */
    double gas_gamma;
/** Location of discontinuity */
    double disLoc;
};

/** Solution of Riemann problem */
class Solution {
  public:
/**
 * @param ncells Number of cells in domain.
 */
    Solution(unsigned ncell)
      : ncell(ncell),
        density(ncell),
        velocity(ncell),
        pressure(ncell),
        internalEnergy(ncell)
    {
    }
    
/** Number of cells */
    unsigned ncell;
/** Density, velocity, pressure and internal-energy */
    std::vector<double> density, velocity, pressure, internalEnergy;
};

void
showInput(const ProblemState& ps) {
  std::cout << "gamma " << ps.gas_gamma << std::endl;

  std::cout << "dl " << ps.dl << std::endl;
  std::cout << "ul " << ps.ul << std::endl;
  std::cout << "pl " << ps.pl << std::endl;

  std::cout << "dr " << ps.dr << std::endl;
  std::cout << "ur " << ps.ur << std::endl;
  std::cout << "pr " << ps.pr << std::endl;
}

void
prefun(const ProblemState& ps, double& F, double& FD,
  double P, double DK, double PK, double CK) {

  double PRATIO, QRT, AK, BK;
  double gas_gamma = ps.gas_gamma;
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

  if (P<=PK)
  {
// rarefaction wave
    PRATIO = P/PK;
    F = G4*CK*(pow(PRATIO,G1) - 1.0);
    FD = (1.0/(DK*CK))*pow(PRATIO,-G2);
  }
  else
  {
// shock wave
    AK  = G5/DK;
    BK  = G6*PK;
    QRT = sqrt(AK/(BK + P));
    F   = (P - PK)*QRT;
    FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT;
  }
}

/**
 * @param ps Problem state
 * @return Guess for pressure in star region
 */
double
guessp(const ProblemState& ps) {
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  double gas_gamma = ps.gas_gamma;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

  double CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ,
    PTL, PTR, QMAX, QUSER, UM;

  QUSER = 2.0;

  CUP  = 0.25*(DL + DR)*(CL + CR);
  PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP;
  PPV  = std::max(0.0, PPV);
  PMIN = std::min(PL,  PR);
  PMAX = std::max(PL,  PR);
  QMAX = PMAX/PMIN;

  if ((QMAX <= QUSER) && ((PMIN<=PPV) && (PPV<=PMAX)))
  {
    //std::cout << "First if" << std::endl;
    PM = PPV;
  }
  else 
  {
    if (PPV<PMIN) 
    {
      //std::cout << "Second if" << std::endl;
      PQ  = pow(PL/PR,G1);
      UM  = (PQ*UL/CL + UR/CR + G4*(PQ - 1.0))/(PQ/CL + 1.0/CR);
      PTL = 1.0 + G7*(UL - UM)/CL;
      PTR = 1.0 + G7*(UM - UR)/CR;
      PM  = 0.5*(PL*pow(PTL,G3) + PR*pow(PTR,G3));
    }
    else 
    {
      //std::cout << "Second else" << std::endl;
      GEL = sqrt((G5/DL)/(G6*PL + PPV));
      GER = sqrt((G5/DR)/(G6*PR + PPV));
      PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER);
    }
  }
  return PM;
}

void
starpu(const ProblemState& ps, double &pm, double& um) {
  double CHANGE, FL, FLD, FR, FRD, P, POLD, 
    PSTART, TOLPRE, U, UDIFF, PSCALE;

  TOLPRE = 1.0e-6;
  unsigned NRITER = 20;

  double gas_gamma = ps.gas_gamma;
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;
  PSCALE = 1.0;

// compute initial guess
  PSTART = guessp(ps);
  std::cout << "Initial guess " << PSTART << std::endl;

  POLD = PSTART;
  UDIFF = UR-UL;

  bool converged = false;
  for (unsigned i=1; i<NRITER; ++i) {
    prefun(ps, FL, FLD, POLD, DL, PL, CL);
    prefun(ps, FR, FRD, POLD, DR, PR, CR);
    P = POLD - (FL + FR + UDIFF)/(FLD + FRD);
    CHANGE = 2.0*fabs((P - POLD)/(P + POLD));
    std::cout << "Iteration " << i << " Change " << CHANGE << std::endl;
    if (CHANGE<=TOLPRE) {
      converged = true;
      break;
    }
    if (P<0.0) P = TOLPRE;
    POLD = P;
  }
  if (converged)
    std::cout << "Newton iteration converged ..." << std::endl;
  else {
    std::cout << "Newton iteration did not converge" << std::endl;
    exit(1);
  }
// compute velocity in star region
  U = 0.5*(UL + UR + FR - FL);
  std::cout << "Converged pressure " << P << std::endl;
  std::cout << "Converged velocity " << U << std::endl;
  pm = P;
  um = U;
}

void
sample(const ProblemState& ps, double PM, double UM, double S,
  double& D, double& U, double& P) {
  double gas_gamma = ps.gas_gamma;
// compute constants related to gamma
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

  double C, CML, CMR, PML, PMR, SHL, SHR, SL, SR, STL, STR;

  //std::cout << "S = " << S << std::endl;
  //std::cout << "UM = " << UM << std::endl;
  if (S<=UM)
  {
// left of contact discontinuity
    //std::cout << "Left of contact .." << std::endl;

    if (PM<=PL)
    {
      //std::cout << "Left rarefaction ...." << std::endl;

// left rarefaction
      SHL = UL-CL;
      if (S<=SHL)
      {
        //std::cout << "Left data state!" << std::endl;

// left data state
        D = DL;
        U = UL;
        P = PL;
      }
      else 
      {
        CML = CL*pow(PM/PL, G1);
        STL = UM-CML;

        if (S>STL)
        {
          //std::cout << "Star left state!" << std::endl;

// Star left state
          D = DL*pow(PM/PL, 1/GAMMA);
          U = UM;
          P = PM;
        } 
        else 
        {
          //std::cout << "Inside left fan!" << std::endl;
// inside left fan
          U = G5*(CL + G7*UL + S);
          C = G5*(CL + G7*(UL - S));
          D = DL*pow(C/CL, G4);
          P = PL*pow(C/CL, G3);
        }
      }
    }
    else
    {
      //std::cout << "Left shock ..." << std::endl;
// left shock
      PML = PM/PL;
      SL = UL - CL*sqrt(G2*PML + G1);

      if (S<=SL)
      {
        //std::cout << "Left data state!" << std::endl;
// point is left data state
        D = DL;
        U = UL;
        P = PL;
      }
      else
      {
        //std::cout << "Star left state!" << std::endl;
// point is star left state
        D = DL*(PML + G6)/(PML*G6 + 1.0);
        U = UM;
        P = PM;
      }
    }
  }
  else
  {
// right of contact discontinuity
    //std::cout << "Right of contact discontinuity ...." << std::endl;
    if (PM>PR)
    {
// right shock
      //std::cout << "Right shock ..." << std::endl;

      PMR = PM/PR;
      SR  = UR + CR*sqrt(G2*PMR + G1);

      if (S>=SR)
      {
        //std::cout << "Right data state!" << std::endl;
// right data state
        D = DR;
        U = UR;
        P = PR;
      }
      else
      {
        //std::cout << "Star right state!" << std::endl;
// right star state
        D = DR*(PMR + G6)/(PMR*G6 + 1.0);
        U = UM;
        P = PM;
      }
    }
    else
    {
      //std::cout << "Right rarefaction ..." << std::endl;
// right rarefaction
      SHR = UR+CR;

      if (S>=SHR)
      {
        //std::cout << "Right data state!" << std::endl;
// right data state
        D = DR;
        U = UR;
        P = PR;
      }
      else
      {
        CMR = CR*pow(PM/PR, G1);
        STR = UM + CMR;

        if (S<=STR)
        {
          //std::cout << "Star right state!" << std::endl;
// star right state
          D = DR*pow(PM/PR, 1.0/GAMMA);
          U = UM;
          P = PM;
        }
        else
        {
          //std::cout << "Inside left fan!" << std::endl;
// inside left fan
          U = G5*(-CR + G7*UR + S);
          C = G5*(CR - G7*(UR - S));
          D = DR*pow(C/CR, G4);
          P = PR*pow(C/CR, G3);
        }
      }
    }
  }
}

void
sampleWithVacuum(const ProblemState& ps, double S, double& D, double& U, double& P) {
  double gas_gamma = ps.gas_gamma;
// compute constants related to gamma
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

// assume vacuum is to the right (JUST FOR NOW)
  double sstar = UL + 2*CL/G8;
  double relVel = UL-CL;
  if (S<=relVel)
  {
// left of rarefaction
    D = DL;
    U = UL;
    P = PL;
  }
  else
  {
    if ((relVel<S) && (S<=sstar))
    {
// inside rarefaction fan
      D = DL*pow(G5 + G6/CL*(UL-S), G4);
      U = G5*(CL + G7*UL + S);
      P = PL*pow(G5 + G6/CL*(UL-S), G3);
    }
    else
    {
// inside vacuum
      D = 0.0;
      U = 0.0;
      P = 0.0;
    }
  }
}

void
writeSolution(const ProblemState& ps, const std::string& outPrefix, const Solution& solution) {
  double xcell;
  double dx = ps.domlen/solution.ncell;

// write density
  std::string dhFileNm(outPrefix);
  dhFileNm += "-density.txt";
  std::ofstream dhFile(dhFileNm.c_str());
  dhFile << "## X [m] Density [kg/m^3]" << std::endl;
  for (unsigned i=0; i<solution.ncell; ++i) {
    xcell = ps.lower + (i+0.5)*dx;
    dhFile << xcell << " " << solution.density[i] << std::endl;
  }

// write pressure
  std::string prFileNm(outPrefix);
  prFileNm += "-pressure.txt";
  std::ofstream prFile(prFileNm.c_str());
  prFile << "## X [m] Pressure [Pa]" << std::endl;
  for (unsigned i=0; i<solution.ncell; ++i) {
    xcell = ps.lower + (i+0.5)*dx;
    prFile << xcell << " " << solution.pressure[i] << std::endl;
  }

// write velocity
  std::string ulFileNm(outPrefix);
  ulFileNm += "-velocity.txt";
  std::ofstream ulFile(ulFileNm.c_str());
  ulFile << "## X [m] Velocity [m/s]" << std::endl;
  for (unsigned i=0; i<solution.ncell; ++i) {
    xcell = ps.lower + (i+0.5)*dx;
    ulFile << xcell << " " << solution.velocity[i] << std::endl;
  }

// write internal energy
  std::string ieFileNm(outPrefix);
  ieFileNm += "-internal-energy.txt";
  std::ofstream ieFile(ieFileNm.c_str());
  ieFile << "## X [m] Internal-energy [J]" << std::endl;
  for (unsigned i=0; i<solution.ncell; ++i) {
    xcell = ps.lower + (i+0.5)*dx;
    ieFile << xcell << " " << solution.internalEnergy[i] << std::endl;
  }
}

void
exactEulerRp(const ProblemState& ps, Solution& solution) {
  double g1, g2, g3, g4, g5, g6, g7, g8;

  double gas_gamma = ps.gas_gamma;
// compute constants related to gamma
  g1 = (gas_gamma - 1)/(2*gas_gamma);
  g2 = (gas_gamma + 1)/(2*gas_gamma);
  g3 = 2*gas_gamma/(gas_gamma - 1);
  g4 = 2/(gas_gamma - 1);
  g5 = 2/(gas_gamma + 1);
  g6 = (gas_gamma - 1)/(gas_gamma + 1);
  g7 = (gas_gamma - 1)/2;
  g8 = gas_gamma - 1;

// check if a vacuum is generated
  if (g4*(ps.cl+ps.cr) <= (ps.ur-ps.ul)) {
    std::cout << "Initial conditions will lead to vacuum. Aborting ..." << std::endl;
    exit(1);
  }

// compute pressure and velocity in "star" region
  double pm, um;
  starpu(ps, pm, um);

// cell spacing
  double dx = ps.domlen/solution.ncell;
// compute solution at each grid point
  for (unsigned i=0; i<solution.ncell; ++i) {
    double xpos = ps.lower + (i+0.5)*dx;
    double s = (xpos-ps.disLoc)/ps.tEnd;
// compute solution at (x,t) = (xpos-ps.disLoc, ps.tEnd)
    double dsol, usol, psol;
    sample(ps, pm, um, s, dsol, usol, psol);
    //std::cout << xpos << " " << dsol << std::endl;
// copy solution into array
    solution.density[i] = dsol;
    solution.velocity[i] = usol;
    solution.pressure[i] = psol;
    solution.internalEnergy[i] = psol/dsol/g8;
  }
}

void
exactEulerRpWithVacuum(const ProblemState& ps, Solution& solution) {

  double gas_gamma = ps.gas_gamma;
  double g8 = gas_gamma - 1;

// cell spacing
  double dx = ps.domlen/solution.ncell;
// compute solution at each grid point
  for (unsigned i=0; i<solution.ncell; ++i) {
    double xpos = ps.lower + (i+0.5)*dx;
    double s = (xpos-ps.disLoc)/ps.tEnd;
// compute solution at (x,t) = (xpos-ps.disLoc, ps.tEnd)
    double dsol, usol, psol;
    sampleWithVacuum(ps, s, dsol, usol, psol);
    //std::cout << xpos << " " << dsol << std::endl;
// copy solution into array
    solution.density[i] = dsol;
    solution.velocity[i] = usol;
    solution.pressure[i] = psol;
    if (dsol < std::numeric_limits<double>::epsilon())
      solution.internalEnergy[i] = 0.0;
    else
      solution.internalEnergy[i] = psol/dsol/g8;
  }
}

int
main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Usage::" << std::endl;
    std::cout << " exacteulerrp <input-file> <output-prefix>" << std::endl;
    std::cout << "  input-file: Name of input file" << std::endl;
    std::cout << "  output-prefix: Prefix for output files" << std::endl;

    exit(1);
  }
// name of input file and output prefix
  std::string inFile(argv[1]);
  std::string outPrefix(argv[2]);

  ProblemState ps;
  unsigned ncell;

// initialize problem state from input file
  NameValuePair nvpair(inFile);
  ps.dl = nvpair.getValue("dl");
  ps.ul = nvpair.getValue("ul");
  ps.pl = nvpair.getValue("pl");
  ps.dr = nvpair.getValue("dr");
  ps.ur = nvpair.getValue("ur");
  ps.pr = nvpair.getValue("pr");
  ps.pscale = 1.0;
  ps.disLoc = nvpair.getValue("location");
  ps.tEnd = nvpair.getValue("tEnd");
  ps.gas_gamma = nvpair.getValue("gamma");
  ncell = (unsigned) nvpair.getValue("ncell");

  if (nvpair.hasValue("lower")) {
    ps.lower = nvpair.getValue("lower");
    ps.upper = nvpair.getValue("upper");
    ps.domlen = ps.upper - ps.lower;
  }
  else {
    ps.domlen = nvpair.getValue("length");
    ps.lower = 0.0;
    ps.upper = ps.domlen;
  }

// compute sound speeds in each region
  if (ps.dl != 0)
    ps.cl = std::sqrt(ps.gas_gamma*ps.pl/ps.dl);
  else
    ps.cl = 0.0;

  if (ps.dr != 0)
    ps.cr = std::sqrt(ps.gas_gamma*ps.pr/ps.dr);
  else
    ps.cr = 0.0;

//  showInput(ps);

// allocate arrays for solution
  Solution solution(ncell);

// compute solution
  if ((ps.dl != 0.0) && (ps.dr != 0.0))
    exactEulerRp(ps, solution);
  else if (ps.dr == 0)
    exactEulerRpWithVacuum(ps, solution);
  else
    exactEulerRpWithVacuum(ps, solution);    

// write solution to text files
  writeSolution(ps, outPrefix, solution);

  return 1;
}
