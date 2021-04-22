/*
 *  Copyright 2021 Chan Beom Park <cbpark@gmail.com>
 *
 *  This code calculates various observables for the e+ e- --> tau+ tau- process
 *  using YAM2 (https://github.com/cbpark/YAM2). The decay topology is
 *
 *  X --> Y1 + Y2 --> vis1(p1) inv1(k1) + vis2(p2) inv2(k2)
 *
 *  where X corresponds to sqrt(s) and Y's are the decaying particle (tau).
 *  vis and inv are visible and invisible particles, respectively.
 */

#include <algorithm>  // std::max, std::minmax
#include <cmath>      // std::hypot
#include <cstdlib>    // std::atof
#include <iostream>
#include <optional>  // std::optional

// ROOT
#include "Math/LorentzVector.h"
// #include "Math/Vector4D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TVector2.h"

// YAM2
#include "YAM2/yam2.h"

using std::cerr;
using std::cout;

/// we use the recommended LorentzVector rather than legacy TLorentzVector.
/// we could set 'using LorentzVector = ROOT::Math::XYZTVector;'
using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;

const auto APPNAME{"bditau"};

/*
 *  Constants
 */
/// the longitudinal momentum of the di-tau system (it's 0 for the CMS,
/// while 7 - 4 = 3 GeV in the LAB frame).
const double PZTOT = 0.0;

/// sqrt(s) of CMS.
const double SQRTS = 10.579;

/// the four-momentum of the di-tau system.
const LorentzVector CMS{0.0, 0.0, PZTOT, SQRTS};

/// zero four-momentum (0, 0, 0, 0) for convenience to calculate M2.
const yam2::FourMomentum ZERO;

/*
 *  Functions for calculating the event variables.
 */
/// E_{miss} = sqrt(ptmiss_x^2 + ptmiss_y^2 + pzmiss^2),
/// analogous to MET in hadron colliders.
double eMiss(const LorentzVector &p1, const LorentzVector &p2,
             const TVector2 &ptmiss) {
    const double pzmiss = PZTOT - p1.pz() - p2.pz();
    // std::hypot(x, y, z) = sqrt(x^2 + y^2 + z^2) since c++17.
    return std::hypot(ptmiss.Px(), ptmiss.Py(), pzmiss);
}

/// the recoil mass.
/// See Eq.(1) in https://arxiv.org/pdf/1506.05992.pdf
double mRecoil(const LorentzVector &p1, const LorentzVector &p2) {
    const auto p_miss = CMS - p1 - p2;
    return p_miss.mass();
}

/// ratio of min(p1, p2) and max(p1, p2).
///
/// Here 'p1' and 'p2' could be either visible or invisible particle momenta.
/// By definition, it is in between 0 and 1.
double momentumRatio(const LorentzVector &p1, const LorentzVector &p2) {
    const auto [pmin, pmax] = std::minmax({p1.P(), p2.P()});
    return pmin / std::max(pmax, 1.0e-10);
}

/// the M2 variable and its solution to the invisible particle momenta.
std::optional<yam2::M2Solution> getM2(const LorentzVector &p1,
                                      const LorentzVector &p2,
                                      const TVector2 &ptmiss,
                                      const yam2::Mass &m_inv) {
    // input kinematic configuration for M2.
    // see the above for the constants (ZERO, SQRTS, PZTOT).
    const auto input = yam2::mkInput({{p1.e(), p1.px(), p1.py(), p1.pz()},
                                      {p2.e(), p2.px(), p2.py(), p2.pz()}},
                                     {ZERO, ZERO}, {ptmiss.Px(), ptmiss.Py()},
                                     m_inv, {}, SQRTS, {PZTOT});
    // the latter arguments are tolerance and maximal evaluation number.
    return yam2::m2Cons(input, 1.0e-6, 1000);
}

/// converts 'yam2::FourMomentum' to 'LorentzVector'.
LorentzVector toLorentzVector(const yam2::FourMomentum &p) {
    return {p.px(), p.py(), p.pz(), p.e()};
}

/*
 *  The main function
 */
int main(int argc, char *argv[]) {
    if (!(argc == 3 || argc == 4)) {
        cout << "usage: ./bin/" << APPNAME
             << " <event.root> <output.root> [mInvisible]\n"
             << "  <event.root>: input ntuple file (required).\n"
             << "  <output.root>: output file to store the result (required).\n"
             << "  [mInvisible]: the input mass for invisible particles "
                "(optional, default = 0)\n";
        return 1;
    }

    // input root file.
    TFile infile(argv[1]);
    cout << APPNAME << ": the input file is " << infile.GetName() << '\n';

    // check the tree.
    auto keys = infile.GetListOfKeys();
    if (keys->GetSize() < 1) {
        cerr << APPNAME << ": the input has no tree.\n";
        infile.Close();
        return 1;
    }

    // get the tree.
    auto tree_name = infile.GetListOfKeys()->At(0)->GetName();
    cout << APPNAME << ": the name of the input tree is " << tree_name << '\n';
    auto event = dynamic_cast<TTree *>(infile.Get(tree_name));

    // the four-momentum of visible particles in the three-prong decay.
    double p1x, p1y, p1z, e1;
    // the three-momentum of visible particles in the one-prong decay.
    double p2x, p2y, p2z;
    // set each particle momentum component from the tree.
    event->SetBranchAddress("tau_3prong_px_CMS", &p1x);
    event->SetBranchAddress("tau_3prong_py_CMS", &p1y);
    event->SetBranchAddress("tau_3prong_pz_CMS", &p1z);
    event->SetBranchAddress("tau_3prong_E_CMS", &e1);
    event->SetBranchAddress("track_1prong_px_CMS", &p2x);
    event->SetBranchAddress("track_1prong_py_CMS", &p2y);
    event->SetBranchAddress("track_1prong_pz_CMS", &p2z);

    // the visible particle momenta.
    LorentzVector p1, p2;
    // the missing transverse momentum.
    TVector2 ptmiss;
    // the MAOS solutions to the invisible momenta.
    LorentzVector k1sol, k2sol;
    // the invisible particle mass.
    yam2::Mass m_inv;
    if (argc == 3) {
        m_inv = yam2::Mass{0.0};
    } else {
        m_inv = yam2::Mass{std::atof(argv[3])};
    }
    cout << APPNAME << ": the invisible mass is " << m_inv.value << '\n';

    // the event variables.
    double e_miss, m_recoil, m2, xi_p, xi_k;

    // ntuple for storing the variables.
    TNtuple vars{"vars", "event variables", "e_miss:m_recoil:m2:xi_p:xi_k"};

    // --------------------------------------------------------------------------
    // event loop
    const auto nentries = event->GetEntries();
    for (auto iev = 0; iev != nentries; ++iev) {
        event->GetEntry(iev);

        p1.SetPxPyPzE(p1x, p1y, p1z, e1);
        // the visible particle in the one-prong decay is assumed to be
        // massless.
        p2.SetPxPyPzE(p2x, p2y, p2z, std::hypot(p2x, p2y, p2z));
        ptmiss.Set(-(p1.Px() + p2.Px()), -(p1.Py() + p2.Py()));

#ifdef DEBUG
        cout << "--- (" << iev + 1 << ")\n";
        cout << "p1: " << p1 << ", p2: " << p2 << '\n';
        cout << "ptmiss: ";
        ptmiss.Print();
#endif

        e_miss = eMiss(p1, p2, ptmiss);

        m_recoil = mRecoil(p1, p2);

        xi_p = momentumRatio(p1, p2);

        const auto m2sol = getM2(p1, p2, ptmiss, m_inv);
        if (!m2sol) {
            cerr << APPNAME << ": failed to find minimum for M2! (" << iev
                 << ")\n";
            m2 = -1.0;
            xi_k = -1.0;
        } else {
            m2 = m2sol.value().m2();
            k1sol = toLorentzVector(m2sol.value().k1());
            k2sol = toLorentzVector(m2sol.value().k2());
            xi_k = momentumRatio(k1sol, k2sol);
        }

#ifdef DEBUG
        cout << "e_miss: " << e_miss << ", m_recoil: " << m_recoil
             << ", m2: " << m2 << ", xi_p: " << xi_p << ", xi_k: " << xi_k
             << '\n';
#endif
        // fill the ntuple
        vars.Fill(e_miss, m_recoil, m2, xi_p, xi_k);
    }
    // event loop ends.
    // --------------------------------------------------------------------------

    infile.Close();
    cout << APPNAME << ": processed " << nentries << " events.\n";

    TFile outfile{argv[2], "recreate"};
    cout << APPNAME << ": the output is stored in " << outfile.GetName()
         << '\n';
    vars.Write();
    outfile.Close();
}
