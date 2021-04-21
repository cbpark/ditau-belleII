/*
 *  Copyright 2021 Chan Beom Park <cbpark@gmail.com>
 */

#include <algorithm>  // std::max, std::minmax
#include <cmath>      // std::hypot
#include <iostream>
#include <optional>  // std::optional

// ROOT
#include "Math/LorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"

// YAM2
#include "YAM2/yam2.h"

using std::cerr;
using std::cout;

/// we use the recommended LorentzVector rather than legacy TLorentzVector.
using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;

const auto appname{"bditau"};

/// the longitudinal momentum of di-tau system (it's 0 for the CMS).
const double PZTOT = 0.0;

/// sqrt(s) of CMS.
const double SQRTS = 10.579;

/// the four-momentum of CMS.
const LorentzVector CMS{0.0, 0.0, PZTOT, SQRTS};

/// the invisible particle mass.
const yam2::Mass MINV{0.0};

/// zero momentum for convenience to calculate M2.
const yam2::FourMomentum ZERO;

/// E_{miss} = sqrt(ptmiss_x^2 + ptmiss_y^2 + pzmiss^2).
double eMiss(const LorentzVector &p1, const LorentzVector &p2,
             const TVector2 &ptmiss);

/// the recoil mass.
double mRecoil(const LorentzVector &p1, const LorentzVector &p2);

/// ratio of min(p1, p2) and max(p1, p2).
///
/// Here 'p1' and 'p2' could be either visible or invisible particle momenta.
/// By definition, it is in between 0 and 1.
double ratioMomentum(const LorentzVector &p1, const LorentzVector &p2);

std::optional<yam2::M2Solution> getM2(const LorentzVector &p1,
                                      const LorentzVector &p2,
                                      const TVector2 &ptmiss);

/// converts 'yam2::FourMomentum' to 'LorentzVector'.
LorentzVector toLorentzVector(const yam2::FourMomentum &p);

// ------------------------------------------------------------------------------
// the main function
int main(int, char *argv[]) {
    // input root file.
    auto infile = TFile(argv[1]);
    cout << appname << ": the input file is " << infile.GetName() << '\n';

    // check the tree.
    auto keys = infile.GetListOfKeys();
    if (keys->GetSize() < 1) {
        cerr << appname << ": the input has no tree.\n";
        infile.Close();
        return 1;
    }

    // get the tree.
    auto tree_name = infile.GetListOfKeys()->At(0)->GetName();
    cout << appname << ": the name of the tree is " << tree_name << '\n';

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

    // the observables.
    double e_miss, m_recoil, m2, xi_p, xi_k;

    // --------------------------------------------------------------------------
    // event loop
    const auto nentries = event->GetEntries();
    for (auto iev = 0; iev != nentries; ++iev) {
        event->GetEntry(iev);

        p1.SetPxPyPzE(p1x, p1y, p1z, e1);
        // the visible particle in the one-prong decay is assumed to be
        // massless. std::hypot(x, y, z) = sqrt(x^2 + y^2 + z^2) since c++17.
        p2.SetPxPyPzE(p2x, p2y, p2z, std::hypot(p2x, p2y, p2z));
        ptmiss.Set(-(p1.Px() + p2.Px()), -(p1.Py() + p2.Py()));

#ifdef DEBUG
        cout << "p1: " << p1 << ", p2: " << p2 << '\n';
        cout << "ptmiss: ";
        ptmiss.Print();
#endif

        e_miss = eMiss(p1, p2, ptmiss);

        m_recoil = mRecoil(p1, p2);

        xi_p = ratioMomentum(p1, p2);

        const auto m2sol = getM2(p1, p2, ptmiss);
        if (!m2sol) {
            cerr << appname << ": failed to find minimum for M2! (" << iev
                 << ")\n";
            m2 = -1.0;
            xi_k = -1.0;
        } else {
            m2 = m2sol.value().m2();
            k1sol = toLorentzVector(m2sol.value().k1());
            k2sol = toLorentzVector(m2sol.value().k2());
            xi_k = ratioMomentum(k1sol, k2sol);
        }
        cout << "e_miss: " << e_miss << ", m_recoil: " << m_recoil
             << ", xi_p: " << xi_p << '\n';
        cout << "m2: " << m2 << ", xi_k: " << xi_k << '\n';
    }
    cout << appname << ": processed " << nentries << " events.\n";
    // event loop ends.
    // --------------------------------------------------------------------------

    infile.Close();
}
// ------------------------------------------------------------------------------

double eMiss(const LorentzVector &p1, const LorentzVector &p2,
             const TVector2 &ptmiss) {
    const double pzmiss = PZTOT - p1.pz() - p2.pz();
    return std::hypot(ptmiss.Px(), ptmiss.Py(), pzmiss);
}

double mRecoil(const LorentzVector &p1, const LorentzVector &p2) {
    const auto p_miss = CMS - p1 - p2;
    return p_miss.mass();
}

double ratioMomentum(const LorentzVector &p1, const LorentzVector &p2) {
    const auto [pmin, pmax] = std::minmax({p1.P(), p2.P()});
    return pmin / std::max(pmax, 1.0e-10);
}

std::optional<yam2::M2Solution> getM2(const LorentzVector &p1,
                                      const LorentzVector &p2,
                                      const TVector2 &ptmiss) {
    // input kinematic configuration for M2.
    // see the above for the constants (ZERO, MINV, SQRTS, PZTOT).
    const auto input = yam2::mkInput({{p1.e(), p1.px(), p1.py(), p1.pz()},
                                      {p2.e(), p2.px(), p2.py(), p2.pz()}},
                                     {ZERO, ZERO}, {ptmiss.Px(), ptmiss.Py()},
                                     MINV, {}, SQRTS, {PZTOT});
    return yam2::m2Cons(input, 1.0e-6, 1000);
}

LorentzVector toLorentzVector(const yam2::FourMomentum &p) {
    return {p.px(), p.py(), p.pz(), p.e()};
}
