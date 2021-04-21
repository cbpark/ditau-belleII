/*
 *  Copyright 2021 Chan Beom Park <cbpark@gmail.com>
 */

#include <Math/LorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector2.h>
#include <cmath>  // std::hypot
#include <iostream>
#include <string>  // std::string

using std::cerr;
using std::cout;

/// we use the recommended LorentzVector rather than legacy TLorentzVector.
/// see https://root.cern/doc/master/classROOT_1_1Math_1_1LorentzVector.html
using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;

const std::string appname = "bditau";

int main(int, char *argv[]) {
    auto infile = TFile(argv[1]);
    cout << appname << ": the input file is " << infile.GetName() << '\n';

    auto keys = infile.GetListOfKeys();
    if (keys->GetSize() < 1) {
        std::cerr << appname << ": the input has no tree.\n";
        infile.Close();
        return 1;
    }

    auto tree_name = infile.GetListOfKeys()->At(0)->GetName();
    cout << appname << ": the name of the tree is " << tree_name << '\n';
    auto event = dynamic_cast<TTree *>(infile.Get(tree_name));

    // the four-momentum of visible particles in the three-prong decay.
    double p1x, p1y, p1z, e1;
    // the three-momentum of visible particles in the one-prong decay.
    double p2x, p2y, p2z;
    event->SetBranchAddress("tau_3prong_px_CMS", &p1x);
    event->SetBranchAddress("tau_3prong_py_CMS", &p1y);
    event->SetBranchAddress("tau_3prong_pz_CMS", &p1z);
    event->SetBranchAddress("tau_3prong_E_CMS", &e1);
    event->SetBranchAddress("track_1prong_px_CMS", &p2x);
    event->SetBranchAddress("track_1prong_py_CMS", &p2y);
    event->SetBranchAddress("track_1prong_pz_CMS", &p2z);

    // --------------------------------------------------------------------------
    // event loop
    const auto nentries = event->GetEntries();
    for (auto iev = 0; iev != nentries; ++iev) {
        event->GetEntry(iev);

        const LorentzVector p1(p1x, p1y, p1z, e1);
        // the visible particle in the one-prong decay is assumed to be
        // massless. std::hypot(x, y, z) = sqrt(x^2 + y^2 + z^2).
        const LorentzVector p2(p2x, p2y, p2z, std::hypot(p2x, p2y, p2z));

        // the missing transverse momentum.
        const TVector2 ptmiss(-(p1.Px() + p2.Px()), -(p1.Py() + p2.Py()));
#ifdef DEBUG
        cout << "p1: " << p1 << "\np2: " << p2 << '\n';
        cout << "ptmiss: ";
        ptmiss.Print();
#endif
    }
    // event loop ends.
    // --------------------------------------------------------------------------

    infile.Close();
}
