// Copyright 2016 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

/*******************************************************************************
 * This file is part of NUISANCE.
 *
 * NUISANCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NUISANCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with NUISANCE. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#include "MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu.h"

#include "TH1D.h"
#include "TH2D.h"
#include "WireCellHelper.h"

#include <cmath>

//********************************************************************
template <distribution_t D, distribution_t... Ds>
MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<D, Ds...>::MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu(
    nuiskey samplekey) {
  //********************************************************************

  // Sample overview ---------------------------------------------------
  std::string descrip =
      "MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu sample. \n"
      "Target: Ar \n"
      "Flux: BNB FHC numu \n"
      "Signal: CC 1Mu1pi0Np WireCell \n"
      "Contact: microboone_info@fnal.gov\n"
      "Reference: https://arxiv.org/abs/2606.23559 \n";

  // Setup common settings
  fSettings = LoadSampleSettings(samplekey);
  fSettings.SetDescription(descrip);
  fSettings.SetXTitle("Bin Number");
  fSettings.SetYTitle("d#sigma (cm^{2}/^{36}/Ar)");
  fSettings.SetAllowedTypes("FULL,DIAG/FREE,SHAPE,FIX/SYSTCOV/STATCOV", "FIX/FULL");
  fSettings.SetEnuRange(0.0, 10.0);
  fSettings.DefineAllowedTargets("Ar");

  // Plot information
  fSettings.SetTitle("MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu");
  fSettings.DefineAllowedSpecies("numu");

  std::string sample_name = fSettings.GetName();

  FinaliseSampleSettings();

  // Scaling Setup ---------------------------------------------------
  // ScaleFactor setup for DiffXSec/cm2/Ar nucleus.
  fScaleFactor = ((GetEventHistogram()->Integral("width") / (fNEvents + 0.0)) * 1E-38 /
                  TotalIntegratedFlux()) *
                 40.0;

  // Setup Histograms.
  CC1Mu1pi0NpHelper<D, Ds...> ana_helper;
  ana_helper.load_measurement();

  // Access bin numbers and widths for each measurement based on input physics values.
  fTable = ana_helper.get_lookuptable();

  // We do not use the template arguments here because this also handles the kAll case.
  fDists = fTable.get_dists();

  // Data histogram.
  fDataHist = static_cast<TH1D *>(ana_helper.get_data());
  fDataHist->SetName(Form("%s_data", sample_name.c_str()));

  ScaleData(1E-36);

  // Measurement covariance.
  fFullCovar = ana_helper.get_cov_m();

  // Convert from 10^-72 to 10^-76, which NUISANCE expects.
  (*fFullCovar) *= 1E4;

  // Set data errors from the covariance matrix.
  // This may not be strictly required, but it suppresses some warnings.
  for (int i = 0; i < fDataHist->GetNbinsX(); ++i) {
    fDataHist->SetBinError(i + 1, std::sqrt((*fFullCovar)(i, i)) * 1E-38);
  }

  // Additional Wiener-SVD Ac smearing matrix.
  fSmearingMatrix = ana_helper.get_ac_m();
  fSaveFine = false;

  // Final setup ---------------------------------------------------
  FinaliseMeasurement();
}

//********************************************************************
template <distribution_t D, distribution_t... Ds>
bool MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<D, Ds...>::isSignal(
    FitEvent *customEvent) {
  const double INVALID = -99999.0;

  const double MuonKEThreshold = 40.0;
  const double ProtonKEThreshold = 60.0;
  const double ChargedPiKEThreshold = 40.0;

  // Observables.
  double MuonKE = INVALID;
  double ProtonKE = INVALID;

  // Counters.
  int NPi0 = 0;
  int NProtonAboveThreshold = 0;
  int NChargedPionAboveThreshold = 0;

  // Muon ------------------------------------------------------------
  FitParticle *Muon = customEvent->GetHMFSParticle(13);
  const bool HasMuon = (Muon != nullptr);

  if (HasMuon) {
    MuonKE = Muon->E() - Muon->M();
  }

  // Final-state particles ------------------------------------------
  const std::vector<int> AllFSIndices = customEvent->GetAllFSParticleIndices(0);

  for (const int idx : AllFSIndices) {
    const int pdg = customEvent->GetParticlePDG(idx);

    // Muon is handled separately.
    if (std::abs(pdg) == 13) {
      continue;
    }

    FitParticle *part = customEvent->GetParticle(idx);
    if (!part) {
      continue;
    }

    const double E = part->E();
    const double M = part->M();

    // Protection against unusual generator records.
    if (E < M) {
      continue;
    }

    const double KE = part->KE();

    if (pdg == 111) {
      ++NPi0;
    } else if (pdg == 2212) {
      if (KE >= ProtonKEThreshold) {
        ++NProtonAboveThreshold;
      }

      // Leading proton by KE.
      if (ProtonKE < 0.0 || KE > ProtonKE) {
        ProtonKE = KE;
      }
    } else if (std::abs(pdg) == 211) {
      if (KE > ChargedPiKEThreshold) {
        ++NChargedPionAboveThreshold;
      }
    }
  }

  // Signal definition -----------------------------------------------
  const bool IsCC1pi0NpSignal =
      HasMuon &&
      MuonKE >= MuonKEThreshold &&
      NPi0 == 1 &&
      NProtonAboveThreshold >= 1 &&
      ProtonKE >= ProtonKEThreshold &&
      NChargedPionAboveThreshold == 0;

  return IsCC1pi0NpSignal;
}

//********************************************************************
template <distribution_t D, distribution_t... Ds>
void MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<D, Ds...>::FillEventVariables(
    FitEvent *customEvent) {
  const double INVALID = -99999.0;

  auto ResetBins = [&]() {
    fXVar = -999;
    for (const auto dist : fDists) {
      fXVars[dist] = -999;
    }
  };

  if (!isSignal(customEvent)) {
    ResetBins();
    return;
  }

  auto SafeCosTheta = [](double pz, double p) {
    if (p <= 0.0) {
      return -5.0;
    }

    return pz / p;
  };

  auto ClampUnit = [](double x) {
    if (x > 1.0) {
      return 1.0;
    }

    if (x < -1.0) {
      return -1.0;
    }

    return x;
  };

  // Initialize observables ------------------------------------------

  // Muon observables.
  double MuonKE = INVALID;
  double CosThetaMu = INVALID;

  // pi0 observables.
  double Pi0KE = INVALID;
  double Pi0CosTheta = INVALID;

  // Leading-proton observables.
  double ProtonKE = INVALID;
  double ProtonCosTheta = INVALID;

  // TKI observables.
  double DeltaPT = INVALID;
  double DeltaAlphaT = INVALID;
  double DeltaPTx = INVALID;
  double DeltaPTy = INVALID;

  // Momentum components.
  double mu_px = INVALID;
  double mu_py = INVALID;
  double mu_pz = INVALID;

  double pi0_px = INVALID;
  double pi0_py = INVALID;
  double pi0_pz = INVALID;

  double p_px = INVALID;
  double p_py = INVALID;
  double p_pz = INVALID;

  // Muon -------------------------------------------------------------
  FitParticle *Muon = customEvent->GetHMFSParticle(13);

  if (Muon) {
    const double mu_E = Muon->E();
    const double mu_M = Muon->M();
    const double mu_p = Muon->p();

    mu_px = Muon->P3()[0];
    mu_py = Muon->P3()[1];
    mu_pz = Muon->P3()[2];

    MuonKE = mu_E - mu_M;
    CosThetaMu = SafeCosTheta(mu_pz, mu_p);
  }

  // Final-state particles ------------------------------------------
  const std::vector<int> AllFSIndices = customEvent->GetAllFSParticleIndices(0);

  for (const int idx : AllFSIndices) {
    const int pdg = customEvent->GetParticlePDG(idx);

    // Skip the muon because it was handled above.
    if (std::abs(pdg) == 13) {
      continue;
    }

    FitParticle *part = customEvent->GetParticle(idx);
    if (!part) {
      continue;
    }

    const double E = part->E();
    const double M = part->M();

    // Protection against unusual generator records.
    if (E < M) {
      continue;
    }

    const double KE = part->KE();
    const double px = part->P3()[0];
    const double py = part->P3()[1];
    const double pz = part->P3()[2];
    const double p = part->p();

    // pi0 observables.
    if (pdg == 111) {
      Pi0KE = KE;
      Pi0CosTheta = SafeCosTheta(pz, p);

      pi0_px = px;
      pi0_py = py;
      pi0_pz = pz;
    }

    // Leading-proton observables.
    else if (pdg == 2212) {
      if (ProtonKE < 0.0 || KE > ProtonKE) {
        ProtonKE = KE;
        ProtonCosTheta = SafeCosTheta(pz, p);

        p_px = px;
        p_py = py;
        p_pz = pz;
      }
    }
  }

  // Calculate TKI variables using the muon, pi0, and leading proton.
  const double delta_pT_x_comp = mu_px + pi0_px + p_px;
  const double delta_pT_y_comp = mu_py + pi0_py + p_py;

  DeltaPT = std::sqrt(delta_pT_x_comp * delta_pT_x_comp +
                      delta_pT_y_comp * delta_pT_y_comp);

  const double muT_mag = std::sqrt(mu_px * mu_px + mu_py * mu_py);

  if (muT_mag > 0.0 && DeltaPT > 0.0) {
    const double numerator = -mu_px * delta_pT_x_comp - mu_py * delta_pT_y_comp;
    const double denominator = muT_mag * DeltaPT;
    const double cos_alphaT = ClampUnit(numerator / denominator);
    const double rad_to_deg = 180.0 / std::acos(-1.0);

    DeltaAlphaT = rad_to_deg * std::acos(cos_alphaT);

    const double muT_hat_x = mu_px / muT_mag;
    const double muT_hat_y = mu_py / muT_mag;

    const double z_cross_muT_hat_x = -muT_hat_y;
    const double z_cross_muT_hat_y = muT_hat_x;

    DeltaPTx = z_cross_muT_hat_x * delta_pT_x_comp +
               z_cross_muT_hat_y * delta_pT_y_comp;

    DeltaPTy = -muT_hat_x * delta_pT_x_comp - muT_hat_y * delta_pT_y_comp;
  }

  int curr_bin = 0;

  // Loop over blocks and store the global bin for each distribution.
  for (const auto dist : fDists) {
    const int nblockbins = fTable.get_nbins(dist);
    int localbin = -1;

    switch (dist) {
      case kCC1pi0NpKEMu:
        localbin = fTable.find_bin(dist, MuonKE);
        break;

      case kCC1pi0NpCosThetaMu:
        localbin = fTable.find_bin(dist, CosThetaMu);
        break;

      case kCC1pi0NpKEpi0:
        localbin = fTable.find_bin(dist, Pi0KE);
        break;

      case kCC1pi0NpCosThetapi0:
        localbin = fTable.find_bin(dist, Pi0CosTheta);
        break;

      case kCC1pi0NpKEp:
        localbin = fTable.find_bin(dist, ProtonKE);
        break;

      case kCC1pi0NpCosThetap:
        localbin = fTable.find_bin(dist, ProtonCosTheta);
        break;

      case kCC1pi0NpdpT:
        localbin = fTable.find_bin(dist, DeltaPT);
        break;

      case kCC1pi0NpdpTx:
        localbin = fTable.find_bin(dist, DeltaPTx);
        break;

      case kCC1pi0NpdpTy:
        localbin = fTable.find_bin(dist, DeltaPTy);
        break;

      case kCC1pi0NpdaT:
        localbin = fTable.find_bin(dist, DeltaAlphaT);
        break;
    }

    // FillHistograms is handled directly because one event can belong to multiple blocks.
    fXVars[dist] = (localbin < 0) ? -999 : localbin + curr_bin;
    curr_bin += nblockbins;
  }
}

//********************************************************************
template <distribution_t D, distribution_t... Ds>
void MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<D, Ds...>::FillHistograms() {
  // Loop over each distribution block and fill it.
  for (const auto dist : fDists) {
    fXVar = fXVars[dist];
    Measurement1D::FillHistograms();
  }
}

//********************************************************************
template <distribution_t D, distribution_t... Ds>
void MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<D, Ds...>::ConvertEventRates() {
  // Standard conversion.
  Measurement1D::ConvertEventRates();

  int curr_bin = 0;

  // Convert each distribution block to a differential cross section.
  for (const auto dist : fDists) {
    const int nblockbins = fTable.get_nbins(dist);

    for (int i = 0; i < nblockbins; ++i) {
      const double bin_width = fTable.get_width(dist, i);
      const int global_bin = curr_bin + i + 1;

      fMCHist->SetBinContent(global_bin, fMCHist->GetBinContent(global_bin) / bin_width);
      fMCHist->SetBinError(global_bin, fMCHist->GetBinError(global_bin) / bin_width);
    }

    curr_bin += nblockbins;
  }

  // Apply Wiener-SVD Ac smearing.
  const int n = fMCHist->GetNbinsX();
  TVectorD v(n);
  TVectorD e(n);

  for (int i = 0; i < n; ++i) {
    v(i) = fMCHist->GetBinContent(i + 1);
    e(i) = std::pow(fMCHist->GetBinError(i + 1), 2);
  }

  TVectorD vs = (*fSmearingMatrix) * v;
  TVectorD es = (*fSmearingMatrix) * e;

  for (int i = 0; i < n; ++i) {
    fMCHist->SetBinContent(i + 1, vs(i));
    fMCHist->SetBinError(i + 1, std::sqrt(es(i)));
  }
}

template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kCC1pi0NpKEMu>;
template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kCC1pi0NpCosThetaMu>;
template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kCC1pi0NpKEpi0>;
template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kCC1pi0NpCosThetapi0>;
template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kCC1pi0NpKEp>;
template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kCC1pi0NpCosThetap>;
template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kCC1pi0NpdpT>;
template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kCC1pi0NpdpTx>;
template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kCC1pi0NpdpTy>;
template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kCC1pi0NpdaT>;
template class MicroBooNE_BNB_NumuCC1pi0Np_2026_XSec_nu<kAllCC1pi0Np>;

