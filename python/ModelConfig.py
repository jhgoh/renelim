import sys
import ROOT
import numpy as np

from Config import ConfigRENE, getFileAndObj


def load_model(config_path="config.yaml", det_idx=0):
    """Build the RooFit model described by the configuration.

    Parameters
    ----------
    config_path : str
        Path to the YAML configuration file.
    det_idx : int
        Index of the detector to use.

    Returns
    -------
    dict
        Dictionary containing workspace, configuration and main objects.
    """
    ###############################################################################
    ## Basic ROOT setup and load custom PDF classes
    ###############################################################################
    ROOT.gROOT.ProcessLineSync(".x src/NuOscIBDPdf.cxx+")
    ROOT.gROOT.ProcessLineSync(".x src/SmearedNuOscIBDPdf.cxx+")

    ws = ROOT.RooWorkspace("ws", "ws")
    config = ConfigRENE(config_path)

    ###############################################################################
    ## Detector informations
    ## We take only one detector for this version.
    ## We also guess range of neutrino energy distribution from the response matrix
    ###############################################################################
    det_names = config.getDetectors("name")
    responses = config.getDetectors("response")
    baselines = [config.getBaselines(name) for name in det_names]

    det_name = det_names[det_idx]
    det_response = responses[det_idx]
    baseline_list = baselines[det_idx]

    f_resp = None
    h_resp = None
    if isinstance(det_response, str):
        f_resp, h_resp = getFileAndObj(det_response)
        ROOT.gROOT.cd()
        nbins_enu = h_resp.GetNbinsX()
        min_enu = h_resp.GetXaxis().GetXmin()
        max_enu = h_resp.GetXaxis().GetXmax()
        nbins_ereco = h_resp.GetNbinsY()
        min_ereco = h_resp.GetYaxis().GetXmin()
        max_ereco = h_resp.GetYaxis().GetXmax()
    elif isinstance(det_response, list):
        respA, respB, respC = det_response
        h_resp = None
        nbins_enu, min_enu, max_enu = 100, 0, 10
        nbins_ereco, min_ereco, max_ereco = 100, 0, 10
    else:
        raise ValueError("detectors[].response has to be string or list")

    reactor_idx = np.argmin(baseline_list)
    baseline = baseline_list[reactor_idx]

    ###############################################################################
    ## Important parameters of interests
    ## Define these variables in advance, to avoid possible buggy behaviours
    ###############################################################################
    v_sin13 = ROOT.RooRealVar("v_sin13", "sin^{2}(#theta_{13})", 0, 0.5)
    v_sin14 = ROOT.RooRealVar("v_sin14", "sin^{2}(#theta_{14})", 0, 0.5)
    v_dm31 = ROOT.RooRealVar("v_dm31", "#Delta^{2}_{31}", 0, 1, unit="eV^{2}")
    v_dm41 = ROOT.RooRealVar("v_dm41", "#Delta^{2}_{41}", 0, 5, unit="eV^{2}")
    ws.Import(v_sin14)
    ws.Import(v_dm41)
    ws.Import(v_sin13)
    ws.Import(v_dm31)
    v_sin13 = ws.var("v_sin13")
    v_sin14 = ws.var("v_sin14")
    v_dm31 = ws.var("v_dm31")
    v_dm41 = ws.var("v_dm41")

    ###############################################################################
    ## Neutrino energy spectrums
    ###############################################################################
    ## (Maybe later) Consider the burn-up effect,
    ## the fuel composition is a subject to be changed in time.
    ## We choose the 6th core only (1st one in the configuration)
    elem_names = config.getReactors("elements")[reactor_idx]["name"]
    elem_fracs = config.getReactors("elements")[reactor_idx]["fraction"]
    formula = "1"
    formula_vars = []
    for i in range(len(elem_names) - 1):
        en = elem_names[i]
        ws.factory(f"v_{en}[{elem_fracs[i]}, 0, 1]")
        ws.var(f"v_{en}").setConstant(True)
        formula += f"-@{i}"
        formula_vars.append(f"v_{en}")
    formula_vars = ",".join(formula_vars)
    ws.factory(f'EXPR::v_{elem_names[-1]}("{formula}", {{ {formula_vars} }})')
    v_elem_fracs = ROOT.RooArgList()
    for i in range(len(elem_names) - 1):
        v_elem_fracs.add(ws.var(f"v_{elem_names[i]}"))
    v_elem_fracs.add(ws.function(f"v_{elem_names[-1]}"))

    ## Load the Neutrino flux model, such as Huber-Mueller, incorporating the fuel compositions
    grps_HM = ROOT.std.vector("TGraph")()
    for en in elem_names:
        _, grp = getFileAndObj(config.get(f"physics.isotope_flux.{en}"))
        ROOT.gROOT.cd()
        grps_HM.push_back(grp.Clone())
        del grp

    ###############################################################################
    ## IBD cross section
    ###############################################################################
    _, grp_xsec = getFileAndObj(config.get("physics.ibd_xsec"))
    ROOT.gROOT.cd()

    ################################################################################
    ## Build the Oscillated neutrino energy spectrum
    ################################################################################
    v_ENu = ROOT.RooRealVar("v_ENu", "Neutrino Energy", min_enu, max_enu, unit="MeV")
    # v_ENu.setBins(nbins_enu)
    v_ENu.setBins(1024)
    ws.Import(v_ENu)
    v_ENu = ws.var("v_ENu")

    _sin13, _sin13err = config.get("physics.oscillation.sin13")
    _dm31, _dm31err = config.get("physics.oscillation.dm31")
    v_sin13.setVal(_sin13)
    v_dm31.setVal(_dm31)

    v_L = ROOT.RooRealVar("v_L", "L", 0, 10000, unit="meter")
    v_L.setVal(baseline)
    ws.Import(v_L)
    v_L = ws.var("v_L")
    v_L.setConstant(True)

    pdf_ENu = ROOT.NuOscIBDPdf(
        "pdf_ENu", "pdf_ENu", v_ENu, v_L, v_sin13, v_dm31, v_sin14, v_dm41, v_elem_fracs, grps_HM, grp_xsec
    )
    ws.Import(pdf_ENu)
    pdf_ENu = ws.pdf("pdf_ENu")

    ################################################################################
    ## Build the convoluted PDF
    ## to do the convolution with conditional variable, the response matrix
    ## has to be transposed.
    ################################################################################
    v_EReco = ROOT.RooRealVar("v_EReco", "Reconstructed Energy", min_ereco, max_ereco, unit="MeV")
    # v_EReco.setBins(nbins_ereco)
    v_EReco.setBins(1024)
    ws.Import(v_EReco)
    v_EReco = ws.var("v_EReco")

    # h_respT = ROOT.TH2D("hRespT", h_resp.GetTitle(),
    #                    h_resp.GetNbinsY(), h_resp.GetYaxis().GetXmin(), h_resp.GetYaxis().GetXmax(),
    #                    h_resp.GetNbinsX(), h_resp.GetXaxis().GetXmin(), h_resp.GetXaxis().GetXmax())
    # for i in range(h_resp.GetNbinsX()):
    #    for j in range(h_resp.GetNbinsY()):
    #        h_respT.SetBinContent(j+1, i+1, h_resp.GetBinContent(i+1, j+1))
    # dh_respT = ROOT.RooDataHist("dhRespT", h_respT.GetTitle(),
    #                            ROOT.RooArgList(v_EReco, v_ENu), h_respT)
    # pdf_respT = ROOT.RooHistPdf("pdf_RespT", "pdf_RespT",
    #                           ROOT.RooArgList(v_EReco, v_ENu), dh_respT)
    # pdf_joint = ROOT.RooProdPdf("pdf_Joint", "Joint pdf",
    #                            ROOT.RooArgList(pdf_ENu, pdf_respT),
    #                            ROOT.RooFit.Conditional(ROOT.RooArgSet(pdf_respT), ROOT.RooArgSet(v_EReco)))
    # pdf_EReco = pdf_joint.createProjection(ROOT.RooArgSet(v_ENu))
    # pdf_EReco.SetName("pdf_EReco")
    # pdf_EReco.SetTitle("PDF of reconstructed energy")
    # fmt: off
    pdf_EReco = ROOT.SmearedNuOscIBDPdf(
        "pdf_EReco", "pdf_EReco", v_EReco, v_ENu, v_L,
        v_sin13, v_dm31, v_sin14, v_dm41,
        v_elem_fracs, grps_HM, grp_xsec, h_resp
    )
    # fmt: on
    ws.Import(pdf_EReco)
    pdf_EReco = ws.pdf("pdf_EReco")

    return {
        "ws": ws,
        "config": config,
        "v_sin14": v_sin14,
        "v_dm41": v_dm41,
        "v_ENu": v_ENu,
        "v_EReco": v_EReco,
        "pdf_ENu": pdf_ENu,
        "pdf_EReco": pdf_EReco,
    }
