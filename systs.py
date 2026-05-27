import awkward as ak
import pandas as pd
import polars as pl
import uproot
from tqdm import tqdm
from enum import IntEnum
import matplotlib.pyplot as plt
import numpy as np
import yaml
from typing import Optional

class Flavor(IntEnum):
    NuE = 12
    NuEBar = -12
    NuMu = 14
    NuMuBar = -14
    NuTau = 16
    NuTauBar = -16
    NC = 0

def flavour_reco() -> pl.Expr:
    return pl.when(pl.col('cvn_nue') > pl.col('cvn_numu')
        ).then(
            pl.when(
                pl.col('cvn_nue') > pl.col('cvn_nc')
            ).then(
                Flavor.NuE
            ).otherwise(
                Flavor.NC
            )
        ).otherwise(
            pl.when(
                pl.col('cvn_numu') > pl.col('cvn_nc')
            ).then(
                Flavor.NuMu
            ).otherwise(
                Flavor.NC
            )
        )

def direc_reco() -> pl.Expr:
    return pl.when(
        (pl.col('npfps') < 3) | (pl.col('recoE_nue') > 1.3) #For high energy events or low number of PFPs, we use the hit direc reco
    ).then(
        pl.col('direc_nc')
    ).otherwise(
        pl.when(
            pl.col('reco_pdg') == Flavor.NuMu
        ).then(
            pl.when(pl.col('direc_numu') != 0).then( #Avoiding cases with zero and moving to the hit direc reco
                pl.col('direc_numu')
            ).otherwise(
                pl.col('direc_nc')
            )
        ).otherwise(
            pl.when(
                pl.col('reco_pdg') == Flavor.NuE
            ).then(
                pl.when(pl.col('direc_nue') != 0).then( #Avoiding cases with zero and moving to the hit direc reco
                    pl.col('direc_nue')
                ).otherwise(
                    pl.col('direc_nc')
                )
            ).otherwise(
                pl.col('direc_nc')
            )
        )
    )

def E_reco() -> pl.Expr:
    return pl.when(
        pl.col('reco_pdg') == Flavor.NuMu
    ).then(
        pl.col('recoE_numu')
    ).otherwise(
        pl.col('recoE_nue')
    )

def Ehad_reco() -> pl.Expr:
    return pl.when(
        pl.col('reco_pdg') == Flavor.NuMu
    ).then(
        pl.col('recoEhad_numu')
    ).otherwise(
        pl.col('recoEhad_nue')
    )

def Elep_reco() -> pl.Expr:
    return pl.when(
        pl.col('reco_pdg') == Flavor.NuMu
    ).then(
        pl.col('recoE_numu') - pl.col('recoEhad_numu')
    ).otherwise(
        pl.col('recoE_nue') - pl.col('recoEhad_nue')
    )

def data_filter() -> pl.Expr:
    return (pl.col('npfps') >= 2) & (pl.col('cvn_nue').is_not_null()) & (pl.col('direc_nc') != 0)


def _to_numpy_1d(values) -> np.ndarray:
    if values is None:
        return np.array([])
    return np.asarray(values)


def _angle_between_vectors(px1: float, py1: float, pz1: float, px2: float, py2: float, pz2: float) -> float:
    mag1 = np.sqrt(px1**2 + py1**2 + pz1**2)
    mag2 = np.sqrt(px2**2 + py2**2 + pz2**2)
    if mag1 <= 0 or mag2 <= 0:
        return np.nan
    cos_theta = (px1 * px2 + py1 * py2 + pz1 * pz2) / (mag1 * mag2)
    return float(np.arccos(np.clip(cos_theta, -1.0, 1.0)))


def match_prefsi_to_postfsi(
    pre_pdg,
    pre_px,
    pre_py,
    pre_pz,
    pre_e,
    post_pdg,
    post_px,
    post_py,
    post_pz,
    post_e,
    rel_energy_weight: float = 0.2,
) -> tuple[list[int], list[float], list[float], int]:
    pre_pdg = _to_numpy_1d(pre_pdg)
    pre_px = _to_numpy_1d(pre_px)
    pre_py = _to_numpy_1d(pre_py)
    pre_pz = _to_numpy_1d(pre_pz)
    pre_e = _to_numpy_1d(pre_e)

    post_pdg = _to_numpy_1d(post_pdg)
    post_px = _to_numpy_1d(post_px)
    post_py = _to_numpy_1d(post_py)
    post_pz = _to_numpy_1d(post_pz)
    post_e = _to_numpy_1d(post_e)

    n_pre = len(pre_pdg)
    match_post_idx = [-1] * n_pre
    delta_theta = [float('nan')] * n_pre
    delta_e = [float('nan')] * n_pre

    if n_pre == 0 or len(post_pdg) == 0:
        return match_post_idx, delta_theta, delta_e, 0

    common_pdgs = np.intersect1d(np.unique(pre_pdg), np.unique(post_pdg))
    for pdg in common_pdgs:
        pre_indices = np.where(pre_pdg == pdg)[0]
        post_indices = np.where(post_pdg == pdg)[0]
        if len(pre_indices) == 0 or len(post_indices) == 0:
            continue

        candidates: list[tuple[float, int, int, float, float]] = []
        for i in pre_indices:
            for j in post_indices:
                dtheta = _angle_between_vectors(pre_px[i], pre_py[i], pre_pz[i], post_px[j], post_py[j], post_pz[j])
                de = float(post_e[j] - pre_e[i])
                rel_de = abs(de) / max(abs(float(pre_e[i])), 1e-6)
                score = (0.0 if np.isnan(dtheta) else dtheta) + rel_energy_weight * rel_de
                candidates.append((score, int(i), int(j), dtheta, de))

        candidates.sort(key=lambda x: x[0])
        used_pre = set()
        used_post = set()
        for _, i, j, dtheta, de in candidates:
            if i in used_pre or j in used_post:
                continue
            match_post_idx[i] = j
            delta_theta[i] = dtheta
            delta_e[i] = de
            used_pre.add(i)
            used_post.add(j)

    n_matched = int(np.sum(np.array(match_post_idx) >= 0))
    return match_post_idx, delta_theta, delta_e, n_matched


class CAFManager:
    def __init__(self, fname:str, genie:str|None = None) -> None:
        self.data = self.read_cafs(fname)
        self.augment_data()
        self.data_filter = data_filter()
        if genie is not None:
            self.add_genie(genie)

    def read_cafs(self, fname:str) -> pl.DataFrame:
        def fix_empty_arrays(data: ak.Array) -> ak.Array:
            return ak.where(ak.num(data) == 0, ak.Array([[-999]] * len(data)), data)

        with uproot.open(fname) as f:
            weights = f['weights'].arrays(library='pd')
            weights['nuPDG'] = ak.flatten(f['cafTree/rec/mc/mc.nu.pdg'].array())
            weights['Ev'] = ak.flatten(f['cafTree/rec/mc/mc.nu.E'].array())
            # weights['isCC'] = ak.flatten(f['cafTree/rec/mc/mc.nu.iscc'].array())
            weights['NuMomY'] = ak.flatten(f['cafTree/rec/mc/mc.nu.momentum.y'].array())
            weights['mode'] = ak.flatten(f['cafTree/rec/mc/mc.nu.mode'].array())

            recoE_numu = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.Enu.lep_calo'].array())
            recoE_nue = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.Enu.e_calo'].array())
            recoEhad_numu = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.Enu.mu_had'].array())
            recoEhad_nue = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.Enu.e_had'].array())
            weights['recoE_numu'] = ak.flatten(recoE_numu)
            weights['recoE_nue'] = ak.flatten(recoE_nue)
            weights['recoEhad_numu'] = ak.flatten(recoEhad_numu)
            weights['recoEhad_nue'] = ak.flatten(recoEhad_nue)


            direc_numu = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.dir.lngtrk.y'].array())
            direc_nue = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.dir.heshw.y'].array())
            direc_nc = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.dir.heshw.y'].array())
            weights['direc_numu'] = -ak.flatten(direc_numu) #Minus sign to have the convention negative=upgoing neutrinos
            weights['direc_nue'] = -ak.flatten(direc_nue) #Minus sign to have the convention negative=upgoing neutrinos
            weights['direc_nc'] = -ak.flatten(direc_nc) #Minus sign to have the convention negative=upgoing neutrinos

            cvn_nue = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.nuhyp.cvn.nue'].array())
            cvn_numu = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.nuhyp.cvn.numu'].array())
            cvn_nc = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.nuhyp.cvn.nc'].array())

            weights['cvn_numu'] = ak.flatten(cvn_numu)
            weights['cvn_nue'] = ak.flatten(cvn_nue)
            weights['cvn_nc'] = ak.flatten(cvn_nc)

            weights['npfps'] = ak.flatten(fix_empty_arrays(f['cafTree/rec/fd/fd.hd.pandora/fd.hd.pandora.npfps'].array()))
        weights['direc_true'] = -weights['NuMomY']/weights['Ev'] #Minus sign to have the convention negative=upgoing neutrinos
        weights['nue_w'] *= weights["xsec"]
        weights['numu_w'] *= weights["xsec"]

        return pl.from_pandas(pd.DataFrame(weights))
    
    def augment_data(self) -> None:
        self.data = self.data.with_columns(
            reco_pdg=flavour_reco()
        ).with_columns(
            reco_direc=direc_reco()
        ).with_columns(
            reco_E=E_reco()
        ).with_columns(
            reco_Ehad=Ehad_reco()
        ).with_columns(
            reco_Elep=Elep_reco()
        )

    def add_genie(self, fname:str) -> None:
        genie = self.load_genie(fname)
        self.add_columns(genie)

    def add_columns(self, data:pl.DataFrame) -> None:
        self.data = pl.concat([self.data, data], how='horizontal')

    def load_genie(self, fname:str) -> pl.DataFrame:
        with uproot.open(fname) as f:
            genie = f["genie_dump"]
            keys = genie.keys()
            loaded = {}
            for key in tqdm(keys):
                loaded[key] = genie[key].array(library="np")
            genie = pl.DataFrame(loaded)
        return genie

    def add_qe_fsi_matching(self, rel_energy_weight: float = 0.2) -> None:
        required_columns = [
            'mode',
            'preFSI_Pdg', 'preFSI_Px', 'preFSI_Py', 'preFSI_Pz', 'preFSI_E',
            'postFSI_Pdg', 'postFSI_Px', 'postFSI_Py', 'postFSI_Pz', 'postFSI_E'
        ]
        missing = [col for col in required_columns if col not in self.data.columns]
        if missing:
            raise ValueError(
                f"Missing columns for QE FSI matching: {missing}. "
                "Ensure CAF and GENIE dump are both loaded in CAFManager."
            )

        qe_mask = (self.data['mode'] == 1).to_list()

        match_indices_all = []
        delta_theta_all = []
        delta_e_all = []
        n_pre_all = []
        n_post_all = []
        n_matched_all = []
        frac_matched_all = []

        iterator = zip(
            qe_mask,
            self.data['preFSI_Pdg'].to_list(),
            self.data['preFSI_Px'].to_list(),
            self.data['preFSI_Py'].to_list(),
            self.data['preFSI_Pz'].to_list(),
            self.data['preFSI_E'].to_list(),
            self.data['postFSI_Pdg'].to_list(),
            self.data['postFSI_Px'].to_list(),
            self.data['postFSI_Py'].to_list(),
            self.data['postFSI_Pz'].to_list(),
            self.data['postFSI_E'].to_list(),
        )

        for is_qe, pre_pdg, pre_px, pre_py, pre_pz, pre_e, post_pdg, post_px, post_py, post_pz, post_e in iterator:
            if not is_qe:
                match_indices_all.append([])
                delta_theta_all.append([])
                delta_e_all.append([])
                n_pre_all.append(0)
                n_post_all.append(0)
                n_matched_all.append(0)
                frac_matched_all.append(0.0)
                continue

            pre_count = int(len(pre_pdg)) if pre_pdg is not None else 0
            post_count = int(len(post_pdg)) if post_pdg is not None else 0
            match_idx, dtheta, de, n_matched = match_prefsi_to_postfsi(
                pre_pdg,
                pre_px,
                pre_py,
                pre_pz,
                pre_e,
                post_pdg,
                post_px,
                post_py,
                post_pz,
                post_e,
                rel_energy_weight=rel_energy_weight,
            )

            match_indices_all.append(match_idx)
            delta_theta_all.append(dtheta)
            delta_e_all.append(de)
            n_pre_all.append(pre_count)
            n_post_all.append(post_count)
            n_matched_all.append(n_matched)
            frac_matched_all.append(float(n_matched / pre_count) if pre_count > 0 else 0.0)

        self.data = self.data.with_columns(
            pl.Series('preFSI_match_post_index', match_indices_all),
            pl.Series('preFSI_delta_theta', delta_theta_all),
            pl.Series('preFSI_delta_E', delta_e_all),
            pl.Series('preFSI_n', n_pre_all),
            pl.Series('postFSI_n', n_post_all),
            pl.Series('prefsi_n_matched', n_matched_all),
            pl.Series('prefsi_frac_matched', frac_matched_all),
        )

    def qe_fsi_match_long_table(self) -> pl.DataFrame:
        required_columns = [
            'mode', 'preFSI_Pdg', 'preFSI_E', 'postFSI_E',
            'preFSI_match_post_index', 'preFSI_delta_theta', 'preFSI_delta_E'
        ]
        missing = [col for col in required_columns if col not in self.data.columns]
        if missing:
            raise ValueError(
                f"Missing columns for long-format FSI table: {missing}. "
                "Run add_qe_fsi_matching() first."
            )

        rows = []
        for event_idx, row in enumerate(self.data.iter_rows(named=True)):
            if row['mode'] != 1:
                continue

            pre_pdgs = _to_numpy_1d(row['preFSI_Pdg'])
            pre_e = _to_numpy_1d(row['preFSI_E'])
            post_e = _to_numpy_1d(row['postFSI_E'])
            match_idx = _to_numpy_1d(row['preFSI_match_post_index'])
            dtheta = _to_numpy_1d(row['preFSI_delta_theta'])
            de = _to_numpy_1d(row['preFSI_delta_E'])

            for i in range(len(pre_pdgs)):
                j = int(match_idx[i]) if i < len(match_idx) else -1
                if j < 0:
                    continue
                if j >= len(post_e):
                    continue

                rows.append({
                    'event_index': event_idx,
                    'pre_index': int(i),
                    'post_index': int(j),
                    'pdg': int(pre_pdgs[i]),
                    'preE': float(pre_e[i]),
                    'postE': float(post_e[j]),
                    'deltaE': float(de[i]),
                    'deltaTheta': float(dtheta[i]),
                })

        if not rows:
            return pl.DataFrame({
                'event_index': [],
                'pre_index': [],
                'post_index': [],
                'pdg': [],
                'preE': [],
                'postE': [],
                'deltaE': [],
                'deltaTheta': [],
            })
        return pl.DataFrame(rows)
    
class Systematic:
    def __init__(self, config:dict) -> None:
        self.config = config
        self.setup_from_config()
    
    def setup_from_config(self) -> None:
        self.isCorrection = self.config.get('isCorrection', False)
        self.cv = self.config.get('centralParamValue', 0.0)
        self.name = self.config.get('prettyName', 'UnnamedSystematic')
        self.id = self.config.get('systParamId', -1)
        self.paramVariations = self.config.get('paramVariations', [])
        if self.isCorrection:
            self.paramVariations = [1]
        self.nshifts = len(self.paramVariations)

    def __repr__(self) -> str:
        return f"Systematic(name={self.name}, id={self.id}, isCorrection={self.isCorrection}, cv={self.cv}, paramVariations={self.paramVariations})"
    
    def _repr_html_(self) -> str:
        correction_badge = "🔧 Correction" if self.isCorrection else "📊 Systematic"
        correction_color = "#e67e22" if self.isCorrection else "#3498db"
        
        # Format parameter variations
        if self.paramVariations:
            var_count = len(self.paramVariations)
            if var_count <= 7:
                var_display = ', '.join(f'<span style="background:#ecf0f1;padding:2px 6px;border-radius:3px;margin:2px;display:inline-block;">{v}</span>' for v in self.paramVariations)
            else:
                first_few = ', '.join(f'<span style="background:#ecf0f1;padding:2px 6px;border-radius:3px;margin:2px;display:inline-block;">{v}</span>' for v in self.paramVariations[:5])
                var_display = f'{first_few} <span style="color:#7f8c8d;">... and {var_count - 5} more</span>'
        else:
            var_display = '<span style="color:#95a5a6;font-style:italic;">None</span>'
        
        html = f"""
        <style>
            .systematic-card {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                border: 2px solid {correction_color};
                border-radius: 10px;
                overflow: hidden;
                margin: 15px 0;
                box-shadow: 0 4px 12px rgba(0,0,0,0.15);
                max-width: 800px;
            }}
            .systematic-header {{
                background: linear-gradient(135deg, {correction_color} 0%, {correction_color}dd 100%);
                color: white;
                padding: 16px 20px;
                display: flex;
                align-items: center;
                justify-content: space-between;
            }}
            .systematic-title {{
                font-size: 20px;
                font-weight: bold;
                margin: 0;
            }}
            .systematic-badge {{
                background-color: rgba(255,255,255,0.2);
                padding: 4px 12px;
                border-radius: 12px;
                font-size: 12px;
                font-weight: 600;
            }}
            .systematic-body {{
                background-color: #ffffff;
                padding: 20px;
            }}
            .systematic-property {{
                display: flex;
                padding: 12px 0;
                border-bottom: 1px solid #ecf0f1;
                align-items: center;
            }}
            .systematic-property:last-child {{
                border-bottom: none;
            }}
            .systematic-label {{
                font-weight: 600;
                color: #2c3e50;
                min-width: 180px;
                display: flex;
                align-items: center;
            }}
            .systematic-label::before {{
                content: "▸";
                margin-right: 8px;
                color: {correction_color};
                font-weight: bold;
            }}
            .systematic-value {{
                color: #34495e;
                flex: 1;
            }}
            .cv-value {{
                background: linear-gradient(135deg, #f39c12 0%, #e67e22 100%);
                color: white;
                padding: 4px 12px;
                border-radius: 6px;
                font-weight: bold;
                display: inline-block;
            }}
            .id-badge {{
                background-color: #34495e;
                color: white;
                padding: 4px 10px;
                border-radius: 6px;
                font-family: 'Courier New', monospace;
                font-weight: bold;
            }}
        </style>
        <div class="systematic-card">
            <div class="systematic-header">
                <div class="systematic-title">⚙️ {self.name}</div>
                <div class="systematic-badge">{correction_badge}</div>
            </div>
            <div class="systematic-body">
                <div class="systematic-property">
                    <div class="systematic-label">Systematic ID</div>
                    <div class="systematic-value"><span class="id-badge">#{self.id}</span></div>
                </div>
                <div class="systematic-property">
                    <div class="systematic-label">Central Value</div>
                    <div class="systematic-value"><span class="cv-value">{self.cv}</span></div>
                </div>
                <div class="systematic-property">
                    <div class="systematic-label">Parameter Variations</div>
                    <div class="systematic-value">{var_display}</div>
                </div>
            </div>
        </div>
        """
        return html

class SystConfig:
    def __init__(self, fname:str) -> None:
        self.load_fcl(fname)
        self.blacklist = [
            'instance_name',
            'tool_options',
            'tool_type',
            'parameter_headers'
        ]

    def load_fcl(self, fname:str) -> None:
        from fhicl_parser import FhiclParser
        parser = FhiclParser()
        with open(fname, 'r') as f:
            config = parser.parse_text(f.read())
        self.config = config

    def get_syst(self, syst_name:str) -> Systematic:
        providers = self.config['generated_systematic_provider_configuration'].keys()
        for provider in providers:
            systs = self.config['generated_systematic_provider_configuration'][provider]
            for syst in systs:
                if syst == syst_name:
                    return Systematic(systs[syst_name])
        raise ValueError(f'Systematic {syst_name} not found in configuration.')
    
    def get_all_systs(self) -> list[Systematic]:
        syst_list = []
        providers = self.config['generated_systematic_provider_configuration'].keys()
        for provider in providers:
            if provider == "syst_providers":
                continue
            systs = self.config['generated_systematic_provider_configuration'][provider]
            for syst in systs:
                if syst in self.blacklist:
                    continue
                syst_list.append(Systematic(systs[syst]))
        return syst_list
    
    def to_oscillation_yaml(self, output_file: Optional[str] = None, 
                           sample_names: Optional[list[str]] = None,
                           error: float = 1.0,
                           flat_prior: bool = False,
                           param_bounds: Optional[list[float]] = None,
                           generated_value: float = 0.0,
                           prefit_value: float = 0.0,
                           spline_mode: Optional[list[int]] = None,
                           spline_type: str = "PerEvent",
                           mcmc_step_scale: float = 0.4,
                           syst_type: str = "Spline",
                           param_group: str = "Xsec") -> str:
        """
        Convert the SystConfig to YAML format for the oscillation fitter.
        
        Parameters:
        -----------
        output_file : str, optional
            Path to save the YAML file. If None, returns the YAML string.
        sample_names : list[str], optional
            List of sample names. Default is ["ND_*", "FD_*", "ATM"]
        error : float, optional
            Error value for the systematic. Default is 1.0
        flat_prior : bool, optional
            Whether to use a flat prior. Default is False
        param_bounds : list[float], optional
            Parameter bounds [min, max]. Default is [-9999.0, 9999.0]
        generated_value : float, optional
            Generated parameter value. Default is 0.0
        prefit_value : float, optional
            Pre-fit parameter value. Default is 0.0
        spline_mode : list[int], optional
            Mode for spline. Default is [0]
        spline_type : str, optional
            Type of spline. Default is "PerEvent"
        mcmc_step_scale : float, optional
            MCMC step scale. Default is 0.4
        syst_type : str, optional
            Type of systematic. Default is "Spline"
        param_group : str, optional
            Parameter group. Default is "Xsec"
            
        Returns:
        --------
        str : YAML formatted string
        """
        # Set defaults
        if sample_names is None:
            sample_names = ["ND_*", "FD_*", "ATM"]
        if param_bounds is None:
            param_bounds = [-9999.0, 9999.0]
        if spline_mode is None:
            spline_mode = [0]
        
        # Get all systematics
        all_systs = self.get_all_systs()
        
        # Build the YAML structure
        yaml_list = []
        
        for syst in all_systs:
            syst_dict = {
                'Systematic': {
                    'SampleNames': sample_names,
                    'Error': error,
                    'FlatPrior': flat_prior,
                    'Names': {
                        'FancyName': syst.name,
                        'ParameterName': syst.name
                    },
                    'ParameterBounds': param_bounds,
                    'ParameterValues': {
                        'Generated': generated_value if not syst.isCorrection else syst.cv,
                        'PreFitValue': prefit_value if not syst.isCorrection else syst.cv
                    },
                    'SplineInformation': {
                        'Mode': spline_mode,
                        'SplineName': syst.name,
                        'Type': spline_type
                    },
                    'StepScale': {
                        'MCMC': mcmc_step_scale
                    },
                    'Type': syst_type,
                    'ParameterGroup': param_group
                }
            }
            yaml_list.append(syst_dict)
        
        # Convert to YAML string
        yaml_str = yaml.dump(yaml_list, default_flow_style=False, sort_keys=False)
        
        # Save to file if specified
        if output_file:
            with open(output_file, 'w') as f:
                f.write(yaml_str)
            print(f"YAML configuration saved to {output_file}")
        
        return yaml_str

    
    def __repr__(self) -> str:
        repr_str = "╔══════════════════════════════════════════════════════════╗\n"
        repr_str += "║           📊 Systematic Configuration                    ║\n"
        repr_str += "╚══════════════════════════════════════════════════════════╝\n\n"
        
        providers = self.config['generated_systematic_provider_configuration'].keys()
        provider_count = 0
        total_systs = 0
        
        for provider in providers:
            if provider == "syst_providers":
                continue
            
            systs = self.config['generated_systematic_provider_configuration'][provider]
            syst_list = [syst for syst in systs if syst not in self.blacklist]
            
            if not syst_list:
                continue
                
            provider_count += 1
            total_systs += len(syst_list)
            
            repr_str += f"┌─ 🔧 Provider: {provider}\n"
            repr_str += f"│  ({len(syst_list)} systematic{'s' if len(syst_list) > 1 else ''})\n"
            repr_str += "│\n"
            
            for i, syst in enumerate(syst_list, 1):
                prefix = "└──" if i == len(syst_list) else "├──"
                repr_str += f"│  {prefix} ✓ {syst}\n"
            
            repr_str += "│\n"
        
        repr_str += "─" * 60 + "\n"
        repr_str += f"📈 Summary: {provider_count} provider{'s' if provider_count != 1 else ''}, "
        repr_str += f"{total_systs} total systematic{'s' if total_systs != 1 else ''}\n"
        
        return repr_str
    
    def _repr_html_(self) -> str:
        html = """
        <style>
            .syst-config {
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                border: 2px solid #3498db;
                border-radius: 8px;
                padding: 0;
                margin: 10px 0;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                overflow: hidden;
            }
            .syst-header {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 15px 20px;
                font-size: 18px;
                font-weight: bold;
            }
            .syst-provider {
                border-bottom: 1px solid #ecf0f1;
                padding: 15px 20px;
                background-color: #f8f9fa;
            }
            .syst-provider:last-child {
                border-bottom: none;
            }
            .provider-name {
                font-weight: bold;
                color: #2c3e50;
                font-size: 16px;
                margin-bottom: 8px;
                display: flex;
                align-items: center;
            }
            .provider-count {
                margin-left: 10px;
                font-size: 12px;
                background-color: #3498db;
                color: white;
                padding: 2px 8px;
                border-radius: 12px;
                font-weight: normal;
            }
            .syst-list {
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(250px, 1fr));
                gap: 8px;
                margin-top: 10px;
            }
            .syst-item {
                background-color: white;
                padding: 8px 12px;
                border-radius: 4px;
                border-left: 3px solid #27ae60;
                font-size: 13px;
                color: #34495e;
                display: flex;
                align-items: center;
            }
            .syst-item::before {
                content: "✓";
                color: #27ae60;
                font-weight: bold;
                margin-right: 8px;
            }
            .syst-footer {
                background-color: #34495e;
                color: white;
                padding: 12px 20px;
                font-size: 14px;
                text-align: center;
            }
            .syst-footer strong {
                color: #3498db;
            }
        </style>
        """
        
        html += '<div class="syst-config">'
        html += '<div class="syst-header">📊 Systematic Configuration</div>'
        
        providers = self.config['generated_systematic_provider_configuration'].keys()
        
        
        provider_count = 0
        total_systs = 0
        
        for provider in providers:
            if provider == "syst_providers":
                continue
            
            systs = self.config['generated_systematic_provider_configuration'][provider]
            syst_list = [syst for syst in systs if syst not in self.blacklist]
            
            if not syst_list:
                continue
                
            provider_count += 1
            total_systs += len(syst_list)
            
            html += '<div class="syst-provider">'
            html += f'<div class="provider-name">🔧 {provider}'
            html += f'<span class="provider-count">{len(syst_list)} systematic{"s" if len(syst_list) > 1 else ""}</span>'
            html += '</div>'
            html += '<div class="syst-list">'
            
            for syst in syst_list:
                html += f'<div class="syst-item">{syst}</div>'
            
            html += '</div></div>'
        
        html += f'<div class="syst-footer">📈 Summary: <strong>{provider_count}</strong> provider{"s" if provider_count != 1 else ""}, '
        html += f'<strong>{total_systs}</strong> total systematic{"s" if total_systs != 1 else ""}</div>'
        html += '</div>'
        
        return html

class VariationManager:
    def __init__(self, fname:str):
        self.fname = fname
        self.variations = {}

    def get_variations(self, syst:Systematic) -> pl.DataFrame:
        if syst.name in self.variations:
            return self.variations[syst.name]
        with uproot.open(self.fname) as f:
            variations = f['SystWeights'][syst.name].array(library='np')
            col_names = [f'{syst.name}_{i}' for i in range(syst.nshifts)]
            variations = pl.from_numpy(variations, schema=col_names)
            variations = variations.select(pl.all().clip(lower_bound=0))  # Ensure no negative weights
        self.variations[syst.name] = variations
        return variations
    
    def get_1sigma_variations(self, syst:Systematic) -> dict[str, pl.Series]:
        variations = self.get_variations(syst)
        if 0 in syst.paramVariations:
            ref_idx = syst.paramVariations.index(0)
        else:
            print(f"No reference index (0) found in paramVariations for syst {syst.name}")
            ref_idx=None
        if 1 in syst.paramVariations:
            p1s_idx = syst.paramVariations.index(1)
        else:
            p1s_idx = ref_idx
            print(f"Warning: No +1 sigma index found in paramVariations for syst {syst.name}, using reference index instead.")
        if -1 in syst.paramVariations:
            m1s_idx = syst.paramVariations.index(-1)
        else:
            m1s_idx = ref_idx
            print(f"Warning: No -1 sigma index found in paramVariations for syst {syst.name}, using reference index instead.")
        if p1s_idx is not None:
            up_1sigma = variations[f'{syst.name}_{p1s_idx}']
        else:
            up_1sigma = pl.Series(name='up_1sigma', values=np.ones(variations.height))
        if m1s_idx is not None:
            down_1sigma = variations[f'{syst.name}_{m1s_idx}']
        else:
            down_1sigma = pl.Series(name='down_1sigma', values=np.ones(variations.height))
        if ref_idx is not None:
            ref = variations[f'{syst.name}_{ref_idx}']
        else:
            ref = pl.Series(name='ref', values=np.ones(variations.height))
        return pl.DataFrame({'up_1sigma': up_1sigma, 'down_1sigma': down_1sigma, 'ref': ref})
            

class SystAnalyzer:
    def __init__(self, caf_manager:CAFManager, systs_config:SystConfig, variation_manager:VariationManager) -> None:
        self.caf_manager = caf_manager
        self.systs_config = systs_config
        self.variation_manager = variation_manager

    def plot_1sigma_band(self, systs:list[str | Systematic] | str | Systematic, bin_var:pl.Expr, bins:np.ndarray, filter:pl.Expr = pl.lit(True), ax=None, ax_ratio=None, weighted:bool = True, oscillated:bool = True, show_ratio:bool = True) -> None:
        # Create figure with subplots if ax is None and show_ratio is True
        if ax is None and show_ratio:
            fig, (ax_main, ax_ratio_internal) = plt.subplots(2, 1, figsize=(10, 8), 
                                                     gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
            plt.sca(ax_main)
            # Use the internally created ratio axis if no external one provided
            if ax_ratio is None:
                ax_ratio = ax_ratio_internal
        elif ax is not None:
            plt.sca(ax)
            ax_main = ax
            # ax_ratio can be passed separately or be None
        else:
            fig = plt.figure()
            ax_main = plt.gca()
            # ax_ratio can be passed separately or be None
            
        sigma_up, sigma_down = None, None

        if type(systs) is not list:
            systs = [systs]

        data = self.caf_manager.data

        if len(systs) == 0:
            raise ValueError("systs list is empty")

        for systematic in systs:
            if type(systematic) is str:
                syst = self.systs_config.get_syst(systematic)
            elif type(systematic) is Systematic:
                syst = systematic
            else:
                raise ValueError("systs must be a list of strings or Systematic objects")
            
            variations = self.variation_manager.get_1sigma_variations(syst)
            local_data = pl.concat([data, variations], how='horizontal')
            local_data = local_data.filter(self.caf_manager.data_filter, filter)

            if weighted and oscillated:
                weights = local_data['final_oscillated_w']
            elif weighted and not oscillated:
                weights = local_data['nue_w'] + local_data['numu_w']
            else:
                weights = np.ones(len(local_data))

            binned_data = local_data.select(bin_var.alias('bin_var'))['bin_var']

            ref, _ = np.histogram(binned_data, bins=bins, weights=local_data['ref']*weights)
            up_1sigma, _ = np.histogram(binned_data, bins=bins, weights=local_data['up_1sigma']*weights)
            down_1sigma, _ = np.histogram(binned_data, bins=bins, weights=local_data['down_1sigma']*weights)

            pos_1sigma = np.maximum(np.maximum(ref - up_1sigma, ref - down_1sigma), 0)
            neg_1sigma = np.minimum(np.minimum(ref - up_1sigma, ref - down_1sigma), 0)

            if sigma_up is None:
                sigma_up = pos_1sigma**2
                sigma_down = neg_1sigma**2
            else:
                sigma_up = sigma_up + pos_1sigma**2
                sigma_down = sigma_down + neg_1sigma**2

        plt.sca(ax_main)
        plt.step(bins[:-1], ref, where='post', label='Nominal')
        # plt.step(edges[:-1], ref + sigma_up**0.5, where='post', label=r'1$\sigma$ up')
        # plt.step(edges[:-1], ref - sigma_down**0.5, where='post', label=r'1$\sigma$ down')
        color = plt.gca().lines[-1].get_color()
        plt.fill_between(bins[:-1], ref - sigma_down**0.5, ref + sigma_up**0.5, step='post', color=color, alpha=0.5, label=r'1$\sigma$ band')
        plt.ylabel('Events')
        
        # Add ratio panel if requested
        if ax_ratio is not None:
            plt.sca(ax_ratio)
            bin_centers = bins[:-1]
            
            # Calculate ratios, handling division by zero
            ratio_up = np.where(ref > 0, (ref + sigma_up**0.5) / ref, 1.0)
            ratio_down = np.where(ref > 0, (ref - sigma_down**0.5) / ref, 1.0)
            
            # Plot reference line and ratio band
            plt.axhline(y=1.0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
            plt.fill_between(bin_centers, ratio_down, ratio_up, step='post', color=color, alpha=0.5)
            
            plt.xlabel(bin_var)
            plt.ylabel('Ratio')
            plt.grid(True, alpha=0.3, linestyle=':')
            
            # Adjust y-limits for better visibility, handling NaN/Inf cases
            # Filter out NaN and Inf values for computing limits
            valid_down = ratio_down[np.isfinite(ratio_down)]
            valid_up = ratio_up[np.isfinite(ratio_up)]
            
            if len(valid_down) > 0 and len(valid_up) > 0:
                y_min = min(0.95, np.min(valid_down) - 0.05)
                y_max = max(1.05, np.max(valid_up) + 0.05)
                # Ensure limits are reasonable (not too extreme)
                y_min = max(y_min, 0.5)  # Don't go below 50%
                y_max = min(y_max, 1.5)  # Don't go above 150%
            else:
                # Default range if no valid data
                y_min, y_max = 0.8, 1.2
            
            # Final safety check
            if np.isfinite(y_min) and np.isfinite(y_max) and y_min < y_max:
                plt.ylim([y_min, y_max])
            else:
                plt.ylim([0.8, 1.2])  # Fallback to default range
            
            # Remove x-tick labels from top plot
            ax_main.set_xticklabels([])
            ax_main.set_xlabel('')
        else:
            plt.xlabel(bin_var)
            
        plt.sca(ax_main)

    def plot_2d_fractional_variation(self, systs:list[str | Systematic]| str | Systematic, bin_var_x:pl.Expr, bins_x:np.ndarray, bin_var_y:pl.Expr, bins_y:np.ndarray, filter:pl.Expr = pl.lit(True), ax=None, weighted:bool = True, oscillated:bool = True, which="up") -> None:
        if ax is not None:
            plt.sca(ax)
        else:
            fig = plt.figure()

        if which not in ["up", "down"]:
            raise ValueError("which must be 'up' or 'down'")
        
        if type(systs) is not list:
            systs = [systs]
        
        sigma_up, sigma_down = None, None
        for syst in systs:

            if type(syst) is str:
                syst = self.systs_config.get_syst(syst)
            elif type(syst) is Systematic:
                pass
            else:
                raise ValueError("systs must be a list of strings or Systematic objects")

            data = self.caf_manager.data

            variations = self.variation_manager.get_1sigma_variations(syst)
            local_data = pl.concat([data, variations], how='horizontal')
            local_data = local_data.filter(self.caf_manager.data_filter, filter)

            if weighted and oscillated:
                weights = local_data['final_oscillated_w']
            elif weighted and not oscillated:
                weights = local_data['nue_w'] + local_data['numu_w']
            else:
                weights = np.ones(len(local_data))

            binned_data_x = local_data.select(bin_var_x.alias('bin_var_x'))['bin_var_x']
            binned_data_y = local_data.select(bin_var_y.alias('bin_var_y'))['bin_var_y']

            ref, _, _ = np.histogram2d(binned_data_x, binned_data_y, bins=[bins_x, bins_y], weights=local_data['ref']*weights)
            up_1sigma, _, _ = np.histogram2d(binned_data_x, binned_data_y, bins=[bins_x, bins_y], weights=local_data['up_1sigma']*weights)
            down_1sigma, _, _ = np.histogram2d(binned_data_x, binned_data_y, bins=[bins_x, bins_y], weights=local_data['down_1sigma']*weights)

            pos_1sigma = np.maximum(np.maximum(ref - up_1sigma, ref - down_1sigma), 0)
            neg_1sigma = np.minimum(np.minimum(ref - up_1sigma, ref - down_1sigma), 0)
            if sigma_up is None:
                sigma_up = pos_1sigma**2
                sigma_down = neg_1sigma**2
            else:
                sigma_up = sigma_up + pos_1sigma**2
                sigma_down = sigma_down + neg_1sigma**2
        frac_up = (sigma_up**0.5/ ref)*100
        frac_down = (sigma_down**0.5/ ref)*100

        if which == "up":
            frac = frac_up
            cmap = 'Reds'
        else:
            frac = -frac_down
            cmap = 'Blues_r'


        X, Y = np.meshgrid(bins_x[:-1], bins_y[:-1], indexing='ij')
        pcm = plt.pcolormesh(X, Y, frac, shading='auto', cmap=cmap)#, vmin=-np.nanmax(np.abs(frac)), vmax=np.nanmax(np.abs(frac)))
        plt.colorbar(pcm, label='Variation (%)')
        plt.xlabel(bin_var_x)
        plt.ylabel(bin_var_y)
        # plt.title(rf'{syst.name} {which}1$\sigma$ Variation')
        # plt.tight_layout()

    def plot_event_rate(self, syst_values:dict[Systematic|str, float], bin_var:pl.Expr, bins:np.ndarray, filter:pl.Expr = pl.lit(True), ax=None, weighted:bool = True, oscillated:bool = True) -> None:
        if ax is not None:
            plt.sca(ax)
        else:
            fig = plt.figure()

        data = self.caf_manager.data
        local_data = data.filter(self.caf_manager.data_filter, filter)

        for syst, value in syst_values.items():
            if type(syst) is str:
                syst = self.systs_config.get_syst(syst)
            elif type(syst) is Systematic:
                pass
            else:
                raise ValueError("systs must be a list of strings or Systematic objects")

            variations = self.variation_manager.get_variations(syst)
            if value not in syst.paramVariations:
                print(f"Warning: value {value} not in paramVariations for syst {syst.name}, using closest available value.")
                closest_value = min(syst.paramVariations, key=lambda x: abs(x - value))
                value = closest_value
            value_idx = syst.paramVariations.index(value)
            weight_col = f'{syst.name}_{value_idx}'
            local_data = pl.concat([local_data, variations.select(pl.col(weight_col).alias(f'{syst.name}_weight'))], how='horizontal')
            local_data = local_data.with_columns(
                (pl.col('nue_w') * pl.col(f'{syst.name}_weight')).alias('nue_w'),
                (pl.col('numu_w') * pl.col(f'{syst.name}_weight')).alias('numu_w'),
            )

        if weighted and oscillated:
            weights = local_data['nue_w']*local_data['osc_from_e_w'] + local_data['numu_w']*local_data['osc_from_mu_w']
        elif weighted and not oscillated:
            weights = local_data['nue_w'] + local_data['numu_w']
        else:
            weights = np.ones(len(local_data))

        binned_data = local_data.select(bin_var.alias('bin_var'))['bin_var']

        hist, _ = np.histogram(binned_data, bins=bins, weights=weights)

        plt.step(bins[:-1], hist, where='post', label='Event Rate')
        plt.xlabel(bin_var)
        plt.ylabel('Events')

class MultiPlot:
    def __init__(self, analyzer:SystAnalyzer, bin_var:pl.Expr, bins:np.ndarray, ax=None) -> None:
        self.colors = []
        self.labels = []
        self.handles = []
        self.bin_var = bin_var
        self.bins = bins
        self.analyzer = analyzer
        
        # Check if ax is a tuple of (main_ax, ratio_ax) or a single axis
        if ax is not None:
            if isinstance(ax, tuple) and len(ax) == 2:
                self.ax_main = ax[0]
                self.ax_ratio = ax[1]
            else:
                self.ax_main = ax
                self.ax_ratio = None
        else:
            _, self.ax_main = plt.subplots()
            self.ax_ratio = None

    def add_syst_band(self, systs, label:str, filter:pl.Expr=pl.lit(True), oscillated:bool = True) -> None:
        # Pass ax_ratio explicitly to plot_1sigma_band
        self.analyzer.plot_1sigma_band(systs, bin_var=self.bin_var, bins=self.bins, filter=filter, 
                                      ax=self.ax_main, ax_ratio=self.ax_ratio, 
                                      oscillated=oscillated, show_ratio=False)
        
        color = self.ax_main.lines[-1].get_color()
        self.colors.append(color)
        self.labels.append(label)
        line = plt.Line2D([0], [0], color=color, alpha=1)
        band = plt.Rectangle((0,0),1,1, color=color, alpha=0.5)
        self.handles.append((line, band))

    def legend(self):
        plt.sca(self.ax_main)
        plt.legend(self.handles, self.labels)


def create_oscillation_yaml_from_systematics(systematics: list[Systematic], 
                                             output_file: Optional[str] = None,
                                             sample_names: Optional[list[str]] = None,
                                             error: float = 1.0,
                                             flat_prior: bool = False,
                                             param_bounds: Optional[list[float]] = None,
                                             generated_value: float = 0.0,
                                             prefit_value: float = 0.0,
                                             spline_mode: Optional[list[int]] = None,
                                             spline_type: str = "PerEvent",
                                             mcmc_step_scale: float = 0.4,
                                             syst_type: str = "Spline",
                                             param_group: str = "Xsec") -> str:
    """
    Create a YAML configuration for the oscillation fitter from a list of Systematic objects.
    
    Parameters:
    -----------
    systematics : list[Systematic]
        List of Systematic objects to convert
    output_file : str, optional
        Path to save the YAML file. If None, returns the YAML string.
    sample_names : list[str], optional
        List of sample names. Default is ["ND_*", "FD_*", "ATM"]
    error : float, optional
        Error value for the systematic. Default is 1.0
    flat_prior : bool, optional
        Whether to use a flat prior. Default is False
    param_bounds : list[float], optional
        Parameter bounds [min, max]. Default is [-9999.0, 9999.0]
    generated_value : float, optional
        Generated parameter value. Default is 0.0
    prefit_value : float, optional
        Pre-fit parameter value. Default is 0.0
    spline_mode : list[int], optional
        Mode for spline. Default is [0]
    spline_type : str, optional
        Type of spline. Default is "PerEvent"
    mcmc_step_scale : float, optional
        MCMC step scale. Default is 0.4
    syst_type : str, optional
        Type of systematic. Default is "Spline"
    param_group : str, optional
        Parameter group. Default is "Xsec"
        
    Returns:
    --------
    str : YAML formatted string
    
    Example:
    --------
    >>> config = SystConfig('systs.fcl')
    >>> systematics = config.get_all_systs()
    >>> yaml_str = create_oscillation_yaml_from_systematics(
    ...     systematics, 
    ...     output_file='oscillation_config.yaml'
    ... )
    """
    # Set defaults
    if sample_names is None:
        sample_names = ["ND_*", "FD_*", "ATM"]
    if param_bounds is None:
        param_bounds = [-9999.0, 9999.0]
    if spline_mode is None:
        spline_mode = [0]
    
    # Build the YAML structure
    yaml_list = []
    
    for syst in systematics:
        syst_dict = {
            'Systematic': {
                'SampleNames': sample_names,
                'Error': error,
                'FlatPrior': flat_prior,
                'Names': {
                    'FancyName': syst.name,
                    'ParameterName': syst.name
                },
                'ParameterBounds': param_bounds,
                'ParameterValues': {
                    'Generated': generated_value if not syst.isCorrection else syst.cv,
                    'PreFitValue': prefit_value if not syst.isCorrection else syst.cv
                },
                'SplineInformation': {
                    'Mode': spline_mode,
                    'SplineName': syst.name,
                    'Type': spline_type
                },
                'StepScale': {
                    'MCMC': mcmc_step_scale
                },
                'Type': syst_type,
                'ParameterGroup': param_group
            }
        }
        yaml_list.append(syst_dict)
    
    # Convert to YAML string
    yaml_str = yaml.dump(yaml_list, default_flow_style=False, sort_keys=False)
    
    # Save to file if specified
    if output_file:
        with open(output_file, 'w') as f:
            f.write(yaml_str)
        print(f"YAML configuration saved to {output_file}")
    
    return yaml_str

def add_headers_to_root_file(systs:SystConfig, fname:str) -> None:
    systs_data = {'name':[], 'cv': [], 'variations': [], 'id': [], 'isCorrection': []}

    for syst in systs.get_all_systs():
        systs_data['name'].append(syst.name)
        systs_data['cv'].append(syst.cv)
        systs_data['variations'].append(syst.paramVariations)
        systs_data['id'].append(syst.id)
        systs_data['isCorrection'].append(syst.isCorrection)

    import pandas as pd
    with uproot.update(fname) as f:
        f['systsHeader'] = pd.DataFrame(systs_data)
    print(f"Added systsHeader to {fname}")

def create_dummy_splines(cafs:CAFManager, fname:str='dummy_splines.root') -> None:
    dummy_spline = cafs.data.with_columns(
        DummySpline = pl.lit([1., 2., 4., 10.])
    ).select("DummySpline")
    dummy_spline_woccqe = cafs.data.with_columns(
        DummySplineNoCCQE = pl.when(
            pl.col('mode') != 0
            ).then(pl.lit([1., 2., 4., 10.])).otherwise(pl.lit([1., 1., 1., 1.])
            )
    ).select("DummySplineNoCCQE")
    all_splines = pl.concat([dummy_spline, dummy_spline_woccqe], how='horizontal')

    systs_data = {
        'name':['DummySpline', 'DummySplineNoCCQE'],
        'cv': [0., 0.],
        'variations': [[0., 1., 2., 3.], [0., 1., 2., 3.]],
        'id': [0, 1],
        'isCorrection': [False, False]}

    import pandas as pd
    with uproot.recreate('dummy_splines.root') as f:
        f['systsHeader'] = pd.DataFrame(systs_data)
        f['SystWeights'] = all_splines.to_pandas()
    print(f"Created dummy splines file: {fname}")
