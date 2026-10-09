#!/usr/bin/env python
# coding: utf-8

# standard library
import argparse
import os
import re
import io
import warnings
from functools import lru_cache

warnings.filterwarnings("ignore")

#numpy, scipy, pandas
import numpy as np
import pandas as pd
import scipy
from scipy.integrate import dblquad

#matplotlib
import matplotlib
import matplotlib.colors as mc
import matplotlib.tri as tri
#get_ipython().run_line_magic('matplotlib', 'inline')
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import patches
from matplotlib import style

#astropy
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.constants as const

#astroquery
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
from astroquery.mast import Observations
from astroquery.mast import Tesscut

#other
from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry
from lightkurve import search_targetpixelfile
from alive_progress import alive_bar

from uncertainties import ufloat
from uncertainties import umath_core as um
from functools import lru_cache


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Run the TESS ON/OFF contamination analysis.")
    parser.add_argument("--target", type=str, default=None,
                        help="Target folder or TOI name, e.g. TOI1011.01 or a path.")
    parser.add_argument("--custom-aperture", type=int, default=0,
                        help="Use a custom circular aperture radius in pixels. 0 means default behavior.")
    parser.add_argument("--sector", type=int, default=None,
                        help="Override the TESS sector used for the TPF/Tesscut query. Default keeps the target metadata behavior.")
    parser.add_argument("--dark", action="store_true",
                        help="Use the dark plotting style.")
    parser.add_argument("--output-dir", type=str, default="on-off",
                        help="Directory where output plots and files are saved.")
    return parser.parse_args(argv)


plt.rcParams.update({'font.size': 20})


DEFAULT_CONFIG = {
    "output_dir": "on-off",
    "dark": 0,
    "custom_aperture": 0,
    "sector": None,
    "search_radius_pixels": 10,
    "pixel_size_arcsec": 21,
    "prf_col_offset": 44,
    "prf_image_size": 13,
    "prf_fits_base_url": "https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/",
}


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def log_stage(stage, message=""):
    """Print a concise progress update for the user while the analysis runs."""
    if message:
        print(f"[progress] {stage}: {message}")
    else:
        print(f"[progress] {stage}")


def safe_target_name(target_arg=None):
    """Return the target folder name or the provided name without trailing slashes."""
    if target_arg is None:
        target_arg = os.getcwd()
    return os.path.basename(os.path.normpath(target_arg))


def load_toi_catalog(csv_path=None):
    """Download and cache the TOI catalog once, then reuse it for repeated runs."""
    if csv_path is None:
        csv_path = os.environ.get(
            "TOI_CATALOG_PATH",
            os.path.join(os.path.expanduser("~"), ".toi-contam", "toi_catalog.csv"),
        )

    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    if not os.path.exists(csv_path):
        import urllib.request
        url_tess = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
        urllib.request.urlretrieve(url_tess, csv_path)
    return pd.read_csv(csv_path, low_memory=False, delimiter=',', index_col=1)


def load_target_metadata(target_name=None, toi_df=None):
    """Return a small dict describing the target and its TOI metadata."""
    if toi_df is None:
        toi_df = load_toi_catalog()

    target_name = safe_target_name(target_name)
    metadata = {
        "name": target_name,
        "toi": False,
        "toi_number": None,
        "ra": None,
        "dec": None,
        "source": None,
        "period": None,
        "duration": None,
        "t0": None,
        "sector_last": None,
    }

    if str(target_name[0:3]) == 'TOI':
        metadata["toi"] = True
        name2 = str(target_name.split("I")[1])
        metadata["toi_number"] = float(name2)
        row = toi_df.loc[metadata["toi_number"]]
        metadata["ra"] = row['RA']
        metadata["dec"] = row['Dec']
        metadata["source"] = row['Source']
        metadata["period"] = row['Period (days)']
        metadata["duration"] = row['Duration (hours)'] / 24
        metadata["t0"] = row['Epoch (BJD)'] - 2457000
        sector_last0 = row['Sectors']
        metadata["sector_last"] = [int(s) for s in re.findall(r'\d+', sector_last0)][-1]

    return metadata


def query_gaia_catalog(coord, radius):
    """Convenience wrapper for the Gaia cone search."""
    Gaia.ROW_LIMIT = -1
    return (Gaia.cone_search(coordinate=coord, radius=radius,
                            table_name='gaiadr3.gaia_source')).get_results()


def query_tic_catalog(coord, radius):
    """Convenience wrapper for the TIC query used by the contamination analysis."""
    return Catalogs.query_region(coord, catalog="TIC", radius=radius)


def prepare_target_context(target_arg=None, toi_df=None):
    """Load the target metadata and convert it to a coordinate + metadata dict."""
    if toi_df is None:
        toi_df = load_toi_catalog()

    target_name = safe_target_name(target_arg)
    metadata = load_target_metadata(target_name, toi_df)
    validate_target_metadata(metadata)
    coord = SkyCoord(
        ra=metadata["ra"], dec=metadata["dec"],
        unit=(u.hourangle, u.degree), frame='icrs'
    ) if metadata["ra"] is not None else None

    return {
        "name": target_name,
        "metadata": metadata,
        "coord": coord,
    }


def build_catalog_tables(coord, query_radius):
    """Query Gaia and TIC around the target and return pandas DataFrames."""
    gaia_table = query_gaia_catalog(coord, query_radius)
    tess_table = query_tic_catalog(coord, query_radius)
    df_gaia = gaia_table.to_pandas()
    df_gaia = df_gaia.dropna(subset=["phot_g_mean_mag"]).reset_index(drop=True)
    df_tess = tess_table.to_pandas()
    return {
        "gaia_table": gaia_table,
        "df_gaia": df_gaia,
        "df_tess": df_tess,
    }


def resolve_aperture_choice(source, custom_aperture=None):
    """Resolve a user-facing aperture choice with a safe fallback."""
    if custom_aperture is not None and custom_aperture > 0:
        return {"source": "mine", "custom_aperture": custom_aperture}
    else:
        return {"source": source, "custom_aperture": None}


def normalize_aperture_radius(aperture_radius):
    """Normalize the custom aperture radius to a positive integer or None."""
    if aperture_radius is None:
        return None
    try:
        radius = int(aperture_radius)
    except (TypeError, ValueError):
        raise ValueError(f"Aperture radius must be an integer, got {aperture_radius!r}.")
    if radius <= 0:
        raise ValueError("Aperture radius must be positive.")
    return radius


def validate_target_metadata(metadata):
    """Validate that the target metadata contains enough information to proceed."""
    if metadata["ra"] is None or metadata["dec"] is None:
        raise ValueError(f"Target metadata is incomplete for {metadata['name']!r}.")


def derive_target_arrays(df_gaia):
    """Return the Gaia and TESS-relevant arrays used throughout the analysis."""
    gmag = df_gaia['phot_g_mean_mag']
    gmag_err = 2.5 * np.log10(1 + 1 / df_gaia['phot_g_mean_flux_over_error'])
    rra = df_gaia['ra']
    ddec = df_gaia['dec']
    dist = df_gaia['dist'] * 3600
    return {
        "gmag": gmag,
        "gmag_err": gmag_err,
        "rra": rra,
        "ddec": ddec,
        "dist": dist,
    }


def compute_depth_info(toi_df, toi_number, toi_flag, fallback_file="depth"):
    """Resolve the transit depth and uncertainty for either TOI or local fallback input."""
    if toi_flag:
        depth = toi_df.loc[toi_number]['Depth (ppm)'] / (10 ** 6)
        e_depth = toi_df.loc[toi_number]['Depth (ppm) err'] / (10 ** 6)
        return {
            "depth": depth,
            "e_depth": e_depth,
            "depth_pdcsap": depth,
            "e_depth_pdcsap": e_depth,
        }

    headings = ['depth', 'com']
    info = pd.read_csv(fallback_file, delimiter=" ", skiprows=0, names=headings, engine='python', index_col=None)
    depth = info['depth'][0]
    return {
        "depth": depth,
        "e_depth": None,
        "depth_pdcsap": depth,
        "e_depth_pdcsap": None,
    }


"""
From "spacetelescope/notebooks/MAST/TESS": retrieve TESS Pixel Response Function images.
These helpers are kept at module scope so the heavy PRF logic is reusable and the main
analysis can stay focused on orchestration.
"""

def path_lookup(ccd, camera, sector):
    """Return the datestring and subdirectory for the selected TESS PRF set."""
    if sector < 1:
        raise ValueError("Sector must be greater than 0.")
    if (camera > 4) | (ccd > 4):
        raise ValueError("Camera or CCD is larger than 4.")

    if sector <= 3:
        add_path = "/start_s0001/"
        datestring = "2018243163600"

        if camera >= 3:
            datestring = "2018243163601"
        elif (camera == 2) & (ccd == 4):
            datestring = "2018243163601"
    else:
        add_path = "/start_s0004/"
        datestring = "2019107181900"
        if (camera == 1) & (ccd >= 2):
            datestring = "2019107181901"
        elif (camera == 2):
            datestring = "2019107181901"
        elif (camera == 3):
            if ccd >= 2:
                datestring = "2019107181902"
            else:
                datestring = "2019107181901"
        elif (camera == 4):
            datestring = "2019107181902"

    return datestring, add_path


def read_one_prf_fits_file(ccd, camera, col, row, path, datestring):
    """Read the interleaved PRF array for one sampled row/col location."""
    fn = "cam%u_ccd%u/tess%13s-prf-%1u-%1u-row%04u-col%04u.fits" % \
        (camera, ccd, datestring, camera, ccd, row, col)
    filepath = os.path.join(path, fn)
    hdulistObj = fits.open(filepath)
    return hdulistObj[0].data


def determineFourClosestPrfLoc(col, row):
    """Return the four nearby PRF sample locations around a target position."""
    posRows = np.array([1, 513, 1025, 1536, 2048])
    posCols = np.array([45, 557, 1069, 1580, 2092])

    difcol = np.abs(posCols - col)
    difrow = np.abs(posRows - row)

    imagePos = []
    for r in posRows[np.argsort(difrow)[0:2]]:
        for c in posCols[np.argsort(difcol)[0:2]]:
            imagePos.append((c, r))

    return imagePos


def getOffsetsFromPixelFractions(col, row):
    """Convert sub-pixel position into the offset used in the interleaved PRF array."""
    gridSize = 9
    colFrac = np.remainder(float(col), 1)
    rowFrac = np.remainder(float(row), 1)
    colOffset = gridSize - np.round(gridSize * colFrac) - 1
    rowOffset = gridSize - np.round(gridSize * rowFrac) - 1
    return int(colOffset), int(rowOffset)


def getRegSampledPrfFitsByOffset(prfArray, colOffset, rowOffset):
    """Extract the regular 13x13 PRF image from the interleaved array."""
    gridSize = 9
    assert colOffset < gridSize
    assert rowOffset < gridSize

    nColOut, nRowOut = prfArray.shape
    nColOut /= float(gridSize)
    nRowOut /= float(gridSize)

    iCol = colOffset + (np.arange(nColOut) * gridSize).astype(int)
    iRow = rowOffset + (np.arange(nRowOut) * gridSize).astype(int)

    tmp = prfArray[iRow, :]
    return tmp[:, iCol]


def interpolatePrf(regPrfArray, col, row, imagePos):
    """Linearly interpolate a PRF from the four surrounding PRF grid points."""
    p11, p21, p12, p22 = regPrfArray
    c0 = imagePos[0][0]
    c1 = imagePos[1][0]
    r0 = imagePos[0][1]
    r1 = imagePos[2][1]

    assert c0 != c1
    assert r0 != r1

    dCol = (col - c0) / (c1 - c0)
    dRow = (row - r0) / (r1 - r0)

    tmp1 = p11 + (p21 - p11) * dCol
    tmp2 = p12 + (p22 - p12) * dCol
    return tmp1 + (tmp2 - tmp1) * dRow


@lru_cache(maxsize=None)
def get_prf_at_col_row_fits(col, row, ccd, camera, sector, path):
    """Return a 13x13 interpolated PRF image evaluated at a given CCD position."""
    col = float(col)
    row = float(row)
    prfImages = []
    datestring, subDirectory = path_lookup(ccd, camera, sector)
    prf_path = path + subDirectory
    colOffset, rowOffset = getOffsetsFromPixelFractions(col, row)
    imagePos = determineFourClosestPrfLoc(col, row)

    for pos in imagePos:
        prfArray = read_one_prf_fits_file(ccd, camera, pos[0], pos[1], prf_path, datestring)
        img = getRegSampledPrfFitsByOffset(prfArray, colOffset, rowOffset)
        prfImages.append(img)

    return interpolatePrf(prfImages, col, row, imagePos)

##################################################################################################
########################################### code start ###########################################
##################################################################################################


def main(argv=None):
    global OUTPUT_DIR, CONFIG, dark

    args = parse_args(argv)
    OUTPUT_DIR = args.output_dir
    createFolder(OUTPUT_DIR)

    sector_for_output = args.sector if args.sector is not None else None
    if sector_for_output is not None:
        OUTPUT_DIR = os.path.join(args.output_dir, f"sector_{sector_for_output}")
        createFolder(OUTPUT_DIR)

    CONFIG = {
        **DEFAULT_CONFIG,
        "dark": int(args.dark),
        "output_dir": args.output_dir,
        "custom_aperture": args.custom_aperture,
        "sector": args.sector,
    }

    dark = CONFIG["dark"]
    if dark == 1:
        style.use('dark_background')
    else:
        style.use('classic')
        
    url = os.getcwd() if args.target is None else args.target

    print("\n" + "=" * 80)
    print("TESS contamination analysis")
    print(f"Target: {safe_target_name(url)}")
    print("=" * 80 + "\n")

    TOI_df = load_toi_catalog()
    target_context = prepare_target_context(url, TOI_df)
    #name = target_context["name"]
    metadata = target_context["metadata"]
    coord = target_context["coord"]

    log_stage("Target metadata ready", f"TOI={metadata['toi']} | source={metadata['source']}")

    TOI = metadata["toi"]
    ra = metadata["ra"]
    dec = metadata["dec"]
    source = metadata["source"]
    #period = metadata["period"]
    #duration = metadata["duration"]
    #T0 = metadata["t0"]
    sector_last = metadata["sector_last"]
    sector_choice = args.sector if args.sector is not None else sector_last
    if args.sector is not None:
        OUTPUT_DIR = os.path.join(args.output_dir, f"sector_{args.sector}")
        createFolder(OUTPUT_DIR)
    toi_number = metadata["toi_number"]
    
    if args.sector is not None:
            available_sectors = []
            if metadata.get("source") is not None:
                available_sectors = [int(s) for s in re.findall(r'\d+', str(TOI_df.loc[toi_number]['Sectors']))]
            if not available_sectors:
                available_sectors = [sector_last]
            if int(args.sector) not in available_sectors:
                raise ValueError(
                    f"Requested sector {args.sector} is not available for target {safe_target_name(url)}. "
                    f"Available sectors: {sorted(available_sectors)}"
                )

    #RAd = coord.ra.deg; Decd = coord.dec.deg

    ############################################################################################
    #Query to find the stars within a circle of radius configured in DEFAULT_CONFIG
    pixel_size = CONFIG["pixel_size_arcsec"] * u.arcsec
    search_radius = CONFIG["search_radius_pixels"]
    ############################################################################################
    query_radius = search_radius * pixel_size

    log_stage("Querying nearby catalogs", f"search radius = {search_radius} TESS pixels")
    catalog_tables = build_catalog_tables(coord, query_radius)
    #r = catalog_tables["gaia_table"]
    df_tess = catalog_tables["df_tess"]
    df_gaia = catalog_tables["df_gaia"]
    log_stage("Catalog query complete", f"Gaia matches = {len(df_gaia)}")

    _COEFFICIENTS = np.array([
        # Camera 1: CCD 1, 2, 3, 4; coefficients a, b, c, d.
        (-0.006325018017572858, 0.0889782745534607, -0.6104603268545619, -20.546580032142504),
        (-0.006300834719139644, 0.09014323129367145, -0.6091236433530977, -20.538130856096615),
        (-0.0063271807682276264, 0.09023001362272956, -0.6095549353074154, -20.540606026864292),
        (-0.005880867575214354, 0.0879901582068267, -0.603535919577106, -20.53782788906888),
        # Camera 2: CCD 1, 2, 3, 4; coefficients a, b, c, d.
        (-0.006840581638144296, 0.09124905006981088, -0.6062625516285151, -20.51858394290852),
        (-0.006739762661601594, 0.08994903273633076, -0.6069026067612524, -20.53266719802872),
        (-0.005845075772658837, 0.08717722091862662, -0.6034964033996696, -20.515385022167607),
        (-0.007751648779293352, 0.09544503915038027, -0.6151103850109424, -20.5097343358954),
        # Camera 3: CCD 1, 2, 3, 4; coefficients a, b, c, d.
        (-0.007168914513397514, 0.09431067790449457, -0.6103959113077466, -20.539797560363578),
        (-0.006771613368979871, 0.09131957273441386, -0.6098902707363071, -20.547479768294163),
        (-0.007264899731961687, 0.09454191554476432, -0.6177355344230424, -20.540218471786318),
        (-0.00952158257773152, 0.10500988473748511, -0.6278561780069848, -20.531827088481958),
        # Camera 4: CCD 1, 2, 3, 4; coefficients a, b, c, d.
        (-0.007470316306022466, 0.09702893320666926, -0.6249704460570471, -20.51461122808701),
        (-0.00836372384377166, 0.10033144705589458, -0.6274354980638552, -20.513047834467667),
        (-0.00810923844177131, 0.09869085435044767, -0.6184637948861488, -20.527549078817017),
        (-0.010647976179780912, 0.10761667299424135, -0.6307008563977301, -20.52866458629858),
    ], dtype=float).reshape(4, 4, 4)


    def _coefficients(cam, ccd):
        if cam not in (1, 2, 3, 4) or ccd not in (1, 2, 3, 4):
            raise ValueError("Camera and CCD must be integers in 1..4")
        return _COEFFICIENTS[int(cam) - 1, int(ccd) - 1]

    def magneiinst(cam, ccd, gnei, bpnei, rpnei):
        """Mantiene il polinomio e il fallback G + zero point originali."""
        a, b, c, d = _coefficients(cam, ccd)
        g, bp, rp = np.broadcast_arrays(gnei, bpnei, rpnei)
        color = bp - rp
        valid = np.isfinite(bp) & np.isfinite(rp) & (bp > 0) & (rp > 0)
        return np.where(valid, a * color**3 + b * color**2 + c * color + d + g, g + d)

    def GetTmagCalib(tinst, cam, ccd):
        """Rimuove lo zero point strumentale come nella versione originale."""
        return np.asarray(tinst) - _coefficients(cam, ccd)[3]

    target_arrays = derive_target_arrays(df_gaia)
    gmag = target_arrays["gmag"]
    #gmag_err = target_arrays["gmag_err"]
    rra = target_arrays["rra"]
    ddec = target_arrays["ddec"]
    dist = target_arrays["dist"]

    ################################################################################################
    #################################### Transit DEPTH (from TIC) ##################################
    ################################################################################################
    depth_info = compute_depth_info(TOI_df, toi_number, TOI)
    depth = depth_info["depth"]
    e_depth = depth_info["e_depth"]
    depth_pdcsap = depth_info["depth_pdcsap"]
    e_depth_pdcsap = depth_info["e_depth_pdcsap"]


    ################################################################################################
    #################################### Target Pixel file (tpf) ###################################
    ################################################################################################
    aperture_choice = resolve_aperture_choice(source, CONFIG["custom_aperture"])
    source = aperture_choice["source"]
    radius_input = normalize_aperture_radius(aperture_choice["custom_aperture"])

    if TOI and source!= 'qlp' and source!='mine':
        search_kwargs = {"sector": sector_choice} if sector_choice is not None else {}
        tpf_first = search_targetpixelfile('TIC '+str(TOI_df.loc[toi_number]['TIC ID']), **search_kwargs)
        tpf = tpf_first[0].download() # uses the default sector unless the user overrides it
        #The aperture extension header contains the same WCS as that in the Pixels extension.
        wcs = WCS(tpf.hdu[2].header)
    else: 
        """catalog_data = Catalogs.query_region(str(RAd)+' '+str(Decd), catalog="TIC", radius = .01)"""  
        tpf = False
        tesscut_kwargs = {"coordinates": coord, "size": 3}
        if sector_choice is not None:
            tesscut_kwargs["sector"] = sector_choice
        hdulist = Tesscut.get_cutouts(**tesscut_kwargs)
        hdu1 = hdulist[0] #hdu
        wcs = WCS(hdu1[2].header) #wcs
        
        
    ###########

    # Pull out the location of the middle of the CCD from the Physics WCS in the header.
    if tpf:
        image_head = tpf.hdu[1].header 
        prime_head = tpf.hdu[0].header
    else:
        n = 0 
        image_head = hdulist[n][1].header
        prime_head = hdulist[n][0].header
        #ap_head = hdulist[n][2].header

    camera = prime_head['CAMERA'] 
    ccd = prime_head['CCD'] 
    sector = prime_head['SECTOR']
    
    log_stage("Sector used", f"{(sector)}")

    Tmag = magneiinst(camera, ccd, df_gaia['phot_g_mean_mag'], df_gaia['phot_bp_mean_mag'], df_gaia['phot_rp_mean_mag'])  #TESS magnitudes
    Tmag_err = Tmag - Tmag + 0.006 

    ###########
        
    #Coordinates
    ra = rra; dec = ddec #Gaia

    #coordinates in pixel (conversion from RA, Dec)
    pix_coords = []; pix_coord = np.zeros([len(ra), 2])
    for i in range(len(ra)):
        RApix = np.ndarray.item(wcs.all_world2pix(ra[i], dec[i], 0)[0])
        Decpix = np.ndarray.item(wcs.all_world2pix(ra[i], dec[i], 0)[1])
        pix_coord[i, 0] = RApix; pix_coord[i, 1] = Decpix
    pix_coords.append(pix_coord)

    ##################################################################################################
    ################################### Aperture comparison images ###################################
    ##################################################################################################
    with plt.style.context('default'):
        
        #Figure: Aperture TESS
        plt.rcParams.update({'font.size': 20})
        #fig, ax = plt.subplots(1, 1, sharey=True, figsize = (15,9))
        fig = plt.figure(figsize = (15,9))

        if tpf:
            tpf_data = tpf.hdu[1].data 
            first_image = tpf_data['FLUX'][1] #data
            # print(tpf.hdu[1].header['CROWDSAP'], 1/(1+cont_ratio))
            # print(tpf.hdu[1].header['FLFRCSAP'])
            s = tpf.flux.shape 
            aperture_mask=tpf.pipeline_mask #mask from TPF
        else:
            first_image = hdu1[1].data['FLUX'][0] #flux
            positions = [(pix_coords[0][0,0], pix_coords[0][0,1])]
            if source == 'qlp':
            
                """ to find optimal aperture radius; from Huang et al. 2020"""

                tic_ids = TOI_df.loc[toi_number]['TIC ID']
                sectors =  TOI_df.loc[toi_number]['Sectors']
                numbers = [int(s) for s in re.findall(r'\d+',sectors)]
                # print(numbers[0])
                try:
                    obsTable = Observations.query_criteria(provenance_name="QLP",
                                           target_name=tic_ids,
                                           sequence_number=numbers[0])
                    data = Observations.get_product_list(obsTable[0])
                except IndexError:
                    obsTable = Observations.query_criteria(provenance_name="QLP",
                                           target_name=tic_ids,
                                           sequence_number=1) #TODO change
                    data = Observations.get_product_list(obsTable[0])
                    
                download_lc = Observations.download_products(data)
                where = 'mastDownload/HLSP/'+str(data[0]['obs_id'])+'/'+str(data[0]['obs_id'])+'.fits'
                lcfile = where
                lchdu = fits.open(lcfile)
                rad_ap = lchdu[1].header['BESTAP'][0]
                # print(rad_ap)
                aperture = CircularAperture(positions, r=int(rad_ap))
            else:
                if radius_input is None:
                    radius_input = int(input("Give me custom aperture radius: "))
                aperture = CircularAperture(positions, r=int(radius_input)) #custom value

        lower = np.percentile(first_image,1) 
        upper = np.percentile(first_image,99) #1%-99% flux

        ax = plt.subplot(111, projection = wcs) #axis
        plt.xlabel('Right Ascension') 
        plt.ylabel('Declination')
        plt.grid(axis = 'both', color = 'white', ls = 'solid')

        plt.xlim(pix_coords[0][0,0]-3,pix_coords[0][0,0]+3) 
        plt.ylim(pix_coords[0][0,1]-3,pix_coords[0][0,1]+3) #4 TESS pixels

        plt.imshow(first_image, origin = 'lower', cmap = plt.cm.viridis,clim=(lower,upper))
        
        ap_color = 'white'

        if tpf:
            #TESS aperture
            for i in range(s[1]):
                for j in range(s[2]):
                    if aperture_mask[i, j]:
                        ax.plot([j-0.5, j+0.5],[i-0.5, i-0.5],
                            transform=ax.get_transform(wcs), color=ap_color, zorder=1)
                        ax.plot([j-0.5, j+0.5],[i+0.5, i+0.5],\
                            transform=ax.get_transform(wcs), color=ap_color, zorder=1)
                        ax.plot([j-0.5, j-0.5],[i-0.5, i+0.5],\
                            transform=ax.get_transform(wcs), color=ap_color, zorder=1)
                        ax.plot([j+0.5, j+0.5],[i-0.5, i+0.5],\
                            transform=ax.get_transform(wcs), color=ap_color, zorder=1)
        else:
            #aperture
            aperture.plot(color=ap_color, lw=2, transform=ax.get_transform(wcs))

        #Gaia stars (Target + contaminants)
        for j in range(len(rra)):
            size = 1024.0 / 2**((gmag[j]-gmag[0])/2)
            if j == 0:
                ax.scatter(rra[j], ddec[j], transform=ax.get_transform('fk5'), s=size, zorder = 4,
                           edgecolor='black', facecolor='yellow', linewidth =1, label = 'Target star',\
                           marker = "$t$")
                ax.scatter(rra[j], ddec[j], transform=ax.get_transform('fk5'), s=size,
                           edgecolor='black', facecolor='red', linewidth =1, zorder = 3, alpha = 0.7)
            else:
                ax.scatter(rra[j], ddec[j], transform=ax.get_transform('fk5'), s=size,
                           edgecolor='black', facecolor='red', linewidth =1, zorder = 4, alpha = 0.7)

        cb = plt.colorbar(label="Flux 1-99% [ADU]") 
        ax = cb.ax 
        text = ax.yaxis.label
        cb.formatter.set_powerlimits((0, 0))
        font = matplotlib.font_manager.FontProperties(family='DejaVu Sans', style='italic', size=25)
        text.set_font_properties(font) 

        plt.scatter([], [], edgecolor='white', facecolor='none', label='TESS \naperture', s = 200, \
                    marker = 'o')

        # Legend
        maglim = max(gmag-gmag[0]) 
        add = 0
        if int(maglim) % 2 != 0:
            add = 1
        maxmag = int(maglim) + add 
        legend_mags = np.linspace(0,maxmag,int((maxmag+2)/3+1))
        fake_sizes = gmag[0] + legend_mags
        for f in fake_sizes:
            size = 1024.0 / 2**((f-gmag[0])/2)
            plt.scatter(0,0,s=size,c='red',alpha=0.8, edgecolor=None,zorder = 10,\
                        label = r'$\Delta m=$ '+str(int(f-gmag[0]))) 

        legend = plt.legend(fancybox=True, framealpha=0.5, loc='best', facecolor = 'black',\
                            numpoints = 1, prop = {'size': 18})
        plt.setp(legend.get_texts(), color='w')
        try:
            plt.savefig(f'{OUTPUT_DIR}/aperture_tess.pdf', format = 'pdf', dpi=fig.dpi, bbox_inches = "tight")
        except ValueError:
            pass

    # The module-level PRF helpers above are the single source of truth.
    # This legacy duplicate block has been removed to keep one implementation.

    # Pull out the location of the middle of the CCD from the Physics WCS in the header.
    if tpf:
        image_head = tpf.hdu[1].header 
        prime_head = tpf.hdu[0].header
    else:
        n = 0 
        image_head = hdulist[n][1].header
        prime_head = hdulist[n][0].header
        #ap_head = hdulist[n][2].header

    col_center = image_head['1CRV4P'] #col
    row_center = image_head['2CRV4P'] #row

    # Define the location for which we want to retrieve the PRF.
    prf_col_offset = CONFIG["prf_col_offset"]

    # This is the directory where MAST stores the prf FITS files.
    path = CONFIG["prf_fits_base_url"]


    ##########################################################################################
    ######################################### k value  #######################################
    ##########################################################################################

    #Squared TESS aperture (from TPF)
    if tpf:
        center = [round(pix_coords[0][0,0]),round(pix_coords[0][0,1])]

        #Center offset between Pixel Response Function (PRF) and Target Pixel File (TPF)
        off_x = center[0] - int(s[1]/2.) #col
        off_y = center[1] - int(s[1]/2.) #row

        #Size difference between TPF (or TessCut) and PRF
        if s[1] % 2 == 0:
            incr = 1 # Even 
        else:
            incr = 0 # Odd
        size_diff = int((13 - s[1] - incr)/2)
        
        #aperture
        ap = []
        for i in range(s[1]):
            for j in range(s[2]):
                if aperture_mask[i, j]: #matrix[row,col]
                    ap.append([j,i]) #very important! CCD col = X, RA; CCD row = Y, Dec   
        ap = np.array(ap) 
        ap = np.array([ap])

        #Mine function, useful to evaluate the k parameter from PRF.
        #k parameter
        
        @lru_cache(maxsize=None)
        def prf_cached(col, row, ccd, camera, sector, path):
            return get_prf_at_col_row_fits(col, row, ccd, camera, sector, path)  

        rel_flux_per_aperture = np.zeros([len(ap),len(rra)])
        err_flux_per_aperture = np.zeros([len(ap),len(rra)])
        err_flux_per_aperture_try = np.zeros([len(ap),len(rra)])
        for kk in range(len(ap)):
            sep_ctr = []
            with alive_bar(len(rra)) as bar:  # alive bar (terminal)
                for i in range(len(rra)):
                    # location of star in pixel space for aperture k
                    mu_x = pix_coords[kk][i, 0] #col (star)
                    mu_y = pix_coords[kk][i, 1] #row (star)
                    # star's flux
                    A = 10**(-Tmag[i]/2.5)
                    A_err = A*((np.log(10)/2.5)*Tmag_err[i])
                    A_err_prov = 10**(-ufloat(Tmag[i], Tmag_err[i])/2.5)
                    #print(A_err - A_err_prov.s)
                    # evaluate flux (from PRF) in each pixel
                    this_flux = 0
                    this_flux_err = 0
                    this_flux_err_try = []
                    this_flux_err_try2 = 0
                    skip = False #to skip unreliable PRF results (pixel outside the 13x13 PRF size)
                    for j in range(len(ap[kk])):
                        this_pixel = ap[kk][j]+size_diff #different size between tpf and prf
                
                        """
                        Evaluation of the total flux falling into the TESS photometric aperture using the TESS 
                        Pixel Response Function. To consider a given star's total flux, we need to 
                        evaluate the flux fraction that falls within the central pixel (star location) 
                        and also the flux fraction that falls within the neighborhood pixels. 
                        Pixel A is the central one, while B, C, and D are the others (clockwise).
                        """
                        #step to find the weight (contribution) of a given pixel
                        posx = round(mu_x) - mu_x #col
                        posy = round(mu_y) - mu_y #row
                        weightx = 0; weighty = 0; weightx2 = 0; weighty2 = 0
                        weightA = 0; weightB = 0; weightC = 0; weightD = 0
                        #col weight
                        if posx > 0:
                            weightx = 1 - posx 
                            weightx2 = posx
                            col_off = 1  #to change the selected pixel (col)
                        else: 
                            weightx = 1 + posx 
                            weightx2 = -posx
                            col_off = -1 #to change the selected pixel (col)
                        #row weight
                        if posy > 0:
                            weighty = 1 - posy 
                            weighty2 = posy
                            row_off = 1  #to change the selected pixel (row)
                        else:
                            weighty = 1 + posy 
                            weighty2 = -posy
                            row_off = -1 #to change the selected pixel (row)
                        #weights (clockwise)
                        weightA = weightx*weighty   #upper right
                        weightB = weightx*weighty2  #lower right
                        weightC = weightx2*weighty2 #lower left
                        weightD = weightx2*weighty  #upper left
                
                        w = [weightA, weightB, weightC, weightD] #weights array
                        pixel_row = [0, row_off, row_off, 0] #row change, array
                        pixel_col = [0, 0, col_off, col_off] #col change, array
                
                        try:
                            this_flux_pix = 0
                            this_flux_pix_err = 0  
                            row = 0; col = 0
                            
                            this_flux_pix_err_try = []
                            this_flux_pix_err_try2 = 0
                            
                            #loop to evaluate the total flux, sum of the four pixel considered
                            for jj in range(4):
                                #+0.5 is necessary because we are using the int() approx.
                                row = int(this_pixel[1]+0.5+pixel_row[jj]-(mu_y-center[1])-off_y)
                                col = int(this_pixel[0]+0.5+pixel_col[jj]-(mu_x-center[0])-off_x)
                                #to skip negative index
                                if row < 0 or col < 0:
                                    skip = True

                                prf = prf_cached(col_center + mu_x + prf_col_offset,
                                                row_center + mu_y, ccd, camera, sector, path)
                
                                #weighted flux (A*weight), multiplied by the PRF. Array with [row][col].
                                prf_val = prf[row,col]
                                this_flux_pix += w[jj]*A*prf_val
                                this_flux_pix_err += (w[jj]*A_err*prf_val)**2
                                                            
                                this_flux_pix_err_try.append(w[jj]*A_err_prov*prf_val) #ufloat
                                                            
                            this_flux_pix_err_try2 = sum(this_flux_pix_err_try)
                                                                                       
                        except IndexError as e:
                            skip = True #to skip external index
                
                        if skip:
                            this_flux += this_flux_pix
                            this_flux_err_try.append(ufloat(0,0))
                            #The PRF provides a flux lower than probably expected in this case...
                        else:
                            this_flux += this_flux_pix
                            this_flux_err += np.sqrt(this_flux_pix_err)
                            this_flux_err_try.append(this_flux_pix_err_try2)                        
                        
                        this_flux_err_try2 = sum(this_flux_err_try)
                        #print(this_flux_err - this_flux_err_try2.s)
                            
                    #separation (arcseconds) from the center of the aperture.
                    sep_ctr.append(np.sqrt((mu_x-center[0]+off_x)**2+(mu_y-center[1]+off_y)**2))
                    rel_flux_per_aperture[kk, i] = this_flux
                    err_flux_per_aperture[kk, i] = this_flux_err
                    err_flux_per_aperture_try[kk, i] = this_flux_err_try2.s
                    
                    bar() 
                    
        #TESS magnitude    
        Tmag_arr = np.array(Tmag) 
        Tmag_err_arr = np.array(Tmag_err) 
        k = rel_flux_per_aperture/(10**(Tmag_arr/(-2.5))) #k parameter

        k_err = 0 #TODO: add the PRF's uncertainty (?)
        
        cont_ratio_prim = np.sum(rel_flux_per_aperture[0,1:])/rel_flux_per_aperture[0][0]
        cont_ratio_sec = (rel_flux_per_aperture[0][0] + np.sum(rel_flux_per_aperture[0,2:]))/rel_flux_per_aperture[0][1]
        dilut_factor = np.sum(rel_flux_per_aperture[0,1:])/np.sum(rel_flux_per_aperture[0,:])
        dilut_factor2 = (rel_flux_per_aperture[0][0] + np.sum(rel_flux_per_aperture[0,2:]))/np.sum(rel_flux_per_aperture[0,:])

        try:
            os.remove(f'{OUTPUT_DIR}/contam.txt')
        except OSError:
            pass

        with open(f'{OUTPUT_DIR}/contam.txt', 'a') as f:
            f.write('dilution_factor = '+str(dilut_factor)+'\n')
            f.write('dilution_factor2 = '+str(dilut_factor2)+'\n')
            f.write('cont ratio secondary = '+str(cont_ratio_sec)+'\n')
            f.write('cont ratio primary = '+str(cont_ratio_prim)+'\n')
            f.close
        
        cont_ratio_prim_err = cont_ratio_prim*np.sqrt((np.sqrt(np.sum(err_flux_per_aperture[0,1:]**2))/np.sum(rel_flux_per_aperture[0,1:]))**2 + 
                            (err_flux_per_aperture[0,0]/rel_flux_per_aperture[0][0])**2)
            
        numero = []                
        for i in range(len(rel_flux_per_aperture[0])-1):
            numero.append(ufloat(rel_flux_per_aperture[0,i+1],err_flux_per_aperture_try[0,i+1]))
                            
        # print(sum(numero))
        #cont_ratio_prim_err_try = sum(numero)/ufloat(rel_flux_per_aperture[0][0],err_flux_per_aperture_try[0][0])
        #print(cont_ratio_prim_err_try.n, cont_ratio_prim_err_try.s)
        #print('cont ratio primary err = ', cont_ratio_prim_err)


        #relative flux
        rel_flux = (k*(10**(-Tmag_arr/2.5)))/(np.sum((10**(-Tmag_arr/2.5))*k))

        #uncertainty in the relative flux
        err_0 = (((np.log(10)/2.5)*Tmag_err_arr))**2 #abs err (instrumental) #+(k_err/k)**2 
        err_1 = (np.sum(err_0*((k*(10**(-Tmag_arr/2.5)))**2)))/((np.sum(k*(10**(-Tmag_arr/2.5))))**2)
        err_2 = np.sqrt(err_0+err_1)
        err_3 = (2.5/np.log(10))*(err_2)

        err_flux = err_2*rel_flux #useful later

        #flux ratio (mag)
        fl = -2.5*np.log10(1-depth/rel_flux)
    #####################################################################################
    #####################################################################################

    #Circular TESS aperture
    if not tpf:
        len_prf = CONFIG["prf_image_size"]
        center = [int(len_prf/2), int(len_prf/2)]
        
        #k parameter (with circular aperture)
        rel_flux_per_aperture = np.zeros([len(aperture),len(rra)])
        for kk in range(len(aperture)):
            sep_ctr = []
            with alive_bar(len(rra)) as bar:  # alive bar (terminal)
                for i in range(len(rra)):
                    # location of star in pixel space for aperture k
                    mu_x = pix_coords[kk][i, 0] #col (star)
                    mu_y = pix_coords[kk][i, 1] #row (star)
                    # star's flux
                    A = 10**(-Tmag[i]/2.5)
                    
                    # evaluate flux (from PRF) in each pixel
                    this_pixel = [center[0]-(mu_x-pix_coords[kk][0,0]), center[1]-(mu_y-pix_coords[kk][0,1])] #col,row
                
                    """
                        Evaluation of the total flux falling into the TESS photometric aperture using the TESS 
                        Pixel Response Function. 
                    """
            
                    #PRF multiplied by the total flux of the star selected
                    this_flux_pix = A*get_prf_at_col_row_fits(col_center + mu_x +prf_col_offset, \
                                                        row_center + mu_y, ccd, camera, sector, path)
                                                    
                    if source == 'qlp':
                        aperture = CircularAperture(this_pixel, r=int(rad_ap))
                    else:
                        """ custom circular aperture, in pixel radius """
                        #input_aperture = input("Give me the aperture in pixel: ")
                        if radius_input is None:
                            radius_input = int(input("Give me custom aperture radius: "))
                        aperture = CircularAperture(this_pixel, r=int(radius_input)) #custom aperture

                    #method "exact", to consider also fractions of pixels inside the circular aperture
                    phot_table = aperture_photometry(this_flux_pix, aperture, method = 'exact')
                    #print(phot_table)
            
                    if not np.isnan(phot_table['aperture_sum']):
                        rel_flux_per_aperture[kk, i] = phot_table['aperture_sum']
                    else:
                        rel_flux_per_aperture[kk,i] = 0 #to avoid "nan" value
                            
                    #separation (arcseconds) from the center of the aperture.
                    sep_ctr.append(np.sqrt((mu_x-positions[0][0])**2+(mu_y-positions[0][1])**2))
                    
                    bar() 
                    
        #TESS magnitude    
        Tmag_arr = np.array(Tmag) 
        Tmag_err_arr = np.array(Tmag_err) 
        k = rel_flux_per_aperture/(10**(Tmag_arr/(-2.5))) #k parameter
        
        cont_ratio_prim = np.sum(rel_flux_per_aperture[0,1:])/rel_flux_per_aperture[0][0]
        cont_ratio_sec = (rel_flux_per_aperture[0][0] + np.sum(rel_flux_per_aperture[0,2:]))/rel_flux_per_aperture[0][1]
        dilut_factor = np.sum(rel_flux_per_aperture[0,1:])/np.sum(rel_flux_per_aperture[0,:])
        dilut_factor2 = (rel_flux_per_aperture[0][0] + np.sum(rel_flux_per_aperture[0,2:]))/np.sum(rel_flux_per_aperture[0,:])

        try:
            os.remove(f'{OUTPUT_DIR}/contam.txt')
        except OSError:
            pass
        
        with open(f'{OUTPUT_DIR}/contam.txt', 'a') as f:
            f.write('dilution_factor = '+str(dilut_factor)+'\n')
            f.write('dilution_factor2 = '+str(dilut_factor2)+'\n')
            f.write('cont ratio secondary = '+str(cont_ratio_sec)+'\n')
            f.write('cont ratio primary = '+str(cont_ratio_prim)+'\n')
            f.close
        
        #relative flux
        rel_flux = (k*(10**(-Tmag_arr/2.5)))/(np.sum((10**(-Tmag_arr/2.5))*k))

        #uncertainty in the relative flux
        err_0 = (((np.log(10)/2.5)*Tmag_err_arr))**2 #abs err (instrumental) #+(k_err/k)**2 
        err_1 = (np.sum(err_0*((k*(10**(-Tmag_arr/2.5)))**2)))/((np.sum(k*(10**(-Tmag_arr/2.5))))**2)
        err_2 = np.sqrt(err_0+err_1)
        err_3 = (2.5/np.log(10))*(err_2)

        err_flux = err_2*rel_flux #useful later

        #flux ratio (mag)
        fl = -2.5*np.log10(1-depth/rel_flux)
        
        
        #Figure: PRF
        with plt.style.context('default'):
            for i in range(2):
                plt.rcParams.update({'font.size': 20})
                pr = get_prf_at_col_row_fits(col_center + pix_coords[0][i, 0] +prf_col_offset, \
                                            row_center + pix_coords[0][i, 1], ccd, camera,\
                                            sector, path)

                fig, ax = plt.subplots(ncols=2, sharey=True, tight_layout=True, figsize=(15, 9))
                ax = plt.subplot(111)
                plt.imshow(pr, cmap = plt.cm.viridis)

                this_pixel = [center[0]-(pix_coords[0][i,0]-pix_coords[0][0,0]), \
                              center[1]-(pix_coords[0][i,1]-pix_coords[0][0,1])]
                if source == 'qlp':
                    aperture = CircularAperture(this_pixel, r=int(rad_ap))
                else:
                    """ circular aperture, custom pixel radius """
                    if radius_input is None:
                        radius_input = int(input("Give me custom aperture radius: "))
                    aperture = CircularAperture(this_pixel, r=int(radius_input)) #custom value
                aperture.plot(color='white', lw=2, label = 'TESS aperture') 
                
                cb = plt.colorbar(label="Flux %")
                plt.legend(loc = 'upper left', framealpha=0.7, scatterpoints=1, labelspacing=1, fontsize = 20)
                plt.xlabel('PRF col') 
                plt.ylabel('PRF row')
                #plt.savefig(f'{OUTPUT_DIR}/prf{str(i+1)}.pdf', format = 'pdf', dpi=fig.dpi, bbox_inches = "tight")
                
    #####################################################################################
    #####################################################################################


    def final_figure(sel, gmag, fl, err_3, rra, dist, cb):
        verify = False
        fig, ax = plt.subplots(1, 1, sharey=True, figsize = (15,9))
        cont_fig = 0
        for i in range(len(rra)):
            norm = plt.Normalize(0, 210)
            plt.axhline(y=0, color = "grey")
        
        if sel == 0:
            plt.scatter(gmag[sel], fl[0][sel], zorder = 4, color = 'red', alpha = 0.5, s = 200)
            plt.errorbar(gmag[sel], fl[0][sel], yerr=err_3[sel], fmt='', ecolor='orange',\
                         elinewidth=3, capsize=10, zorder = 2, color = 'red', alpha = 0.5, markersize = 10)
        else:
            for i in range(1, len(rra)):
                if not np.isnan(fl[0][i]):
                    plt.errorbar(gmag[i], fl[0][i], yerr=err_3[i], fmt='o', ecolor='orange',\
                                  elinewidth=2, capsize=0, zorder = 4, color = 'red', alpha = 0.5, markersize = 8)
                    cont_fig += 1
                    ann = ax.annotate(cont_fig, (gmag[i], fl[0][i]), annotation_clip=False)
                    ann.set_clip_box(ax.bbox)

        ax = cb.ax; text = ax.yaxis.label
        font = matplotlib.font_manager.FontProperties(family='DejaVu Sans', style='italic', size=20)
        text.set_font_properties(font)

        plt.scatter([], [], c='red', alpha=0.6, label='Expected On-Off', s = 35)
        plt.legend(loc = 'best', scatterpoints=1, labelspacing=1, fontsize = 18)
        
        plt.xlabel('Gaia magnitude', fontsize = 20)
        plt.ylabel('Magnitude variation (On-Off)', fontsize = 20)
        
        if sel == 0:
            plt.savefig(f'{OUTPUT_DIR}/onoff_target.pdf', format = 'pdf', dpi=fig.dpi, bbox_inches = "tight")
        else:
            plt.savefig(f'{OUTPUT_DIR}/onoff_comparison.pdf', format = 'pdf', dpi=fig.dpi, bbox_inches = "tight")

    def save_result_tables(gmag, rra, ddec, dist, sep_ctr, k, fl, err_3, df_gaia):
        onoff_teo = []; err_onoff_teo = []; Gaia_id = []; sep_c = []; k_val = []
        for i in range(len(rra)):
            sep_c.append(sep_ctr[i]*21)
            k_val.append(k[0][i])
            onoff_teo.append(fl[0][i])
            err_onoff_teo.append(err_3[i])
            Gaia_id.append('Gaia DR3 '+str(df_gaia['source_id'][i]))

        df_fin = pd.DataFrame({"GaiaID":Gaia_id, "on-off_teo": onoff_teo, "e_on-off_teo": err_onoff_teo, \
                               "Gaia_mag": gmag, "sep_target": dist, "sep_aperture's_center": sep_c, \
                               "k_parameter": k_val, "Gaia_RA":rra, "Gaia_Dec":ddec})
        df_fin2 = df_fin.to_string(formatters={
            'on-off_teo': '{:,.5f}'.format,
            'err_on-off_teo': '{:,.5f}'.format,
            'sep_target': '{:,.2f}'.format,
            "sep_aperture's_center": '{:,.2f}'.format,
            "k_parameter": '{:,.7f}'.format
        })
        df_fin = pd.read_csv(io.StringIO(df_fin2), sep = "\s+")
        df_fin.to_csv(f'{OUTPUT_DIR}/On-Off_Photometry.csv', index = False, header=True, sep = ',')

    def save_contaminant_tables(depth, rel_flux, err_flux, gmag, dist, k, fl, sep_ctr, df_gaia, rra):
        depth2 = depth/rel_flux
        e_depth2 = depth2*(err_flux/rel_flux)

        # print("depth =",depth, ", undiluted depth =",depth2[0][0])

        first = True
        dep = []; depmag = []; dis = []; kpar = []; gg = []; gaga = []; dist_min = 0; sep_ctr2 = []
        for i in range(len(rra)):
            if depth2[0][i] < 1:
                dep.append(depth2[0][i])
                depmag.append(-2.5*np.log10(1-depth2[0][i]))
                dis.append(dist[i])
                kpar.append(k[0][i])
                gg.append(gmag[i])
                gaga.append('Gaia DR3 '+str(df_gaia['source_id'][i]))
                sep_ctr2.append(sep_ctr[i]*21)
                # print('depth (flux) = ',depth2[0][i],', depth (mag) = ',-2.5*np.log10(1-depth2[0][i]),\
                #       ', sep target (arcsec) = ', dist[i],', k = ', k[0][i],', g mag = ', gmag[i])

                if first:
                    if i != 0:
                        # print('distance of the first possible contaminant analyzed = ', dist[i])
                        dist_min = dist[i]
                        first = False

        df2 = pd.DataFrame({"Gaia_ID": gaga, "depth_(flux)":dep, "depth_(mag)":depmag, "sep_target":dis,"sep_center": sep_ctr2, \
                            "k_value": kpar, "Gaia_mag": gg})
        df22 = df2.to_string(formatters={
            'depth_(flux)': '{:,.5f}'.format,
            'depth_(mag)': '{:,.5f}'.format,
            'sep_target': '{:,.2f}'.format,
            "sep_center": '{:,.2f}'.format
        })
        df2 = pd.read_csv(io.StringIO(df22), sep = "\s+")
        df2.to_csv(f'{OUTPUT_DIR}/Useful_results.csv', index = False, header = True, sep = ',')

        try:
            os.remove(f'{OUTPUT_DIR}/First_contam.txt')
        except OSError:
            pass

        try:
            os.remove(f'{OUTPUT_DIR}/Last_contam.txt')
        except OSError:
            pass

        with open(f'{OUTPUT_DIR}/First_contam.txt', 'a') as f:
            f.write(str(dist_min)+'(arcseconds)')
            f.close

        with open(f'{OUTPUT_DIR}/Last_contam.txt', 'a') as f:
            f.write(str(dis[-1])+'(arcseconds)')
            f.close

        return depth2, e_depth2, dist_min

    def compute_contamination_summary():
        """Orchestrate the final plotting and summary export for the contamination analysis."""
        final_figure(0, gmag, fl, err_3, rra, dist, cb)
        final_figure(1, gmag, fl, err_3, rra, dist, cb)
        save_result_tables(gmag, rra, ddec, dist, sep_ctr, k, fl, err_3, df_gaia)
        depth2, e_depth2, dist_min = save_contaminant_tables(depth, rel_flux, err_flux, gmag, dist, k, fl, sep_ctr, df_gaia, rra)
        return depth2, e_depth2, dist_min


    log_stage("Computing contamination summary", "writing the contamination tables and final figures")
    #################### contaminant limit.. when depth (Delta F/F) -> 1
    depth2, e_depth2, dist_min = compute_contamination_summary()

    print("\n" + "=" * 80)
    print("Analysis complete")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 80)

    ################################################################################################################
    ############################ planetary radius evaluation (after dilution estimation) ###########################
    ################################################################################################################

    RSUN = const.R_sun.cgs.value 
    REARTH = const.R_earth.cgs.value
    factor = RSUN/REARTH

    s_rad = df_tess['rad'][0] #tic rad
    e_s_rad = df_tess['e_rad'][0] #tic err rad

    #new planetary radius (undiluted, teoric)
    rp = (np.sqrt(depth2[0][0])*s_rad)*(RSUN/REARTH)
    e_rp = rp*(np.sqrt((e_s_rad/s_rad)**2+((0.5*e_depth2[0][0])/(depth2[0][0]**(-3/2)))**2))

    rp_pdcsap = (np.sqrt(depth_pdcsap)*s_rad)*(RSUN/REARTH)
    e_rp_pdcsap = rp*(np.sqrt((e_s_rad/s_rad)**2+((0.5*e_depth_pdcsap)/(depth_pdcsap**(-3/2)))**2))

    #diluted planetary radius (obs)
    rp_obs = (np.sqrt(depth)*s_rad)*(RSUN/REARTH)
    e_rp_obs = rp*(np.sqrt((e_s_rad/s_rad)**2+((0.5*e_depth)/(depth**(-3/2)))**2))


    try:
        os.remove(f'{OUTPUT_DIR}/New_planetary_radius.txt')
    except OSError:
        pass

    try:
        with open(f'{OUTPUT_DIR}/New_planetary_radius.txt', 'a') as f:
            f.write(str(round(rp,3)) +' '+ str(round(e_rp,3))+' (Earth radii) undiluted, custom-estimated\n')
            f.write(str(round(rp_pdcsap,3)) +' '+ str(round(e_rp_pdcsap,3))+' (Earth radii) undiluted, pdcsap\n')
            f.write(str(round(rp_obs,3)) +' '+ str(round(e_rp_obs,3))+' (Earth radii) diluted, sap')
            f.close
    except NameError:
        pass

    ################################################################################################################
    ########################## Tess magnitude for a possible contaminating EB (unresolved) #########################
    ################################################################################################################

    #data
    lum = df_tess['lum'][0] #target luminosity (L_sun)
    lum_eb_log = [-3.57,-2.52,-2.03,-1.64,-1.39,-1.1,-0.76,-0.34,-0.17]; """Pecaut & Mamajek (2013, ApJS, 208, 9)"""
    lum_eb = 10**(np.array(lum_eb_log))
    lum_ratio = lum_eb/lum

    #arbitrary EB radii (solar units)
    r_st = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] 
    r_st = np.array(r_st)
    ratio = r_st/s_rad #radii ratio

    #fl_mag = -2.5*np.log10(1-np.array(r_st)**2) #depth mag
    fl_mag = (-2.5*np.log10(1-((ratio**2)*(1-(lum_ratio*(ratio**(-2))))))) #depth mag
    #added contribution from eb luminosity


    for i in range(2):
        if i == 1:
            fl_mag = r_st
            #r_st = np.sqrt(1-10**(fl_mag/(-2.5)))
            r_st = np.sqrt(((1-10**(fl_mag/(-2.5)))+lum_ratio)*(s_rad**2))
        summa = np.sum((10**(-Tmag_arr/2.5))*k) - k[0][0]*(10**(-Tmag_arr[0]/2.5))

        rel = depth/(1-10**(-0.4*fl_mag))
        numerat = summa*rel
        denom = k[0][0]*(1-rel)
        Teb = -2.5*np.log10(numerat/denom) #tess magnitude for the EB
        # print('Tess magnitude EB: ',Teb)

        s_rad_array = []; earth_radii = []
        for j in range(len(r_st)):
            s_rad_array.append(s_rad)
            #earth_radii.append(r_st[j]*s_rad*factor)
            earth_radii.append(r_st[j]*factor)

        #df_eb = pd.DataFrame({"eb TESS magnitude":Teb,"depth [mag]":fl_mag, "eb/target ratio":r_st, "Target radius [sun]":s_rad_array, "eb radius [earth radii]":earth_radii})
        df_eb = pd.DataFrame({"eb TESS magnitude":Teb,"depth [mag]":fl_mag, "eb radius":r_st, "Target radius [sun]":s_rad_array, "eb radius [earth radii]":earth_radii})
        if i == 1:
            df_eb.to_csv(f'{OUTPUT_DIR}/eb_prevision2b.csv', index = False, header = True, sep = ',')
        else:
            df_eb.to_csv(f'{OUTPUT_DIR}/eb_previsionb.csv', index = False, header = True, sep = ',')
            
            
"""
Info: not used for now, but could be useful in the future to perform a light-curve inspection of the target and contaminant stars.
            
    ################################################################################################################
    ############################### TARGET & CONTAMINANT STARS Light-curve inspection ##############################
    ################################################################################################################

    def LC(ID, flux_tot, time, aperture, duration, period, T0):
        
        plt.rcParams.update({'font.size': 20})
        
        phot_table = []
        for i in range(len(flux_tot)):
            phot_table.append(aperture_photometry(flux_tot[i], aperture, method = 'exact')['aperture_sum'][0])       

        #Remove outliers +-5sigma separation
        clipped_flux = slide_clip(time, phot_table, window_length=5*duration,
                                  low=5, high=5, 
                                  method='mad',  # mad or std
                                  center='median'  # median or mean
                                  )

        #fig = plt.figure(figsize = (15,9))
        #plt.scatter(time, clipped_flux, s = 10, color='black')
        
        #Mask the transit signal
        mask = transit_mask(time=time, period=period, duration=duration, T0=T0)
            
        #Flatten the light-curve
        flatten_lc, trend_lc = flatten(time, clipped_flux, window_length=5*duration, \
                                       edge_cutoff=0.5,      # length (in units of time) to be cut off each edge.
                                       break_tolerance=0.5,  # Split into segments at breaks longer than that
                                       return_trend=True,    # Return trend and flattened light curve
                                       robust = True, cval=5.0, method='biweight', mask = mask)
                              
        #fig = plt.figure(figsize = (15,9))
        ##plt.scatter(time, flatten_lc, s = 10, color='black')
        #plt.scatter(time[~mask], flatten_lc[~mask], s = 10, color='black')
        #plt.scatter(time[mask], flatten_lc[mask], s = 10, color='orange')
        
        #Phase-folded light-curve
        time_folded = (time - T0 + period/2) % period - period/2
        idx = np.where((time_folded > -1.5*duration) & (time_folded < 1.5*duration))
        
        plt.figure(figsize=(14, 6))
        plt.scatter(time_folded[idx[:]], flatten_lc[idx[:]])
        plt.xlabel("Phase")
        plt.ylabel("Normalized flux")
        plt.savefig('on-off/tess_lc_GaiaDR3'+str(ID)+'.pdf', format = 'pdf', dpi=fig.dpi, bbox_inches = "tight")
        
    def star_LC(ID, RA, Dec, sector, source, duration, period, T0):
        
        #star coordiantes
        coord = SkyCoord(ra=RA, dec=Dec, unit=(u.degree, u.degree), frame='icrs')
        
        #Creation of a pixel file, using Tesscut
        size = 3; #TODO: decide what area I want to show/use
        hdulist = Tesscut.get_cutouts(coordinates = coord, size = size, sector = sector)
        hdu1 = hdulist[0] #hdu
        
        #wcs coordinate system
        wcs = WCS(hdu1[2].header,relax=True) #wcs
        
        ##bad pixels removal
        if source == 'qlp':
            bad_bits = np.array([1,3,4,5,6,8,11,12,13]) #QLP
        else:
            bad_bits = np.array([1,2,3,4,5,6,8,10,12]) #PDCSAP
        
        quality = hdu1[1].data['QUALITY']

        value = 0
        for v in bad_bits:
            value = value + 2**(v-1)
        bad_data = np.bitwise_and(quality, value) >= 1 
        ##
        
        #flux and time without bad pixels
        flux_tot = hdu1[1].data['FLUX'][~bad_data] #flux
        time = hdu1[1].data['TIME'][~bad_data] #time
        first_image = hdu1[1].data['FLUX'][0] #flux
        
        radius_input = 0.5; #TODO: decide what radius I want to use, to perform aperture photometry
        #radius_input = input("Give me custom aperture radius: ")
        
        #Figure: TESS-cutout aperture around selected star
        with plt.style.context('default'):
            plt.rcParams.update({'font.size': 20})
            fig, ax = plt.subplots(ncols=2, sharey=True, figsize = (15,9))

            ap_color = 'black'

            #Coordinates transformation from deg to pixels
            RApix = wcs.all_world2pix(RA, Dec, 0)[0].item()
            Decpix = wcs.all_world2pix(RA, Dec, 0)[1].item()
            positions = [(RApix, Decpix)]
        
            #circular aperture
            aperture = CircularAperture(positions, r=radius_input)

            #flux range
            lower = np.percentile(first_image,1) 
            upper = np.percentile(first_image,99)

            #axis
            ax = plt.subplot(111, projection = wcs) 
            
            #labels and grid
            plt.grid(axis = 'both', color = 'white', ls = 'solid')

            #star selected
            ax.scatter(RA, Dec, transform=ax.get_transform('fk5'), s=500, zorder = 4, edgecolor='black',\
                       facecolor='yellow', linewidth =0.5, label = 'Chosen star', marker = "*")
            ax.scatter(RA, Dec, transform=ax.get_transform('fk5'), s=1000, zorder = 3, edgecolor='black',\
                       facecolor='red', linewidth =1, alpha = 0.7)

            #photometric aperture chosen
            aperture.plot(color=ap_color, lw=2, transform=ax.get_transform(wcs))
            im = plt.imshow(first_image, origin = 'lower', cmap = plt.cm.viridis,clim=(lower,upper))   
            
            ax = plt.gca()

            # We need to set the ticklabel position explicitly
            ax.coords['ra'].set_ticklabel_position('l')
            ax.coords['ra'].set_axislabel_position('l')
            ax.coords['ra'].set_axislabel('Right Ascension')
            ax.coords['dec'].set_ticklabel_position('b')
            ax.coords['dec'].set_axislabel_position('b')
            ax.coords['dec'].set_axislabel('Declination')     

            #colorbar
            cb = plt.colorbar(label="Flux 1-99% [ADU]") 
            ax = cb.ax 
            text = ax.yaxis.label
            cb.formatter.set_powerlimits((0, 0))
            font = matplotlib.font_manager.FontProperties(family='DejaVu Sans', style='italic', size=25)
            text.set_font_properties(font) 
        
            #legend
            legend = plt.legend(fancybox=True, framealpha=0.8, loc='best',\
                                numpoints = 1, prop = {'size': 18})
        
            plt.savefig('on-off/tess_cutout_GaiaDR3'+str(ID)+'.pdf', format = 'pdf', dpi=fig.dpi, bbox_inches = "tight")
        
            LC(ID, flux_tot, time, aperture, duration, period, T0)


    from wotan import slide_clip
    from wotan import transit_mask, flatten

    for i in range(len(rra)):
        if not np.isnan(fl[0][i]):
            ID = df_gaia['source_id'][i]
            star_LC(ID, rra[i], ddec[i], sector_choice, source, duration, period, T0)
"""  


if __name__ == "__main__":
    try:
        main()
    except ValueError as exc:
        print(f"Error: {exc}")
        raise SystemExit(1)
