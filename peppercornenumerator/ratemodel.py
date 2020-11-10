#
#  peppercornenumerator/ratemodel.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

import math

def polymer_link_length(before, after):
    """ Length estimates in nucleotides for linkers connecting two domains.

    Args:
        before (Loop): A loop object for which we measure the length.
        after (Loop): A loop object for which we measure the length.

    Returns:
        [float]: The length of the shorter loop.
    """
    L_stem = 2.0 / 0.43 # rough equivalent number of single-stranded nucleotides to span a stem
    L_before = float('inf') if before.is_open else \
            1 + before.bases + before.stems + L_stem * before.stems 
    L_after = float('inf') if after.is_open else \
            1 + after.bases + after.stems + L_stem * after.stems
    if after.is_open and before.is_open:
        raise ValueError("Computing polymer lengths for disconnected complexes!")
    return min(L_before, L_after)

def polymer_link_rate(hllen, ha = 1e6, hb = -2.5, kmax = 33_000):
    """ Unimolecular hairpin closing rate, as a function of hairpin loop length. 

    References for default values: 
        * Bonnet et al. (1998)
        * Kuznetsov et al. (2008)
        * Nayak et al. (2012)
        * Tsukanov et al. (2013)

    Args:
        hllen (flt): Hairpin loop length. The number of unpaired bases in a
            hairpin loop, or an arbitrary linker length in terms of an
            equivalently long chain of unpaired bases.
        ha (flt, optional): The rate for closing the first base-pair. 
            Defaults to 1_000_000 [/s].
        hb (flt, optional): Exponential factor to relate hllen to the
            probability of forming the initial contact. Defaults to -2.5.
        kmax (flt, optional): Fastest hairpin closing rate. Defaults to 33_000 [/s]. 
            (E.g. hairpins shorter than 4 nt cannot close any faster than this.)

    Returns:
        [flt]: The rate for closing the initial base-pair of a hairpin.
    """
    return min(ha * hllen**hb, kmax)

def unimolecular_binding_rate(length, before, after):
    """
    Rate constant formula for unimolecular binding of a domain of the given length.
    Could be zippering, hairpin closing, bubble closing, bulge formation, multiloop formation,
    depending on the flanking loops, which may be open or closed.
    """
    if not before.is_open and before.stems == 1 and before.bases == 0:
        if not after.is_open and after.stems == 1 and after.bases == 0:
            return 1e4  # bubble closing rate from Altan-Bonnet 2003
        # zippering from Wetmur&Davidson 1968, Gueron&Leroy 1995, Srinivas et al 2013, low end
        return (1e6) / length
    if not after.is_open and after.stems == 1 and after.bases == 0:
        # zippering from Wetmur&Davidson 1968, Gueron&Leroy 1995, Srinivas et al 2013, low end
        return (1e6) / length
    # bulge closing assumed to be similar to faster of two hairpin closings
    L = polymer_link_length(before, after)
    # hairpin closing adapted (simplified) from data in Bonnet et al 1998
    return polymer_link_rate(L)

def bimolecular_binding_rate(length, k_nuc = 3e5):
    """ Rate constant formula for bimolecular association (binding).

    Hybridization is determined by well-aligned initial contacts, therefore 
    dependent linearly on length

    Literature: 
        Wetmur 1976 review,
        Srinivas et al. 2013 AEL model, 
        Zhang & Winfree 2009.

    Args:
        length (int): The length of a domain.
        k_nuc (flt, optional): The hybridization rate per nucleotide. 
            Defaults to 3*10^5.

    Returns:
        [float]: The rate constant for bimolecular binding.
    """
    return length * k_nuc

def opening_rate(length, dG_bp = -1.7, dG_assoc = 1.9, kelvin = 298.15, dissoc = True):
    """ Rate constant for opening a duplex of a given `length`.

    k_open = k_bind * exp( (length * dG_bp + dG_assoc) / RT )

    Args: 
        dG_bp (flt): Energy contribution of a single WC base pair. 
            Defaults to -1.7 kcal/mol.
        dG_assoc (flt): Association energy penalty.
            Defaults to 1.9 kcal/mol.
        kelvin (flt): Temperature in K. Defaults to 298.15 (25 C).

    Returns:
        [float]: The rate constant.
    """
    # R = 0.001987 kcal/mol/K
    # T = (273.15 + 25) K
    RT = 0.001987 * kelvin

    if True or dissoc: # stay with simple model!
        k_bind21 = bimolecular_binding_rate(length)
        dG = (length * dG_bp) + dG_assoc
        rate = k_bind21 * math.exp(dG/RT)
    else:
        k_bind11 = 1e6 * length # using zipping rate, not bubble-closing
        dG = length * dG_bp
        rate = k_bind11 * math.exp(dG/RT)
    return rate

def branch_3way_remote_rate(length, before, after):
    """ Rate constant for a 3-way branch-migration reaction.

    Args: 
        length (int): the length of the domain
        before (:obj:Loop): the loop object on the initiation site.
        after (:obj:Loop): the other (remote) loop object.

    Returns:
        k1 = P(link) * k_init / l
    """
    # step = 0.1e-3  # 100us branch migration step time from Srinivas et al 2013 (not relevant)
    # k_init = k_uni exp(-dGsp / RT) with k_uni = 7.5e7, dGsp = 7.3 kcal/mol,
    # T = 273.15 + 25 K, R = 0.001987 kcal/mol/K

    # Experimentally derived time for a single step of 3-way branch migration.
    t_step = 1e-4  # sec
    k_step = 1/t_step

    if after._parts[0] is None: 
        # Initialtion penalty if there is *no* danling end on the displaced side.
        t_init = 3.0e-3 
        k_init = 1/t_init
    else :
        # Initialtion penalty if there is a danling end on the displaced side,
        # in which case it is just a regular displacement step.
        t_init = t_step
        k_init = k_step

    k1 = None # The probability of a successfull branch migration scales with 1/L
    k2 = None # The time spent on branch migration scales with 1/L^2

    # "standard" 3-way bm initiation (plus "after" being closed)
    if not before.is_open and before.stems == 1 and before.bases == 0:
        k1 = k_init / length
        #k2 = k_init / length + k_step / length**2
    else :
        # consider a slowdown analogous to Genot et al 2011 (remote) supp info derivation
        # bulge closing assumed to be similar to faster of two hairpin closings
        L = polymer_link_length(before, after)

        # how much slower than our (admittedly slow) zippering rate is this?
        ratio = polymer_link_rate(L) / 1e6

        # we slow down initiation and thus success probability (note: ratio < 1/30)
        k1 = ratio * k_init / length
        #k2 = ratio * k_init / length + k_step / length**2

    return k1

def branch_4way_remote_rate(length, before, after):
    """ Rate constant for a 4-way branch-migration reaction.
    """
    # rates recalculated from Nadine Dabby, PhD Thesis, 2013, based on
    # assuming last 6 bp dissociate faster than 4way bm step
    open_step = 107   # sec, = 1/k_first  (this is for open 4-way case only)
    # sec, = 1/k_bm     (this is used for initiating closed 4-way case;
    # consistent with Panyutin&Hsieh 1993)
    closed_step = 3.0

    # open_step = 200 # fudge !
    if not before.is_open and not after.is_open:
        init = closed_step
        if before.bases == 0 and before.stems == 1 and after.bases == 0 and after.stems == 1:
            return 1.0 / init / length   # closed & ready-to-rock-and-roll 4 way initiation
    if before.is_open:
        init = open_step
        if after.bases == 0 and after.stems == 1:
            # we care about probability of taking this path, not actual time
            return 1.0 / init / length
    if after.is_open:
        init = open_step
        if before.bases == 0 and before.stems == 1:
            # we care about probability of taking this path, not actual time
            return 1.0 / init / length

    # consider a slowdown analogous to Genot et al 2011 (remote) supp info derivation
    # bulge closing assumed to be similar to faster of two hairpin closings
    L = polymer_link_length(before, after)
    # how much slower than our (admittedly slow) zippering rate is this?
    ratio = polymer_link_rate(L) / (1e6)
    # we slow down initiation and thus success probability (note: ratio < 1/30)
    return ratio / init / length

