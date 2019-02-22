#
#  output.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/21/10.
#  modifications by Stefan Badelt.
#

from peppercornenumerator import __version__
from peppercornenumerator.utils import natural_sort, convert_units
from peppercornenumerator.utils import PeppercornUsageError
from peppercornenumerator.condense import PepperCondensation

def format_units(cplx, unit = 'M'):
    ini = cplx._concentration[0]
    num = float(cplx._concentration[1])
    old_unit = cplx._concentration[2]
    num = convert_units(num, old_unit, unit)

    return ' '.join([ini, str(num), unit])

# DEPRECATED
def write_kernel(enumerator, fh = None, detailed = True, condensed = False, 
        composite = None, molarity = 'M', time = 's'):
    
    print("# Deprecated function: use peppercornenumerator.output.write_pil instead of write_kernel")
    return write_pil(enumerator, fh, detailed, condensed, composite, molarity, time)

def write_pil(enumerator, fh = None, detailed = True, condensed = False, 
        composite = None, molarity = 'M', time = 's'):
    """Write the contents of :obj:`Enumerator()` into a *.pil file.

    Args:
      enumerator (:obj:`Enumerator()`): The enumertor object.
      fh (filehandle): The filehandle to which the output is written to. 
        If fh is None, then output is written to a string and returned.
      detailed (bool, optional): Print detailed CRN. Defaults to True.
      condensed (bool, optional): Print condensed CRN. Defaults to False.
      composite (dict, optional): A dictionary of composite domains (or
        strands) to be included into the *.pil file.
    """

    out = []
    def output_string(string):
        if fh is None:
            out.append(string)
        else :
            fh.write(string)

    output_string("# File generated by peppercorn-{}\n".format(__version__))

    # Print domains
    seen = set()
    output_string("\n# Domains ({}) \n".format(len(enumerator.domains)))
    for dom in natural_sort(enumerator.domains):
        if dom.is_complement and not dom.nucleotides: 
            dom = ~dom
        if dom not in seen :
            if dom.nucleotides:
                output_string("sequence {:s} = {} : {:d}\n".format(
                    dom.name, dom.nucleotides, len(dom)))
            else :
                output_string("length {:s} = {:d}\n".format(dom.name, len(dom)))
            seen.add(dom)

    if composite:
        output_string("\n# Strands or composite domains ({}) \n".format(len(composite)))
        for comp in composite:
            if comp[-1]=='*' and comp[:-1] in composite:
                continue
            output_string("sup-sequence {} = {} : {:d}\n".format(comp, ' '.join(
                map(str, composite[comp])), sum(map(len,composite[comp]))))

    # Print resting complexes
    output_string("\n# Resting complexes ({}) \n".format(len(
        enumerator.resting_complexes)))
    for cplx in natural_sort(enumerator.resting_complexes):
        if cplx._concentration :
            output_string("{:s} = {:s} @{}\n".format(cplx.name, 
                cplx.kernel_string, format_units(cplx, unit = molarity)))
        else:
            output_string("{:s} = {:s}\n".format(cplx.name, cplx.kernel_string))
 
    if condensed :
        # Print resting macrostates
        output_string("\n# Resting macrostates ({}) \n".format(
            len(enumerator.resting_macrostates)))
        for resting in natural_sort(enumerator.resting_macrostates):
            output_string("macrostate {:s} = [{}]\n".format(
                resting, ', '.join(map(str,resting.complexes))))

        # Print reactions
        output_string("\n# Condensed reactions ({}) \n".format(len(enumerator.condensed_reactions)))
        for rxn in natural_sort(enumerator.condensed_reactions):
            output_string("reaction {:s}\n".format(rxn.full_string(molarity, time)))

    if detailed :
        # Print transient complexes
        output_string("\n# Transient complexes ({}) \n".format(len(enumerator.transient_complexes)))
        for cplx in natural_sort(enumerator.transient_complexes):
            if cplx._concentration :
                output_string("{:s} = {:s} @{}\n".format(cplx.name, cplx.kernel_string, 
                    format_units(cplx, unit = molarity)))
            else:
                output_string("{:s} = {:s}\n".format(cplx.name, cplx.kernel_string))

        # Print reactions
        output_string("\n# Detailed reactions ({}) \n".format(len(enumerator.reactions)))
        for rxn in natural_sort(enumerator.reactions):
            output_string("reaction {:s}\n".format(rxn.full_string(molarity, time)))

    return ''.join(out)

def write_crn(enumerator, crn, condensed = False, molarity = 'M', time = 's'):
    crn.write("# File generated by peppercorn-{}\n".format(__version__))

    if condensed :
        # Print reactions
        crn.write("\n# Condensed reactions: concentration = {}, time = {}\n".format(
            molarity, time))
        for reaction in natural_sort(enumerator.condensed_reactions):
            units = [molarity] * (reaction.arity[0] - 1) + [time]
            rate = reaction.rateformat(units)
            crn.write("{} [k = {}]\n".format(reaction, rate.constant))
    else :
        # Print reactions
        crn.write("\n# Detailed reactions: concentration = {}, time = {}\n".format(
            molarity, time))
        for reaction in natural_sort(enumerator.reactions):
            units = [molarity] * (reaction.arity[0] - 1) + [time]
            rate = reaction.rateformat(units)
            crn.write("{} [k = {}]\n".format(reaction, rate.constant))

