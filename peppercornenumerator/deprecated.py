
class oPepperComplex(ComplexA):
    """
    Peppercorn complex object. 

    Overwrites some functions with new names, adds some convenient stuff..
    """

    @property
    def concentration(self):
        if self._concentration is not None:
            return self._concentration
        return None

    @concentration.setter
    def concentration(self, trip):
        if trip is None:
            self._concentration = None
        else:
            (mode, value, unit) = trip
            assert isinstance(value, (int, float))
            self._concentration = PepperComplex.CONCENTRATION(mode, value, unit)

    def concentrationformat(self, out):
        mod = self._concentration.mode
        val = self._concentration.value
        uni = self._concentration.unit
        val = convert_units(val, uni, out) 
        return PepperComplex.CONCENTRATION(mod, val, out)

    def full_string(self):
        return "Complex(%s): %s %s" % (
            self.name, str(self.strands), str(self.structure))

    def triple(self, *loc):
        # overwrite standard func
        return (self.get_domain(loc), self.get_paired_loc(loc), loc)

    @property
    def pk_domains(self):
        pd = []
        for (x,y) in self.enclosed_domains:
            pd.append((self.get_domain((x,y)), x, y))
        return pd

class oPepperReaction(ReactionS):
    @property
    def reverse_reaction(self):
        return self._reverse_reaction

    @reverse_reaction.setter
    def reverse_reaction(self, rxn):
        def rev_rtype(rtype, arity):
            """Returns the reaction type of a corresponding reverse reaction. """
            if rtype == 'open' and arity == (1,1):
                return 'bind11'
            elif rtype == 'open' and arity == (1,2):
                return 'bind21'
            elif rtype == 'bind11' or rtype == 'bind21':
                return 'open'
            else:
                return rtype

        if rxn is not False:
            assert rxn.rtype == rev_rtype(self.rtype, self.arity)
        self._reverse_reaction = rxn

    def full_string(self, molarity='M', time='s'):
        """Prints the reaction in PIL format.
        Reaction objects *always* specify rate in /M and /s.  """

        if self.rate :
            newunits = [molarity] * (self.arity[0] - 1) + [time]
            newrate = self.rateformat(newunits)
            rate = newrate.constant
            assert newunits == newrate.units
            units = ''.join(map('/{}'.format, newrate.units))
        else :
            rate = float('nan')
            units = ''

        if self.rtype :
            return '[{:14s} = {:12g} {:4s} ] {} -> {}'.format(self.rtype, rate, units,
                    " + ".join(map(str, self.reactants)), " + ".join(map(str, self.products)))
        else :
            return '[{:12g} {:4s} ] {} -> {}'.format(rate, units,
                    " + ".join(map(str, self.reactants)), " + ".join(map(str, self.products)))

