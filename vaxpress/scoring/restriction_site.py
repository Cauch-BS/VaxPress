#
# VaxPress
#
# Copyright 2023 Seoul National University
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# “Software”), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


from . import ScoringFunction


class RestrictionSiteFitness(ScoringFunction):
    # Prefix in command line arguments, e.g. "--resite-weight"
    name = "resite"
    description = "Restriction Site"
    priority = 100

    # Definitions for command line arguments.
    arguments = [
        (
            "weight",
            dict(
                type=float,
                default=0.0,  # penalize for the violations
                help="scoring weight for the presence of the restriction site "
                "(default: -100.0)",
            ),
        ),
        (
            "restriction-sites",
            dict(
                default=None,
                type=str,
                help="restriction site sequence containing only A, C, G, T or "
                "U (default: None)",
            ),
        ),
        (
            "original-sites",
            dict(
                default=0.0,
                type=str,
                help="number of restriction sites in the original sequence ",
            ),
        ),
    ]

    # Argument "_length_cds" is always passed to the constructor. Leave it as
    # it is even if you don't use it.
    def __init__(self, weight, restriction_sites, original_sites, _length_cds):
        self.weight = weight
        self.restriction_sites = (
            restriction_sites.split(",") if restriction_sites else []
        )
        self.original_sites = original_sites
        if not self.restriction_sites:
            self.weight = 0.0

        for i, site in enumerate(self.restriction_sites):
            if not site:
                raise ValueError("restriction site cannot be empty")
            new_site = site.upper().replace("T", "U")
            if set(new_site) - set("ACGU"):
                raise ValueError("restriction site contains invalid characters")
            self.restriction_sites[i] = new_site

    def score(self, seqs):
        metrics = []
        scores = []
        for seq in seqs:
            found = -self.original_sites
            for site in self.restriction_sites:
                found += seq.count(site)
            metrics.append(found)
            scores.append(found * self.weight)
        return {"resite": scores}, {"resite": metrics}
