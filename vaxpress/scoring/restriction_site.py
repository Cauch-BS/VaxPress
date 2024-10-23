# Example code that optimizes a sequence to have a specific restriction enzyme
# recognition site.

from vaxpress.scoring import ScoringFunction


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
                default=-100.0,  # penalize for the violations
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
    ]

    # Argument "_length_cds" is always passed to the constructor. Leave it as
    # it is even if you don't use it.
    def __init__(self, weight, restriction_sites, _length_cds):
        self.weight = weight
        self.restriction_sites = restriction_sites.split(",")
        if self.restriction_sites is None:
            raise EOFError  # disable this scoring function

        for i, site in enumerate(self.restriction_sites):
            if not site:
                raise ValueError("restriction site cannot be empty")
            new_site = site.upper().replace("T", "U")
            if set(new_site) - set("ACGU"):
                raise ValueError("restriction site contains invalid characters")
            self.restriction_sites[i] = new_site

    def score(self, seqs):
        has_site = []
        scores = []
        for seq in seqs:
            found = -len(self.restriction_sites)
            for site in self.restriction_sites:
                if site in seq:
                    found += 1
            has_site.append(found)
            scores.append(found * self.weight)
        return {"resite": scores}, {"resite": has_site}
