Examples
========

Runnable scripts demonstrating the conceptual breakthroughs behind
**tbee**, from Bloch's band theory through the Thouless quantum pump --
see :doc:`/history` for the full chronology each script illustrates.

See also the narrative walkthrough: :doc:`/tutorial`.

Each script in this gallery is self-contained and can be run directly with
``python examples/<section>/<script>.py``. Every script also carries an
RST module docstring as its title/description and uses ``# %%`` markers
to split narrative text from code, which is exactly what Sphinx-Gallery
renders into the pages below -- the script *is* the source of truth for
what you see, not a copy of it. Every numeric claim a script's narrative
makes is checked with an ``assert`` right there in the code -- nothing is
asserted in the docs that isn't also verified in code.

Sections
--------

- **tight_binding** -- the generic real-space/reciprocal-space
  machinery: a graphene flake and its Bloch band structure, including
  Wallace's 1947 linear (Dirac) dispersion near the K point.
- **magnetic_field** -- the Peierls substitution: an Aharonov-Bohm ring's
  flux-periodic spectrum, and the fractal Hofstadter butterfly swept
  continuously in flux.
- **disorder** -- Anderson localization: the Inverse Participation Ratio
  vs. disorder strength, and an extended state next to a localized one.
- **topology** -- topological band theory: the SSH model's bulk-boundary
  correspondence, the Haldane model's Chern-number phase transition and
  Berry curvature, zigzag graphene and Kane-Mele helical edge states on
  ribbons, and the Thouless quantum pump's quantized charge transport.
- **flat_bands** -- the kagome and Lieb lattices' exactly flat bands, the
  natural home for strong-correlation physics via Lieb's theorem.
- **dynamics** -- real-time wavepacket propagation: Bloch oscillations and
  the Wannier-Stark ladder under a uniform tilt.
