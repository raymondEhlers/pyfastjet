# Notes on jet finding examples

I included a few examples with different levels of complexity. Note that these samples are biased towards my
background (heavy ions, beginning w/ EIC), so my use patterns probably aren't entirely representative. But
they're probably a good start.

And I'm of course happy to answer any questions if things are unclear.

## Examples

- Basic jet finding using the anti-kt algorithm.
- Jet background subtraction (using `useConstituentSubtraction = false`)
- A more complicated example, subtracting background contributions from the jet constituents (using
`useConstituentSubtraction = true`)

Note that for the background subtraction examples to be meaningful, you should pass in a substantial number of
particles. I've tested the code, and it seems to run fine, but I haven't tested the outputs yet in any detail.
If that's helpful, I can revisit it later.

## Stray notes

### Naming of "fastjet-core"

I see that you include the `fastjet` repo as `fastjet-core` in your repo. Just as a heads up in case you're
not aware: there is already a `fjcore` package from the fastjet authors, which includes some basic fastjet
classes, but is missing some functionality. In principle, it's all fine, but it may be confusing.

### Centauro jet finding algorithm for the EIC (and others)

Adding support for the Centauro jet finding algorithm would be helpful given that EIC activities are ramping
up quickly this year. Assuming fastjet contrib 1.045 is available (this is the current version), it's as
simple as:

```c++
fastjet::contrib::CentauroPlugin centauroPlugin(jetR);
fastjet::JetDefinition jetFinder(centauroPlugin);
```

Note that I _think_ this also requires a boost into a new reference frame, but this can be taken care of
elsewhere (or via additional arguments to the Centauro plugin, but I wouldn't worry overly much about
that for now). This is low priority overall, but I think it should be fairly easy to support, at least in a
minimal capacity.


