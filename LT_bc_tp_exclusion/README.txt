These are some old scripts used to clean up the lineage tracking data.
There were clear cases of sequencing library cross-contamination in 
some of this data, so I went through and looked at all the different
environment-specific barcodes in each dataset. Ideally these barcodes
would have told us exactly what to exclude in every case, but it turns
out there were some cases where unexpected environment barcodes really
were present in evolutions where we wouldn't have expected them (possibly
due to some cross contamination during library creation?), so the problem
was a little more tricky. Here, I manually looked through trajectories
and picked out barcodes that had unrealistic (e.g. jumping up to high freq.
for a single timepoint) trajectories. These barcodes are excluded here and
several timepoints with particularly bad issues are also excluded. This
is all a bit manual and ad-hoc but it provides a much more realistic 
picture of lineage dynamics.
