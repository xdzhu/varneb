# E089: same-reference, matched-strain PO-to-T comparison

Core conclusion: the two ordinary-converged training paths retain similar
T-triplet displacement coverage at their sampled peaks, but their atomic
gradient coverage differs substantially under the registered strain change;
compact geometry alone is not an energetic-response or accuracy certificate.

Archetype: quantitative grid. Backend: existing Python/Matplotlib workflow.
Target: evidence-bound JCTC working manuscript, not a final journal-spec audit.
Final size: 183 by 202 mm; editable PDF/SVG plus PNG/TIFF and scalar CSV.

Panel map: (a,b) epsilon0 and +0.01 energy profiles, each relative to its own
matched relaxed PO+ endpoint; (c,d) common unstrained-T mass-chart displacement
fractions; (e,f) dual atomic-gradient fractions in that same reference chart.
The columns use individual normalized generalized arcs, not pointwise matched
geometries. Vertical guides mark highest sampled images, not certified TSs.

Evidence: actual terminal receipts, nine portable native SCFs per path,
unchanged original T/Cmma references, ranks3/4/4 and all four registered Cmma
frames. The two PO+ SCFs differ physically with strain; their energy zeros
must not be advertised as an identical cached endpoint. Keep absolute E,
gradient norms, omitted norms, raw stresses and Green strain in scalar sources.

Statistics: eighteen image records including reused endpoints; zero independent
replicates, no sampling distribution or confidence intervals. Four Cmma frames
are retained registrations, not statistical uncertainty. No fitted MEP or
interpolation. Zero-denominator fractions remain undefined and are not drawn.

Style: unbold larger (a) labels aligned across columns; aligned x/y axes,
visible top/right spines and ticks; no panel titles. Larger axis/tick/legend
text, white semi-transparent framed legends outside curves and panel labels.
Data provenance and image-integrity checks precede export. Review risks:
gradient percentages at small absolute norm, unmatched rank, incomplete G2,
ordinary convergence versus true stationarity, trained strain contrast versus
independent prediction. None is hidden by the figure. No new DFT, changed
input, G3 selection, holdout access or TS/model-accuracy claim.
