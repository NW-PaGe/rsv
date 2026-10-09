# CHANGELOG

We use this CHANGELOG to document breaking changes, bug fixes, and config value changes that affect the WA build. Changes to the Nexstrain build as a whole are logged by the Nextstrain team in the repo root: CHANGELOG.md.

## 2026

* 09 Oct 2026: Updated config to use augur subsampling instead of augur filter in response to edits to upstream build. Under this update, the sampling scheme has modest changes, making use of 4 sample groups instead of 2:
  * The prior approach created 2 sample subsets (defined via separate rules): 1 WA focused & 1 global. The WA subset was always filtered at 12Y, while the global subset was filtered based on the build resolution defined in the config.
  * This keeps the WA vs global sampling approach but also adds in the time-based subsampling approach used in the primary Nextstrain build
    * The time-based sampling includes a recent set (6Y-present) and a set of A.D- (subtype A) or B.D- (subtype B) only samples from 6Y - 12Y ago. With this change, uncommon and extinct historic clades will now no longer show up in the historic data in the WA-focused tree.

* 07 Aug 2026: Increased min_length threshold for genome from 10,000 to 14,5000. Results in:
  *  RSV A: loss of 52 seqeunces, which includes a loss of 35 WA sequences
  *  RSV B: forgot to log

* 06 Aug 2026: Changed divergence units from mutations-per-site to mutations.