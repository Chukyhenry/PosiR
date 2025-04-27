*   This is a resubmission of version 0.1.1.

## R CMD check results

*   Duration: [Your local check duration, e.g., 14.8s]
*   0 errors ✔ | 0 warnings ✔ | 1 NOTE ✖ (local check)

## NOTES

### Local Check NOTE:
*   checking for future file timestamps ... NOTE
    *   unable to verify current time
*   Comment: This seems related to the checking environment and not the package itself.

### Win-builder NOTE:
*   checking CRAN incoming feasibility ... NOTE
    *   Possibly misspelled words in DESCRIPTION: Kuchibhotla (11:37), PoSI (9:59), al (11:52), et (11:49)
*   Comment: These words are correctly spelled:
    *   'Kuchibhotla' is a proper name (author).
    *   'PoSI' is the acronym for Post-Selection Inference but has been replaced.
    *   'et' and 'al' are standard citation abbreviations ('et al.') but has been replaced.
*   Added the required LICENSE file to resolve the ‘Invalid license file pointers’ WARNING.
