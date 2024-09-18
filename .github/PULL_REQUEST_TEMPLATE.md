# [BUG/FEAT/DOC/*etc.*] **Description**

> **Please review the [contributing guidelines](.github/CONTRIBUTING.md) before continuing.**

> 1-2 sentence summary of the change and (if applicable) which issue this relates to.

## Justification

> Often, the justification of bugs is that "they exist", so if your pull request fixes a single bug in a single commit, this is not necessary.

> Feel free to create any sections necessary to fully describe the pull request.

## Changes to GITM outputs

> Detail if GITM outputs will be changed and, if so, how. For example, would all runs be changed or only those with specific settings enabled?

---

## Notes

>> This section should be removed before submitting the pull request.

- You must the changes made with this pull request. Please detail the configuration you tested with and whether this change affects the outputs of GITM.
- All changes to code will be tested in the following order. If any step fails, the pull request cannot be accepted:
    1. Code is checked against the [code standards](CONTRIBUTING.md#formatting-code). If `fprettify` recommends any changes the test fails, you you should run the format check before submitting.
    2. The GITM manual is built.
    3. GITM is compiled with `gfortran10` and then run using the UAM.in files within [`srcTests/`](../srcTests/). This ensures no changes are made which unintentionally break something else.
- Pull requests should be made against the `develop` branch of GITMCode/GITM, not `main`.

In summary, ensure you have:

> - Tested your changes & compared the outputs to the most recent version on the `develop` branch. Detail here if any outputs differ & why this is necessary.
> - Updated any relevant `UAM.in` files if you are adding an optional feature
> - Updated the documentation, if necessary
> - Commented you code in hard-to-understand areas
> - Clean up unnecessary comments, debugging steps, etc. within your code
> - Ensure your code is formatted correctly. Use the included Python script for assistance.
> - Pull requests should target the develop branch. Not main/master!
