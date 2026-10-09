# Uploading the NSC browser to GitHub

Requested destination: https://github.com/AustinJin19/NSC_cross_species_genome_browser

The destination repository is private. Browser authentication is available; command-line Git authentication must be configured separately if using Git locally.

The source archive contains browser code, tests, requirements, and documentation. It excludes experimental data, genome files, generated indexes/databases, installed dependencies, and local outputs. The live working project is not a Git repository.

## Recommended workflow

1. Obtain the existing repository URL and clone it into a separate directory, preserving the live browser workspace.
2. Create a branch such as `add-nsc-browser`.
3. Extract the source archive into a temporary directory and compare its contents against the clone. Copy the source/documentation into the agreed location. Do not overwrite an existing README, license, or ignore rules without reviewing the differences.
4. Review `git status`, `git diff`, and staged files. The root `.gitignore` is a strict allowlist for this standalone browser layout; merge it carefully if the destination repository has other code.
5. Commit and push the branch, then open a pull request against the repository's default branch.

Example commands after copying and reviewing the files:

```bash
git switch -c add-nsc-browser
git add README.md docs cross_species_visual_tool/web_app .gitignore
git diff --cached --stat
git diff --cached --check
git commit -m "Add NSC comparative genome browser"
git push -u origin add-nsc-browser
```

Select a software license explicitly; none has been assigned automatically. Review data redistribution and branding permissions independently. Use a data repository or documented external downloads for genomic datasets rather than including the local 14 GB Hi-C file in Git.

GitHub stores this code; it does not run the Python backend. A clone requires the separately supplied data layout described in `CODE_SUMMARY.md`.
