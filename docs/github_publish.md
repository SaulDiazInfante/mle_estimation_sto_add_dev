# GitHub Publish Guide

This project is now structured for GitHub-based development, CI, and release automation.

## Local git setup

Run these commands from the repository root:

```bash
git init -b main
make setup-git-hooks
git add .
git commit -m "Initial project import"
```

## Create the GitHub repository

Create an empty repository in GitHub's web UI, then connect it locally:

```bash
git remote add origin https://github.com/<your-user-or-org>/<repo-name>.git
git push -u origin main
```

If you prefer SSH:

```bash
git remote add origin git@github.com:<your-user-or-org>/<repo-name>.git
git push -u origin main
```

## Recommended GitHub settings

- Set the default branch to `main`.
- Protect `main` and require the `CI / build-and-test` status check before merge.
- Enable GitHub Pages and set the source to GitHub Actions so the `Docs` workflow can publish the Doxygen site.
- If the repository is private, make sure your GitHub plan supports Pages for private repositories. GitHub Free supports Pages for public repositories only.
- Require pull requests for changes to `main`.
- Require at least one approving review if you are collaborating with others.
- Enable GitHub Actions for the repository.
- Enable Dependabot version updates if you want automated maintenance for Actions.

## Release workflow

After the repository is pushed, create a version tag to publish a release artifact:

```bash
git tag v0.1.0
git push origin v0.1.0
```

That tag triggers the release workflow, which builds the project, runs tests, packages the Linux executable, and attaches it to a GitHub Release.

## TDD workflow on GitHub

- Add or update tests in the same branch as the code change.
- Open the pull request only after `make test` passes locally.
- Use the smoke and unit test output in CI as the merge gate.
- Keep issue descriptions focused on observable behavior so tests can be written first.

## Large file policy

- Keep generated outputs and local datasets out of git.
- Run `make check-large-files` before large refactors or imports.
- The local pre-commit hook and CI both enforce the tracked-file size policy.

## One open decision

Choose and add a project license before publishing if you want the repository to be clearly reusable by others.
