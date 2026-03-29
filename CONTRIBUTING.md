# Contributing to STMotif

Thank you for considering contributing to STMotif! Here are some
guidelines to help you get started.

## Reporting bugs

Please open an [issue](https://github.com/heraldoborges/STMotif/issues)
with:

- A minimal reproducible example
- The output of `sessionInfo()`
- A description of what you expected vs. what happened

## Suggesting features

Feature requests are welcome. Please open an issue describing:

- The problem you're trying to solve
- Your proposed solution (if any)
- Relevant use cases

## Pull requests

1. Fork the repository and create a branch from `main`.
2. Run `devtools::document()` if you changed any roxygen2 comments.
3. Run `devtools::test()` and make sure all tests pass.
4. Run `devtools::check()` and make sure there are no errors or warnings.
5. Update `NEWS.md` with a description of your changes.
6. Submit the pull request.

## Code style

- Use `<-` for assignment.
- Use `TRUE`/`FALSE`, never `T`/`F`.
- Qualify external package calls with `pkg::fun()` or use `@importFrom`.
- Document all exported functions with roxygen2.
- Add tests for new functionality.

## Code of Conduct

Please note that this project follows a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating,
you agree to abide by its terms.
