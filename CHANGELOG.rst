# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.2.1] - 2023-06-29

### Changed
- More docstring about the code from `pymatgen` code for the `POTCAR` metadata parser.
- Changed from `node` 12 to 16 for the Github Actions.

## [3.2.0] - 2023-06-05

### Added
- Added a way to parse generalized k-point grid, where the generators are typically used.
- `POTCAR` metadata parser from `pymatgen` to not having to install `pymatgen` as dependency.

### Changed
- Bugfix in the partial density of states parsing from `DOSCAR`.
- Improved detection of a successful VASP start (not full execution).
- Moved to a `toml` description of the package.
- Changed contact info for maintainer.
- Dependencies of the `tests` and `pre-commit` extras follow versions of `aiida-vasp` and `aiida-core`. Except for this release where we use a more recent `pylint>=2.15`.

### Removed
- Future and past dependency. Old backwards compatibility for Python 2 that was dormant in the code.

## [3.1.0] - 2022-05-27

### Added
- Added NBANDS to `run_status`.
- Added entry `bandocc` to the stream parser to be notified if the topmost band is occupied.

## [3.0.0] - 2022-05-24

### Added
- Added CHGCAR parser which dumps data in numpy arrays.
- Added EIGENVAL and DOSCAR parser.

### Changed
- Removed symmetry output. What remains is the number of space group operations, original cell type and symmetrized cell type.
- Fixed parsing of wildcard or N/A containing entries.
- Various other bugfixes.

## [2.0.1] - 2021-01-23

### Changed
- Fixed the check for the truncated xml file that contained a bug.

## [2.0.0] - 2021-01-04

### Added
- Posibility to return multiple total energy types.
- Posibility to parse the total energies for each electronic step. Since the total energy array for ionic and electronic steps is staggered,, we flatten it and supply an additional array with key `electronic_steps` where each entry indicates how many electronic steps was performed for each ionic step. Using this, the staggered array can be rebuilt.

### Changed
- The return of `get_energies` is now a dict in order to be able to return multiple energy types.
- `final` key for the ionic steps was changed to `last`.
- For static runs, there is now only one ionic entry that is returned.

## [1.2.0] - 2020-12-29

### Added
- Included additional regex errors for the stream parser
- Posibility to write cartesian coordinates to the POSCAR
- Added parsing of additional total energies, all electronic steps and the final which should be similar to the final electronic step without corrections.
- Parsing of number of electronic steps.
- Parsing of timing data.
- Parsing of VASP version.

## [1.1.2] - 2020-10-28

### Changed
- Bugfix in the history flag for the stream parser.

### Added
- Property function for the stream parser to return a bool if the xml file was truncated.

## [1.1.1] - 2020-09-07

### Changed
- Fixed docstring

## [1.1.0] - 2020-08-25

### Added
- Parser for the standard stream of VASP.

## [1.0.0] - 2020-06-10

Considered the first stable release.

### Added
- Magnetization from [@JPchico](https://github.com/JPchico)
- Enabled GitHub Actions
