# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
