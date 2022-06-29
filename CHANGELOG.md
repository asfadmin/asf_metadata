# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [PEP 440](https://www.python.org/dev/peps/pep-0440/)
and uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!--

------
## Example template!!

## [version](https://github.com/asfadmin/asf_metadata/compare/vOLD...vNEW)

### Added:
-

### Changed:
-

### Fixed:
- 

### Removed:
-
------

-->

------

## [0.0.5](https://github.com/asfadmin/asf_metadata/compare/v0.0.4...v0.0.5)

### Fixed

- Include `templates` dir to package, by setting `include_package_data=True` in `setup.py`.

------

## [0.0.4](https://github.com/asfadmin/asf_metadata/compare/v0.0.3...v0.0.4)

### Fixed

- Let users import `asf_metadata.metadata` by adding `__init__.py` to that directory.

------

## [0.0.3](https://github.com/asfadmin/asf_metadata/compare/v0.0.2...v0.0.3)

### Changed

- Pipeline changes for pypi.

------

## [0.0.2](https://github.com/asfadmin/asf_metadata/compare/v0.0.1...v0.0.2)

### Added

- Sentinel Burst map logic.
- Added pypi release pipeline.

### Fixed

- Moved minimal package requirements to `setup.py`, from `requirements.in`.
- Moved `generate_iso_template`, `iso_metadata_xml2json`, and `generate_iso_metadata` methods to import directly under `asf_metadata`, to import easier.

------

## [0.0.1](https://github.com/asfadmin/asf_metadata/compare/v0.0.0...v0.0.1)

### Added

- Initial Release! First deployment to PyPi.
