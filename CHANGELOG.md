# ChangeLog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- New ROI detection using k-means with a larger number of clusters than 2. This allows a finer detection of the sample-related region.

## [1.0.4.1] - 2018-10-08
### Removed
- Removed dependency from 'autothresholdr' package. Now Otsu is performed using
  the function threshold(x, 'auto') from 'imager'.

## [1.0.4] - 2018-10-06
### Fixed
- Fixed a bug in the function .match.mz.array.
- Improved the comments in the function .match.mz.array.
