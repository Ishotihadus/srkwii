# srkwii /serikawai:/

Personal toolbox for MATLAB.

[Serika Hakozaki](https://www.project-imas.com/wiki/Serika_Hakozaki).

## Installation

```matlab
addpath /path/to/srkwii
```

## Requirements

MATLAB, as new version as possible.

Some functions require

- Parallel Computing Toolbox
- Signal Processing Toolbox

## Usage

### SPTK

[SPTK](http://sp-tk.sourceforge.net/) wrapper.

- `srkwii.sptk.setsptkpath(path)`: set SPTK binary directory. If not set, `/opt/SPTK-3.11-double/bin` is used.
- `srkwii.sptk.getsptkpath`
    - `srkwii.sptk.getsptkpath()`
    - `srkwii.sptk.getsptkpath(command)`

## References

- WORLD analysis / synthesis
    - [M. Morise, F. Yokomori, and K. Ozawa, “WORLD: a vocoder-based high-quality speech synthesis system for real-time applications,” IEICE transactions on information and systems, vol. E99-D, no. 7, pp. 1877–1884, 2016.](https://search.ieice.org/bin/summary.php?id=e99-d_7_1877)
    - D4C: [M. Morise, “D4C, a band-aperiodicity estimator for high-quality speech synthesis,” Speech Communication, vol. 84, pp. 57–65, Nov. 2016.](https://www.sciencedirect.com/science/article/pii/S0167639316300413)
- Non-negative matrix factorization (NMF)
    - [D.D. Lee and H.S. Seung, “Algorithms for non-negative matrix factorization,” in Proc. Advances in Neural Information Processing Systems 13 (NIPS 2000), pp. 556–562, 2001.](https://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization)
    - Beta divergence: [M. Nakano, H. Kameoka, J.L. Roux, Y. Kitano, N. Ono, and S. Sagayama, “Convergence-guaranteed multiplicative algorithms for nonnegative matrix factorization with β-divergence”, in Proc. 2010 IEEE international workshop on machine learning for signal processing, pp. 283–288, 2010.](https://ieeexplore.ieee.org/abstract/document/5589233)
    - L1/2 quasi-norm sparsity penalty: [C. Joder, F. Weninger, D. Virette, and B. Schuller, “A comparative study on sparsity penalties for NMF-based speech separation: beyond Lp-norms,” in Proc. 2013 IEEE international conference on acoustics, speech and signal processing (ICASSP 2013), pp. 858–862, 2013.](https://ieeexplore.ieee.org/abstract/document/6637770)
    - Sparse KL divergence: [A. Cichocki, R. Zdunek, and S. Amari, “New algorithms for non-negative matrix factorization in applications to blind source separation,” in Proc. 2006 IEEE International conference on acoustics speech and signal processing (ICASSP 2006), pp. V-621–624, 2006.](https://ieeexplore.ieee.org/abstract/document/1661352)
    - Cepstral regularization / Log-euclidean distance: [H. Kameoka, T. Higuchi, M. Tanaka, and L. Li, “Nonnegative matrix factorization with basis clustering using cepstral distance regularization,” IEEE/ACM transactions on audio, speech, and language processing, vol. 26, no. 6, 2018.](https://ieeexplore.ieee.org/abstract/document/8264769)

## License

This software is distributed under MIT License.
