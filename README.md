ditau-belleII
=============

## Dependencies

* [NLopt](https://nlopt.readthedocs.io): See [NLopt installation](https://nlopt.readthedocs.io/en/latest/NLopt_Installation/). For some Linux distributions, it can be installed by using the system package manager. For example,

```
# For Arch Linux/Manjaro:
sudo pacman -S nlopt

# For CentOS/Fedora:
sudo dnf install NLopt-devel

# For Debian/Ubuntu:
sudo apt-get install libnlopt-cxx-dev

# For openSUSE:
sudo zypper install nlopt
```

For Ubuntu prior to 20.04 LTS (Focal Fossa), install `libnlopt-dev`. (Check whether there exists `/usr/include/nlopt.hpp`.) In CentOS, the [EPEL](https://fedoraproject.org/wiki/EPEL) repository must be installed and enabled. In macOS, you can install NLopt using [Homebrew](https://brew.sh/).

``` no-hightlight
brew install nlopt
```

* [YAM2](https://github.com/cbpark/YAM2): See [How to build](https://github.com/cbpark/YAM2/blob/master/README.md). We assume that the installation path for YAM2 is `/usr/local`:

``` no-hightlight
cd YAM2
make
sudo DESTDIR=/usr/local make install
```

* [ROOT](https://root.cern/): See [Installing ROOT](https://root.cern/install/). Make sure that `root-config` is in the `PATH`:

``` no-highlight
$ root-config --version
6.22/08
```

## Building and running

Before building, see [`Makefile`](./Makefile) to check whether the paths are all correct. Then, run

``` no-hightlight
make
```

If the build is successful, an executable `bditau` will be created in the `bin` directory. The output will be stored in a ntuple (`vars`).

``` no-highlight
$ ./bin/bditau
usage: ./bin/bditau <event.root> <output.root> [mInvisible]
  <event.root>: input ntuple file (required).
  <output.root>: output file to store the result (required).
  [mInvisible]: the input mass for invisible particles (optional, default = 0)

$ ./bin/bditau event.root output.root
bditau: the input file is event.root
bditau: the name of the input tree is tau3x1
bditau: the invisible mass is 0
bditau: processed 10 events.
bditau: the output is stored in output.root

$ echo 'TFile f("output.root"); tree = (TTree *) f.Get("vars"); tree->Print()' | root
   ------------------------------------------------------------------
  | Welcome to ROOT 6.22/08                        https://root.cern |
  | (c) 1995-2020, The ROOT Team; conception: R. Brun, F. Rademakers |
  | Built for linuxx8664gcc on Mar 10 2021, 14:20:04                 |
  | From tags/v6-22-08@v6-22-08                                      |
  | Try '.help', '.demo', '.license', '.credits', '.quit'/'.q'       |
   ------------------------------------------------------------------

******************************************************************************
*Tree    :vars      : event variables                                        *
*Entries :       10 : Total =            3758 bytes  File  Size =        988 *
*        :          : Tree compression factor =   1.00                       *
******************************************************************************
*Br    0 :e_miss    : Float_t                                                *
*Entries :       10 : Total  Size=        689 bytes  One basket in memory    *
*Baskets :        0 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br    1 :m_recoil  : Float_t                                                *
*Entries :       10 : Total  Size=        701 bytes  One basket in memory    *
*Baskets :        0 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br    2 :m2        : Float_t                                                *
*Entries :       10 : Total  Size=        665 bytes  One basket in memory    *
*Baskets :        0 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br    3 :xi_p      : Float_t                                                *
*Entries :       10 : Total  Size=        677 bytes  One basket in memory    *
*Baskets :        0 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br    4 :xi_k      : Float_t                                                *
*Entries :       10 : Total  Size=        677 bytes  One basket in memory    *
*Baskets :        0 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
```

## Debugging

Turn on the `DEBUG` symbol:

``` no-highlight
$ make clean && CXXFLAGS=-DDEBUG make
```
