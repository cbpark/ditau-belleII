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
